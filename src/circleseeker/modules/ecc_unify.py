"""Merge confirmed and inferred eccDNA tables with redundancy detection.

This module handles:
1. Detection of overlaps between confirmed and inferred eccDNA
2. Removal of redundant inferred entries
3. Merging all tables into a unified output
4. Renumbering of final eccDNA IDs
"""

from __future__ import annotations

import re
import json
from pathlib import Path
from typing import Any, Optional

import pandas as pd
from circleseeker.utils.logging import get_logger


def parse_region(r: str) -> tuple[str, int, int]:
    """Parse region string into chromosome, start, end.

    Handles both single regions (chr1:100-200) and first segment of chimeric regions.
    """
    # If it's a chimeric region string, take the first segment
    if ";" in str(r):
        r = str(r).split(";")[0]

    m = re.match(r"([^:]+):(\d+)-(\d+)", str(r))
    if not m:
        raise ValueError(f"Bad region: {r}")
    return m.group(1), int(m.group(2)), int(m.group(3))


def parse_chimeric_regions(regions_str: str) -> list[tuple[str, int, int]]:
    """Parse chimeric region string into list of (chr, start, end) tuples.

    Example: "chr1:100-200;chr2:300-400" -> [("chr1", 100, 200), ("chr2", 300, 400)]
    """
    segments = []
    for segment in str(regions_str).split(";"):
        try:
            chr_name, coords = segment.strip().split(":")
            start, end = coords.split("-")
            segments.append((chr_name, int(start), int(end)))
        except (ValueError, AttributeError, IndexError):
            continue
    return segments


def reciprocal_overlap_ok(
    a_start: int, a_end: int, b_start: int, b_end: int, thr: float = 0.99, tol: int = 10
) -> bool:
    """Check if two regions have reciprocal overlap above threshold.

    Args:
        a_start, a_end: First region coordinates
        b_start, b_end: Second region coordinates
        thr: Minimum reciprocal overlap fraction (default 0.99)
        tol: Tolerance in bp for boundary matching (default 10)

    Returns:
        True if regions have sufficient reciprocal overlap
    """
    s = max(a_start, b_start)
    e = min(a_end, b_end)
    if s >= e:
        return False

    ov = e - s
    la = max(1, a_end - a_start)
    lb = max(1, b_end - b_start)

    # Calculate reciprocal overlap fractions
    ov_frac_a = ov / la
    ov_frac_b = ov / lb

    # Without tolerance, both fractions must exceed threshold
    if tol == 0:
        return (ov_frac_a >= thr) and (ov_frac_b >= thr)

    # With tolerance, check if coordinates are close enough
    # to be considered the same region
    start_diff = abs(a_start - b_start)
    end_diff = abs(a_end - b_end)

    # If boundaries are within tolerance, consider it a match
    if start_diff <= tol and end_diff <= tol:
        return True

    # Otherwise apply standard reciprocal overlap test
    return (ov_frac_a >= thr) and (ov_frac_b >= thr)


def build_chr_index(
    df: pd.DataFrame, type_filter: Optional[str] = None
) -> dict[str, list[tuple[int, int, int]]]:
    """Build chromosome-indexed structure for efficient overlap queries.

    Uses vectorized string operations for better performance with large DataFrames.

    Args:
        df: DataFrame with 'Regions' column
        type_filter: Optional eccDNA_type to filter for (e.g., 'UeccDNA')
    """
    idx: dict[str, list[tuple[int, int, int]]] = {}

    # Filter by type if specified
    if type_filter and "eccDNA_type" in df.columns:
        df = df[df["eccDNA_type"] == type_filter]

    if df.empty:
        return idx

    # Get regions column (try both cases)
    regions_col = "Regions" if "Regions" in df.columns else "regions"
    if regions_col not in df.columns:
        return idx

    # Reset index to ensure contiguous integer index for proper mapping
    df = df.reset_index(drop=True)

    # Get regions and handle chimeric (take first segment)
    regions = df[regions_col].astype(str).str.split(";").str[0]

    # Vectorized regex extraction: chr:start-end
    pattern = r"([^:]+):(\d+)-(\d+)"
    extracted = regions.str.extract(pattern)
    extracted.columns = ["chr", "start", "end"]

    # Filter valid rows (all three columns present)
    valid_mask = extracted.notna().all(axis=1)
    valid_df = extracted[valid_mask].copy()
    valid_df["start"] = pd.to_numeric(valid_df["start"], errors="coerce").astype("Int64")
    valid_df["end"] = pd.to_numeric(valid_df["end"], errors="coerce").astype("Int64")

    # Build index from valid entries
    for i, (chr_val, start, end) in valid_df.iterrows():
        if pd.notna(chr_val) and pd.notna(start) and pd.notna(end):
            idx.setdefault(chr_val, []).append((int(start), int(end), i))

    # Sort by start position for each chromosome
    for ch in idx:
        idx[ch].sort(key=lambda t: t[0])

    return idx


def find_redundant_simple(
    inferred_df: pd.DataFrame, confirmed_df: pd.DataFrame, thr: float = 0.99, tol: int = 10
) -> set[str]:
    """Find inferred simple eccDNA that overlap with confirmed UeccDNA.

    Args:
        inferred_df: DataFrame with inferred simple eccDNA
        confirmed_df: DataFrame with all confirmed eccDNA (will filter for UeccDNA)
        thr: Reciprocal overlap threshold
        tol: Coordinate tolerance

    Returns:
        Set of redundant eccDNA IDs
    """
    # Build index from confirmed UeccDNA regions only
    idx = build_chr_index(confirmed_df, type_filter="UeccDNA")
    redundant_ids = set()

    # Check each inferred entry
    for _, row in inferred_df.iterrows():
        # Try to get coordinates from various possible column names
        chr_val = row.get("chr", row.get("Chr", None))
        start_val = row.get("start0", row.get("Start0", None))
        end_val = row.get("end0", row.get("End0", None))

        if pd.notna(chr_val) and pd.notna(start_val) and pd.notna(end_val):
            ch = str(chr_val).strip()
            if not ch:
                continue
            try:
                s = int(start_val)
                e = int(end_val)
            except (TypeError, ValueError):
                continue
        else:
            # Try parsing from regions column if coordinates not available
            regions = row.get("regions", row.get("Regions", None))
            if regions and pd.notna(regions):
                try:
                    ch, s, e = parse_region(regions)
                except (ValueError, TypeError, IndexError):
                    continue
            else:
                continue

        # Check against confirmed entries on same chromosome
        cand = idx.get(ch, [])
        hit = False

        for ts, te, _ in cand:
            if ts > e + tol:
                break  # No more possible overlaps
            if te < s - tol:
                continue  # Not overlapping yet

            if reciprocal_overlap_ok(s, e, ts, te, thr, tol):
                hit = True
                break

        if hit:
            redundant_ids.add(row["eccDNA_id"])

    return redundant_ids


def find_redundant_chimeric(
    inferred_df: pd.DataFrame,
    confirmed_df: pd.DataFrame,
    thr: float = 0.99,
    tol: int = 10,
    method: str = "overlap",
) -> set[str]:
    """Find inferred chimeric eccDNA that match confirmed CeccDNA.

    Supports two matching strategies:
    - 'exact': Exact string matching (legacy method)
    - 'overlap': Segment-wise reciprocal overlap (recommended)

    Args:
        inferred_df: DataFrame with inferred chimeric eccDNA
        confirmed_df: DataFrame with all confirmed eccDNA (will filter for CeccDNA)
        thr: Reciprocal overlap threshold (default 0.99)
        tol: Coordinate tolerance in bp (default 10)
        method: Matching method - 'exact' or 'overlap' (default 'overlap')

    Returns:
        Set of redundant eccDNA IDs
    """
    if method == "exact":
        return find_redundant_chimeric_exact(inferred_df, confirmed_df)
    elif method == "overlap":
        return find_redundant_chimeric_overlap(inferred_df, confirmed_df, thr, tol)
    else:
        raise ValueError(f"Unknown method '{method}'. Use 'exact' or 'overlap'.")


def find_redundant_chimeric_exact(
    inferred_df: pd.DataFrame, confirmed_df: pd.DataFrame
) -> set[str]:
    """Find redundant chimeric eccDNA using exact string matching (legacy).

    This is the original implementation preserved for backward compatibility.
    """
    # Handle empty DataFrames
    if inferred_df.empty or confirmed_df.empty:
        return set()

    # Build set of confirmed CeccDNA region strings
    confirmed_c = confirmed_df[confirmed_df["eccDNA_type"] == "CeccDNA"]
    confirmed_regions = set()

    for _, row in confirmed_c.iterrows():
        regions = row.get("Regions", row.get("regions", ""))
        if regions and pd.notna(regions):
            normalized = regions.replace(" ", "")
            confirmed_regions.add(normalized)

    # Check each inferred chimeric entry
    redundant_ids = set()

    for ecc_id, group in inferred_df.groupby("eccDNA_id"):
        group = group.sort_values("seg_index")
        segments = []

        for _, row in group.iterrows():
            chr_val = row.get("chr", "")
            start_val = row.get("start0", "")
            end_val = row.get("end0", "")

            if chr_val and pd.notna(start_val) and pd.notna(end_val):
                segments.append(f"{chr_val}:{int(start_val)}-{int(end_val)}")

        if segments:
            inferred_region = ";".join(segments)
            if inferred_region in confirmed_regions:
                redundant_ids.add(ecc_id)

    return redundant_ids


def find_redundant_chimeric_overlap(
    inferred_df: pd.DataFrame, confirmed_df: pd.DataFrame, thr: float = 0.99, tol: int = 10
) -> set[str]:
    """Find redundant chimeric eccDNA using segment-wise reciprocal overlap.

    This method compares each segment individually using the same strategy
    as simple eccDNA overlap detection.
    """
    # Handle empty DataFrames
    if inferred_df.empty or confirmed_df.empty:
        return set()

    # Filter confirmed CeccDNA and parse their segments
    confirmed_c = confirmed_df[confirmed_df["eccDNA_type"] == "CeccDNA"]
    confirmed_segs_list = []

    for idx, row in confirmed_c.iterrows():
        regions = row.get("Regions", row.get("regions", ""))
        if regions and pd.notna(regions):
            segs = parse_chimeric_regions(regions)
            if segs:
                confirmed_segs_list.append((idx, segs))

    # Check each inferred chimeric eccDNA
    redundant_ids = set()

    for ecc_id, group in inferred_df.groupby("eccDNA_id"):
        group = group.sort_values("seg_index")
        inferred_segs = []

        for _, row in group.iterrows():
            chr_val = row.get("chr", "")
            start_val = row.get("start0", "")
            end_val = row.get("end0", "")

            if chr_val and pd.notna(start_val) and pd.notna(end_val):
                inferred_segs.append((str(chr_val), int(start_val), int(end_val)))

        if not inferred_segs:
            continue

        for conf_idx, conf_segs in confirmed_segs_list:
            if _segments_match(inferred_segs, conf_segs, thr, tol):
                redundant_ids.add(ecc_id)
                break

    return redundant_ids


def _segments_match(
    segs_a: list[tuple[str, int, int]], segs_b: list[tuple[str, int, int]], thr: float, tol: int
) -> bool:
    """Check if two lists of segments match using segment-wise reciprocal overlap."""
    if len(segs_a) != len(segs_b):
        return False

    if not segs_a:
        return False

    def ordered_match(a: list[tuple[str, int, int]], b: list[tuple[str, int, int]]) -> bool:
        for seg_a, seg_b in zip(a, b):
            chr_a, start_a, end_a = seg_a
            chr_b, start_b, end_b = seg_b

            if chr_a != chr_b:
                return False

            if not reciprocal_overlap_ok(start_a, end_a, start_b, end_b, thr, tol):
                return False

        return True

    if ordered_match(segs_a, segs_b):
        return True

    seg_count = len(segs_a)
    if seg_count < 2:
        return False

    for offset in range(1, seg_count):
        rotated = segs_b[offset:] + segs_b[:offset]
        if ordered_match(segs_a, rotated):
            return True

    return False


def prepare_inferred_simple(df: pd.DataFrame, redundant_ids: set[str]) -> pd.DataFrame:
    """Prepare inferred simple table for merging.

    Converts from inferred format to standard format and removes redundant entries.
    """
    # Filter out redundant entries
    filtered = df[~df["eccDNA_id"].isin(redundant_ids)].copy()

    if filtered.empty:
        return pd.DataFrame()

    # Build regions column from coordinates
    filtered["Regions"] = filtered.apply(
        lambda r: f"{r['chr']}:{int(r['start0'])}-{int(r['end0'])}", axis=1
    )

    # Map to standard columns
    result = pd.DataFrame(
        {
            "eccDNA_id": filtered["eccDNA_id"],
            "Regions": filtered["Regions"],
            "Strand": filtered.get("strand", "+"),
            "Length": filtered["length"],
            "eccDNA_type": "UeccDNA",  # Simple inferred maps to UeccDNA
            "State": "Inferred",
            "Seg_total": 1,
            "Hit_count": 1,
        }
    )

    return result


def prepare_inferred_chimeric(df: pd.DataFrame, redundant_ids: set[str]) -> pd.DataFrame:
    """Prepare inferred chimeric table for merging.

    Converts from segment-based format to unified format and removes redundant entries.
    """
    # Filter out redundant entries
    filtered = df[~df["eccDNA_id"].isin(redundant_ids)].copy()

    if filtered.empty:
        return pd.DataFrame()

    # Group by eccDNA_id and build unified entries
    result_rows = []

    for ecc_id, group in filtered.groupby("eccDNA_id"):
        # Sort by segment index
        group = group.sort_values("seg_index")

        # Build regions and strands
        regions = []
        strands = []

        for _, row in group.iterrows():
            regions.append(f"{row['chr']}:{int(row['start0'])}-{int(row['end0'])}")
            strands.append(row.get("strand", "+"))

        # Take first row for metadata
        first_row = group.iloc[0]

        result_rows.append(
            {
                "eccDNA_id": ecc_id,
                "Regions": ";".join(regions),
                "Strand": ";".join(strands),
                "Length": first_row["length"],
                "eccDNA_type": "CeccDNA",  # Chimeric inferred maps to CeccDNA
                "State": "Inferred",
                "Seg_total": first_row["seg_total"],
                "Hit_count": 1,
            }
        )

    return pd.DataFrame(result_rows)


def renumber_eccdna(df: pd.DataFrame) -> pd.DataFrame:
    """Continue numbering for inferred eccDNA after confirmed ones.

    Keeps original confirmed IDs, adds new sequential IDs for inferred.

    Args:
        df: Merged DataFrame with all eccDNA

    Returns:
        DataFrame with eccDNA_id (final) and original_id columns
    """
    df = df.copy()

    # Save original IDs
    df["original_id"] = df["eccDNA_id"]

    # Handle empty DataFrame
    if len(df) == 0:
        return df

    # Sort by type and state
    type_order = {"UeccDNA": 0, "MeccDNA": 1, "CeccDNA": 2}
    state_order = {"Confirmed": 0, "Inferred": 1}

    df["_type_order"] = df["eccDNA_type"].map(type_order)
    df["_state_order"] = df["State"].map(state_order)

    # Sort
    df = df.sort_values(["_type_order", "_state_order"]).reset_index(drop=True)

    # Assign IDs - keep confirmed, continue numbering for inferred
    final_ids = []

    for ecc_type in ["UeccDNA", "MeccDNA", "CeccDNA"]:
        pattern = re.compile(rf"^{re.escape(ecc_type)}(\d+)$")
        # Get confirmed entries - keep their original IDs
        confirmed_df = df[(df["eccDNA_type"] == ecc_type) & (df["State"] == "Confirmed")]

        valid_numbers = []
        for orig_id in confirmed_df["original_id"].astype(str):
            match = pattern.match(orig_id)
            if match:
                try:
                    valid_numbers.append(int(match.group(1)))
                except ValueError:
                    continue

        max_num = max(valid_numbers) if valid_numbers else 0
        next_num = max_num
        seen_numbers: set[int] = set()

        for _, row in confirmed_df.iterrows():
            orig_id = str(row["original_id"])
            match = pattern.match(orig_id)
            if match:
                num = int(match.group(1))
                if num not in seen_numbers:
                    final_ids.append(orig_id)
                    seen_numbers.add(num)
                    continue

            next_num += 1
            final_ids.append(f"{ecc_type}{next_num}")

        # Get inferred entries - continue numbering from confirmed
        inferred_df = df[(df["eccDNA_type"] == ecc_type) & (df["State"] == "Inferred")]
        if len(inferred_df) > 0:
            # Continue numbering for inferred
            for i in range(len(inferred_df)):
                next_num += 1
                final_ids.append(f"{ecc_type}{next_num}")

    # Set the final eccDNA_id
    df["eccDNA_id"] = final_ids

    # Clean up temporary columns
    df = df.drop(columns=["_type_order", "_state_order"])

    # Reorder columns - eccDNA_id first, original_id second, then rest
    cols = ["eccDNA_id", "original_id"] + [
        c for c in df.columns if c not in ["eccDNA_id", "original_id"]
    ]
    return df[cols]


def generate_overlap_stats_json(
    confirmed_df: pd.DataFrame,
    inferred_simple: pd.DataFrame | None,
    inferred_chimeric: pd.DataFrame | None,
    simple_redundant: set[str],
    chimeric_redundant: set[str],
    output_file: Path | str | None = None,
) -> dict:
    """Generate overlap statistics in JSON format for downstream processing.

    Args:
        confirmed_df: Confirmed eccDNA DataFrame
        inferred_simple: Inferred simple eccDNA DataFrame
        inferred_chimeric: Inferred chimeric eccDNA DataFrame
        simple_redundant: Set of redundant simple eccDNA IDs
        chimeric_redundant: Set of redundant chimeric eccDNA IDs
        output_file: Optional JSON file path to save statistics

    Returns:
        Dictionary with overlap statistics
    """
    import json

    stats: dict[str, Any] = {
        "confirmed": {
            "UeccDNA": len(confirmed_df[confirmed_df["eccDNA_type"] == "UeccDNA"]),
            "MeccDNA": len(confirmed_df[confirmed_df["eccDNA_type"] == "MeccDNA"]),
            "CeccDNA": len(confirmed_df[confirmed_df["eccDNA_type"] == "CeccDNA"]),
            "total": len(confirmed_df),
        },
        "inferred_simple": {
            "total": 0,
            "overlapping": 0,
            "overlap_percentage": 0.0,
            "non_redundant": 0,
            "non_redundant_percentage": 0.0,
            "redundant_ids": [],
        },
        "inferred_chimeric": {
            "total": 0,
            "overlapping": 0,
            "overlap_percentage": 0.0,
            "non_redundant": 0,
            "non_redundant_percentage": 0.0,
            "redundant_ids": [],
        },
        "summary": {
            "total_confirmed": len(confirmed_df),
            "total_inferred": 0,
            "total_overlapping": 0,
            "total_non_redundant": 0,
            "final_total": len(confirmed_df),
        },
    }

    # Process simple inferred
    if inferred_simple is not None and not inferred_simple.empty:
        total_simple = len(inferred_simple)
        redundant_simple = len(simple_redundant)
        non_redundant_simple = total_simple - redundant_simple

        stats["inferred_simple"] = {
            "total": total_simple,
            "overlapping": redundant_simple,
            "overlap_percentage": (
                round(redundant_simple / total_simple * 100, 2) if total_simple > 0 else 0
            ),
            "non_redundant": non_redundant_simple,
            "non_redundant_percentage": (
                round(non_redundant_simple / total_simple * 100, 2) if total_simple > 0 else 0
            ),
            "redundant_ids": sorted(list(simple_redundant)),
        }

    # Process chimeric inferred
    if inferred_chimeric is not None and not inferred_chimeric.empty:
        unique_chimeric = inferred_chimeric["eccDNA_id"].nunique()
        redundant_chimeric = len(chimeric_redundant)
        non_redundant_chimeric = unique_chimeric - redundant_chimeric

        stats["inferred_chimeric"] = {
            "total": unique_chimeric,
            "overlapping": redundant_chimeric,
            "overlap_percentage": (
                round(redundant_chimeric / unique_chimeric * 100, 2) if unique_chimeric > 0 else 0
            ),
            "non_redundant": non_redundant_chimeric,
            "non_redundant_percentage": (
                round(non_redundant_chimeric / unique_chimeric * 100, 2)
                if unique_chimeric > 0
                else 0
            ),
            "redundant_ids": sorted(list(chimeric_redundant)),
        }

    # Update summary
    total_inferred = stats["inferred_simple"]["total"] + stats["inferred_chimeric"]["total"]
    total_overlapping = (
        stats["inferred_simple"]["overlapping"] + stats["inferred_chimeric"]["overlapping"]
    )
    total_non_redundant = (
        stats["inferred_simple"]["non_redundant"] + stats["inferred_chimeric"]["non_redundant"]
    )

    stats["summary"] = {
        "total_confirmed": len(confirmed_df),
        "total_inferred": total_inferred,
        "total_overlapping": total_overlapping,
        "total_non_redundant": total_non_redundant,
        "final_total": len(confirmed_df) + total_non_redundant,
    }

    # Add formatted summary strings for easy display
    simple = stats["inferred_simple"]
    chimeric = stats["inferred_chimeric"]
    stats["formatted_summary"] = {
        "inferred_simple": (
            f"Inferred UeccDNA: {simple['total']} sequences, "
            f"{simple['overlapping']} ({simple['overlap_percentage']:.2f}%) "
            f"overlap with confirmed UeccDNA"
        ),
        "inferred_chimeric": (
            f"Inferred CeccDNA: {chimeric['total']} sequences, "
            f"{chimeric['overlapping']} ({chimeric['overlap_percentage']:.2f}%) "
            f"overlap with confirmed CeccDNA"
        ),
    }

    # Save to JSON file if specified
    if output_file:
        with open(output_file, "w") as f:
            json.dump(stats, f, indent=2)

    return stats


def generate_overlap_report(
    confirmed_df: pd.DataFrame,
    inferred_simple: pd.DataFrame | None,
    inferred_chimeric: pd.DataFrame | None,
    simple_redundant: set[str],
    chimeric_redundant: set[str],
    output_file: Path | str | None = None,
) -> str:
    """Generate detailed overlap statistics report.

    Args:
        confirmed_df: Confirmed eccDNA DataFrame
        inferred_simple: Inferred simple eccDNA DataFrame
        inferred_chimeric: Inferred chimeric eccDNA DataFrame
        simple_redundant: Set of redundant simple eccDNA IDs
        chimeric_redundant: Set of redundant chimeric eccDNA IDs
        output_file: Optional file path to save the report

    Returns:
        Report text
    """
    report_lines = []
    report_lines.append("=" * 80)
    report_lines.append("eccDNA Overlap Analysis Report")
    report_lines.append("=" * 80)
    report_lines.append("")

    # Confirmed eccDNA summary
    report_lines.append("## Confirmed eccDNA Summary")
    report_lines.append("-" * 40)
    for ecc_type in ["UeccDNA", "MeccDNA", "CeccDNA"]:
        count = len(confirmed_df[confirmed_df["eccDNA_type"] == ecc_type])
        report_lines.append(f"  {ecc_type}: {count} sequences")
    report_lines.append(f"  Total: {len(confirmed_df)} sequences")
    report_lines.append("")

    # Inferred Simple (UeccDNA) overlap analysis
    report_lines.append("## Inferred Simple eccDNA (→ UeccDNA) Analysis")
    report_lines.append("-" * 40)

    if inferred_simple is not None and not inferred_simple.empty:
        total_simple = len(inferred_simple)
        redundant_simple = len(simple_redundant)
        non_redundant_simple = total_simple - redundant_simple
        overlap_pct = (redundant_simple / total_simple * 100) if total_simple > 0 else 0

        report_lines.append(f"  Total inferred: {total_simple} sequences")
        report_lines.append(
            f"  Overlapping with confirmed: {redundant_simple} ({overlap_pct:.2f}%)"
        )
        report_lines.append(
            f"  Non-redundant (to be added): {non_redundant_simple} ({100-overlap_pct:.2f}%)"
        )

        if redundant_simple > 0:
            report_lines.append("")
            report_lines.append("  Redundant IDs:")
            for rid in sorted(simple_redundant):
                # Find the matching confirmed eccDNA
                row = inferred_simple[inferred_simple["eccDNA_id"] == rid].iloc[0]
                chr_val = row.get("chr", "")
                start_val = row.get("start0", "")
                end_val = row.get("end0", "")
                report_lines.append(f"    - {rid}: {chr_val}:{start_val}-{end_val}")
    else:
        report_lines.append("  No inferred simple eccDNA provided")

    report_lines.append("")

    # Inferred Chimeric (CeccDNA) overlap analysis
    report_lines.append("## Inferred Chimeric eccDNA (→ CeccDNA) Analysis")
    report_lines.append("-" * 40)

    if inferred_chimeric is not None and not inferred_chimeric.empty:
        unique_chimeric = inferred_chimeric["eccDNA_id"].nunique()
        redundant_chimeric = len(chimeric_redundant)
        non_redundant_chimeric = unique_chimeric - redundant_chimeric
        overlap_pct = (redundant_chimeric / unique_chimeric * 100) if unique_chimeric > 0 else 0

        report_lines.append(f"  Total inferred: {unique_chimeric} sequences")
        report_lines.append(
            f"  Overlapping with confirmed: {redundant_chimeric} ({overlap_pct:.2f}%)"
        )
        report_lines.append(
            f"  Non-redundant (to be added): {non_redundant_chimeric} ({100-overlap_pct:.2f}%)"
        )

        if redundant_chimeric > 0:
            report_lines.append("")
            report_lines.append("  Redundant IDs:")
            for rid in sorted(chimeric_redundant):
                report_lines.append(f"    - {rid}")
    else:
        report_lines.append("  No inferred chimeric eccDNA provided")

    report_lines.append("")

    # Summary
    report_lines.append("## Summary")
    report_lines.append("-" * 40)

    simple_added = 0
    chimeric_added = 0

    if inferred_simple is not None and not inferred_simple.empty:
        simple_added = len(inferred_simple) - len(simple_redundant)

    if inferred_chimeric is not None and not inferred_chimeric.empty:
        chimeric_added = inferred_chimeric["eccDNA_id"].nunique() - len(chimeric_redundant)

    report_lines.append(f"  Total confirmed eccDNA: {len(confirmed_df)}")
    report_lines.append(f"  Non-redundant inferred to be added: {simple_added + chimeric_added}")
    report_lines.append(f"    - Simple (UeccDNA): {simple_added}")
    report_lines.append(f"    - Chimeric (CeccDNA): {chimeric_added}")
    report_lines.append(
        f"  Final total eccDNA: {len(confirmed_df) + simple_added + chimeric_added}"
    )

    report_lines.append("")
    report_lines.append("=" * 80)

    report_text = "\n".join(report_lines)

    # Save to file if specified
    if output_file:
        with open(output_file, "w") as f:
            f.write(report_text)

    return report_text


def generate_overlap_statistics(
    confirmed_df: pd.DataFrame,
    inferred_simple: pd.DataFrame | None,
    inferred_chimeric: pd.DataFrame | None,
    simple_redundant: set[str],
    chimeric_redundant: set[str],
    output_file: Path | str | None = None,
) -> dict:
    """Generate structured overlap statistics for reporting.

    Args:
        confirmed_df: Confirmed eccDNA DataFrame
        inferred_simple: Inferred simple eccDNA DataFrame
        inferred_chimeric: Inferred chimeric eccDNA DataFrame
        simple_redundant: Set of redundant simple eccDNA IDs
        chimeric_redundant: Set of redundant chimeric eccDNA IDs
        output_file: Optional JSON file path to save statistics

    Returns:
        Dictionary containing overlap statistics
    """
    stats: dict[str, Any] = {
        "confirmed": {
            "UeccDNA": len(confirmed_df[confirmed_df["eccDNA_type"] == "UeccDNA"]),
            "MeccDNA": len(confirmed_df[confirmed_df["eccDNA_type"] == "MeccDNA"]),
            "CeccDNA": len(confirmed_df[confirmed_df["eccDNA_type"] == "CeccDNA"]),
            "total": len(confirmed_df),
        },
        "inferred_simple": {
            "total": 0,
            "overlapping": 0,
            "overlapping_pct": 0.0,
            "non_redundant": 0,
            "non_redundant_pct": 0.0,
            "redundant_ids": [],
        },
        "inferred_chimeric": {
            "total": 0,
            "overlapping": 0,
            "overlapping_pct": 0.0,
            "non_redundant": 0,
            "non_redundant_pct": 0.0,
            "redundant_ids": [],
        },
        "summary": {
            "total_confirmed": len(confirmed_df),
            "total_inferred": 0,
            "total_overlapping": 0,
            "total_non_redundant": 0,
            "final_total": len(confirmed_df),
        },
    }

    # Process inferred simple
    if inferred_simple is not None and not inferred_simple.empty:
        total_simple = len(inferred_simple)
        redundant_simple = len(simple_redundant)
        non_redundant_simple = total_simple - redundant_simple
        overlap_pct = (redundant_simple / total_simple * 100) if total_simple > 0 else 0

        stats["inferred_simple"] = {
            "total": total_simple,
            "overlapping": redundant_simple,
            "overlapping_pct": round(overlap_pct, 2),
            "non_redundant": non_redundant_simple,
            "non_redundant_pct": round(100 - overlap_pct, 2),
            "redundant_ids": sorted(list(simple_redundant)),
        }

    # Process inferred chimeric
    if inferred_chimeric is not None and not inferred_chimeric.empty:
        unique_chimeric = inferred_chimeric["eccDNA_id"].nunique()
        redundant_chimeric = len(chimeric_redundant)
        non_redundant_chimeric = unique_chimeric - redundant_chimeric
        overlap_pct = (redundant_chimeric / unique_chimeric * 100) if unique_chimeric > 0 else 0

        stats["inferred_chimeric"] = {
            "total": unique_chimeric,
            "overlapping": redundant_chimeric,
            "overlapping_pct": round(overlap_pct, 2),
            "non_redundant": non_redundant_chimeric,
            "non_redundant_pct": round(100 - overlap_pct, 2),
            "redundant_ids": sorted(list(chimeric_redundant)),
        }

    # Update summary
    stats["summary"]["total_inferred"] = (
        stats["inferred_simple"]["total"] + stats["inferred_chimeric"]["total"]
    )
    stats["summary"]["total_overlapping"] = (
        stats["inferred_simple"]["overlapping"] + stats["inferred_chimeric"]["overlapping"]
    )
    stats["summary"]["total_non_redundant"] = (
        stats["inferred_simple"]["non_redundant"] + stats["inferred_chimeric"]["non_redundant"]
    )
    stats["summary"]["final_total"] = (
        stats["summary"]["total_confirmed"] + stats["summary"]["total_non_redundant"]
    )

    # Save to JSON file if specified
    if output_file:
        with open(output_file, "w") as f:
            json.dump(stats, f, indent=2)

    return stats


def merge_eccdna_tables(
    confirmed_file: Path | str | pd.DataFrame,
    inferred_simple: Path | str | pd.DataFrame | None = None,
    inferred_chimeric: Path | str | pd.DataFrame | None = None,
    overlap_threshold: float = 0.99,
    tolerance: int = 10,
    renumber: bool = True,
    overlap_report_file: Path | str | None = None,
    overlap_stats_json: Path | str | None = None,
) -> tuple[pd.DataFrame, str, dict]:
    """Merge confirmed and inferred eccDNA tables with redundancy removal.

    Args:
        confirmed_file: Path to confirmed eccDNA CSV file (contains all types)
        inferred_simple: Path/DataFrame of inferred simple eccDNA or None
        inferred_chimeric: Path/DataFrame of inferred chimeric eccDNA or None
        overlap_threshold: Reciprocal overlap threshold for redundancy
        tolerance: Coordinate tolerance in bp
        renumber: Whether to renumber final IDs
        overlap_report_file: Optional file path to save overlap report
        overlap_stats_json: Optional file path to save overlap statistics as JSON
                           (if None, will default to output_prefix_overlap_stats.json)

    Returns:
        Tuple of (Merged DataFrame, Overlap report text, Overlap statistics dict)
    """
    # Load confirmed table
    if isinstance(confirmed_file, (str, Path)):
        confirmed_df = pd.read_csv(confirmed_file)
        # Generate default overlap_stats_json path if not specified
        if overlap_stats_json is None:
            confirmed_path = Path(confirmed_file)
            overlap_stats_json = confirmed_path.parent / f"{confirmed_path.stem}_overlap_stats.json"
    else:
        confirmed_df = confirmed_file
        # If the caller passes a DataFrame and doesn't provide an output path,
        # avoid writing overlap stats into the current working directory.

    # Load inferred tables if needed
    simple_df = None
    chimeric_df = None
    simple_redundant = set()
    chimeric_redundant = set()

    if inferred_simple is not None:
        if isinstance(inferred_simple, (str, Path)):
            simple_df = pd.read_csv(inferred_simple)
        else:
            simple_df = inferred_simple

    if inferred_chimeric is not None:
        if isinstance(inferred_chimeric, (str, Path)):
            chimeric_df = pd.read_csv(inferred_chimeric)
        else:
            chimeric_df = inferred_chimeric

    # Find redundant entries before merging
    if simple_df is not None and not simple_df.empty:
        simple_redundant = find_redundant_simple(
            simple_df, confirmed_df, overlap_threshold, tolerance
        )

    if chimeric_df is not None and not chimeric_df.empty:
        chimeric_redundant = find_redundant_chimeric(
            chimeric_df, confirmed_df, overlap_threshold, tolerance
        )

    # Generate overlap statistics (JSON)
    overlap_stats = generate_overlap_statistics(
        confirmed_df,
        simple_df,
        chimeric_df,
        simple_redundant,
        chimeric_redundant,
        overlap_stats_json,
    )

    # Generate text report
    overlap_report = generate_overlap_report(
        confirmed_df,
        simple_df,
        chimeric_df,
        simple_redundant,
        chimeric_redundant,
        overlap_report_file,
    )

    # Start with confirmed entries
    logger = get_logger("ecc_unify")
    tables = [confirmed_df]

    # Add non-redundant inferred simple
    if simple_df is not None and not simple_df.empty:
        if len(simple_redundant) < len(simple_df):
            prepared_simple = prepare_inferred_simple(simple_df, simple_redundant)
            if not prepared_simple.empty:
                tables.append(prepared_simple)
                logger.info(f"Added {len(prepared_simple)} non-redundant simple inferred eccDNA")
                logger.info(f"Removed {len(simple_redundant)} redundant simple inferred eccDNA")

    # Add non-redundant inferred chimeric
    if chimeric_df is not None and not chimeric_df.empty:
        unique_ids = chimeric_df["eccDNA_id"].nunique()
        if len(chimeric_redundant) < unique_ids:
            prepared_chimeric = prepare_inferred_chimeric(chimeric_df, chimeric_redundant)
            if not prepared_chimeric.empty:
                tables.append(prepared_chimeric)
                logger.info(
                    f"Added {len(prepared_chimeric)} non-redundant chimeric inferred eccDNA"
                )
                logger.info(f"Removed {len(chimeric_redundant)} redundant chimeric inferred eccDNA")

    # Merge all tables
    if not tables:
        raise ValueError("No data to merge")

    merged = pd.concat(tables, ignore_index=True)

    # Renumber if requested
    if renumber:
        merged = renumber_eccdna(merged)

    # Print summary
    logger = get_logger("ecc_unify")
    logger.info("\nFinal eccDNA counts:")
    for ecc_type in ["UeccDNA", "MeccDNA", "CeccDNA"]:
        type_df = merged[merged["eccDNA_type"] == ecc_type]
        confirmed = len(type_df[type_df["State"] == "Confirmed"])
        inferred = len(type_df[type_df["State"] == "Inferred"])
        logger.info(
            f"  {ecc_type}: {len(type_df)} total ({confirmed} confirmed, {inferred} inferred)"
        )

    # Notify about JSON stats file
    if overlap_stats_json:
        logger.info(f"\nOverlap statistics saved to: {overlap_stats_json}")

    return merged, overlap_report, overlap_stats


def main() -> None:
    """Command-line interface."""
    import argparse

    parser = argparse.ArgumentParser(description="Merge eccDNA tables with redundancy detection")
    parser.add_argument(
        "-c", "--confirmed", required=True, help="Confirmed eccDNA table (CSV with all types)"
    )
    parser.add_argument("-s", "--inferred-simple", help="Inferred simple eccDNA table (CSV)")
    parser.add_argument("-i", "--inferred-chimeric", help="Inferred chimeric eccDNA table (CSV)")
    parser.add_argument("-o", "--output", required=True, help="Output merged table (CSV)")
    parser.add_argument("--overlap-report", help="Output overlap analysis report file (text)")
    parser.add_argument("--overlap-stats", help="Output overlap statistics file (JSON)")
    parser.add_argument(
        "--overlap-threshold",
        type=float,
        default=0.99,
        help="Reciprocal overlap threshold (default: 0.99)",
    )
    parser.add_argument(
        "--tolerance", type=int, default=10, help="Coordinate tolerance in bp (default: 10)"
    )
    parser.add_argument("--no-renumber", action="store_true", help="Don't renumber eccDNA IDs")

    args = parser.parse_args()

    # Merge tables
    result, overlap_report, overlap_stats = merge_eccdna_tables(
        args.confirmed,
        args.inferred_simple,
        args.inferred_chimeric,
        args.overlap_threshold,
        args.tolerance,
        renumber=not args.no_renumber,
        overlap_report_file=args.overlap_report,
        overlap_stats_json=args.overlap_stats,
    )

    # Save result
    logger = get_logger("ecc_unify")
    result.to_csv(args.output, index=False)
    logger.info(f"\nSaved {len(result)} eccDNA entries to {args.output}")

    # Print overlap report if not saved to file
    if not args.overlap_report:
        logger.info("\n" + overlap_report)


if __name__ == "__main__":
    main()
