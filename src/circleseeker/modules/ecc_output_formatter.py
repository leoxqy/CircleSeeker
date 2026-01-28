"""
ECC Output Formatter - Generate standardized eccDNA output files.

This module transforms internal eccDNA data into the standardized output format:
- eccDNA_summary.csv: One row per eccDNA with key metrics
- eccDNA_regions.csv: All genomic regions (sources, candidates, segments)
- eccDNA_reads.csv: Read-level support information
- Type-specific BED and FASTA files

Output directory structure:
    output/
    ├── eccDNA_summary.csv
    ├── eccDNA_regions.csv
    ├── eccDNA_reads.csv
    ├── eccDNA_all.fasta
    ├── Uecc/
    │   ├── uecc.bed
    │   └── uecc.fasta
    ├── Mecc/
    │   ├── mecc_sites.bed
    │   └── mecc.fasta
    ├── Cecc/
    │   ├── cecc_segments.bed
    │   ├── cecc_junctions.bedpe
    │   └── cecc.fasta
    ├── report.html
    └── summary.txt
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Optional

import pandas as pd
from circleseeker.utils.logging import get_logger


# ================== ID Formatting ==================

def format_eccdna_id(ecc_type: str, number: int, width: int = 4) -> str:
    """Format eccDNA ID with zero-padded number.

    Args:
        ecc_type: Type prefix (Uecc, Mecc, Cecc)
        number: Sequence number
        width: Zero-padding width (default 4)

    Returns:
        Formatted ID like "UeccDNA0001"
    """
    type_map = {
        "Uecc": "UeccDNA",
        "Mecc": "MeccDNA",
        "Cecc": "CeccDNA",
        "UeccDNA": "UeccDNA",
        "MeccDNA": "MeccDNA",
        "CeccDNA": "CeccDNA",
    }
    prefix = type_map.get(ecc_type, ecc_type)
    return f"{prefix}{number:0{width}d}"


def renumber_with_padding(df: pd.DataFrame, width: int = 4) -> pd.DataFrame:
    """Renumber eccDNA IDs with zero-padding, Confirmed before Inferred.

    Args:
        df: DataFrame with eccDNA_type, State columns
        width: Zero-padding width

    Returns:
        DataFrame with updated eccDNA_id and _old_eccDNA_id columns
    """
    df = df.copy()

    # Preserve original eccDNA_id for mapping
    # Note: unified.csv may already have 'original_id' column (e.g., ICeccDNA1)
    # We need to map both old eccDNA_id AND original_id to new IDs
    df["_old_eccDNA_id"] = df["eccDNA_id"]

    # Sort by type and state
    type_order = {"UeccDNA": 0, "MeccDNA": 1, "CeccDNA": 2}
    state_order = {"Confirmed": 0, "Inferred": 1}

    df["_type_order"] = df["eccDNA_type"].map(type_order).fillna(3)
    df["_state_order"] = df["State"].map(state_order).fillna(2)
    df = df.sort_values(["_type_order", "_state_order"]).reset_index(drop=True)

    # Assign new IDs – iterate in canonical order, then append any
    # unexpected types so len(new_ids) always matches len(df).
    known_types = ["UeccDNA", "MeccDNA", "CeccDNA"]
    remaining_types = [t for t in df["eccDNA_type"].unique() if t not in known_types]

    new_ids: list[str] = []
    for ecc_type in known_types + remaining_types:
        count = (df["eccDNA_type"] == ecc_type).sum()
        new_ids.extend(format_eccdna_id(ecc_type, i, width) for i in range(1, count + 1))

    df["eccDNA_id"] = new_ids
    df = df.drop(columns=["_type_order", "_state_order"])

    return df


# ================== Location Formatting ==================

def format_location(chr_val: str, start: int, end: int, strand: str) -> str:
    """Format a single location string.

    Returns:
        Location string like "Chr1:100-200(+)"
    """
    return f"{chr_val}:{start}-{end}({strand})"


def format_mecc_location(regions_df: pd.DataFrame) -> str:
    """Format MeccDNA location with | separator for candidate sites.

    Args:
        regions_df: DataFrame with chr, start, end, strand for one MeccDNA

    Returns:
        Location string like "Chr1:100-200(+)|Chr2:300-400(+)"
    """
    locs = (
        regions_df["chr"].astype(str) + ":" +
        regions_df["start"].astype(str) + "-" +
        regions_df["end"].astype(str) + "(" +
        regions_df["strand"].astype(str) + ")"
    )
    return "|".join(locs)


def format_cecc_location(regions_df: pd.DataFrame) -> str:
    """Format CeccDNA location with ; separator for segments.

    Args:
        regions_df: DataFrame with chr, start, end, strand for one CeccDNA

    Returns:
        Location string like "Chr1:100-200(+);Chr2:300-400(+)"
    """
    # Sort by segment index
    if "region_idx" in regions_df.columns:
        regions_df = regions_df.sort_values("region_idx")

    locs = (
        regions_df["chr"].astype(str) + ":" +
        regions_df["start"].astype(str) + "-" +
        regions_df["end"].astype(str) + "(" +
        regions_df["strand"].astype(str) + ")"
    )
    return ";".join(locs)


# ================== Summary Table Generation ==================

def generate_summary_table(
    unified_df: pd.DataFrame,
    regions_df: pd.DataFrame,
) -> pd.DataFrame:
    """Generate eccDNA_summary.csv from unified data.

    Args:
        unified_df: Unified eccDNA DataFrame
        regions_df: Regions DataFrame with all genomic locations

    Returns:
        Summary DataFrame with one row per eccDNA
    """
    logger = get_logger("ecc_output_formatter")

    summary_rows = []

    for _, row in unified_df.iterrows():
        ecc_id = row["eccDNA_id"]
        ecc_type = row["eccDNA_type"]
        state = row.get("State", "Confirmed")
        length = row.get("Length", row.get("length", 0))

        # Get regions for this eccDNA
        ecc_regions = regions_df[regions_df["eccDNA_id"] == ecc_id]

        # Determine chr, start, end, strand and location based on type
        if ecc_type == "UeccDNA":
            # Single source
            if len(ecc_regions) > 0:
                reg = ecc_regions.iloc[0]
                chr_val = reg["chr"]
                start = reg["start"]
                end = reg["end"]
                strand = reg["strand"]
                location = format_location(chr_val, start, end, strand)
            else:
                chr_val, start, end, strand = ".", ".", ".", "."
                location = ""
            segment_count = 1

        elif ecc_type == "MeccDNA":
            # Multiple candidate sites - use primary for single columns
            primary = ecc_regions[ecc_regions["role"] == "primary"]
            if len(primary) > 0:
                reg = primary.iloc[0]
                chr_val = reg["chr"]
                start = reg["start"]
                end = reg["end"]
                strand = reg["strand"]
            elif len(ecc_regions) > 0:
                reg = ecc_regions.iloc[0]
                chr_val = reg["chr"]
                start = reg["start"]
                end = reg["end"]
                strand = reg["strand"]
            else:
                chr_val, start, end, strand = ".", ".", ".", "."

            # Location shows all candidate sites with |
            location = format_mecc_location(ecc_regions)
            segment_count = 1  # MeccDNA is still one segment, multiple candidates

        elif ecc_type == "CeccDNA":
            # Multiple segments
            chr_val = "multi"
            start = "."
            end = "."
            strand = "."
            location = format_cecc_location(ecc_regions)
            segment_count = len(ecc_regions)

        else:
            chr_val, start, end, strand = ".", ".", ".", "."
            location = ""
            segment_count = 1

        # Extract metrics
        read_count = row.get("reads_count", row.get("read_count", 1))
        copy_number = row.get("copy_number", row.get("repeat_number", 1))
        confidence = row.get("confidence_score", row.get("confidence", ""))

        summary_rows.append({
            "eccDNA_id": ecc_id,
            "type": ecc_type.replace("DNA", ""),  # Uecc, Mecc, Cecc
            "state": state,
            "chr": chr_val,
            "start": start,
            "end": end,
            "strand": strand,
            "length": length,
            "location": location,
            "segment_count": segment_count,
            "read_count": read_count,
            "copy_number": copy_number,
            "confidence": confidence,
        })

    summary_df = pd.DataFrame(summary_rows)
    logger.info(f"Generated summary table with {len(summary_df)} entries")

    return summary_df


# ================== Regions Table Generation ==================

def generate_regions_table(
    uecc_df: Optional[pd.DataFrame] = None,
    mecc_df: Optional[pd.DataFrame] = None,
    cecc_df: Optional[pd.DataFrame] = None,
    inferred_simple_df: Optional[pd.DataFrame] = None,
    inferred_chimeric_df: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """Generate eccDNA_regions.csv from type-specific data.

    Args:
        uecc_df: UeccDNA core CSV data (Confirmed)
        mecc_df: MeccDNA sites CSV data (Confirmed)
        cecc_df: CeccDNA segments CSV data (Confirmed)
        inferred_simple_df: Inferred simple eccDNA CSV
        inferred_chimeric_df: Inferred chimeric eccDNA CSV

    Returns:
        Unified regions DataFrame
    """
    logger = get_logger("ecc_output_formatter")

    regions_rows = []

    def _resolve_col(df: pd.DataFrame, primary: str, fallback: str, default: object = "") -> pd.Series:
        """Resolve column name with fallback."""
        if primary in df.columns:
            return df[primary]
        if fallback in df.columns:
            return df[fallback]
        return pd.Series(default, index=df.index)

    # Process UeccDNA (Confirmed) — vectorized
    if uecc_df is not None and not uecc_df.empty:
        uecc_block = pd.DataFrame({
            "eccDNA_id": uecc_df["eccDNA_id"].values,
            "region_idx": 1,
            "chr": _resolve_col(uecc_df, "chr", "eChr", "").values,
            "start": _resolve_col(uecc_df, "start0", "eStart0", 0).values,
            "end": _resolve_col(uecc_df, "end0", "eEnd0", 0).values,
            "strand": _resolve_col(uecc_df, "strand", "eStrand", "+").values,
            "length": _resolve_col(uecc_df, "length", "eLength", 0).values,
            "role": "source",
        })
        regions_rows.extend(uecc_block.to_dict("records"))

    # Process MeccDNA (Confirmed)
    if mecc_df is not None and not mecc_df.empty:
        for ecc_id, group in mecc_df.groupby("eccDNA_id"):
            # Use itertuples for performance
            for i, row in enumerate(group.itertuples(index=False), start=1):
                # Determine role: primary for first/best site, candidate for others
                is_primary = (i == 1) or getattr(row, "is_primary", False)
                role = "primary" if is_primary else "candidate"

                regions_rows.append({
                    "eccDNA_id": ecc_id,
                    "region_idx": i,
                    "chr": getattr(row, "chr", ""),
                    "start": getattr(row, "start0", 0),
                    "end": getattr(row, "end0", 0),
                    "strand": getattr(row, "strand", "+"),
                    "length": getattr(row, "length", 0),
                    "role": role,
                })

    # Process CeccDNA (Confirmed)
    if cecc_df is not None and not cecc_df.empty:
        for ecc_id, group in cecc_df.groupby("eccDNA_id"):
            group = group.sort_values("seg_index") if "seg_index" in group.columns else group
            n_segs = len(group)

            # Use itertuples for performance
            for i, row in enumerate(group.itertuples(index=False), start=1):
                # Determine role based on position
                if i == 1:
                    role = "head"
                elif i == n_segs:
                    role = "tail"
                else:
                    role = "middle"

                regions_rows.append({
                    "eccDNA_id": ecc_id,
                    "region_idx": i,
                    "chr": getattr(row, "chr", ""),
                    "start": getattr(row, "start0", 0),
                    "end": getattr(row, "end0", 0),
                    "strand": getattr(row, "strand", "+"),
                    "length": getattr(row, "length", getattr(row, "seg_length", 0)),
                    "role": role,
                })

    # Process Inferred Simple (UeccDNA) — vectorized
    if inferred_simple_df is not None and not inferred_simple_df.empty:
        inf_block = pd.DataFrame({
            "eccDNA_id": inferred_simple_df["eccDNA_id"].values,
            "region_idx": 1,
            "chr": _resolve_col(inferred_simple_df, "chr", "Chr", "").values,
            "start": _resolve_col(inferred_simple_df, "start0", "Start0", 0).values,
            "end": _resolve_col(inferred_simple_df, "end0", "End0", 0).values,
            "strand": _resolve_col(inferred_simple_df, "strand", "Strand", "+").values,
            "length": _resolve_col(inferred_simple_df, "length", "Length", 0).values,
            "role": "source",
        })
        regions_rows.extend(inf_block.to_dict("records"))

    # Process Inferred Chimeric (CeccDNA)
    if inferred_chimeric_df is not None and not inferred_chimeric_df.empty:
        for ecc_id, group in inferred_chimeric_df.groupby("eccDNA_id"):
            # Sort by seg_index if available
            if "seg_index" in group.columns:
                group = group.sort_values("seg_index")
            n_segs = len(group)

            for i, (_, row) in enumerate(group.iterrows(), start=1):
                # Use junction_role if available, otherwise determine by position
                role = row.get("junction_role", "")
                if not role or pd.isna(role):
                    if i == 1:
                        role = "head"
                    elif i == n_segs:
                        role = "tail"
                    else:
                        role = "middle"

                regions_rows.append({
                    "eccDNA_id": ecc_id,
                    "region_idx": i,
                    "chr": row.get("chr", ""),
                    "start": row.get("start0", 0),
                    "end": row.get("end0", 0),
                    "strand": row.get("strand", "+"),
                    "length": row.get("length", 0),
                    "role": role,
                })

    regions_df = pd.DataFrame(
        regions_rows,
        columns=["eccDNA_id", "region_idx", "chr", "start", "end", "strand", "length", "role"],
    )
    logger.info(f"Generated regions table with {len(regions_df)} entries")

    return regions_df


# ================== Reads Table Generation ==================

def generate_reads_table(
    uecc_df: Optional[pd.DataFrame] = None,
    mecc_df: Optional[pd.DataFrame] = None,
    cecc_df: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """Generate eccDNA_reads.csv with read-level information.

    Args:
        uecc_df: UeccDNA core CSV data with read_name column
        mecc_df: MeccDNA sites CSV data
        cecc_df: CeccDNA segments CSV data

    Returns:
        Reads DataFrame with one row per eccDNA-read pair
    """
    logger = get_logger("ecc_output_formatter")

    reads_rows = []

    def parse_reads(read_str: str, ecc_id: str, row: pd.Series) -> list[dict[str, Any]]:
        """Parse semicolon-separated read names and create rows."""
        if pd.isna(read_str) or not read_str:
            return []

        rows: list[dict[str, Any]] = []
        for read_name in str(read_str).split(";"):
            read_name = read_name.strip()
            if read_name:
                rows.append({
                    "eccDNA_id": ecc_id,
                    "read_name": read_name,
                    "copy_number": row.get("copy_number", row.get("repeat_number", 1)),
                    "match_degree": row.get("match_degree", ""),
                    "identity": row.get("identity_best", row.get("identity", "")),
                    "mapq": row.get("mapq_best", row.get("mapq", "")),
                })
        return rows

    # Process each type
    for df in [uecc_df, mecc_df, cecc_df]:
        if df is None or df.empty:
            continue

        # Get unique eccDNA entries (for Mecc/Cecc, group first)
        if "read_name" in df.columns:
            for _, row in df.drop_duplicates("eccDNA_id").iterrows():
                ecc_id = row["eccDNA_id"]
                read_str = row.get("read_name", "")
                reads_rows.extend(parse_reads(read_str, ecc_id, row))

    reads_df = pd.DataFrame(reads_rows)
    if not reads_df.empty:
        reads_df = reads_df.drop_duplicates(["eccDNA_id", "read_name"])

    logger.info(f"Generated reads table with {len(reads_df)} entries")

    return reads_df


# ================== BED File Generation ==================

def generate_uecc_bed(regions_df: pd.DataFrame, output_path: Path) -> None:
    """Generate uecc.bed file.

    Format: chr, start, end, name, score (read_count), strand
    """
    uecc_regions = regions_df[regions_df["role"] == "source"].copy()

    if uecc_regions.empty:
        return

    bed_df = pd.DataFrame({
        "chr": uecc_regions["chr"].values,
        "start": uecc_regions["start"].astype(int).values,
        "end": uecc_regions["end"].astype(int).values,
        "name": uecc_regions["eccDNA_id"].values,
        "score": 1,
        "strand": uecc_regions["strand"].values,
    })
    bed_df.to_csv(output_path, sep="\t", index=False, header=False)


def generate_mecc_bed(regions_df: pd.DataFrame, output_path: Path) -> None:
    """Generate mecc_sites.bed file.

    Format: chr, start, end, name (with role), score, strand
    """
    mecc_regions = regions_df[regions_df["role"].isin(["primary", "candidate"])].copy()

    if mecc_regions.empty:
        return

    bed_df = pd.DataFrame({
        "chr": mecc_regions["chr"].values,
        "start": mecc_regions["start"].astype(int).values,
        "end": mecc_regions["end"].astype(int).values,
        "name": mecc_regions["eccDNA_id"].astype(str) + "|" + mecc_regions["role"].astype(str),
        "score": 1,
        "strand": mecc_regions["strand"].values,
    })
    bed_df.to_csv(output_path, sep="\t", index=False, header=False)


def generate_cecc_bed(regions_df: pd.DataFrame, output_path: Path) -> None:
    """Generate cecc_segments.bed file.

    Format: chr, start, end, name (with seg info), score, strand
    """
    cecc_regions = regions_df[regions_df["role"].isin(["head", "middle", "tail"])].copy()

    if cecc_regions.empty:
        return

    cecc_regions = cecc_regions.sort_values(["eccDNA_id", "region_idx"])
    seg_counts = cecc_regions.groupby("eccDNA_id")["region_idx"].transform("count")

    bed_df = pd.DataFrame({
        "chr": cecc_regions["chr"].values,
        "start": cecc_regions["start"].astype(int).values,
        "end": cecc_regions["end"].astype(int).values,
        "name": (
            cecc_regions["eccDNA_id"].astype(str) + "|seg" +
            cecc_regions["region_idx"].astype(str) + "/" +
            seg_counts.astype(str) + "|" +
            cecc_regions["role"].astype(str)
        ),
        "score": 1,
        "strand": cecc_regions["strand"].values,
    })
    bed_df.to_csv(output_path, sep="\t", index=False, header=False)


def generate_cecc_bedpe(regions_df: pd.DataFrame, output_path: Path) -> None:
    """Generate cecc_junctions.bedpe file.

    Format: chr1, start1, end1, chr2, start2, end2, name, score, strand1, strand2
    """
    cecc_regions = regions_df[regions_df["role"].isin(["head", "middle", "tail"])].copy()

    if cecc_regions.empty:
        return

    bedpe_rows = []
    for ecc_id, group in cecc_regions.groupby("eccDNA_id"):
        group = group.sort_values("region_idx")
        segs = list(group.itertuples())

        # Create junctions between consecutive segments
        for i in range(len(segs) - 1):
            seg1 = segs[i]
            seg2 = segs[i + 1]

            bedpe_rows.append({
                "chr1": seg1.chr,
                "start1": int(seg1.end) - 1,
                "end1": int(seg1.end),
                "chr2": seg2.chr,
                "start2": int(seg2.start),
                "end2": int(seg2.start) + 1,
                "name": f"{ecc_id}|{seg1.region_idx}->{seg2.region_idx}",
                "score": 1,
                "strand1": seg1.strand,
                "strand2": seg2.strand,
            })

    bedpe_df = pd.DataFrame(bedpe_rows)
    bedpe_df.to_csv(output_path, sep="\t", index=False, header=False)


# ================== FASTA File Generation ==================

def generate_fasta_files(
    sequences: dict[str, str],
    output_dir: Path,
    summary_df: pd.DataFrame,
) -> None:
    """Generate FASTA files for all eccDNA and by type.

    FASTA header format:
        >{eccDNA_id}|{location}|length={length}|state={state}

    Args:
        sequences: Dict mapping eccDNA_id to sequence
        output_dir: Output directory
        summary_df: Summary DataFrame to determine types
    """
    logger = get_logger("ecc_output_formatter")

    # Create type subdirectories
    (output_dir / "Uecc").mkdir(parents=True, exist_ok=True)
    (output_dir / "Mecc").mkdir(parents=True, exist_ok=True)
    (output_dir / "Cecc").mkdir(parents=True, exist_ok=True)

    # Build header info from summary_df (vectorized)
    info_cols = ["location", "length", "state", "type"]
    for col in info_cols:
        if col not in summary_df.columns:
            summary_df[col] = ""
    header_info = summary_df.set_index("eccDNA_id")[info_cols].to_dict("index")

    def format_header(ecc_id: str) -> str:
        """Format FASTA header with metadata."""
        info = header_info.get(ecc_id, {})
        location = info.get("location", "")
        length = info.get("length", "")
        state = info.get("state", "")
        return f">{ecc_id}|{location}|length={length}|state={state}"

    # Group sequences by type
    type_seqs: dict[str, dict[str, str]] = {"Uecc": {}, "Mecc": {}, "Cecc": {}}

    for ecc_id in sequences:
        info = header_info.get(ecc_id, {})
        ecc_type = info.get("type", "")
        if ecc_type in type_seqs:
            type_seqs[ecc_type][ecc_id] = sequences[ecc_id]

    # Write all sequences
    all_fasta_path = output_dir / "eccDNA_all.fasta"
    with open(all_fasta_path, "w") as f:
        for ecc_id, seq in sequences.items():
            f.write(f"{format_header(ecc_id)}\n{seq}\n")
    logger.info(f"Wrote {len(sequences)} sequences to {all_fasta_path}")

    # Write type-specific FASTA files
    for ecc_type, seqs in type_seqs.items():
        if seqs:
            fasta_path = output_dir / ecc_type / f"{ecc_type.lower()}.fasta"
            with open(fasta_path, "w") as f:
                for ecc_id, seq in seqs.items():
                    f.write(f"{format_header(ecc_id)}\n{seq}\n")
            logger.info(f"Wrote {len(seqs)} {ecc_type} sequences to {fasta_path}")


# ================== Main Output Function ==================

def format_output(
    unified_csv: Path,
    uecc_core_csv: Optional[Path] = None,
    mecc_sites_csv: Optional[Path] = None,
    cecc_segments_csv: Optional[Path] = None,
    uecc_fasta: Optional[Path] = None,
    mecc_fasta: Optional[Path] = None,
    cecc_fasta: Optional[Path] = None,
    output_dir: Path = Path("output"),
    id_width: int = 4,
) -> None:
    """Main function to generate all output files.

    Args:
        unified_csv: Path to unified eccDNA CSV
        uecc_core_csv: Path to UeccDNA core CSV
        mecc_sites_csv: Path to MeccDNA sites CSV
        cecc_segments_csv: Path to CeccDNA segments CSV
        uecc_fasta: Path to UeccDNA FASTA
        mecc_fasta: Path to MeccDNA FASTA
        cecc_fasta: Path to CeccDNA FASTA
        output_dir: Output directory
        id_width: Zero-padding width for IDs
    """
    logger = get_logger("ecc_output_formatter")

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    unified_df = pd.read_csv(unified_csv)

    uecc_df = pd.read_csv(uecc_core_csv) if uecc_core_csv and uecc_core_csv.exists() else None
    mecc_df = pd.read_csv(mecc_sites_csv) if mecc_sites_csv and mecc_sites_csv.exists() else None
    cecc_df = pd.read_csv(cecc_segments_csv) if cecc_segments_csv and cecc_segments_csv.exists() else None

    # Renumber IDs with padding
    unified_df = renumber_with_padding(unified_df, id_width)

    # Build ID mapping for updating other tables
    # renumber_with_padding() stores pre-renumber IDs in '_old_eccDNA_id';
    # also honour 'original_id' if it already exists (e.g. ICeccDNA1).
    id_map = dict(zip(unified_df["_old_eccDNA_id"], unified_df["eccDNA_id"]))
    if "original_id" in unified_df.columns:
        id_map.update(zip(unified_df["original_id"], unified_df["eccDNA_id"]))

    # Update IDs in type-specific tables
    if uecc_df is not None:
        uecc_df["eccDNA_id"] = uecc_df["eccDNA_id"].map(id_map).fillna(uecc_df["eccDNA_id"])
    if mecc_df is not None:
        mecc_df["eccDNA_id"] = mecc_df["eccDNA_id"].map(id_map).fillna(mecc_df["eccDNA_id"])
    if cecc_df is not None:
        cecc_df["eccDNA_id"] = cecc_df["eccDNA_id"].map(id_map).fillna(cecc_df["eccDNA_id"])

    # Generate regions table
    regions_df = generate_regions_table(uecc_df, mecc_df, cecc_df)

    # Generate summary table
    summary_df = generate_summary_table(unified_df, regions_df)

    # Generate reads table
    reads_df = generate_reads_table(uecc_df, mecc_df, cecc_df)

    # Save CSV files
    summary_df.to_csv(output_dir / "eccDNA_summary.csv", index=False)
    regions_df.to_csv(output_dir / "eccDNA_regions.csv", index=False)
    reads_df.to_csv(output_dir / "eccDNA_reads.csv", index=False)

    logger.info(f"Saved eccDNA_summary.csv ({len(summary_df)} rows)")
    logger.info(f"Saved eccDNA_regions.csv ({len(regions_df)} rows)")
    logger.info(f"Saved eccDNA_reads.csv ({len(reads_df)} rows)")

    # Create type subdirectories
    (output_dir / "Uecc").mkdir(exist_ok=True)
    (output_dir / "Mecc").mkdir(exist_ok=True)
    (output_dir / "Cecc").mkdir(exist_ok=True)

    # Generate BED files
    generate_uecc_bed(regions_df, output_dir / "Uecc" / "uecc.bed")
    generate_mecc_bed(regions_df, output_dir / "Mecc" / "mecc_sites.bed")
    generate_cecc_bed(regions_df, output_dir / "Cecc" / "cecc_segments.bed")
    generate_cecc_bedpe(regions_df, output_dir / "Cecc" / "cecc_junctions.bedpe")

    # Load and merge FASTA sequences
    sequences: dict[str, str] = {}

    def load_fasta(fasta_path: Path) -> dict[str, str]:
        """Load FASTA file into dict."""
        seqs: dict[str, str] = {}
        if not fasta_path.exists():
            return seqs

        current_id: Optional[str] = None
        current_seq: list[str] = []

        with open(fasta_path) as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_id:
                        seqs[current_id] = "".join(current_seq)
                    current_id = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)

            if current_id:
                seqs[current_id] = "".join(current_seq)

        return seqs

    # Load sequences and update IDs
    for fasta_path in [uecc_fasta, mecc_fasta, cecc_fasta]:
        if fasta_path and fasta_path.exists():
            for old_id, seq in load_fasta(fasta_path).items():
                new_id = id_map.get(old_id, old_id)
                sequences[new_id] = seq

    # Generate FASTA files
    generate_fasta_files(sequences, output_dir, summary_df)

    logger.info(f"Output files generated in {output_dir}")


# ================== CLI ==================

def main() -> None:
    """Command-line interface."""
    import argparse

    parser = argparse.ArgumentParser(description="Format eccDNA output files")
    parser.add_argument("-u", "--unified-csv", required=True, help="Unified eccDNA CSV")
    parser.add_argument("--uecc-core", help="UeccDNA core CSV")
    parser.add_argument("--mecc-sites", help="MeccDNA sites CSV")
    parser.add_argument("--cecc-segments", help="CeccDNA segments CSV")
    parser.add_argument("--uecc-fasta", help="UeccDNA FASTA")
    parser.add_argument("--mecc-fasta", help="MeccDNA FASTA")
    parser.add_argument("--cecc-fasta", help="CeccDNA FASTA")
    parser.add_argument("-o", "--output-dir", default="output", help="Output directory")
    parser.add_argument("--id-width", type=int, default=4, help="ID zero-padding width")

    args = parser.parse_args()

    format_output(
        unified_csv=Path(args.unified_csv),
        uecc_core_csv=Path(args.uecc_core) if args.uecc_core else None,
        mecc_sites_csv=Path(args.mecc_sites) if args.mecc_sites else None,
        cecc_segments_csv=Path(args.cecc_segments) if args.cecc_segments else None,
        uecc_fasta=Path(args.uecc_fasta) if args.uecc_fasta else None,
        mecc_fasta=Path(args.mecc_fasta) if args.mecc_fasta else None,
        cecc_fasta=Path(args.cecc_fasta) if args.cecc_fasta else None,
        output_dir=Path(args.output_dir),
        id_width=args.id_width,
    )


if __name__ == "__main__":
    main()
