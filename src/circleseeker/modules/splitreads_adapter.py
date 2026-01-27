"""SplitReads-Core output adapter for CircleSeeker pipeline.

This module converts SplitReads-Core eccDNA_final.txt output into a format
compatible with the iecc_curator module (Step 12). The output format follows
a standard structure that iecc_curator expects for inferred eccDNA processing.
"""

from __future__ import annotations

from pathlib import Path
import re
from typing import Any

import pandas as pd

from circleseeker.utils.logging import get_logger


def parse_splitreads_region(region_str: str) -> tuple[str, int, int, str]:
    """Parse SplitReads-Core merge_region format: chr:start-end_strand

    Args:
        region_str: Region string like "chr2:47242017-47243061_-"

    Returns:
        Tuple of (chromosome, start, end, strand)

    Example:
        >>> parse_splitreads_region("chr2:47242017-47243061_-")
        ('chr2', 47242017, 47243061, '-')
    """
    if not isinstance(region_str, str) or not region_str.strip():
        raise ValueError(f"Invalid SplitReads region: {region_str}")

    region_str = region_str.strip()

    # Split strand
    if region_str.endswith("_+") or region_str.endswith("_-"):
        region, strand = region_str.rsplit("_", 1)
    else:
        region = region_str
        strand = "+"  # Default to positive strand

    # Parse chr:start-end
    if ":" not in region:
        raise ValueError(f"Invalid SplitReads region (missing ':'): {region_str}")
    chrom, coords = region.split(":", 1)
    if not chrom:
        raise ValueError(f"Invalid SplitReads region (empty chromosome): {region_str}")
    if "-" not in coords:
        raise ValueError(f"Invalid SplitReads region (missing '-'): {region_str}")
    start_str, end_str = coords.split("-", 1)
    try:
        start = int(start_str)
        end = int(end_str)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid SplitReads region coordinates: {region_str}") from exc

    return chrom, start, end, strand


def parse_splitreads_regions(region_str: str) -> list[tuple[str, int, int, str]]:
    """Parse SplitReads-Core merge_region field that may contain multiple segments."""
    if pd.isna(region_str):
        raise ValueError(f"Unable to parse SplitReads region: {region_str}")
    region_str = str(region_str).strip()
    if not region_str:
        raise ValueError(f"Unable to parse SplitReads region: {region_str}")

    parts = [part.strip() for part in re.split(r"[;,|]+|\s+", region_str) if part.strip()]
    if not parts:
        raise ValueError(f"Unable to parse SplitReads region: {region_str}")
    return [parse_splitreads_region(part) for part in parts]


def convert_splitreads_to_overview_format(
    splitreads_file: Path | str,
    *,
    prob_present_default: float = 0.95,
    prob_artifact_default: float = 0.0,
) -> pd.DataFrame:
    """Convert SplitReads-Core eccDNA_final.txt to overview.tsv format.

    Args:
        splitreads_file: Path to SplitReads-Core eccDNA_final.txt output file

    Returns:
        DataFrame in overview.tsv format with columns:
        - circle_id
        - regions
        - circle_length
        - segment_count (num_segments)
        - num_split_reads
        - prob_present (placeholder default; SplitReads-Core does not emit probabilities)
        - prob_artifact (placeholder default; SplitReads-Core does not emit artifact scores)
        - af_nanopore (calculated from coverage)
    """
    logger = get_logger(__name__)

    if not 0.0 <= float(prob_present_default) <= 1.0:
        raise ValueError(f"prob_present_default must be in [0,1], got {prob_present_default}")
    if not 0.0 <= float(prob_artifact_default) <= 1.0:
        raise ValueError(f"prob_artifact_default must be in [0,1], got {prob_artifact_default}")

    splitreads_file = Path(splitreads_file)
    if not splitreads_file.exists():
        raise FileNotFoundError(f"SplitReads-Core output file not found: {splitreads_file}")

    # Read SplitReads-Core output
    df = pd.read_csv(splitreads_file, sep="\t")

    logger.info(f"Read {len(df)} eccDNA entries from SplitReads-Core output")

    # Convert to overview format
    converted_rows = []

    for idx, row in df.iterrows():
        region_str = row["merge_region"]

        # Parse region(s); SplitReads-Core uses semicolon/comma to separate segments
        try:
            regions = parse_splitreads_regions(region_str)
        except ValueError as e:
            logger.warning(
                "Skipping row %s: failed to parse region %r - %s",
                idx, region_str, e
            )
            continue

        reported_segments = row.get("num_region")
        reported_count = None
        if pd.notna(reported_segments):
            try:
                reported_count = int(reported_segments)
            except (TypeError, ValueError):
                reported_count = None

        if reported_count is not None and len(regions) != reported_count:
            logger.warning(
                "eccDNA %s reports num_region=%s but parsed %s segments. "
                "Original merge_region: %r",
                row.get("id"),
                row.get("num_region"),
                len(regions),
                region_str,
            )

        # Overview format expects semicolon-separated regions without strand info
        regions_formatted = ";".join(f"{chrom}:{start}-{end}" for chrom, start, end, _ in regions)

        # Calculate approximate abundance from coverage
        # SplitReads-Core coverage = totalbase / merge_len
        # We'll use this as a proxy for hifi_abundance
        coverage = row.get("coverage", 0)
        hifi_abundance = min(coverage / 100.0, 1.0)  # Normalize to 0-1 range

        converted_row = {
            "circle_id": row.get("id", f"splitreads_{idx}"),
            "regions": regions_formatted,
            "circle_length": row.get("merge_len", 0),
            "segment_count": len(regions),
            "num_split_reads": row.get("numreads", 0),
            # SplitReads-Core doesn't provide probability scores; keep explicit placeholders.
            "prob_present": float(prob_present_default),
            "prob_artifact": float(prob_artifact_default),
            "af_nanopore": hifi_abundance,
            # Additional SplitReads-Core metadata
            "splitreads_ctc": row.get("ctc", True),
            "splitreads_totalbase": row.get("totalbase", 0),
            "splitreads_coverage": coverage,
        }

        converted_rows.append(converted_row)

    result_df = pd.DataFrame(converted_rows)
    logger.info(f"Converted {len(result_df)} eccDNA entries to overview format")

    return result_df


def write_overview_tsv(
    splitreads_file: Path | str,
    output_tsv: Path | str,
    *,
    prob_present_default: float = 0.95,
    prob_artifact_default: float = 0.0,
) -> Path:
    """Convert SplitReads-Core output and write as overview TSV.

    Args:
        splitreads_file: Path to SplitReads-Core eccDNA_final.txt
        output_tsv: Path to output TSV file

    Returns:
        Path to the created TSV file
    """
    logger = get_logger(__name__)

    df = convert_splitreads_to_overview_format(
        splitreads_file,
        prob_present_default=prob_present_default,
        prob_artifact_default=prob_artifact_default,
    )

    output_tsv = Path(output_tsv)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    # Write as TSV (tab-separated)
    df.to_csv(output_tsv, sep="\t", index=False)

    logger.info(f"Wrote overview TSV to {output_tsv}")

    return output_tsv


def _parse_args() -> Any:
    """Parse CLI arguments for direct script execution."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Convert SplitReads-Core output to overview TSV format"
    )
    parser.add_argument("-i", "--input", required=True, help="Input SplitReads-Core eccDNA_final.txt file")
    parser.add_argument(
        "-o", "--output", required=True, help="Output TSV file (overview format)"
    )
    return parser.parse_args()


def main() -> None:
    """CLI entry point."""
    args = _parse_args()
    write_overview_tsv(args.input, args.output)


if __name__ == "__main__":
    main()
