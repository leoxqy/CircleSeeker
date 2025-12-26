"""Cresil output adapter for CircleSeeker pipeline.

This module converts Cresil eccDNA_final.txt output into a format compatible
with the iecc_curator module (Step 12), allowing seamless integration of Cresil
as an alternative to Cyrcular for circular DNA detection.
"""

from __future__ import annotations

from pathlib import Path
import re
import pandas as pd

from circleseeker.utils.logging import get_logger


def parse_cresil_region(region_str: str) -> tuple[str, int, int, str]:
    """Parse Cresil merge_region format: chr:start-end_strand

    Args:
        region_str: Region string like "chr2:47242017-47243061_-"

    Returns:
        Tuple of (chromosome, start, end, strand)

    Example:
        >>> parse_cresil_region("chr2:47242017-47243061_-")
        ('chr2', 47242017, 47243061, '-')
    """
    # Split strand
    if region_str.endswith("_+") or region_str.endswith("_-"):
        region, strand = region_str.rsplit("_", 1)
    else:
        region = region_str
        strand = "+"  # Default to positive strand

    # Parse chr:start-end
    chrom, coords = region.split(":")
    start_str, end_str = coords.split("-")
    start, end = int(start_str), int(end_str)

    return chrom, start, end, strand


def parse_cresil_regions(region_str: str) -> list[tuple[str, int, int, str]]:
    """Parse Cresil merge_region field that may contain multiple segments."""
    parts = [part.strip() for part in re.split(r"[;,]+|\s+", region_str) if part.strip()]
    if not parts:
        raise ValueError(f"Unable to parse Cresil region: {region_str}")
    return [parse_cresil_region(part) for part in parts]


def convert_cresil_to_cyrcular_format(cresil_file: Path | str) -> pd.DataFrame:
    """Convert Cresil eccDNA_final.txt to Cyrcular overview.tsv compatible format.

    Args:
        cresil_file: Path to Cresil eccDNA_final.txt output file

    Returns:
        DataFrame in Cyrcular overview.tsv format with columns:
        - circle_id
        - regions
        - circle_length
        - segment_count (num_segments)
        - num_split_reads
        - prob_present (set to 0.95 as default)
        - prob_artifact (set to 0.0 as default)
        - af_nanopore (calculated from coverage)
    """
    logger = get_logger(__name__)

    cresil_file = Path(cresil_file)
    if not cresil_file.exists():
        raise FileNotFoundError(f"Cresil output file not found: {cresil_file}")

    # Read Cresil output
    df = pd.read_csv(cresil_file, sep="\t")

    logger.info(f"Read {len(df)} eccDNA entries from Cresil output")

    # Convert to Cyrcular format
    converted_rows = []

    for _, row in df.iterrows():
        region_str = row["merge_region"]

        # Parse region(s); Cresil uses semicolon/comma to separate segments
        regions = parse_cresil_regions(region_str)

        if len(regions) != row.get("num_region", len(regions)):
            logger.warning(
                "eccDNA %s reports num_region=%s but parsed %s segments",
                row.get("id"),
                row.get("num_region"),
                len(regions),
            )

        # Cyrcular format expects semicolon-separated regions without strand info
        regions_formatted = ";".join(f"{chrom}:{start}-{end}" for chrom, start, end, _ in regions)

        # Calculate approximate abundance from coverage
        # Cresil coverage = totalbase / merge_len
        # We'll use this as a proxy for af_nanopore
        coverage = row.get("coverage", 0)
        af_nanopore = min(coverage / 100.0, 1.0)  # Normalize to 0-1 range

        converted_row = {
            "circle_id": row["id"],
            "regions": regions_formatted,
            "circle_length": row["merge_len"],
            "segment_count": row["num_region"],
            "num_split_reads": row.get("numreads", 0),
            # Cresil doesn't provide probability scores, use defaults
            "prob_present": 0.95,  # High confidence default
            "prob_artifact": 0.0,  # Assume not artifact
            "af_nanopore": af_nanopore,
            # Additional Cresil-specific metadata
            "cresil_ctc": row.get("ctc", True),
            "cresil_totalbase": row.get("totalbase", 0),
            "cresil_coverage": coverage,
        }

        converted_rows.append(converted_row)

    result_df = pd.DataFrame(converted_rows)
    logger.info(f"Converted {len(result_df)} eccDNA entries to Cyrcular-compatible format")

    return result_df


def write_cyrcular_compatible_tsv(cresil_file: Path | str, output_tsv: Path | str) -> Path:
    """Convert Cresil output and write as Cyrcular-compatible TSV.

    Args:
        cresil_file: Path to Cresil eccDNA_final.txt
        output_tsv: Path to output TSV file

    Returns:
        Path to the created TSV file
    """
    logger = get_logger(__name__)

    df = convert_cresil_to_cyrcular_format(cresil_file)

    output_tsv = Path(output_tsv)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    # Write as TSV (tab-separated)
    df.to_csv(output_tsv, sep="\t", index=False)

    logger.info(f"Wrote Cyrcular-compatible TSV to {output_tsv}")

    return output_tsv


def _parse_args():
    """Parse CLI arguments for direct script execution."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Convert Cresil output to Cyrcular-compatible format"
    )
    parser.add_argument("-i", "--input", required=True, help="Input Cresil eccDNA_final.txt file")
    parser.add_argument(
        "-o", "--output", required=True, help="Output TSV file (Cyrcular-compatible format)"
    )
    return parser.parse_args()


def main():
    """CLI entry point."""
    args = _parse_args()
    write_cyrcular_compatible_tsv(args.input, args.output)


if __name__ == "__main__":
    main()
