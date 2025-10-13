"""Curate and split Cyrcular overview TSV into simple/chimeric eccDNA tables.

This module processes Cyrcular output and generates:
1. CSV tables for simple and chimeric eccDNA candidates
2. Optional FASTA sequences extracted from reference genome
   - Simple eccDNA: single sequence from reference
   - Chimeric eccDNA: concatenated sequences from multiple segments
"""

from __future__ import annotations

import logging
from pathlib import Path
import re
from typing import Tuple, Optional

import pandas as pd
from circleseeker.utils.logging import get_logger

try:
    import pysam
    PYSAM_AVAILABLE = True
except ImportError:
    PYSAM_AVAILABLE = False
    pysam = None


def select_best_duplicate(group: pd.DataFrame) -> pd.DataFrame:
    group = group.copy()

    if "prob_present" in group.columns:
        prob_col = "prob_present"
    elif "prob_joint_event" in group.columns:
        prob_col = "prob_joint_event"
    else:
        prob_col = None

    if prob_col:
        max_splits = group["num_split_reads"].max()
        if max_splits > 0:
            group["score"] = (
                group[prob_col] * 0.5
                + (group["num_split_reads"] / max_splits) * 0.3
                + (1 - group.get("prob_artifact", 0)) * 0.2
            )
        else:
            group["score"] = group[prob_col] * 0.7 + (1 - group.get("prob_artifact", 0)) * 0.3
    else:
        max_splits = group["num_split_reads"].max()
        if max_splits > 0:
            group["score"] = (
                (group["num_split_reads"] / max_splits) * 0.6
                + (1 - group.get("prob_artifact", 0)) * 0.4
            )
        else:
            group["score"] = 1 - group.get("prob_artifact", 0)

    return group.nlargest(1, "score")


def curate_ecc_tables(input_tsv: Path | str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    df = pd.read_csv(input_tsv, sep="\t")
    df = df[df["circle_length"] >= 100]

    df_list = [select_best_duplicate(g) for _, g in df.groupby("regions")]
    df_dedup = pd.concat(df_list, ignore_index=True) if df_list else pd.DataFrame()

    simple_circles = []
    chimeric_circles = []
    simple_idx = 1
    chimeric_idx = 1

    for _, row in df_dedup.iterrows():
        if pd.isna(row.get("regions")) or pd.isna(row.get("segment_count")):
            continue
        regions_field = str(row["regions"]).strip()
        regions_list = [part.strip() for part in re.split(r"[;,]+", regions_field) if part.strip()]
        if not regions_list:
            continue

        if int(row["segment_count"]) == 1:
            region = regions_list[0]
            if ":" not in region or "-" not in region:
                continue
            chrom, coords = region.split(":")
            start, end = coords.split("-")
            start_pos, end_pos = int(start), int(end)

            # Determine strand based on original coordinate order
            if start_pos <= end_pos:
                strand = "+"
                final_start, final_end = start_pos, end_pos
            else:
                strand = "-"
                final_start, final_end = end_pos, start_pos

            simple_circles.append(
                {
                    "eccDNA_id": f"IUeccDNA{simple_idx}",
                    "chr": chrom,
                    "start0": final_start - 1,
                    "end0": final_end,
                    "strand": strand,
                    "length": row["circle_length"],
                    "eccdna_type": "Uecc",
                    "state": "Inferred",
                    "num_split_reads": row.get("num_split_reads", 0),
                    "prob_present": row.get("prob_present", row.get("prob_joint_event", 0)),
                    "prob_artifact": row.get("prob_artifact", 0),
                    "hifi_abundance": row.get("af_nanopore", 0),
                }
            )
            simple_idx += 1
        else:
            total_segments = len(regions_list)
            is_valid = all(":" in r and "-" in r for r in regions_list)
            if not is_valid:
                continue
            for seg_idx, region in enumerate(regions_list, 1):
                chrom, coords = region.strip().split(":")
                start, end = coords.split("-")
                start_pos, end_pos = int(start), int(end)

                # Determine strand based on original coordinate order
                if start_pos <= end_pos:
                    strand = "+"
                    final_start, final_end = start_pos, end_pos
                else:
                    strand = "-"
                    final_start, final_end = end_pos, start_pos

                junction_role = (
                    "head" if seg_idx == 1 else "tail" if seg_idx == total_segments else "middle"
                )
                chimeric_circles.append(
                    {
                        "eccDNA_id": f"ICeccDNA{chimeric_idx}",
                        "chr": chrom,
                        "start0": final_start - 1,
                        "end0": final_end,
                        "strand": strand,
                        "length": row["circle_length"],
                        "eccdna_type": "Cecc",
                        "state": "Inferred",
                        "seg_index": seg_idx,
                        "seg_total": total_segments,
                        "junction_role": junction_role,
                        "read_count": 1,
                        "num_split_reads": row.get("num_split_reads", 0),
                        "prob_present": row.get("prob_present", row.get("prob_joint_event", 0)),
                        "prob_artifact": row.get("prob_artifact", 0),
                        "hifi_abundance": row.get("af_nanopore", 0),
                    }
                )
            chimeric_idx += 1

    simple_df = pd.DataFrame(simple_circles)
    chimeric_df = pd.DataFrame(chimeric_circles)
    return simple_df, chimeric_df


def generate_fasta_sequences(
    simple_df: pd.DataFrame,
    chimeric_df: pd.DataFrame,
    reference_fasta: Path | str,
    output_prefix: Path | str
) -> Tuple[Optional[Path], Optional[Path]]:
    """
    Generate FASTA files for inferred eccDNA sequences.

    Args:
        simple_df: Simple eccDNA DataFrame
        chimeric_df: Chimeric eccDNA DataFrame
        reference_fasta: Path to reference genome FASTA
        output_prefix: Output file prefix

    Returns:
        Tuple of (simple_fasta_path, chimeric_fasta_path) or (None, None) if pysam unavailable
    """
    logger = get_logger(__name__)

    if not PYSAM_AVAILABLE:
        logger.warning("pysam not available. Skipping FASTA generation.")
        return None, None

    reference_fasta = Path(reference_fasta)
    if not reference_fasta.exists():
        logger.warning(f"Reference FASTA not found: {reference_fasta}")
        return None, None
        
    output_prefix = Path(output_prefix)
    simple_fasta = output_prefix.parent / f"{output_prefix.name}_UeccDNA_I.fasta"
    chimeric_fasta = output_prefix.parent / f"{output_prefix.name}_CeccDNA_I.fasta"
    
    try:
        # Open reference genome
        fasta = pysam.FastaFile(str(reference_fasta))
        
        # Generate simple eccDNA FASTA
        simple_fasta_path = None
        if not simple_df.empty:
            with open(simple_fasta, 'w') as f:
                for _, row in simple_df.iterrows():
                    chr_name = str(row['chr'])
                    start = int(row['start0'])
                    end = int(row['end0'])
                    eccDNA_id = row['eccDNA_id']
                    length = int(row['length'])
                    
                    try:
                        # Extract sequence (pysam uses 0-based coordinates)
                        sequence = fasta.fetch(chr_name, start, end)
                        if sequence:
                            header = f">{eccDNA_id}|{chr_name}:{start}-{end}|length={length}|type=simple"
                            f.write(header + "\n")
                            # Write sequence in 60-character lines
                            for i in range(0, len(sequence), 60):
                                f.write(sequence[i:i+60] + "\n")
                        else:
                            logger.warning(f"Could not extract sequence for {eccDNA_id}")
                    except Exception as e:
                        logger.warning(f"Error extracting sequence for {eccDNA_id}: {e}")
            simple_fasta_path = simple_fasta
                        
        # Generate chimeric eccDNA FASTA (concatenated segments)
        chimeric_fasta_path = None
        if not chimeric_df.empty:
            with open(chimeric_fasta, 'w') as f:
                # Group by eccDNA_id and process each chimeric circle
                for eccDNA_id, group in chimeric_df.groupby('eccDNA_id'):
                    # Sort by seg_index to maintain order
                    group = group.sort_values('seg_index')
                    
                    segments = []
                    total_length = 0
                    seg_info = []
                    
                    for _, row in group.iterrows():
                        chr_name = str(row['chr'])
                        start = int(row['start0'])
                        end = int(row['end0'])
                        seg_idx = int(row['seg_index'])
                        
                        try:
                            # Extract sequence for this segment
                            sequence = fasta.fetch(chr_name, start, end)
                            if sequence:
                                segments.append(sequence)
                                total_length += len(sequence)
                                seg_info.append(f"seg{seg_idx}:{chr_name}:{start}-{end}")
                            else:
                                logger.warning(f"Could not extract segment {seg_idx} for {eccDNA_id}")
                        except Exception as e:
                            logger.warning(f"Error extracting segment {seg_idx} for {eccDNA_id}: {e}")
                    
                    if segments:
                        # Concatenate all segments
                        full_sequence = "".join(segments)
                        seg_count = len(segments)
                        header = f">{eccDNA_id}|segments={seg_count}|length={total_length}|type=chimeric|{';'.join(seg_info)}"
                        f.write(header + "\n")
                        # Write sequence in 60-character lines
                        for i in range(0, len(full_sequence), 60):
                            f.write(full_sequence[i:i+60] + "\n")
            chimeric_fasta_path = chimeric_fasta
            
        fasta.close()
        return simple_fasta_path, chimeric_fasta_path

    except Exception as e:
        logger.error(f"Error generating FASTA sequences: {e}")
        return None, None


def write_curated_tables_with_fasta(
    simple_df: pd.DataFrame, 
    chimeric_df: pd.DataFrame, 
    output_prefix: Path | str,
    reference_fasta: Optional[Path | str] = None,
    organize_files: bool = True
) -> tuple[Path, Path, Optional[Path], Optional[Path]]:
    """
    Write curated tables and optionally generate FASTA sequences.
    
    Args:
        simple_df: Simple eccDNA DataFrame
        chimeric_df: Chimeric eccDNA DataFrame
        output_prefix: Output file prefix
        reference_fasta: Optional path to reference genome for FASTA generation
        organize_files: If True, organize outputs into Inferred_eccDNA folder
        
    Returns:
        Tuple of (simple_csv, chimeric_csv, simple_fasta, chimeric_fasta)
    """
    output_prefix = Path(output_prefix)
    output_dir = output_prefix.parent
    
    # Create Inferred_eccDNA folder if organize_files is True
    if organize_files:
        inferred_dir = output_dir / f"{output_prefix.name}_Inferred_eccDNA"
        inferred_dir.mkdir(parents=True, exist_ok=True)
        # Update paths to write into the organized folder
        simple_output = inferred_dir / f"{output_prefix.name}_simple.csv"
        chimeric_output = inferred_dir / f"{output_prefix.name}_chimeric.csv"
    else:
        # Use original behavior
        simple_output = output_dir / f"{output_prefix.name}_simple.csv"
        chimeric_output = output_dir / f"{output_prefix.name}_chimeric.csv"

    # Write CSV files
    if not simple_df.empty:
        simple_df.to_csv(simple_output, sep=",", index=False)
    if not chimeric_df.empty:
        chimeric_df.to_csv(chimeric_output, sep=",", index=False)
    
    # Generate FASTA files if reference provided
    simple_fasta, chimeric_fasta = None, None
    if reference_fasta:
        # For FASTA generation, we need to update the output_prefix to point to the new directory
        if organize_files:
            fasta_prefix = inferred_dir / output_prefix.name
        else:
            fasta_prefix = output_prefix
            
        simple_fasta, chimeric_fasta = generate_fasta_sequences(
            simple_df, chimeric_df, reference_fasta, fasta_prefix
        )
        
    return simple_output, chimeric_output, simple_fasta, chimeric_fasta


def write_curated_tables(
    simple_df: pd.DataFrame, 
    chimeric_df: pd.DataFrame, 
    output_prefix: Path | str
) -> tuple[Path, Path]:
    """
    Write curated tables (backward compatible version).
    
    Args:
        simple_df: Simple eccDNA DataFrame
        chimeric_df: Chimeric eccDNA DataFrame
        output_prefix: Output file prefix
        
    Returns:
        Tuple of (simple_csv, chimeric_csv)
    """
    # For backward compatibility, don't organize files
    results = write_curated_tables_with_fasta(simple_df, chimeric_df, output_prefix, organize_files=False)
    return results[0], results[1]  # Return only CSV paths for backward compatibility


# Backward-compatible helper name
def process_eccDNA(
    input_tsv: Path | str,
    output_prefix: Path | str,
    reference_fasta: Optional[Path | str] = None,
    organize_files: bool = True
) -> None:
    """
    Process eccDNA data with optional FASTA generation.

    Args:
        input_tsv: Input Cyrcular overview TSV file
        output_prefix: Output file prefix
        reference_fasta: Optional path to reference genome for FASTA generation
        organize_files: If True, organize outputs into Inferred_eccDNA folder
    """
    logger = get_logger(__name__)
    simple_df, chimeric_df = curate_ecc_tables(input_tsv)

    if reference_fasta:
        # Use full version with FASTA generation
        tsv_files = write_curated_tables_with_fasta(
            simple_df, chimeric_df, output_prefix, reference_fasta, organize_files
        )
        logger.info("Generated files:")
        if organize_files:
            logger.info(f"  Output folder: {Path(output_prefix).name}_Inferred_eccDNA/")
        logger.info(f"  Simple CSV: {tsv_files[0].name}")
        logger.info(f"  Chimeric CSV: {tsv_files[1].name}")
        if tsv_files[2]:
            logger.info(f"  Simple FASTA: {tsv_files[2].name}")
        if tsv_files[3]:
            logger.info(f"  Chimeric FASTA: {tsv_files[3].name}")
    else:
        # Without reference, still organize files if requested
        if organize_files:
            tsv_files = write_curated_tables_with_fasta(
                simple_df, chimeric_df, output_prefix, organize_files=organize_files
            )
            logger.info("Generated files:")
            logger.info(f"  Output folder: {Path(output_prefix).name}_Inferred_eccDNA/")
            logger.info(f"  Simple CSV: {tsv_files[0].name}")
            logger.info(f"  Chimeric CSV: {tsv_files[1].name}")
        else:
            # Use backward compatible version
            write_curated_tables(simple_df, chimeric_df, output_prefix)


def _parse_args():
    """Parse CLI arguments for direct script execution."""
    import argparse
    parser = argparse.ArgumentParser(
        description="Curate Cyrcular overview TSV into inferred eccDNA tables and optional FASTA"
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input Cyrcular overview TSV file"
    )
    parser.add_argument(
        "-o", "--output-prefix",
        required=True,
        help="Output prefix (no extension)"
    )
    parser.add_argument(
        "-r", "--reference-fasta",
        required=False,
        default=None,
        help="Optional: Reference FASTA for sequence generation"
    )
    parser.add_argument(
        "--no-organize",
        action="store_true",
        help="Write directly to output dir without creating subfolder"
    )
    return parser.parse_args()


def main():
    args = _parse_args()
    inp = Path(args.input)
    out_prefix = Path(args.output_prefix)
    ref = args.reference_fasta
    organize = not args.no_organize

    process_eccDNA(inp, out_prefix, reference_fasta=ref, organize_files=organize)


if __name__ == "__main__":
    main()
