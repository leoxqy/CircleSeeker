#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ECC Packager — Generate standardized eccDNA output files.

This module transforms internal CircleSeeker data into the standardized output format:

Target output layout
--------------------
<out_dir>/
├── eccDNA_summary.csv          # Main table (one row per eccDNA)
├── eccDNA_regions.csv          # All genomic regions
├── eccDNA_reads.csv            # Read-level support
├── eccDNA_all.fasta            # All sequences
│
├── Uecc/
│   ├── uecc.bed
│   └── uecc.fasta
│
├── Mecc/
│   ├── mecc_sites.bed
│   └── mecc.fasta
│
├── Cecc/
│   ├── cecc_segments.bed
│   ├── cecc_junctions.bedpe
│   └── cecc.fasta
│
├── report.html
└── summary.txt

Key features:
- eccDNA IDs use zero-padded format: UeccDNA0001, MeccDNA0001, CeccDNA0001
- Confirmed entries numbered first, Inferred continue the sequence
- Location format:
  - Uecc: Chr1:100-200(+)
  - Mecc: Chr1:100-200(+)|Chr2:300-400(+)  (| = OR, candidate sites)
  - Cecc: Chr1:100-200(+);Chr2:300-400(+)  (; = AND, segments)
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path
from typing import Optional

import pandas as pd

from circleseeker.utils.logging import get_logger
from circleseeker.modules.ecc_output_formatter import (
    renumber_with_padding,
    generate_summary_table,
    generate_regions_table,
    generate_reads_table,
    generate_uecc_bed,
    generate_mecc_bed,
    generate_cecc_bed,
    generate_cecc_bedpe,
    generate_fasta_files,
)


# ---------------------------- CLI ---------------------------- #


def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Package CircleSeeker outputs into standardized format.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    p.add_argument(
        "-s",
        "--sample-name",
        required=True,
        dest="sample_name",
        help="Sample name for output folder and file prefixes.",
    )

    # Input files
    p.add_argument("--unified-csv", required=True, help="Path to unified eccDNA CSV.")
    p.add_argument("--uecc-dir", help="Directory containing confirmed Uecc files.")
    p.add_argument("--mecc-dir", help="Directory containing confirmed Mecc files.")
    p.add_argument("--cecc-dir", help="Directory containing confirmed Cecc files.")
    p.add_argument("--inferred-dir", help="Directory containing inferred eccDNA files.")

    # Optional report files
    p.add_argument("-H", "--html", help="Path to existing HTML report.")
    p.add_argument("-T", "--text", help="Path to existing TXT summary.")

    # Output
    p.add_argument(
        "-o",
        "--out-dir",
        default=".",
        help="Base output directory.",
    )

    # Options
    p.add_argument("--id-width", type=int, default=4, help="Zero-padding width for IDs.")
    p.add_argument("--dry-run", action="store_true", help="Print actions without executing.")
    p.add_argument("--overwrite", action="store_true", help="Overwrite existing files.")
    p.add_argument("-v", "--verbose", action="store_true", help="Verbose logging.")

    return p


# ------------------------- Utilities ------------------------- #


def log(msg: str, *, verbose: bool = True) -> None:
    if verbose:
        logger = get_logger("ecc_packager")
        logger.info(msg)


def ensure_dir(path: Path, *, dry: bool, verbose: bool) -> None:
    if path.exists():
        return
    log(f"[mkdir] {path}", verbose=verbose)
    if not dry:
        path.mkdir(parents=True, exist_ok=True)


def copy_file(src: Optional[Path], dst: Path, *, overwrite: bool, dry: bool, verbose: bool) -> None:
    if src is None:
        return
    if not src.exists():
        log(f"[warn] Missing source: {src}", verbose=True)
        return
    if dst.exists() and not overwrite:
        log(f"[skip] Exists: {dst}", verbose=verbose)
        return
    ensure_dir(dst.parent, dry=dry, verbose=verbose)
    log(f"[copy] {src} -> {dst}", verbose=verbose)
    if not dry:
        shutil.copy2(src, dst)


def first_glob(dirpath: Optional[Path], pattern: str) -> Optional[Path]:
    if dirpath is None:
        return None
    if not dirpath.exists():
        return None
    matches = sorted(dirpath.glob(pattern))
    return matches[0] if matches else None


def extract_base_id(header: str) -> str:
    """Extract base eccDNA ID from FASTA header.

    Examples:
        UeccDNA1|Chr1:100-200(+)|length=1000 -> UeccDNA1
        MeccDNA_000001.1_1 -> MeccDNA1
        CeccDNA_000001 -> CeccDNA1
    """
    import re

    # Remove > prefix if present
    if header.startswith(">"):
        header = header[1:]

    # Get first part (before space or |)
    first_part = header.split()[0].split("|")[0]

    # Extract type and number
    # Pattern: (UeccDNA|MeccDNA|CeccDNA)(_?)(\d+)
    match = re.match(r"(UeccDNA|MeccDNA|CeccDNA)_?(\d+)", first_part)
    if match:
        return f"{match.group(1)}{int(match.group(2))}"

    return first_part


def load_fasta(fasta_path: Path) -> dict[str, str]:
    """Load FASTA file into dict, extracting base IDs."""
    seqs: dict[str, str] = {}
    if not fasta_path or not fasta_path.exists():
        return seqs

    current_id: Optional[str] = None
    current_seq: list[str] = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id:
                    seqs[current_id] = "".join(current_seq)
                # Extract base ID from full header
                current_id = extract_base_id(line)
                current_seq = []
            else:
                current_seq.append(line)

        if current_id:
            seqs[current_id] = "".join(current_seq)

    return seqs


# ---------------------- Core packaging ----------------------- #


def run(args: argparse.Namespace) -> int:
    logger = get_logger("ecc_packager")
    sample = args.sample_name
    verbose = bool(args.verbose)
    dry = args.dry_run

    # Parse paths
    unified_csv = Path(args.unified_csv) if args.unified_csv else None
    uecc_dir = Path(args.uecc_dir) if args.uecc_dir else None
    mecc_dir = Path(args.mecc_dir) if args.mecc_dir else None
    cecc_dir = Path(args.cecc_dir) if args.cecc_dir else None
    inferred_dir = Path(args.inferred_dir) if args.inferred_dir else None
    html_in = Path(args.html) if args.html else None
    txt_in = Path(args.text) if args.text else None

    # Output directory
    out_dir = Path(args.out_dir)

    if dry:
        logger.info("[dry-run] Would create output structure in: %s", out_dir)
        return 0

    # Create output directories
    ensure_dir(out_dir, dry=dry, verbose=verbose)
    ensure_dir(out_dir / "Uecc", dry=dry, verbose=verbose)
    ensure_dir(out_dir / "Mecc", dry=dry, verbose=verbose)
    ensure_dir(out_dir / "Cecc", dry=dry, verbose=verbose)

    # Load unified CSV
    if not unified_csv or not unified_csv.exists():
        logger.error("Unified CSV not found: %s", unified_csv)
        return 1

    unified_df = pd.read_csv(unified_csv)
    logger.info("Loaded unified CSV with %d entries", len(unified_df))

    # Renumber IDs with padding
    unified_df = renumber_with_padding(unified_df, args.id_width)

    # Build ID mapping: map both old eccDNA_id and original_id (for Inferred) to new IDs
    id_map = dict(zip(unified_df["_old_eccDNA_id"], unified_df["eccDNA_id"]))

    # Also map original_id (e.g., ICeccDNA1) to new IDs for Inferred entries
    if "original_id" in unified_df.columns:
        orig_id_map = dict(zip(unified_df["original_id"], unified_df["eccDNA_id"]))
        id_map.update(orig_id_map)

    # Load type-specific CSVs (Confirmed)
    uecc_core_csv = first_glob(uecc_dir, "*_UeccDNA.core.csv")
    mecc_sites_csv = first_glob(mecc_dir, "*_MeccSites.core.csv")
    cecc_segments_csv = first_glob(cecc_dir, "*_CeccSegments.core.csv")

    uecc_df = pd.read_csv(uecc_core_csv) if uecc_core_csv else None
    mecc_df = pd.read_csv(mecc_sites_csv) if mecc_sites_csv else None
    cecc_df = pd.read_csv(cecc_segments_csv) if cecc_segments_csv else None

    # Load Inferred CSVs
    inferred_simple_csv = first_glob(inferred_dir, "*_simple.csv") if inferred_dir else None
    inferred_chimeric_csv = first_glob(inferred_dir, "*_chimeric.csv") if inferred_dir else None

    inferred_simple_df = pd.read_csv(inferred_simple_csv) if inferred_simple_csv else None
    inferred_chimeric_df = pd.read_csv(inferred_chimeric_csv) if inferred_chimeric_csv else None

    # Update IDs in type-specific tables
    if uecc_df is not None:
        uecc_df["eccDNA_id"] = uecc_df["eccDNA_id"].map(id_map).fillna(uecc_df["eccDNA_id"])
    if mecc_df is not None:
        mecc_df["eccDNA_id"] = mecc_df["eccDNA_id"].map(id_map).fillna(mecc_df["eccDNA_id"])
    if cecc_df is not None:
        cecc_df["eccDNA_id"] = cecc_df["eccDNA_id"].map(id_map).fillna(cecc_df["eccDNA_id"])
    if inferred_simple_df is not None:
        inferred_simple_df["eccDNA_id"] = (
            inferred_simple_df["eccDNA_id"].map(id_map).fillna(inferred_simple_df["eccDNA_id"])
        )
    if inferred_chimeric_df is not None:
        inferred_chimeric_df["eccDNA_id"] = (
            inferred_chimeric_df["eccDNA_id"].map(id_map).fillna(inferred_chimeric_df["eccDNA_id"])
        )

    # Drop inferred entries that were not kept in the unified table.
    # These are typically redundant inferred calls that overlap confirmed eccDNA and were removed
    # during ecc_unify, but still exist in the raw inferred tables.
    valid_ids = set(unified_df["eccDNA_id"].astype(str))
    if inferred_simple_df is not None and not inferred_simple_df.empty:
        inferred_simple_df = inferred_simple_df[inferred_simple_df["eccDNA_id"].isin(valid_ids)].copy()
        if inferred_simple_df.empty:
            inferred_simple_df = None
    if inferred_chimeric_df is not None and not inferred_chimeric_df.empty:
        inferred_chimeric_df = inferred_chimeric_df[
            inferred_chimeric_df["eccDNA_id"].isin(valid_ids)
        ].copy()
        if inferred_chimeric_df.empty:
            inferred_chimeric_df = None

    # Generate regions table (include inferred data)
    regions_df = generate_regions_table(
        uecc_df, mecc_df, cecc_df,
        inferred_simple_df=inferred_simple_df,
        inferred_chimeric_df=inferred_chimeric_df,
    )

    # Generate summary table
    summary_df = generate_summary_table(unified_df, regions_df)

    # Generate reads table
    reads_df = generate_reads_table(uecc_df, mecc_df, cecc_df)

    # Save CSV files
    summary_df.to_csv(out_dir / "eccDNA_summary.csv", index=False)
    regions_df.to_csv(out_dir / "eccDNA_regions.csv", index=False)
    reads_df.to_csv(out_dir / "eccDNA_reads.csv", index=False)

    logger.info("Saved eccDNA_summary.csv (%d rows)", len(summary_df))
    logger.info("Saved eccDNA_regions.csv (%d rows)", len(regions_df))
    logger.info("Saved eccDNA_reads.csv (%d rows)", len(reads_df))

    # Generate BED files
    generate_uecc_bed(regions_df, out_dir / "Uecc" / "uecc.bed")
    generate_mecc_bed(regions_df, out_dir / "Mecc" / "mecc_sites.bed")
    generate_cecc_bed(regions_df, out_dir / "Cecc" / "cecc_segments.bed")
    generate_cecc_bedpe(regions_df, out_dir / "Cecc" / "cecc_junctions.bedpe")

    # Load and process FASTA files
    sequences = {}

    # Load confirmed FASTA files
    uecc_fasta = first_glob(uecc_dir, "*_UeccDNA_C.fasta")
    mecc_fasta = first_glob(mecc_dir, "*_MeccDNA_C.fasta")
    cecc_fasta = first_glob(cecc_dir, "*_CeccDNA_C.fasta")

    for fasta_path in [uecc_fasta, mecc_fasta, cecc_fasta]:
        if fasta_path:
            for old_id, seq in load_fasta(fasta_path).items():
                new_id = id_map.get(old_id, old_id)
                sequences[new_id] = seq

    # Load inferred FASTA files
    if inferred_dir:
        inferred_uecc_fasta = first_glob(inferred_dir, "*_UeccDNA_I.fasta")
        inferred_cecc_fasta = first_glob(inferred_dir, "*_CeccDNA_I.fasta")

        for fasta_path in [inferred_uecc_fasta, inferred_cecc_fasta]:
            if fasta_path:
                for old_id, seq in load_fasta(fasta_path).items():
                    new_id = id_map.get(old_id, old_id)
                    sequences[new_id] = seq

    # Generate FASTA files
    generate_fasta_files(sequences, out_dir, summary_df)

    # Copy report files
    if html_in:
        copy_file(html_in, out_dir / "report.html", overwrite=args.overwrite, dry=dry, verbose=verbose)
    if txt_in:
        copy_file(txt_in, out_dir / "summary.txt", overwrite=args.overwrite, dry=dry, verbose=verbose)

    logger.info("[done] Output folder: %s", out_dir)
    return 0


# --------------------------- Entrypoint ------------------------ #
if __name__ == "__main__":
    args = build_argparser().parse_args()
    raise SystemExit(run(args))
