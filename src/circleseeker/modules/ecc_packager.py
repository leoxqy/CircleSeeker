#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ECC Packager — copy & rename only (explicit 4 source dirs)

What it does
------------
- Takes **four source directories** explicitly: --uecc-dir, --mecc-dir, --cecc-dir, --inferred-dir
- Takes **three single files**: --merged-csv, --html, --text (all required)
- Copies files into a unified output structure under <out>/<sample_name>/ with standardized names
- **No** new report generation; only renaming & copying what you already have

Target output layout
--------------------
<out_dir>/<sample_name>/
├── <sample_name>_Confirmed_UeccDNA/
│   ├── <sample_name>_UeccDNA_C.fasta
│   ├── <sample_name>_UeccDNA.bed
│   └── <sample_name>_UeccDNA.core.csv
│
├── <sample_name>_Confirmed_MeccDNA/
│   ├── <sample_name>_MeccDNA_C.fasta
│   ├── <sample_name>_MeccSites.bed
│   ├── <sample_name>_MeccBestSite.bed
│   └── <sample_name>_MeccSites.core.csv
│
├── <sample_name>_Confirmed_CeccDNA/
│   ├── <sample_name>_CeccDNA_C.fasta
│   ├── <sample_name>_CeccSegments.bed
│   ├── <sample_name>_CeccJunctions.bedpe
│   └── <sample_name>_CeccSegments.core.csv
│
├── <sample_name>_Inferred_eccDNA/
│   ├── <sample_name>_UeccDNA_I.csv
│   ├── <sample_name>_UeccDNA_I.fasta
│   ├── <sample_name>_chimeric.csv
│   └── <sample_name>_CeccDNA_I.fasta
│
├── <sample_name>_merged_output.csv
├── <sample_name>_report.html
└── <sample_name>_summary.txt

Usage example (matches your paths)
----------------------------------
python ecc_packager.py \
  -s sample_name \
  --uecc-dir sampleecc_C/sample_Uecc_C \
  --mecc-dir sampleecc_C/sample_Mecc_C \
  --cecc-dir sampleecc_C/sample_Cecc_C \
  --inferred-dir sample_Inferred_eccDNA \
  -m sample_merged_output.csv \
  -H sample_report.html \
  -T sample_summary.txt \
  -o ./out \
  -v

Notes
-----
- CLI 使用时依然要求提供四个目录和三个文件；若通过管线调用而缺少其中任意项，会记录警告并跳过对应拷贝。
- Missing individual files inside a source dir will be warned and skipped.
- Use --overwrite to replace existing destination files; --dry-run to preview.
"""

from __future__ import annotations
import argparse
from pathlib import Path
import shutil
from typing import Optional
from circleseeker.utils.logging import get_logger

# ---------------------------- CLI ---------------------------- #


def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Copy and rename CircleSeeker outputs into a standard layout (no generation).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    p.add_argument(
        "-s",
        "--sample-name",
        required=True,
        dest="sample_name",
        help="Sample name for output folder and file prefixes.",
    )

    # Explicit source directories (all required)
    p.add_argument("--uecc-dir", required=True, help="Directory containing confirmed Uecc files.")
    p.add_argument("--mecc-dir", required=True, help="Directory containing confirmed Mecc files.")
    p.add_argument("--cecc-dir", required=True, help="Directory containing confirmed Cecc files.")
    p.add_argument(
        "--inferred-dir", required=True, help="Directory containing inferred eccDNA files."
    )

    # Required single-file inputs
    p.add_argument("-m", "--merged-csv", required=True, help="Path to merged eccDNA CSV.")
    p.add_argument("-H", "--html", required=True, help="Path to existing HTML report.")
    p.add_argument("-T", "--text", required=True, help="Path to existing TXT summary.")

    # Output
    p.add_argument(
        "-o",
        "--out-dir",
        default=".",
        help="Base output directory (the <sample_name> folder will be created inside).",
    )

    # Behavior
    p.add_argument("--dry-run", action="store_true", help="Print actions without copying.")
    p.add_argument("--overwrite", action="store_true", help="Overwrite existing destination files.")
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
        log(f"[warn] Missing source: <None> -> {dst}", verbose=True)
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
    matches = sorted(dirpath.glob(pattern))
    return matches[0] if matches else None


# ---------------------- Core packaging ----------------------- #


def run(args) -> int:
    sample = args.sample_name

    # Sources
    uecc_dir = Path(args.uecc_dir) if args.uecc_dir else None
    mecc_dir = Path(args.mecc_dir) if args.mecc_dir else None
    cecc_dir = Path(args.cecc_dir) if args.cecc_dir else None
    inferred_dir = Path(args.inferred_dir) if args.inferred_dir else None

    merged_csv = Path(args.merged_csv) if args.merged_csv else None
    html_in = Path(args.html) if args.html else None
    txt_in = Path(args.text) if args.text else None

    # Destinations
    base = Path(args.out_dir) / sample
    U_dir = base / f"{sample}_Confirmed_UeccDNA"
    M_dir = base / f"{sample}_Confirmed_MeccDNA"
    C_dir = base / f"{sample}_Confirmed_CeccDNA"
    I_dir = base / f"{sample}_Inferred_eccDNA"

    dest_html = base / f"{sample}_report.html"
    dest_txt = base / f"{sample}_summary.txt"
    dest_merged = base / f"{sample}_merged_output.csv"

    # Create folders
    for d in [base, U_dir, M_dir, C_dir, I_dir]:
        ensure_dir(d, dry=args.dry_run, verbose=bool(args.verbose))

    # Copy Confirmed U
    copy_file(
        first_glob(uecc_dir, "*_UeccDNA_C.fasta"),
        U_dir / f"{sample}_UeccDNA_C.fasta",
        overwrite=args.overwrite,
        dry=args.dry_run,
        verbose=bool(args.verbose),
    )
    copy_file(
        first_glob(uecc_dir, "*_UeccDNA.bed"),
        U_dir / f"{sample}_UeccDNA.bed",
        overwrite=args.overwrite,
        dry=args.dry_run,
        verbose=bool(args.verbose),
    )
    copy_file(
        first_glob(uecc_dir, "*_UeccDNA.core.csv"),
        U_dir / f"{sample}_UeccDNA.core.csv",
        overwrite=args.overwrite,
        dry=args.dry_run,
        verbose=bool(args.verbose),
    )

    # Copy Confirmed M
    copy_file(
        first_glob(mecc_dir, "*_MeccDNA_C.fasta"),
        M_dir / f"{sample}_MeccDNA_C.fasta",
        overwrite=args.overwrite,
        dry=args.dry_run,
        verbose=bool(args.verbose),
    )
    copy_file(
        first_glob(mecc_dir, "*_MeccSites.bed"),
        M_dir / f"{sample}_MeccSites.bed",
        overwrite=args.overwrite,
        dry=args.dry_run,
        verbose=bool(args.verbose),
    )
    copy_file(
        first_glob(mecc_dir, "*_MeccBestSite.bed"),
        M_dir / f"{sample}_MeccBestSite.bed",
        overwrite=args.overwrite,
        dry=args.dry_run,
        verbose=bool(args.verbose),
    )
    copy_file(
        first_glob(mecc_dir, "*_MeccSites.core.csv"),
        M_dir / f"{sample}_MeccSites.core.csv",
        overwrite=args.overwrite,
        dry=args.dry_run,
        verbose=bool(args.verbose),
    )

    # Copy Confirmed C
    copy_file(
        first_glob(cecc_dir, "*_CeccDNA_C.fasta"),
        C_dir / f"{sample}_CeccDNA_C.fasta",
        overwrite=args.overwrite,
        dry=args.dry_run,
        verbose=bool(args.verbose),
    )
    copy_file(
        first_glob(cecc_dir, "*_CeccSegments.bed"),
        C_dir / f"{sample}_CeccSegments.bed",
        overwrite=args.overwrite,
        dry=args.dry_run,
        verbose=bool(args.verbose),
    )
    copy_file(
        first_glob(cecc_dir, "*_CeccJunctions.bedpe"),
        C_dir / f"{sample}_CeccJunctions.bedpe",
        overwrite=args.overwrite,
        dry=args.dry_run,
        verbose=bool(args.verbose),
    )
    copy_file(
        first_glob(cecc_dir, "*_CeccSegments.core.csv"),
        C_dir / f"{sample}_CeccSegments.core.csv",
        overwrite=args.overwrite,
        dry=args.dry_run,
        verbose=bool(args.verbose),
    )

    # Copy Inferred
    copy_file(
        first_glob(inferred_dir, "*_simple.csv"),
        I_dir / f"{sample}_UeccDNA_I.csv",
        overwrite=args.overwrite,
        dry=args.dry_run,
        verbose=bool(args.verbose),
    )
    copy_file(
        first_glob(inferred_dir, "*_UeccDNA_I.fasta"),
        I_dir / f"{sample}_UeccDNA_I.fasta",
        overwrite=args.overwrite,
        dry=args.dry_run,
        verbose=bool(args.verbose),
    )
    copy_file(
        first_glob(inferred_dir, "*_chimeric.csv"),
        I_dir / f"{sample}_chimeric.csv",
        overwrite=args.overwrite,
        dry=args.dry_run,
        verbose=bool(args.verbose),
    )
    copy_file(
        first_glob(inferred_dir, "*_CeccDNA_I.fasta"),
        I_dir / f"{sample}_CeccDNA_I.fasta",
        overwrite=args.overwrite,
        dry=args.dry_run,
        verbose=bool(args.verbose),
    )

    # Copy required single-file inputs
    copy_file(
        merged_csv,
        dest_merged,
        overwrite=args.overwrite,
        dry=args.dry_run,
        verbose=bool(args.verbose),
    )
    copy_file(
        html_in, dest_html, overwrite=args.overwrite, dry=args.dry_run, verbose=bool(args.verbose)
    )
    copy_file(
        txt_in, dest_txt, overwrite=args.overwrite, dry=args.dry_run, verbose=bool(args.verbose)
    )

    log(f"[done] Output folder: {base}", verbose=True)
    return 0


# --------------------------- Entrypoint ------------------------ #
if __name__ == "__main__":
    args = build_argparser().parse_args()
    raise SystemExit(run(args))
