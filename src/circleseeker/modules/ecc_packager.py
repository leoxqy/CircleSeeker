#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ECC Packager — Generate standardized eccDNA output files.

This module transforms internal CircleSeeker data into the standardized output format:

Target output layout
--------------------
<out_dir>/
├── {prefix}_eccDNA_summary.csv   # Main table (one row per eccDNA)
├── {prefix}_eccDNA_regions.csv   # All genomic regions
├── {prefix}_eccDNA_reads.csv     # Read-level support
├── {prefix}_eccDNA_all.fasta     # All sequences
│
├── UeccDNA/
│   ├── {prefix}_uecc.bed
│   └── {prefix}_uecc.fasta
│
├── MeccDNA/
│   ├── {prefix}_mecc_sites.bed
│   └── {prefix}_mecc.fasta
│
├── CeccDNA/
│   ├── {prefix}_cecc_segments.bed
│   ├── {prefix}_cecc_junctions.bedpe
│   └── {prefix}_cecc.fasta
│
├── {prefix}_report.html
└── {prefix}_summary.txt

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
import math
import re
import shutil
from pathlib import Path
from typing import Optional

import numpy as np
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


# ── MeccDNA ONT LR model coefficients (trained on ara rep1 10X) ────
# Features derived at packaging time from summary, reads, and regions tables.
MECC_ONT_FEATURE_NAMES = [
    "mean_mapq",           # mean mapping quality of supporting reads
    "n_candidate_sites",   # number of candidate genomic sites
    "confidence",          # pipeline confidence score
    "mean_match_degree",   # mean match degree of reads
    "copy_number",         # estimated copy number
    "log_length",          # log10(eccDNA length)
]
_MECC_ONT_MODEL_COEFS = np.array([
    -0.7251,   # mean_mapq          (low mapq -> multi-copy -> TP)
    -0.9180,   # n_candidate_sites  (fewer sites -> more confident -> TP)
    -0.3572,   # confidence         (inverted: high conf often FP in ONT)
    +0.8749,   # mean_match_degree  (high match -> TP)
    +0.2493,   # copy_number        (high CN -> TP)
    +0.4574,   # log_length         (longer -> TP)
])
_MECC_ONT_MODEL_INTERCEPT = 2.7477
_MECC_ONT_SCALER_MEAN = np.array([
    13.1908, 3.4677, 0.3386, 98.4415, 6.4218, 2.807,
])
_MECC_ONT_SCALER_STD = np.array([
    17.8667, 1.7339, 0.3409, 1.8683, 4.386, 0.2644,
])
_MECC_ONT_DEFAULT_THRESHOLD = 0.81


def _sigmoid(x: float) -> float:
    return 1.0 / (1.0 + math.exp(-x)) if x > -500 else 0.0


def score_mecc_ont_lr(
    summary_df: pd.DataFrame,
    reads_df: pd.DataFrame,
    regions_df: pd.DataFrame,
    threshold: float = _MECC_ONT_DEFAULT_THRESHOLD,
) -> pd.DataFrame:
    """Apply MeccDNA ONT LR model to score and filter Mecc candidates.

    Low-scoring Mecc entries are demoted to Uecc (not deleted).
    Features are derived from summary, reads, and regions tables.

    Args:
        summary_df: Summary table (one row per eccDNA).
        reads_df: Read-level table with identity, mapq, match_degree columns.
        regions_df: Regions table with one row per genomic site.
        threshold: LR probability threshold; Mecc with score < threshold
            are demoted to Uecc.

    Returns:
        Updated summary_df with demoted entries re-typed as Uecc.
    """
    logger = get_logger("ecc_packager")

    if summary_df.empty or threshold <= 0.0:
        return summary_df

    mecc_mask = summary_df["type"] == "Mecc"
    n_mecc = int(mecc_mask.sum())
    if n_mecc == 0:
        return summary_df

    summary_df = summary_df.copy()

    # --- Aggregate read-level features for each Mecc eccDNA ---
    mecc_ids = set(summary_df.loc[mecc_mask, "eccDNA_id"])

    read_features: dict[str, dict[str, float]] = {}
    if not reads_df.empty:
        mecc_reads = reads_df[reads_df["eccDNA_id"].isin(mecc_ids)].copy()
        if not mecc_reads.empty:
            for col in ("mapq", "match_degree"):
                mecc_reads[col] = pd.to_numeric(mecc_reads[col], errors="coerce")
            agg = mecc_reads.groupby("eccDNA_id").agg(
                mean_mapq=("mapq", "mean"),
                mean_match_degree=("match_degree", "mean"),
            )
            read_features = agg.to_dict("index")

    # --- Count candidate sites from regions table ---
    site_counts: dict[str, int] = {}
    if not regions_df.empty:
        mecc_regions = regions_df[regions_df["eccDNA_id"].isin(mecc_ids)]
        if not mecc_regions.empty:
            site_counts = mecc_regions.groupby("eccDNA_id").size().to_dict()

    # --- Build feature matrix for Mecc rows ---
    mecc_idx = summary_df.index[mecc_mask]
    n = len(mecc_idx)
    features = np.zeros((n, len(MECC_ONT_FEATURE_NAMES)))

    for j, idx in enumerate(mecc_idx):
        row = summary_df.loc[idx]
        ecc_id = row["eccDNA_id"]
        rf = read_features.get(ecc_id, {})

        length_val = float(row.get("length", 1) or 1)
        features[j, 0] = rf.get("mean_mapq", 0.0)                # mean_mapq
        features[j, 1] = float(site_counts.get(ecc_id, 2))       # n_candidate_sites
        features[j, 2] = float(row.get("confidence", 0) or 0)    # confidence
        features[j, 3] = rf.get("mean_match_degree", 0.0)        # mean_match_degree
        features[j, 4] = float(row.get("copy_number", 1) or 1)   # copy_number
        features[j, 5] = np.log10(max(length_val, 1.0))          # log_length

    # --- Standardize and score ---
    normed = (features - _MECC_ONT_SCALER_MEAN) / np.maximum(_MECC_ONT_SCALER_STD, 1e-6)
    logits = normed @ _MECC_ONT_MODEL_COEFS + _MECC_ONT_MODEL_INTERCEPT
    scores = np.array([_sigmoid(float(l)) for l in logits])

    # --- Demote low-score Mecc → Uecc ---
    demoted_count = 0
    for j, idx in enumerate(mecc_idx):
        if scores[j] < threshold:
            summary_df.at[idx, "type"] = "Uecc"
            demoted_count += 1

    if demoted_count:
        logger.info(
            "MeccDNA ONT LR filter: demoted %d/%d low-score Mecc -> Uecc (threshold=%.2f)",
            demoted_count, n_mecc, threshold,
        )
    else:
        logger.debug("MeccDNA ONT LR filter: all %d Mecc passed (threshold=%.2f)", n_mecc, threshold)

    return summary_df


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


# -------------------- CeccDNA QC filter ---------------------- #


def _filter_cecc_qc(
    unified_df: pd.DataFrame,
    cecc_df: Optional[pd.DataFrame],
    inferred_chimeric_df: Optional[pd.DataFrame],
    *,
    min_confirmed_length: int = 1400,
    max_inferred_segments: int = 5,
    logger=None,
) -> tuple[pd.DataFrame, Optional[pd.DataFrame], Optional[pd.DataFrame]]:
    """Apply quality-control filters to CeccDNA entries.

    Addresses two empirically observed FP modes in ONT CeccBuild:

    1. **Small Confirmed Cecc from NUMT/homology artifacts** — Short reads
       from MeccDNA or background that happen to split-align across
       chromosomes due to nuclear-organellar homology (NUMTs/NUPTs).
       These produce Confirmed CeccDNA with very short total length.
       Filter: remove Confirmed CeccDNA with ``length < min_confirmed_length``.

    2. **Mega-cluster Inferred Cecc** — The SplitReads inference engine
       occasionally aggregates many unrelated cross-chr signals into a
       single chimeric call with an implausibly large number of fragments.
       Filter: remove Inferred CeccDNA with ``segment_count > max_inferred_segments``.

    Args:
        unified_df: Unified eccDNA table (one row per eccDNA).
        cecc_df: CeccDNA segments core CSV (may be None).
        inferred_chimeric_df: Inferred chimeric CSV (may be None).
        min_confirmed_length: Minimum total length for Confirmed CeccDNA.
        max_inferred_segments: Maximum segment count for Inferred CeccDNA.
        logger: Optional logger.

    Returns:
        (unified_df, cecc_df, inferred_chimeric_df) with FP entries removed.
    """
    if logger is None:
        logger = get_logger("ecc_packager")

    removed_ids: set[str] = set()

    # --- Rule 1: Confirmed CeccDNA minimum length ---
    if min_confirmed_length > 0:
        mask_confirmed_cecc = (
            (unified_df.get("eccDNA_type", unified_df.get("type", "")) == "CeccDNA")
            & (unified_df.get("State", unified_df.get("state", "")) == "Confirmed")
        )
        if mask_confirmed_cecc.any():
            length_col = "Length" if "Length" in unified_df.columns else "length"
            if length_col in unified_df.columns:
                too_short = mask_confirmed_cecc & (
                    pd.to_numeric(unified_df[length_col], errors="coerce").fillna(0)
                    < min_confirmed_length
                )
                short_ids = set(unified_df.loc[too_short, "eccDNA_id"].astype(str))
                if short_ids:
                    removed_ids |= short_ids
                    logger.info(
                        "CeccDNA QC: removed %d Confirmed CeccDNA with length < %d: %s",
                        len(short_ids),
                        min_confirmed_length,
                        ", ".join(sorted(short_ids)),
                    )

    # --- Rule 2: Inferred CeccDNA maximum segments ---
    if max_inferred_segments > 0:
        mask_inferred_cecc = (
            (unified_df.get("eccDNA_type", unified_df.get("type", "")) == "CeccDNA")
            & (unified_df.get("State", unified_df.get("state", "")) == "Inferred")
        )
        if mask_inferred_cecc.any():
            # seg_total may be in unified_df or need to be inferred from cecc_df
            seg_col = None
            for candidate in ("Seg_total", "seg_total", "segment_count", "n_segments"):
                if candidate in unified_df.columns:
                    seg_col = candidate
                    break

            if seg_col is not None:
                too_many = mask_inferred_cecc & (
                    pd.to_numeric(unified_df[seg_col], errors="coerce").fillna(0)
                    > max_inferred_segments
                )
                mega_ids = set(unified_df.loc[too_many, "eccDNA_id"].astype(str))
            elif inferred_chimeric_df is not None and "eccDNA_id" in inferred_chimeric_df.columns:
                # Compute seg_total from inferred chimeric table
                seg_counts = inferred_chimeric_df.groupby("eccDNA_id").size()
                mega_ids_series = seg_counts[seg_counts > max_inferred_segments].index
                inferred_ids = set(unified_df.loc[mask_inferred_cecc, "eccDNA_id"].astype(str))
                mega_ids = set(str(i) for i in mega_ids_series) & inferred_ids
            else:
                mega_ids = set()

            if mega_ids:
                removed_ids |= mega_ids
                logger.info(
                    "CeccDNA QC: removed %d Inferred CeccDNA with segments > %d: %s",
                    len(mega_ids),
                    max_inferred_segments,
                    ", ".join(sorted(mega_ids)),
                )

    # --- Apply removal ---
    if not removed_ids:
        return unified_df, cecc_df, inferred_chimeric_df

    n_before = len(unified_df)
    unified_df = unified_df[~unified_df["eccDNA_id"].astype(str).isin(removed_ids)].copy()
    n_after = len(unified_df)
    logger.info(
        "CeccDNA QC total: removed %d entries from unified table (%d -> %d)",
        n_before - n_after,
        n_before,
        n_after,
    )

    if cecc_df is not None and "eccDNA_id" in cecc_df.columns:
        cecc_df = cecc_df[~cecc_df["eccDNA_id"].astype(str).isin(removed_ids)].copy()

    if inferred_chimeric_df is not None and "eccDNA_id" in inferred_chimeric_df.columns:
        inferred_chimeric_df = inferred_chimeric_df[
            ~inferred_chimeric_df["eccDNA_id"].astype(str).isin(removed_ids)
        ].copy()

    return unified_df, cecc_df, inferred_chimeric_df


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
    ensure_dir(out_dir / "UeccDNA", dry=dry, verbose=verbose)
    ensure_dir(out_dir / "MeccDNA", dry=dry, verbose=verbose)
    ensure_dir(out_dir / "CeccDNA", dry=dry, verbose=verbose)

    # Load unified CSV
    if not unified_csv or not unified_csv.exists():
        logger.error("Unified CSV not found: %s", unified_csv)
        return 1

    unified_df = pd.read_csv(
        unified_csv,
        dtype={"low_mapq": "boolean", "low_identity": "boolean"},
    )
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

    # ── CeccDNA quality-control filter ──
    unified_df, cecc_df, inferred_chimeric_df = _filter_cecc_qc(
        unified_df,
        cecc_df,
        inferred_chimeric_df,
        min_confirmed_length=1400,
        max_inferred_segments=5,
        logger=logger,
    )

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

    # ── MeccDNA ONT LR quality-control filter ──
    # Score Mecc entries using read-level features (mapq, match_degree) and
    # region-level features (candidate site count).  Low-scoring Mecc are
    # demoted to Uecc.  The model was trained on ONT data, so only apply when
    # the reads table contains the ONT-specific `eccDNA_copy_number` column
    # (present in ONT pipeline output but not NGS).  NGS already runs its
    # own Mecc LR filter in ngs_circle_detect.
    _is_ont_reads = (
        not reads_df.empty
        and "eccDNA_copy_number" in reads_df.columns
        and "mapq" in reads_df.columns
        and "match_degree" in reads_df.columns
    )
    if _is_ont_reads and (summary_df["type"] == "Mecc").any():
        summary_df = score_mecc_ont_lr(
            summary_df, reads_df, regions_df,
            threshold=_MECC_ONT_DEFAULT_THRESHOLD,
        )

    # ── Coordinate-level NMS dedup for Cecc and Mecc ──
    # At higher coverage, the same eccDNA is detected multiple times with
    # slightly shifted coordinates.  cd-hit catches identical sequences but
    # misses coordinate-shifted duplicates.  Apply NMS: sort by read_count
    # descending, remove entries whose fragments overlap a kept entry.
    #
    # Overlap logic differs by type:
    #   Mecc: ANY fragment overlap > 50% → duplicate (single-region eccDNA)
    #   Cecc: ALL fragments must overlap → duplicate (different CeccDNA can
    #         legitimately share one fragment but differ on others)
    def _parse_frags(row):
        frags = []
        loc = str(row.get("location", ""))
        for m in re.finditer(r"(\w+):(\d+)-(\d+)", loc):
            frags.append((m.group(1), int(m.group(2)), int(m.group(3))))
        if not frags:
            c = str(row.get("chr", ""))
            if c and c != "multi":
                try:
                    frags = [(c, int(row["start"]), int(row["end"]))]
                except (ValueError, TypeError):
                    pass
        return frags

    def _frag_overlap(f1, f2):
        """Overlap fraction between two fragments (0 if diff chr)."""
        c1, s1, e1 = f1
        c2, s2, e2 = f2
        if c1 != c2:
            return 0.0
        ovl = max(0, min(e1, e2) - max(s1, s2))
        min_len = min(e1 - s1, e2 - s2)
        return ovl / min_len if min_len > 0 else 0.0

    for etype in ("Cecc", "Mecc"):
        mask = summary_df["type"] == etype
        if mask.sum() <= 1:
            continue
        sub = summary_df[mask].copy()
        det_regions = [(idx, _parse_frags(row), int(row.get("read_count", 1)))
                       for idx, row in sub.iterrows()]
        det_regions.sort(key=lambda x: -x[2])

        keep_idx = []
        for idx, frags, rc in det_regions:
            if not frags:
                keep_idx.append((idx, frags, rc))
                continue
            is_dup = False
            for k_idx, k_frags, _ in keep_idx:
                if not k_frags:
                    continue
                if etype == "Cecc":
                    # ALL fragments of the new entry must overlap a kept fragment
                    n_matched = 0
                    for ci, si, ei in frags:
                        for cj, sj, ej in k_frags:
                            if _frag_overlap((ci, si, ei), (cj, sj, ej)) > 0.5:
                                n_matched += 1
                                break
                    is_dup = (n_matched == len(frags))
                else:
                    # Mecc: ANY fragment overlap is enough
                    for ci, si, ei in frags:
                        for cj, sj, ej in k_frags:
                            if _frag_overlap((ci, si, ei), (cj, sj, ej)) > 0.5:
                                is_dup = True
                                break
                        if is_dup:
                            break
                if is_dup:
                    break
            if not is_dup:
                keep_idx.append((idx, frags, rc))

        n_removed = mask.sum() - len(keep_idx)
        if n_removed > 0:
            remove_idx = set(sub.index) - {k[0] for k in keep_idx}
            summary_df = summary_df.drop(index=remove_idx).reset_index(drop=True)
            logger.info(
                f"{etype} NMS dedup: removed {n_removed} coordinate-overlapping "
                f"detections, kept {len(keep_idx)}"
            )

    # Save CSV files
    summary_df.to_csv(out_dir / f"{sample}_eccDNA_summary.csv", index=False)
    regions_df.to_csv(out_dir / f"{sample}_eccDNA_regions.csv", index=False)
    reads_df.to_csv(out_dir / f"{sample}_eccDNA_reads.csv", index=False)

    logger.info("Saved %s_eccDNA_summary.csv (%d rows)", sample, len(summary_df))
    logger.info("Saved %s_eccDNA_regions.csv (%d rows)", sample, len(regions_df))
    logger.info("Saved %s_eccDNA_reads.csv (%d rows)", sample, len(reads_df))

    # Generate BED files
    generate_uecc_bed(regions_df, out_dir / "UeccDNA" / f"{sample}_uecc.bed")
    generate_mecc_bed(regions_df, out_dir / "MeccDNA" / f"{sample}_mecc_sites.bed")
    generate_cecc_bed(regions_df, out_dir / "CeccDNA" / f"{sample}_cecc_segments.bed")
    generate_cecc_bedpe(regions_df, out_dir / "CeccDNA" / f"{sample}_cecc_junctions.bedpe")

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
    generate_fasta_files(sequences, out_dir, summary_df, prefix=sample)

    # Copy report files
    if html_in:
        copy_file(html_in, out_dir / f"{sample}_report.html", overwrite=args.overwrite, dry=dry, verbose=verbose)
    if txt_in:
        copy_file(txt_in, out_dir / f"{sample}_summary.txt", overwrite=args.overwrite, dry=dry, verbose=verbose)

    logger.info("[done] Output folder: %s", out_dir)
    return 0


# --------------------------- Entrypoint ------------------------ #
if __name__ == "__main__":
    args = build_argparser().parse_args()
    raise SystemExit(run(args))
