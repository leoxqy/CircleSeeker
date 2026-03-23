"""
NGS Circle Detection — Detect eccDNA from short-read breakpoint evidence.

Uses back-jump split reads as the primary eccDNA signal: when a read spans
an eccDNA circle junction, the reference position "jumps backward" along
the read (from the circle's end back to its start).  This physical signal
distinguishes eccDNA junctions from deletions (forward-jump) and other SVs.

A logistic regression model combines multiple evidence features:
  - n_backjump:      back-jump split reads (strongest positive signal)
  - n_fwdjump:       forward-jump split reads (SV noise, negative)
  - n_disc_outward:  outward-facing discordant pairs (negative in practice)
  - n_disc_inward:   inward-facing discordant pairs (negative)
  - n_clip_left/right: soft-clipped reads at boundaries
  - log_size:        log10 of candidate region size (large = more likely FP)

Model coefficients were derived from cross-genome validation
(trained on A.thaliana, validated on human; F1=0.915 on held-out genome).
"""

from __future__ import annotations

import logging
import math
import re
import subprocess
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import pysam

from circleseeker.utils.logging import get_logger


# ── Fixed model coefficients (cross-genome validated) ──────────────────
# Trained on ara_UMC_5200 rep1 30X, tested on human_U_10000 rep1 30X.
# Features are StandardScaler-normalized before applying.

_MODEL_COEFS = np.array([
    +1.630,   # n_backjump
    -0.451,   # n_fwdjump
    -1.216,   # n_disc_outward
    -0.699,   # n_disc_inward
    +0.360,   # n_clip_left
    +0.124,   # n_clip_right
    -0.539,   # log_size
])
_MODEL_INTERCEPT = 0.0  # balanced class_weight shifts this

# Normalization parameters (mean, std) from training set
_SCALER_MEAN = np.array([13.06, 0.003, 0.002, 5.96, 2.49, 2.55, 3.03])
_SCALER_STD = np.array([12.50, 0.08, 0.06, 6.80, 3.20, 3.30, 0.55])

FEATURE_NAMES = [
    "n_backjump", "n_fwdjump", "n_disc_outward", "n_disc_inward",
    "n_clip_left", "n_clip_right", "log_size",
]

# ── CeccDNA LR model coefficients (trained on ara rep1, validated rep2/3) ──
_CECC_MODEL_COEFS = np.array([
    +0.0693,   # n_junctions
    +0.0260,   # mean_reads_per_junction
    -0.0221,   # min_reads_per_junction
    +0.0319,   # split_disc_ratio
    +0.1033,   # log_total_reads
    -0.1057,   # n_fragments
    -0.1057,   # n_chromosomes
    +0.2093,   # log_total_length
    +0.1299,   # max_min_frag_ratio
    -0.0999,   # strand_closure
    +0.0833,   # n_strand_flips
    +0.1721,   # frag_overlap_uecc_ratio
    +0.0376,   # all_frags_overlap
    +0.0260,   # mean_bj_at_frags
    +0.1258,   # mean_disc_at_frags
    -0.0188,   # mean_junction_mapq
    +0.0289,   # min_junction_mapq
])
_CECC_MODEL_INTERCEPT = 0.4381
_CECC_SCALER_MEAN = np.array([
    4.5000, 18.1130, 6.9342, 0.4095, 1.8198,
    2.1053, 2.1053, 3.4735, 3.2979, 0.6842,
    4.0526, 0.4846, 0.0921, 8.9572, 9.8958,
    50.3421, 11.2500,
])
_CECC_SCALER_STD = np.array([
    1.7206, 9.7342, 5.6530, 0.0799, 0.2961,
    0.3832, 0.3832, 0.2412, 2.1262, 0.4648,
    1.6535, 0.2313, 0.2892, 6.2337, 7.6093,
    3.0677, 2.4822,
])

CECC_FEATURE_NAMES = [
    "n_junctions", "mean_reads_per_junction", "min_reads_per_junction",
    "split_disc_ratio", "log_total_reads",
    "n_fragments", "n_chromosomes", "log_total_length", "max_min_frag_ratio",
    "strand_closure", "n_strand_flips",
    "frag_overlap_uecc_ratio", "all_frags_overlap",
    "mean_bj_at_frags", "mean_disc_at_frags",
    "mean_junction_mapq", "min_junction_mapq",
]

# ── MeccDNA LR model coefficients (trained on ara rep1-3) ──────────
MECC_FEATURE_NAMES = [
    "n_backjump", "n_backjump_low", "bj_ratio", "n_reads", "log_length",
    "rpk", "n_discordant", "disc_ratio", "n_softclip", "score",
    "n_genomic_loci",
]
_MECC_MODEL_COEFS = np.array([
    -0.1216,   # n_backjump
    +0.0368,   # n_backjump_low
    -0.1664,   # bj_ratio
    -0.0619,   # n_reads
    -0.0785,   # log_length
    +0.0162,   # rpk
    -0.0619,   # n_discordant
    -0.0765,   # disc_ratio
    -0.0177,   # n_softclip
    -0.0597,   # score
    -0.1005,   # n_genomic_loci
])
_MECC_MODEL_INTERCEPT = 0.2021
_MECC_SCALER_MEAN = np.array([
    1.1944, 5.4023, 0.1158, 6.9921, 2.7810,
    11.5876, 0.5695, 0.0358, 2.4994, 0.5361,
    4.2588,
])
_MECC_SCALER_STD = np.array([
    3.5275, 4.3459, 0.2427, 6.3784, 0.2668,
    9.9181, 3.0086, 0.1462, 6.9089, 0.1015,
    4.4005,
])


def _sigmoid(x: float) -> float:
    return 1.0 / (1.0 + math.exp(-x)) if x > -500 else 0.0


# ── Configuration ──────────────────────────────────────────────────────

@dataclass
class NGSDetectConfig:
    min_mapq: int = 10
    min_clip_len: int = 20
    min_split_len: int = 30
    max_normal_insert: int = 1000
    bin_size: int = 1000          # breakpoint clustering bin size
    min_evidence: int = 2         # minimum total evidence per candidate
    min_circle_size: int = 100
    max_circle_size: int = 100000
    score_threshold: float = 0.3  # default threshold (0.3 → F1≈0.915)
    cecc_score_threshold: float = 0.43  # CeccDNA LR threshold (val Prec=86%)
    mecc_score_threshold: float = 0.32  # MeccDNA LR threshold
    threads: int = 4
    minimap2_preset: str = "sr"


# ── Main class ─────────────────────────────────────────────────────────

class NGSCircleDetect:
    """End-to-end NGS eccDNA detection: FASTQ → BAM → circles."""

    def __init__(
        self,
        reference: Path,
        config: Optional[NGSDetectConfig] = None,
        logger: Optional[logging.Logger] = None,
    ):
        self.reference = Path(reference)
        self.config = config or NGSDetectConfig()
        self.logger = logger or get_logger(self.__class__.__name__)

    # ── Public API ─────────────────────────────────────────────────────

    def run_from_fastq(
        self,
        r1_fastq: Path,
        r2_fastq: Path,
        output_dir: Path,
        prefix: str = "ngs",
    ) -> pd.DataFrame:
        """Full pipeline: align + detect + score + classify."""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        bam = output_dir / f"{prefix}_sorted.bam"
        if not bam.exists():
            self._align(r1_fastq, r2_fastq, bam)

        return self.run_from_bam(bam, output_dir, prefix)

    def run_from_bam(
        self,
        bam_path: Path,
        output_dir: Path,
        prefix: str = "ngs",
    ) -> pd.DataFrame:
        """Detect eccDNA from an existing sorted BAM.

        Pipeline order: CeccDNA-first to prevent intra-chr signal theft.
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Step 1: Extract evidence (intra-chr + cross-chr)
        regions, crosschr, crosschr_readnames = self._extract_evidence(bam_path)
        self.logger.info(
            f"NGS detect: {len(regions)} intra-chr regions, "
            f"{len(crosschr)} cross-chr junctions, "
            f"{len(crosschr_readnames)} cross-chr reads"
        )

        # Step 2: CeccDNA detection FIRST (before Uecc/Mecc)
        cecc_df, cecc_fragments = self._detect_cecc_early(crosschr)
        self.logger.info(f"NGS detect: {len(cecc_df)} CeccDNA circles (raw)")

        # Step 3: Build exclusion zones from ALL CeccDNA fragments (liberal)
        exclusion_zones = self._build_exclusion_zones(cecc_fragments)

        # Step 4-6: Uecc/Mecc pipeline (with CeccDNA exclusion)
        um_df = pd.DataFrame()
        if regions:
            scored, cecc_rescued = self._score_regions(
                regions, exclusion_zones=exclusion_zones
            )
            if cecc_rescued:
                # Add rescued DISC-dominant regions as Cecc fragments
                rescued_df = pd.DataFrame(cecc_rescued)
                cecc_idx = len(cecc_df) + 1
                rescued_df["eccDNA_id"] = [
                    f"NGS_Cecc_frag{i:05d}"
                    for i in range(cecc_idx, cecc_idx + len(rescued_df))
                ]
                cecc_df = pd.concat([cecc_df, rescued_df], ignore_index=True)
                self.logger.info(
                    f"NGS detect: rescued {len(cecc_rescued)} DISC-dominant "
                    f"regions as Cecc fragments"
                )
            # Apply CeccDNA LR scoring
            if not cecc_df.empty:
                cecc_df = self._compute_cecc_features(cecc_df, regions)
                cecc_df = self._score_cecc_lr(cecc_df)
            kept = self._dedup(scored)
            self.logger.info(f"NGS detect: {len(kept)} after dedup")
            results = self._filter_and_classify(kept)
            self.logger.info(f"NGS detect: {len(results)} Uecc/Mecc circles")
            if results:
                um_df = pd.DataFrame(results)
                um_df = self.reclassify_mecc(um_df)
                um_df = self._score_mecc_lr(um_df)
                um_df = self._dedup_mecc_by_sequence(um_df)

        # Step 7: Reclassify Uecc in CeccDNA exclusion zones
        # Disabled: low precision reclassification hurts CeccDNA Prec target (>=90%)

        # Step 8: Filter out Cecc_frag ghost detections
        # Rescued fragments with no junctions and no fragment structure
        # are DISC-dominant single regions, not valid standalone CeccDNA.
        if not cecc_df.empty:
            before_n = len(cecc_df)
            has_junctions = cecc_df.get("n_junctions", pd.Series(dtype=float)).fillna(0) > 0
            has_fragments = cecc_df.get("fragments", pd.Series(dtype=object)).fillna("").astype(str).str.len() > 0
            cecc_keep = has_junctions | has_fragments
            n_ghost = int((~cecc_keep).sum())
            if n_ghost:
                cecc_df = cecc_df[cecc_keep].reset_index(drop=True)
                self.logger.info(
                    f"NGS detect: filtered {n_ghost} Cecc_frag ghost detections "
                    f"(no junctions, no fragments), kept {len(cecc_df)}/{before_n}"
                )

        # Step 9: Filter Mecc ghost large circles
        # Mecc detections with 0 high-MAPQ backjump and very low reads density
        # are noise from discordant reads in repetitive regions.
        # Use n_backjump (high-MAPQ only) because all Mecc have n_backjump_low > 0
        # by definition (that's how they become Mecc in _filter_and_classify).
        if not um_df.empty:
            mecc_mask = um_df["eccdna_type"] == "Mecc"
            n_bj_hi = um_df["n_backjump"].fillna(0)
            length_kb = um_df["length"].fillna(1) / 1000.0
            reads_per_kb = um_df["n_reads"].fillna(0) / length_kb.clip(lower=0.001)
            ghost_mecc = mecc_mask & (n_bj_hi == 0) & (reads_per_kb < 3.0)
            n_ghost_mecc = int(ghost_mecc.sum())
            if n_ghost_mecc:
                um_df = um_df[~ghost_mecc].reset_index(drop=True)
                self.logger.info(
                    f"NGS detect: filtered {n_ghost_mecc} Mecc ghost detections "
                    f"(0 high-MAPQ backjumps, reads/kb < 3.0)"
                )

        # Merge results
        dfs = [d for d in [um_df, cecc_df] if not d.empty]
        df = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()

        if not df.empty:
            out_csv = output_dir / f"{prefix}_ngs_circles.csv"
            df.to_csv(out_csv, index=False)

        return df

    # ── Alignment ──────────────────────────────────────────────────────

    def _align(self, r1: Path, r2: Path, output_bam: Path) -> None:
        self.logger.info("NGS detect: aligning with minimap2 -sr")
        cmd = (
            f"minimap2 -ax {self.config.minimap2_preset} "
            f"-t {self.config.threads} "
            f"{self.reference} {r1} {r2} 2>/dev/null | "
            f"samtools sort -@ {self.config.threads} -o {output_bam}"
        )
        subprocess.run(cmd, shell=True, timeout=7200)
        subprocess.run(f"samtools index {output_bam}", shell=True, timeout=300)

    # ── Evidence extraction ────────────────────────────────────────────

    def _extract_evidence(self, bam_path: Path) -> tuple[dict, list]:
        """Extract per-region evidence from BAM.

        Returns:
            (regions, crosschr_junctions, crosschr_readnames)
            - regions: dict of intra-chromosomal evidence (for Uecc/Mecc)
            - crosschr_junctions: list of cross-chr junctions, each as
              (chr_a, pos_a, strand_a, chr_b, pos_b, strand_b, read_name, evidence_type)
            - crosschr_readnames: set of read names involved in cross-chr junctions
        """
        cfg = self.config
        bam = pysam.AlignmentFile(str(bam_path), "rb")

        regions = defaultdict(lambda: {
            "bj": set(), "bj_low": set(), "fj": set(),
            "do": set(), "di": set(),
            "cl": set(), "cr": set(),
            "positions": [],
        })
        crosschr_junctions = []

        for read in bam.fetch():
            if read.is_unmapped or read.is_secondary or read.is_duplicate:
                continue

            chrom = read.reference_name
            rstart, rend = read.reference_start, read.reference_end

            # ── Split reads (all MAPQ, stratified internally) ──
            if read.has_tag("SA"):
                self._process_split(
                    read, chrom, rstart, rend, regions, crosschr_junctions
                )

            # ── Discordant pairs (require min MAPQ) ──
            if (
                read.mapping_quality >= cfg.min_mapq
                and read.is_paired
                and not read.mate_is_unmapped
                and read.is_read1
            ):
                self._process_discordant(
                    read, chrom, rstart, rend, regions, crosschr_junctions
                )

            # ── Soft-clips ──
            if read.mapping_quality >= cfg.min_mapq and read.cigartuples:
                self._process_softclip(read, chrom, rstart, rend, regions)

        bam.close()
        crosschr_readnames = {item[6] for item in crosschr_junctions}
        return dict(regions), crosschr_junctions, crosschr_readnames

    def _process_split(self, read, chrom, rstart, rend, regions, crosschr):
        cfg = self.config
        for entry in read.get_tag("SA").strip(";").split(";"):
            parts = entry.split(",")
            if len(parts) < 6:
                continue
            sa_chr, sa_pos = parts[0], int(parts[1]) - 1
            sa_mapq = int(parts[4])

            sa_rev = parts[2] == "-"
            cigar_ops = re.findall(r"(\d+)([MIDNSHP=X])", parts[3])
            sa_ref_len = sum(int(n) for n, op in cigar_ops if op in "MDN=X")

            if sa_chr != chrom:
                # Cross-chromosome split read → CeccDNA evidence
                # Require min MAPQ for both primary and SA
                effective_mapq = min(read.mapping_quality, sa_mapq)
                if effective_mapq < cfg.min_mapq:
                    continue
                pri_bp = rend if not read.is_reverse else rstart
                sa_bp = sa_pos if not sa_rev else sa_pos + sa_ref_len
                pri_strand = "-" if read.is_reverse else "+"
                sa_strand = "-" if sa_rev else "+"
                crosschr.append((chrom, pri_bp, pri_strand, sa_chr, sa_bp, sa_strand, read.query_name, "split", effective_mapq))
                continue

            sa_leading_clip = int(cigar_ops[0][0]) if cigar_ops[0][1] in "SH" else 0

            pri_qstart = read.query_alignment_start
            pri_rev = read.is_reverse

            if pri_qstart < sa_leading_clip:
                first_r = rstart if not pri_rev else rend
                second_r = sa_pos if not sa_rev else sa_pos + sa_ref_len
            else:
                first_r = sa_pos if not sa_rev else sa_pos + sa_ref_len
                second_r = rstart if not pri_rev else rend

            left = min(rstart, sa_pos)
            right = max(rend, sa_pos + sa_ref_len)
            size = right - left
            if size < cfg.min_circle_size or size > cfg.max_circle_size:
                continue

            key = (chrom, left // cfg.bin_size, right // cfg.bin_size)
            if second_r < first_r:
                # Back-jump: stratify by min(primary, SA) MAPQ
                effective_mapq = min(read.mapping_quality, sa_mapq)
                if effective_mapq >= cfg.min_mapq:
                    regions[key]["bj"].add(read.query_name)
                else:
                    regions[key]["bj_low"].add(read.query_name)
            else:
                regions[key]["fj"].add(read.query_name)
            regions[key]["positions"].append((left, right))

    def _process_discordant(self, read, chrom, rstart, rend, regions, crosschr):
        cfg = self.config
        if chrom != read.next_reference_name:
            # Cross-chromosome discordant pair → CeccDNA evidence
            mate_chr = read.next_reference_name
            pri_pos = rend if not read.is_reverse else rstart
            mate_pos = read.next_reference_start
            pri_strand = "-" if read.is_reverse else "+"
            mate_strand = "-" if read.mate_is_reverse else "+"
            crosschr.append((chrom, pri_pos, pri_strand, mate_chr, mate_pos, mate_strand, read.query_name, "disc", read.mapping_quality))
            return
        mate_pos = read.next_reference_start
        left = min(rstart, mate_pos)
        right = max(rend, mate_pos + 150)
        size = right - left
        if size < cfg.min_circle_size or size > cfg.max_circle_size:
            return

        key = (chrom, left // cfg.bin_size, right // cfg.bin_size)
        if read.is_reverse == read.mate_is_reverse:
            regions[key]["do"].add(read.query_name)
        elif abs(read.template_length) > cfg.max_normal_insert:
            regions[key]["di"].add(read.query_name)

    def _process_softclip(self, read, chrom, rstart, rend, regions):
        cfg = self.config
        if read.cigartuples[0][0] == 4 and read.cigartuples[0][1] >= cfg.min_clip_len:
            key = (chrom, rstart // cfg.bin_size, rstart // cfg.bin_size)
            regions[key]["cl"].add(read.query_name)
        if read.cigartuples[-1][0] == 4 and read.cigartuples[-1][1] >= cfg.min_clip_len:
            key = (chrom, rend // cfg.bin_size, rend // cfg.bin_size)
            regions[key]["cr"].add(read.query_name)

    # ── Exclusion zone helpers ────────────────────────────────────────

    @staticmethod
    def _region_in_exclusion(
        chrom: str, start: int, end: int, exclusion_zones: dict,
    ) -> bool:
        """Check if a region overlaps any CeccDNA exclusion zone.

        Uses max of both overlap fractions so that a small exclusion zone
        fully within a large region still triggers (and vice versa).
        """
        for zone_s, zone_e in exclusion_zones.get(chrom, []):
            if zone_s >= end:
                break
            if zone_e <= start:
                continue
            ov = min(zone_e, end) - max(zone_s, start)
            ov_region = ov / max(end - start, 1)
            ov_zone = ov / max(zone_e - zone_s, 1)
            if max(ov_region, ov_zone) > 0.3:
                return True
        return False

    # ── Scoring ────────────────────────────────────────────────────────

    def _score_regions(
        self, regions: dict, exclusion_zones: Optional[dict] = None,
    ) -> tuple[list[dict], list[dict]]:
        scored = []
        cecc_rescued = []  # DISC-dominant regions in exclusion zones → Cecc fragments
        # LR coefficients without log_size penalty
        coefs_no_size = _MODEL_COEFS.copy()
        coefs_no_size[6] = 0.0

        for key, ev in regions.items():
            total = (len(ev["bj"]) + len(ev.get("bj_low", set()))
                     + len(ev["fj"]) + len(ev["do"]) + len(ev["di"]))
            if total < self.config.min_evidence:
                continue

            chrom = key[0]
            if ev["positions"]:
                starts = [p[0] for p in ev["positions"]]
                ends = [p[1] for p in ev["positions"]]
                med_s, med_e = int(np.median(starts)), int(np.median(ends))
            else:
                med_s = key[1] * self.config.bin_size
                med_e = key[2] * self.config.bin_size

            size = med_e - med_s
            if size < self.config.min_circle_size:
                continue

            # For scoring, combine high+low MAPQ back-jumps as total signal
            total_bj = len(ev["bj"]) + len(ev.get("bj_low", set()))
            n_fj = len(ev["fj"])

            # Check exclusion zones: rescue DISC-dominant regions as Cecc fragments
            if exclusion_zones and self._region_in_exclusion(
                chrom, med_s, med_e, exclusion_zones
            ):
                n_disc = len(ev["do"]) + len(ev["di"])
                if total_bj <= n_disc and n_disc >= 12:
                    # DISC-dominant with solid evidence → CeccDNA fragment
                    n_total_resc = len(ev["bj"] | ev.get("bj_low", set()) | ev["fj"] | ev["do"] | ev["di"])
                    cecc_rescued.append({
                        "chr": chrom, "start": med_s, "end": med_e,
                        "length": size,
                        "eccdna_type": "Cecc",
                        "state": "Confirmed",
                        "score": round(min(1.0, n_disc / 20.0), 4),
                        "n_backjump": len(ev["bj"]),
                        "n_backjump_low": len(ev.get("bj_low", set())),
                        "n_discordant": n_disc,
                        "n_softclip": len(ev["cl"]) + len(ev["cr"]),
                        "n_reads": n_total_resc,
                        # LR features (defaults for rescued fragments)
                        "n_junctions": 0,
                        "n_fragments": 1,
                        "n_chromosomes": 1,
                        "mean_reads_per_junction": 0,
                        "min_reads_per_junction": 0,
                        "split_disc_ratio": 0,
                        "log_total_reads": round(float(np.log10(max(n_total_resc, 1))), 4),
                        "max_min_frag_ratio": 1.0,
                        "log_total_length": round(float(np.log10(max(size, 1))), 4),
                        "strand_closure": 0,
                        "n_strand_flips": 0,
                        "mean_junction_mapq": 0,
                        "min_junction_mapq": 0,
                    })
                    continue

            features = np.array([
                total_bj, n_fj,
                len(ev["do"]), len(ev["di"]),
                len(ev["cl"]), len(ev["cr"]),
                np.log10(max(size, 1)),
            ])

            # Hybrid scoring: LR (without size penalty) + backjump dominance
            normed = (features - _SCALER_MEAN) / np.maximum(_SCALER_STD, 1e-6)
            logit = float(np.dot(normed, coefs_no_size) + _MODEL_INTERCEPT)
            score = _sigmoid(logit)

            # Backjump dominance rescue: if backjumps clearly dominate
            # forward-jumps, the candidate is likely real regardless of LR score
            if total_bj >= 2 and n_fj == 0:
                score = max(score, 0.5)
            elif total_bj >= 1 and total_bj > n_fj * 2:
                score = max(score, 0.35)

            scored.append({
                "chr": chrom, "start": med_s, "end": med_e, "size": size,
                "score": score,
                "n_backjump": len(ev["bj"]),
                "n_backjump_low": len(ev.get("bj_low", set())),
                "n_fwdjump": n_fj,
                "n_disc_outward": len(ev["do"]), "n_disc_inward": len(ev["di"]),
                "n_clip_left": len(ev["cl"]), "n_clip_right": len(ev["cr"]),
                "n_reads": len(ev["bj"] | ev.get("bj_low", set()) | ev["fj"] | ev["do"] | ev["di"]),
            })

        return scored, cecc_rescued

    # ── Dedup ──────────────────────────────────────────────────────────

    @staticmethod
    def _dedup(candidates: list[dict]) -> list[dict]:
        def ofrac(s1, e1, s2, e2):
            o = max(0, min(e1, e2) - max(s1, s2))
            return min(o / max(e1 - s1, 1), o / max(e2 - s2, 1))

        candidates.sort(key=lambda x: -x["score"])
        kept = []
        for cand in candidates:
            if not any(
                k["chr"] == cand["chr"]
                and ofrac(cand["start"], cand["end"], k["start"], k["end"]) > 0.5
                for k in kept
            ):
                kept.append(cand)
        return kept

    # ── Classify ───────────────────────────────────────────────────────

    def _filter_and_classify(self, candidates: list[dict]) -> list[dict]:
        """Classify circles as UeccDNA or MeccDNA based on MAPQ stratification.

        - High-MAPQ back-jumps (bj) → UeccDNA (unique mapping)
        - Low-MAPQ back-jumps (bj_low) → MeccDNA (multi-mapping = multi-locus)
        - Regions with BOTH → classified by dominant signal
        """
        results = []
        uecc_idx = 1
        mecc_idx = 1

        for cand in candidates:
            if cand["score"] < self.config.score_threshold:
                continue

            n_bj_hi = cand.get("n_backjump", 0)
            n_bj_lo = cand.get("n_backjump_low", 0)

            # Classify by dominant MAPQ signal
            if n_bj_hi > 0 and n_bj_hi >= n_bj_lo:
                etype = "Uecc"
                eid = f"NGS_Uecc{uecc_idx:05d}"
                uecc_idx += 1
            elif n_bj_lo > 0:
                etype = "Mecc"
                eid = f"NGS_Mecc{mecc_idx:05d}"
                mecc_idx += 1
            else:
                continue

            results.append({
                "eccDNA_id": eid,
                "chr": cand["chr"],
                "start": cand["start"],
                "end": cand["end"],
                "length": cand["size"],
                "eccdna_type": etype,
                "state": "Confirmed",
                "score": round(cand["score"], 4),
                "n_backjump": n_bj_hi,
                "n_backjump_low": n_bj_lo,
                "n_discordant": cand.get("n_disc_outward", 0) + cand.get("n_disc_inward", 0),
                "n_softclip": cand.get("n_clip_left", 0) + cand.get("n_clip_right", 0),
                "n_reads": cand.get("n_reads", 0),
            })

        return results

    # ── MeccDNA reclassification ──────────────────────────────────────

    def reclassify_mecc(
        self,
        circles_df: pd.DataFrame,
        min_loci: int = 2,
        min_query_cov: float = 0.7,
    ) -> pd.DataFrame:
        """Reclassify UeccDNA as MeccDNA by re-aligning circle sequences.

        For each detected circle, extract its sequence from the reference
        and align back to the full genome.  If the sequence has full-length
        hits at ≥2 distinct genomic locations, the circle is from a
        repetitive/multi-copy region and is reclassified as MeccDNA.

        Args:
            circles_df: DataFrame from run_from_bam (UeccDNA detections)
            min_loci: Minimum distinct loci for MeccDNA (default 2)
            min_query_cov: Minimum query coverage for a "full-length" hit

        Returns:
            Updated DataFrame with eccdna_type reclassified where applicable.
        """
        if circles_df.empty:
            return circles_df

        import subprocess
        import tempfile

        ref_fasta = pysam.FastaFile(str(self.reference))
        tmpdir = Path(tempfile.mkdtemp())

        # Write circle sequences
        circle_fa = tmpdir / "circles.fasta"
        with open(circle_fa, "w") as f:
            for _, row in circles_df.iterrows():
                try:
                    seq = ref_fasta.fetch(
                        str(row["chr"]), int(row["start"]), int(row["end"])
                    )
                    f.write(f">{row['eccDNA_id']}\n{seq}\n")
                except (ValueError, KeyError):
                    continue
        ref_fasta.close()

        # Align circle sequences to reference (allow many secondaries)
        paf_out = tmpdir / "circles_vs_ref.paf"
        subprocess.run(
            f"minimap2 -cx asm5 --secondary=yes -N 20 -p 0.1 "
            f"-t {self.config.threads} "
            f"{self.reference} {circle_fa} > {paf_out} 2>/dev/null",
            shell=True,
            timeout=300,
        )

        # Parse: count distinct full-length loci per circle
        circle_loci: dict[str, list[tuple]] = defaultdict(list)
        with open(paf_out) as f:
            for line in f:
                cols = line.rstrip().split("\t")
                eid = cols[0]
                qlen = int(cols[1])
                qs, qe = int(cols[2]), int(cols[3])
                tname = cols[5]
                ts, te = int(cols[7]), int(cols[8])

                if (qe - qs) / qlen < min_query_cov:
                    continue
                circle_loci[eid].append((tname, ts, te))

        # Deduplicate loci and reclassify
        def _ofrac(s1, e1, s2, e2):
            ov = max(0, min(e1, e2) - max(s1, s2))
            return min(ov / max(e1 - s1, 1), ov / max(e2 - s2, 1))

        mecc_ids = set()
        mecc_loci_count = {}
        for eid, loci in circle_loci.items():
            unique = []
            for loc in sorted(loci, key=lambda x: (-abs(x[2]-x[1]))):
                if not any(
                    loc[0] == u[0] and _ofrac(loc[1], loc[2], u[1], u[2]) > 0.5
                    for u in unique
                ):
                    unique.append(loc)
            if len(unique) >= min_loci:
                mecc_ids.add(eid)
                mecc_loci_count[eid] = len(unique)

        # Update DataFrame — bidirectional reclassification
        circles_df = circles_df.copy()
        all_aligned = set(circle_loci.keys())

        # ≥2 loci → Mecc, BUT require backjump evidence for borderline cases
        # (n_loci == 2 without any backjump is likely sequence similarity noise)
        confident_mecc_ids = set()
        for eid in mecc_ids:
            row = circles_df.loc[circles_df["eccDNA_id"] == eid]
            if row.empty:
                continue
            r = row.iloc[0]
            n_loci = mecc_loci_count.get(eid, 0)
            n_bj_hi = int(r.get("n_backjump", 0))
            n_bj_lo = int(r.get("n_backjump_low", 0))
            # >=3 loci: confident Mecc (strong multi-copy evidence)
            # ==2 loci: only Mecc if low-MAPQ backjump dominates.
            #   If high-MAPQ backjumps dominate (n_bj_hi > n_bj_lo),
            #   the circle maps uniquely → likely Uecc in paralogous region.
            if n_loci >= 3:
                confident_mecc_ids.add(eid)
            elif n_loci == 2 and n_bj_lo > n_bj_hi:
                confident_mecc_ids.add(eid)

        mask_to_mecc = circles_df["eccDNA_id"].isin(confident_mecc_ids)
        circles_df.loc[mask_to_mecc, "eccdna_type"] = "Mecc"
        circles_df.loc[mask_to_mecc, "n_genomic_loci"] = (
            circles_df.loc[mask_to_mecc, "eccDNA_id"].map(mecc_loci_count)
        )

        # 1 locus or unconfident 2-locus → Uecc
        demote_ids = (all_aligned - confident_mecc_ids) | (mecc_ids - confident_mecc_ids)
        mask_to_uecc = (
            (circles_df["eccDNA_id"].isin(demote_ids))
            & (circles_df["eccdna_type"] == "Mecc")
        )
        circles_df.loc[mask_to_uecc, "eccdna_type"] = "Uecc"
        circles_df.loc[mask_to_uecc, "n_genomic_loci"] = 1

        n_demoted_2loci = len(mecc_ids - confident_mecc_ids)
        if n_demoted_2loci:
            self.logger.info(
                f"MeccDNA reclassification: demoted {n_demoted_2loci} "
                f"2-loci detections without backjump → Uecc"
            )

        n_to_mecc = mask_to_mecc.sum()
        n_to_uecc = mask_to_uecc.sum()
        self.logger.info(
            f"MeccDNA reclassification: {n_to_mecc} → Mecc (≥{min_loci} loci), "
            f"{n_to_uecc} → Uecc (1 locus)"
        )

        # Cleanup
        import shutil
        shutil.rmtree(tmpdir, ignore_errors=True)

        return circles_df

    # ── CeccDNA post-filter ─────────────────────────────────────────────

    def _filter_cecc_by_intrachr(
        self, cecc_df: pd.DataFrame, regions: dict,
        min_bj_for_overlap: int = 3,
    ) -> pd.DataFrame:
        """Remove CeccDNA whose ALL fragments overlap strong intra-chr BJ regions.

        Multi-copy (Mecc) regions create false cross-chr junctions from
        multi-mapping reads. These regions also have strong intra-chr back-jump
        signals. Real CeccDNA fragments should NOT have independent back-jump
        evidence at each fragment location.
        """
        # Build lookup of regions with strong back-jump evidence
        bj_regions: dict[str, list] = defaultdict(list)
        for key, ev in regions.items():
            n_bj = len(ev["bj"]) + len(ev.get("bj_low", set()))
            if n_bj >= min_bj_for_overlap and ev["positions"]:
                chrom = key[0]
                starts = [p[0] for p in ev["positions"]]
                ends = [p[1] for p in ev["positions"]]
                bj_regions[chrom].append(
                    (int(np.median(starts)), int(np.median(ends)))
                )
        for c in bj_regions:
            bj_regions[c].sort()

        if not bj_regions:
            return cecc_df

        def _frag_has_bj(fc, fs, fe):
            for s, e in bj_regions.get(fc, []):
                if s >= fe:
                    break
                if e <= fs:
                    continue
                ov = min(e, fe) - max(s, fs)
                if ov / max(fe - fs, 1) > 0.3:
                    return True
            return False

        keep_mask = []
        for _, row in cecc_df.iterrows():
            frags_str = row.get("fragments", "")
            if pd.isna(frags_str) or not frags_str:
                # Single-region entry (e.g. rescued fragment) — already vetted
                keep_mask.append(True)
                continue
            else:
                frag_list = []
                for frag in str(frags_str).split(";"):
                    parts = frag.split(":")
                    if len(parts) == 2:
                        coords = parts[1].split("-")
                        if len(coords) == 2:
                            frag_list.append((parts[0], int(coords[0]), int(coords[1])))
            if not frag_list:
                keep_mask.append(True)
                continue

            # If ALL fragments overlap strong BJ regions → likely Mecc/Uecc, not CeccDNA
            all_overlap = all(_frag_has_bj(fc, fs, fe) for fc, fs, fe in frag_list)
            keep_mask.append(not all_overlap)

        n_filtered = sum(1 for k in keep_mask if not k)
        if n_filtered:
            self.logger.info(
                f"CeccDNA intrachr filter: removed {n_filtered} with "
                f"strong BJ evidence at all fragments"
            )
        return cecc_df[keep_mask].reset_index(drop=True)

    # ── CeccDNA LR feature computation ────────────────────────────────

    def _compute_cecc_features(
        self, cecc_df: pd.DataFrame, regions: dict,
        min_bj_for_overlap: int = 3,
    ) -> pd.DataFrame:
        """Compute overlap and intra-chr evidence features for CeccDNA."""
        # Build per-chromosome evidence lookup from intra-chr regions
        evidence_lookup: dict[str, list] = defaultdict(list)
        for key, ev in regions.items():
            chrom = key[0]
            n_bj = len(ev["bj"]) + len(ev.get("bj_low", set()))
            n_disc = len(ev["do"]) + len(ev["di"])
            if ev["positions"]:
                starts = [p[0] for p in ev["positions"]]
                ends = [p[1] for p in ev["positions"]]
                med_s, med_e = int(np.median(starts)), int(np.median(ends))
            else:
                med_s = key[1] * self.config.bin_size
                med_e = key[2] * self.config.bin_size
            evidence_lookup[chrom].append((med_s, med_e, n_bj, n_disc))
        for c in evidence_lookup:
            evidence_lookup[c].sort()

        overlap_ratios = []
        all_overlaps = []
        mean_bjs = []
        mean_discs = []

        for _, row in cecc_df.iterrows():
            frags_str = row.get("fragments", "")
            if pd.isna(frags_str) or not frags_str:
                frags = [(row["chr"], int(row["start"]), int(row["end"]))]
            else:
                frags = []
                for frag in str(frags_str).split(";"):
                    parts = frag.split(":")
                    if len(parts) == 2:
                        coords = parts[1].split("-")
                        if len(coords) == 2:
                            frags.append((parts[0], int(coords[0]), int(coords[1])))
            if not frags:
                frags = [(row["chr"], int(row["start"]), int(row["end"]))]

            n_bj_overlap = 0
            bj_vals = []
            disc_vals = []
            for fc, fs, fe in frags:
                max_bj = 0
                max_disc = 0
                has_bj = False
                for rs, re, rn_bj, rn_disc in evidence_lookup.get(fc, []):
                    if rs >= fe:
                        break
                    if re <= fs:
                        continue
                    ov = min(re, fe) - max(rs, fs)
                    if ov / max(fe - fs, 1) > 0.3:
                        if rn_bj >= min_bj_for_overlap:
                            has_bj = True
                        max_bj = max(max_bj, rn_bj)
                        max_disc = max(max_disc, rn_disc)
                if has_bj:
                    n_bj_overlap += 1
                bj_vals.append(max_bj)
                disc_vals.append(max_disc)

            n_frags = len(frags)
            overlap_ratios.append(round(n_bj_overlap / max(n_frags, 1), 4))
            all_overlaps.append(
                1 if n_bj_overlap == n_frags and n_frags > 0 else 0
            )
            mean_bjs.append(
                round(float(np.mean(bj_vals)), 3) if bj_vals else 0.0
            )
            mean_discs.append(
                round(float(np.mean(disc_vals)), 3) if disc_vals else 0.0
            )

        cecc_df = cecc_df.copy()
        cecc_df["frag_overlap_uecc_ratio"] = overlap_ratios
        cecc_df["all_frags_overlap"] = all_overlaps
        cecc_df["mean_bj_at_frags"] = mean_bjs
        cecc_df["mean_disc_at_frags"] = mean_discs
        return cecc_df

    def _score_cecc_lr(self, cecc_df: pd.DataFrame) -> pd.DataFrame:
        """Apply CeccDNA LR model to score and filter candidates."""
        if cecc_df.empty:
            return cecc_df

        threshold = self.config.cecc_score_threshold
        if threshold <= 0.0:
            return cecc_df  # LR filter disabled

        features = np.zeros((len(cecc_df), len(CECC_FEATURE_NAMES)))
        for i, col in enumerate(CECC_FEATURE_NAMES):
            if col in cecc_df.columns:
                features[:, i] = cecc_df[col].fillna(0).values

        normed = (features - _CECC_SCALER_MEAN) / np.maximum(_CECC_SCALER_STD, 1e-6)
        logits = normed @ _CECC_MODEL_COEFS + _CECC_MODEL_INTERCEPT
        scores = np.array([_sigmoid(float(l)) for l in logits])

        cecc_df = cecc_df.copy()
        cecc_df["score_cecc_lr"] = np.round(scores, 4)

        mask = scores >= threshold
        n_filtered = int((~mask).sum())
        if n_filtered:
            self.logger.info(
                f"CeccDNA LR filter: removed {n_filtered}/{len(cecc_df)}, "
                f"kept {int(mask.sum())} (threshold={threshold})"
            )

        return cecc_df[mask].reset_index(drop=True)

    def _score_mecc_lr(self, um_df: pd.DataFrame) -> pd.DataFrame:
        """Apply MeccDNA LR model to score and filter Mecc candidates.

        Only filters Mecc detections; Uecc pass through unchanged.
        Derived features (bj_ratio, log_length, rpk, disc_ratio) are
        computed on the fly from existing columns.
        """
        if um_df.empty:
            return um_df

        threshold = self.config.mecc_score_threshold
        if threshold <= 0.0:
            return um_df

        mecc_mask = um_df["eccdna_type"] == "Mecc"
        if not mecc_mask.any():
            return um_df

        um_df = um_df.copy()
        # Compute derived features for ALL rows (cheap, avoids index issues)
        n_bj = um_df["n_backjump"].fillna(0)
        n_bj_lo = um_df["n_backjump_low"].fillna(0)
        um_df["bj_ratio"] = n_bj / (n_bj + n_bj_lo).clip(lower=1)
        um_df["log_length"] = np.log10(um_df["length"].fillna(1).clip(lower=1))
        um_df["rpk"] = um_df["n_reads"].fillna(0) / (um_df["length"].fillna(1) / 1000).clip(lower=0.001)
        um_df["disc_ratio"] = um_df["n_discordant"].fillna(0) / um_df["n_reads"].fillna(1).clip(lower=1)

        # Extract features for Mecc rows
        mecc_df = um_df[mecc_mask]
        features = np.zeros((len(mecc_df), len(MECC_FEATURE_NAMES)))
        for i, col in enumerate(MECC_FEATURE_NAMES):
            if col in mecc_df.columns:
                features[:, i] = mecc_df[col].fillna(0).values

        normed = (features - _MECC_SCALER_MEAN) / np.maximum(_MECC_SCALER_STD, 1e-6)
        logits = normed @ _MECC_MODEL_COEFS + _MECC_MODEL_INTERCEPT
        scores = np.array([_sigmoid(float(l)) for l in logits])

        um_df.loc[mecc_mask, "score_mecc_lr"] = np.round(scores, 4)

        # Demote low-score Mecc to Uecc (not delete — they may be real Uecc)
        low_score = mecc_mask & (um_df["score_mecc_lr"].fillna(1.0) < threshold)
        n_demoted = int(low_score.sum())
        if n_demoted:
            um_df.loc[low_score, "eccdna_type"] = "Uecc"
            self.logger.info(
                f"MeccDNA LR filter: demoted {n_demoted}/{int(mecc_mask.sum())} "
                f"low-score Mecc → Uecc (threshold={threshold})"
            )

        return um_df

    def _dedup_mecc_by_sequence(
        self, um_df: pd.DataFrame, similarity: float = 0.9
    ) -> pd.DataFrame:
        """Deduplicate MeccDNA by sequence similarity using cd-hit-est.

        Multiple Mecc detections from different genomic loci of the same
        multi-copy eccDNA are clustered by sequence similarity. Within each
        cluster, the representative with the highest evidence is kept;
        others are dropped (they are redundant copies, not distinct eccDNA).

        Uecc detections pass through unchanged.
        """
        if um_df.empty:
            return um_df

        mecc_mask = um_df["eccdna_type"] == "Mecc"
        if mecc_mask.sum() <= 1:
            return um_df

        import subprocess
        import tempfile
        import shutil

        if not shutil.which("cd-hit-est"):
            self.logger.warning("cd-hit-est not found, skipping Mecc dedup")
            return um_df

        mecc_df = um_df[mecc_mask].copy()
        ref_fasta = pysam.FastaFile(str(self.reference))
        tmpdir = Path(tempfile.mkdtemp())

        # Write Mecc sequences
        fa_path = tmpdir / "mecc_seqs.fasta"
        valid_ids = set()
        with open(fa_path, "w") as f:
            for idx, row in mecc_df.iterrows():
                try:
                    seq = ref_fasta.fetch(
                        str(row["chr"]), int(row["start"]), int(row["end"])
                    )
                    if len(seq) >= 50:
                        f.write(f">{idx}\n{seq}\n")
                        valid_ids.add(idx)
                except (ValueError, KeyError):
                    continue
        ref_fasta.close()

        if len(valid_ids) <= 1:
            shutil.rmtree(tmpdir, ignore_errors=True)
            return um_df

        # Run cd-hit-est
        out_path = tmpdir / "mecc_clustered"
        cmd = (
            f"cd-hit-est -i {fa_path} -o {out_path} "
            f"-c {similarity} -n 8 -d 0 -M 0 -T {self.config.threads} "
            f"-g 1 -s 0.8 -aS 0.8 > /dev/null 2>&1"
        )
        try:
            subprocess.run(cmd, shell=True, timeout=120)
        except (subprocess.TimeoutExpired, Exception) as e:
            self.logger.warning(f"cd-hit-est failed: {e}, skipping Mecc dedup")
            shutil.rmtree(tmpdir, ignore_errors=True)
            return um_df

        # Parse clusters
        clstr_path = Path(f"{out_path}.clstr")
        if not clstr_path.exists():
            shutil.rmtree(tmpdir, ignore_errors=True)
            return um_df

        clusters = []  # list of lists of DataFrame indices
        current_cluster = []
        with open(clstr_path) as f:
            for line in f:
                if line.startswith(">Cluster"):
                    if current_cluster:
                        clusters.append(current_cluster)
                    current_cluster = []
                else:
                    # Parse: "0\t500nt, >12345... *" or "0\t500nt, >12345... at 95.00%"
                    parts = line.strip().split(">")
                    if len(parts) >= 2:
                        idx_str = parts[1].split("...")[0].strip()
                        try:
                            current_cluster.append(int(idx_str))
                        except ValueError:
                            continue
            if current_cluster:
                clusters.append(current_cluster)

        # For each cluster with >1 member, keep the best representative
        drop_indices = set()
        for cluster in clusters:
            if len(cluster) <= 1:
                continue
            # Score: n_backjump_low + n_reads (prefer most evidence)
            best_idx = None
            best_score = -1
            for idx in cluster:
                if idx not in mecc_df.index:
                    continue
                row = mecc_df.loc[idx]
                s = float(row.get("n_backjump_low", 0)) + float(row.get("n_reads", 0))
                if s > best_score:
                    best_score = s
                    best_idx = idx
            for idx in cluster:
                if idx != best_idx and idx in mecc_df.index:
                    drop_indices.add(idx)

        shutil.rmtree(tmpdir, ignore_errors=True)

        if drop_indices:
            n_before = int(mecc_mask.sum())
            um_df = um_df.drop(index=drop_indices).reset_index(drop=True)
            n_after = (um_df["eccdna_type"] == "Mecc").sum()
            self.logger.info(
                f"MeccDNA sequence dedup (cd-hit-est {similarity:.0%}): "
                f"{n_before} → {n_after} Mecc "
                f"({n_before - n_after} redundant copies removed)"
            )

        return um_df

    # ── CeccDNA detection (early — runs before Uecc/Mecc) ──────────────

    def _detect_cecc_early(
        self,
        crosschr: list,
        min_junction_reads: int = 3,
        min_total_reads: int = 8,
        min_output_reads: int = 15,
        min_output_junctions: int = 2,
        cluster_dist: int = 500,
        max_fragment_size: int = 100000,
        max_total_size: int = 500000,
    ) -> tuple[pd.DataFrame, list]:
        """Detect CeccDNA early from cross-chr split reads and discordant pairs.

        Runs BEFORE Uecc/Mecc to claim cross-chr signals first.
        Uses two-tier thresholds: liberal for exclusion zone building,
        strict for reported CeccDNA output.

        Returns:
            (cecc_df, cecc_fragment_intervals) where cecc_fragment_intervals
            is a list of (chr, start, end) tuples for building exclusion zones.
        """
        if not crosschr:
            return pd.DataFrame(), []

        # Step 1: Bin and cluster junctions (8-element tuples)
        junction_bins = defaultdict(lambda: {
            "reads": set(), "split_reads": set(), "disc_reads": set(),
            "positions_a": [], "positions_b": [],
            "strands_a": [], "strands_b": [],
            "mapqs": [],
        })

        for item in crosschr:
            chr_a, pos_a, strand_a, chr_b, pos_b, strand_b, rname, etype = item[:8]
            mapq = item[8] if len(item) > 8 else 0
            if chr_a > chr_b:
                chr_a, pos_a, strand_a, chr_b, pos_b, strand_b = \
                    chr_b, pos_b, strand_b, chr_a, pos_a, strand_a
            bin_a = pos_a // cluster_dist
            bin_b = pos_b // cluster_dist
            key = (chr_a, bin_a, chr_b, bin_b)
            junction_bins[key]["reads"].add(rname)
            if etype == "split":
                junction_bins[key]["split_reads"].add(rname)
            else:
                junction_bins[key]["disc_reads"].add(rname)
            junction_bins[key]["mapqs"].append(mapq)
            junction_bins[key]["positions_a"].append(pos_a)
            junction_bins[key]["positions_b"].append(pos_b)
            junction_bins[key]["strands_a"].append(strand_a)
            junction_bins[key]["strands_b"].append(strand_b)

        # Step 2: Filter junctions with enough support
        junctions = []
        for key, ev in junction_bins.items():
            n_reads = len(ev["reads"])
            if n_reads < min_junction_reads:
                continue
            chr_a, _, chr_b, _ = key
            med_a = int(np.median(ev["positions_a"]))
            med_b = int(np.median(ev["positions_b"]))
            # Dominant strand per side
            dominant_strand_a = max(
                set(ev["strands_a"]), key=ev["strands_a"].count
            )
            dominant_strand_b = max(
                set(ev["strands_b"]), key=ev["strands_b"].count
            )
            mapqs = ev["mapqs"] if ev["mapqs"] else [0]
            junctions.append({
                "chr_a": chr_a, "pos_a": med_a, "strand_a": dominant_strand_a,
                "chr_b": chr_b, "pos_b": med_b, "strand_b": dominant_strand_b,
                "n_reads": n_reads,
                "n_split": len(ev["split_reads"]),
                "n_disc": len(ev["disc_reads"]),
                "mean_mapq": float(np.mean(mapqs)),
                "min_mapq": min(mapqs),
            })

        if not junctions:
            return pd.DataFrame(), []

        self.logger.info(
            f"CeccDNA early: {len(junctions)} cross-chr junctions "
            f"(≥{min_junction_reads} reads)"
        )

        # Step 2b: Collect ALL junction breakpoint regions for exclusion zones
        # Each cross-chr junction indicates a potential CeccDNA fragment boundary
        all_cecc_fragments = []
        for j in junctions:
            all_cecc_fragments.append(
                (j["chr_a"], max(0, j["pos_a"] - 500), j["pos_a"] + 500)
            )
            all_cecc_fragments.append(
                (j["chr_b"], max(0, j["pos_b"] - 500), j["pos_b"] + 500)
            )

        # Step 3: Build junction graph using union-find
        parent = {}

        def find(x):
            while parent.get(x, x) != x:
                parent[x] = parent.get(parent[x], parent[x])
                x = parent[x]
            return x

        def union(a, b):
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb

        for j in junctions:
            node_a = (j["chr_a"], j["pos_a"] // cluster_dist)
            node_b = (j["chr_b"], j["pos_b"] // cluster_dist)
            union(node_a, node_b)

        # Step 4: Group junctions by connected component
        components = defaultdict(list)
        for j in junctions:
            node_a = (j["chr_a"], j["pos_a"] // cluster_dist)
            root = find(node_a)
            components[root].append(j)

        # Step 5: Build CeccDNA candidates with strand closure validation
        results = []
        all_cecc_fragments = []
        cecc_idx = 1

        for comp_junctions in components.values():
            fragments = defaultdict(list)
            total_reads = 0
            for j in comp_junctions:
                fragments[j["chr_a"]].append(j["pos_a"])
                fragments[j["chr_b"]].append(j["pos_b"])
                total_reads += j["n_reads"]

            chroms = set(fragments.keys())
            if len(chroms) < 2:
                continue

            # Build fragment list with padding
            frag_list = []
            oversized = False
            for chrom, positions in fragments.items():
                frag_start = min(positions) - 200
                frag_end = max(positions) + 200
                if frag_end - frag_start > max_fragment_size:
                    oversized = True
                    break
                frag_list.append((chrom, max(0, frag_start), frag_end))

            if oversized:
                continue

            frag_list.sort()
            total_length = sum(e - s for _, s, e in frag_list)
            if total_length > max_total_size:
                continue

            # Collect fragments for exclusion zones BEFORE read threshold
            # (liberal: even low-evidence components mark exclusion zones)
            all_cecc_fragments.extend(frag_list)

            if total_reads < min_total_reads:
                continue

            # Strand closure: even number of strand flips = valid circle
            n_strand_flips = sum(
                1 for j in comp_junctions if j["strand_a"] != j["strand_b"]
            )
            strand_closure = (n_strand_flips % 2 == 0)
            n_junctions = len(comp_junctions)

            # Strict filter for CeccDNA output
            if total_reads < min_output_reads:
                continue
            if n_junctions < min_output_junctions:
                continue

            rep_chr, rep_start, rep_end = frag_list[0]
            frag_str = ";".join(f"{c}:{s}-{e}" for c, s, e in frag_list)

            # Score with soft strand closure penalty
            base_score = min(1.0, total_reads / 20.0)
            if not strand_closure:
                base_score *= 0.7

            # Compute LR features
            junction_reads = [j["n_reads"] for j in comp_junctions]
            mean_rpj = float(np.mean(junction_reads))
            min_rpj = min(junction_reads)
            n_split_total = sum(j.get("n_split", 0) for j in comp_junctions)
            n_disc_total = sum(j.get("n_disc", 0) for j in comp_junctions)
            sdr = n_split_total / max(n_split_total + n_disc_total, 1)
            log_treads = float(np.log10(max(total_reads, 1)))
            frag_lengths = [e - s for _, s, e in frag_list]
            mmfr = max(frag_lengths) / max(min(frag_lengths), 1)
            log_tlen = float(np.log10(max(total_length, 1)))
            # MAPQ features
            junc_mapqs = [j.get("mean_mapq", 0) for j in comp_junctions]
            junc_min_mapqs = [j.get("min_mapq", 0) for j in comp_junctions]
            mean_junc_mapq = float(np.mean(junc_mapqs))
            min_junc_mapq = min(junc_min_mapqs)

            eid = f"NGS_Cecc{cecc_idx:05d}"
            cecc_idx += 1

            results.append({
                "eccDNA_id": eid,
                "chr": rep_chr,
                "start": rep_start,
                "end": rep_end,
                "length": total_length,
                "eccdna_type": "Cecc",
                "state": "Confirmed",
                "score": round(base_score, 4),
                "n_backjump": 0,
                "n_backjump_low": 0,
                "n_discordant": total_reads,
                "n_softclip": 0,
                "n_reads": total_reads,
                "n_fragments": len(frag_list),
                "n_chromosomes": len(chroms),
                "n_junctions": n_junctions,
                "fragments": frag_str,
                "mean_reads_per_junction": round(mean_rpj, 3),
                "min_reads_per_junction": min_rpj,
                "split_disc_ratio": round(sdr, 4),
                "log_total_reads": round(log_treads, 4),
                "max_min_frag_ratio": round(mmfr, 3),
                "log_total_length": round(log_tlen, 4),
                "strand_closure": int(strand_closure),
                "n_strand_flips": n_strand_flips,
                "mean_junction_mapq": round(mean_junc_mapq, 1),
                "min_junction_mapq": min_junc_mapq,
            })

        return pd.DataFrame(results), all_cecc_fragments

    # ── Exclusion zone builder ────────────────────────────────────────

    @staticmethod
    def _build_exclusion_zones(
        cecc_fragments: list, padding: int = 200,
    ) -> dict:
        """Merge CeccDNA fragment intervals with padding into exclusion zones.

        Returns:
            dict of chr -> list of (start, end) merged intervals.
        """
        if not cecc_fragments:
            return {}

        by_chr = defaultdict(list)
        for chrom, start, end in cecc_fragments:
            by_chr[chrom].append((max(0, start - padding), end + padding))

        merged = {}
        for chrom, intervals in by_chr.items():
            intervals.sort()
            result = [intervals[0]]
            for s, e in intervals[1:]:
                if s <= result[-1][1]:
                    result[-1] = (result[-1][0], max(result[-1][1], e))
                else:
                    result.append((s, e))
            merged[chrom] = result

        return merged
