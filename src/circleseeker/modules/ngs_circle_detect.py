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
        """Detect eccDNA from an existing sorted BAM."""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Step 1: Extract evidence
        regions = self._extract_evidence(bam_path)
        self.logger.info(f"NGS detect: {len(regions)} candidate regions")

        if not regions:
            return pd.DataFrame()

        # Step 2: Build feature matrix and score
        scored = self._score_regions(regions)

        # Step 3: Deduplicate overlapping candidates
        kept = self._dedup(scored)
        self.logger.info(f"NGS detect: {len(kept)} after dedup")

        # Step 4: Filter by score and classify
        results = self._filter_and_classify(kept)
        self.logger.info(f"NGS detect: {len(results)} final circles")

        if not results:
            return pd.DataFrame()

        df = pd.DataFrame(results)

        # Step 5: Reclassify multi-locus circles as MeccDNA
        df = self.reclassify_mecc(df)

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

    def _extract_evidence(self, bam_path: Path) -> dict:
        """Extract per-region evidence from BAM.

        Reads are NOT filtered by MAPQ. Instead, back-jump reads are
        stored in separate buckets by MAPQ:
          - "bj": MAPQ ≥ mapq_threshold (UeccDNA signal)
          - "bj_low": MAPQ < mapq_threshold (MeccDNA signal)
        This allows downstream classification of UeccDNA vs MeccDNA
        based on mapping ambiguity.
        """
        cfg = self.config
        bam = pysam.AlignmentFile(str(bam_path), "rb")

        regions = defaultdict(lambda: {
            "bj": set(), "bj_low": set(), "fj": set(),
            "do": set(), "di": set(),
            "cl": set(), "cr": set(),
            "positions": [],
        })

        for read in bam.fetch():
            if read.is_unmapped or read.is_secondary or read.is_duplicate:
                continue

            chrom = read.reference_name
            rstart, rend = read.reference_start, read.reference_end

            # ── Split reads (all MAPQ, stratified internally) ──
            if read.has_tag("SA"):
                self._process_split(read, chrom, rstart, rend, regions)

            # ── Discordant pairs (require min MAPQ) ──
            if (
                read.mapping_quality >= cfg.min_mapq
                and read.is_paired
                and not read.mate_is_unmapped
                and read.is_read1
            ):
                self._process_discordant(read, chrom, rstart, rend, regions)

            # ── Soft-clips ──
            if read.mapping_quality >= cfg.min_mapq and read.cigartuples:
                self._process_softclip(read, chrom, rstart, rend, regions)

        bam.close()
        return dict(regions)

    def _process_split(self, read, chrom, rstart, rend, regions):
        cfg = self.config
        for entry in read.get_tag("SA").strip(";").split(";"):
            parts = entry.split(",")
            if len(parts) < 6:
                continue
            sa_chr, sa_pos = parts[0], int(parts[1]) - 1
            sa_mapq = int(parts[4])
            if sa_chr != chrom or sa_mapq < cfg.min_mapq:
                continue

            sa_rev = parts[2] == "-"
            cigar_ops = re.findall(r"(\d+)([MIDNSHP=X])", parts[3])
            sa_ref_len = sum(int(n) for n, op in cigar_ops if op in "MDN=X")
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
                # Back-jump: stratify by MAPQ
                if read.mapping_quality >= cfg.min_mapq:
                    regions[key]["bj"].add(read.query_name)
                else:
                    regions[key]["bj_low"].add(read.query_name)
            else:
                regions[key]["fj"].add(read.query_name)
            regions[key]["positions"].append((left, right))

    def _process_discordant(self, read, chrom, rstart, rend, regions):
        cfg = self.config
        if chrom != read.next_reference_name:
            return  # cross-chr handled separately in future
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

    # ── Scoring ────────────────────────────────────────────────────────

    def _score_regions(self, regions: dict) -> list[dict]:
        scored = []
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
            features = np.array([
                total_bj, len(ev["fj"]),
                len(ev["do"]), len(ev["di"]),
                len(ev["cl"]), len(ev["cr"]),
                np.log10(max(size, 1)),
            ])

            # Normalize and score
            normed = (features - _SCALER_MEAN) / np.maximum(_SCALER_STD, 1e-6)
            logit = float(np.dot(normed, _MODEL_COEFS) + _MODEL_INTERCEPT)
            score = _sigmoid(logit)

            scored.append({
                "chr": chrom, "start": med_s, "end": med_e, "size": size,
                "score": score,
                "n_backjump": len(ev["bj"]),
                "n_backjump_low": len(ev.get("bj_low", set())),
                "n_fwdjump": len(ev["fj"]),
                "n_disc_outward": len(ev["do"]), "n_disc_inward": len(ev["di"]),
                "n_clip_left": len(ev["cl"]), "n_clip_right": len(ev["cr"]),
                "n_reads": len(ev["bj"] | ev.get("bj_low", set()) | ev["fj"] | ev["do"] | ev["di"]),
            })

        return scored

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

        # Update DataFrame
        circles_df = circles_df.copy()
        mask = circles_df["eccDNA_id"].isin(mecc_ids)
        circles_df.loc[mask, "eccdna_type"] = "Mecc"
        circles_df.loc[mask, "n_genomic_loci"] = circles_df.loc[mask, "eccDNA_id"].map(
            mecc_loci_count
        )

        n_mecc = mask.sum()
        self.logger.info(
            f"MeccDNA reclassification: {n_mecc}/{len(circles_df)} circles "
            f"reclassified as MeccDNA (≥{min_loci} genomic loci)"
        )

        # Cleanup
        import shutil
        shutil.rmtree(tmpdir, ignore_errors=True)

        return circles_df
