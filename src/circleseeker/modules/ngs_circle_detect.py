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

        # Step 1: Extract evidence (intra-chr + cross-chr)
        regions, crosschr = self._extract_evidence(bam_path)
        self.logger.info(
            f"NGS detect: {len(regions)} intra-chr regions, "
            f"{len(crosschr)} cross-chr junctions"
        )

        # Step 2-5: Uecc/Mecc pipeline
        um_df = pd.DataFrame()
        if regions:
            scored = self._score_regions(regions)
            kept = self._dedup(scored)
            self.logger.info(f"NGS detect: {len(kept)} after dedup")
            results = self._filter_and_classify(kept)
            self.logger.info(f"NGS detect: {len(results)} Uecc/Mecc circles")
            if results:
                um_df = pd.DataFrame(results)
                um_df = self.reclassify_mecc(um_df)

        # Step 6: CeccDNA detection from cross-chr evidence
        # Run after Uecc/Mecc so we can filter CeccDNA overlapping known circles
        cecc_df = self._detect_cecc(crosschr, existing_circles=um_df)
        self.logger.info(f"NGS detect: {len(cecc_df)} CeccDNA circles")

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
            (regions, crosschr_junctions)
            - regions: dict of intra-chromosomal evidence (for Uecc/Mecc)
            - crosschr_junctions: list of cross-chr split read junctions
              for CeccDNA detection, each as (chr_a, pos_a, chr_b, pos_b, read_name)
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
        return dict(regions), crosschr_junctions

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
                # Breakpoint: end of primary clip → start of SA
                pri_bp = rend if not read.is_reverse else rstart
                sa_bp = sa_pos if not sa_rev else sa_pos + sa_ref_len
                crosschr.append((chrom, pri_bp, sa_chr, sa_bp, read.query_name))
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
            crosschr.append((chrom, pri_pos, mate_chr, mate_pos, read.query_name))
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

    # ── Scoring ────────────────────────────────────────────────────────

    def _score_regions(self, regions: dict) -> list[dict]:
        scored = []
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

        # Update DataFrame — bidirectional reclassification
        circles_df = circles_df.copy()
        all_aligned = set(circle_loci.keys())

        # ≥2 loci → Mecc (regardless of initial classification)
        mask_to_mecc = circles_df["eccDNA_id"].isin(mecc_ids)
        circles_df.loc[mask_to_mecc, "eccdna_type"] = "Mecc"
        circles_df.loc[mask_to_mecc, "n_genomic_loci"] = (
            circles_df.loc[mask_to_mecc, "eccDNA_id"].map(mecc_loci_count)
        )

        # 1 locus → Uecc (fix Mecc that were misclassified by MAPQ)
        single_locus_ids = all_aligned - mecc_ids
        mask_to_uecc = (
            (circles_df["eccDNA_id"].isin(single_locus_ids))
            & (circles_df["eccdna_type"] == "Mecc")
        )
        circles_df.loc[mask_to_uecc, "eccdna_type"] = "Uecc"
        circles_df.loc[mask_to_uecc, "n_genomic_loci"] = 1

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

    # ── CeccDNA detection ─────────────────────────────────────────────

    def _detect_cecc(
        self,
        crosschr: list,
        existing_circles: Optional[pd.DataFrame] = None,
        min_junction_reads: int = 5,
        min_total_reads: int = 15,
        cluster_dist: int = 500,
        max_fragment_size: int = 100000,
        max_total_size: int = 500000,
    ) -> pd.DataFrame:
        """Detect CeccDNA from cross-chromosome split reads and discordant pairs.

        Strategy:
        1. Cluster cross-chr junctions by (chr_a, chr_b, pos_a_bin, pos_b_bin)
        2. Filter clusters with enough supporting reads
        3. Build a junction graph: nodes = chromosomal fragments, edges = junctions
        4. Connected components with ≥2 chromosomes = CeccDNA candidates
        5. Filter by fragment size and total size to remove over-clustered SVs
        """
        if not crosschr:
            return pd.DataFrame()

        # Step 1: Bin and cluster junctions
        junction_bins = defaultdict(lambda: {
            "reads": set(), "positions_a": [], "positions_b": [],
        })

        for chr_a, pos_a, chr_b, pos_b, rname in crosschr:
            if chr_a > chr_b:
                chr_a, pos_a, chr_b, pos_b = chr_b, pos_b, chr_a, pos_a
            bin_a = pos_a // cluster_dist
            bin_b = pos_b // cluster_dist
            key = (chr_a, bin_a, chr_b, bin_b)
            junction_bins[key]["reads"].add(rname)
            junction_bins[key]["positions_a"].append(pos_a)
            junction_bins[key]["positions_b"].append(pos_b)

        # Step 2: Filter junctions with enough support
        junctions = []
        for key, ev in junction_bins.items():
            n_reads = len(ev["reads"])
            if n_reads < min_junction_reads:
                continue
            chr_a, _, chr_b, _ = key
            med_a = int(np.median(ev["positions_a"]))
            med_b = int(np.median(ev["positions_b"]))
            junctions.append({
                "chr_a": chr_a, "pos_a": med_a,
                "chr_b": chr_b, "pos_b": med_b,
                "n_reads": n_reads,
            })

        if not junctions:
            return pd.DataFrame()

        self.logger.info(
            f"CeccDNA: {len(junctions)} cross-chr junctions "
            f"(≥{min_junction_reads} reads)"
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

        # Build existing circle lookup for filtering CeccDNA FP
        existing_by_chr = defaultdict(list)
        if existing_circles is not None and not existing_circles.empty:
            for _, row in existing_circles.iterrows():
                existing_by_chr[row["chr"]].append(
                    (int(row["start"]), int(row["end"]))
                )
            for c in existing_by_chr:
                existing_by_chr[c].sort()

        def _frag_overlaps_existing(chrom, start, end):
            """Check if a fragment overlaps any existing Uecc/Mecc detection."""
            for s, e in existing_by_chr.get(chrom, []):
                if s >= end:
                    break
                if e <= start:
                    continue
                ov = min(e, end) - max(s, start)
                if ov / max(end - start, 1) > 0.3:
                    return True
            return False

        # Step 5: Build CeccDNA candidates with filtering
        results = []
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

            if total_reads < min_total_reads:
                continue

            # Build fragment list with padding
            frag_list = []
            oversized = False
            for chrom, positions in fragments.items():
                frag_start = min(positions) - 200
                frag_end = max(positions) + 200
                frag_size = frag_end - frag_start
                if frag_size > max_fragment_size:
                    oversized = True
                    break
                frag_list.append((chrom, max(0, frag_start), frag_end))

            if oversized:
                continue

            frag_list.sort()
            total_length = sum(e - s for _, s, e in frag_list)
            if total_length > max_total_size:
                continue

            # Filter: if ALL fragments overlap existing Uecc/Mecc, this is
            # likely a repetitive region (MeccDNA/NUMT) not a real CeccDNA
            if existing_by_chr and all(
                _frag_overlaps_existing(c, s, e) for c, s, e in frag_list
            ):
                continue

            n_junctions = len(comp_junctions)
            rep_chr, rep_start, rep_end = frag_list[0]
            frag_str = ";".join(f"{c}:{s}-{e}" for c, s, e in frag_list)

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
                "score": round(min(1.0, total_reads / 20.0), 4),
                "n_backjump": 0,
                "n_backjump_low": 0,
                "n_discordant": total_reads,
                "n_softclip": 0,
                "n_reads": total_reads,
                "n_fragments": len(frag_list),
                "n_chromosomes": len(chroms),
                "n_junctions": n_junctions,
                "fragments": frag_str,
            })

        return pd.DataFrame(results)
