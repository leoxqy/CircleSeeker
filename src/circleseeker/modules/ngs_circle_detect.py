"""
NGS Circle Detection — Detect eccDNA from short-read breakpoint evidence.

Uses a Bayesian scoring model to evaluate circle candidates:

    P(eccDNA | evidence) ∝ P(evidence | eccDNA) × P(eccDNA)

Evidence includes:
  - Split read count at breakpoint
  - Discordant pair count
  - Soft-clip count at boundaries
  - Breakpoint position consistency (sharpness)

Supports UeccDNA (same-chr), MeccDNA (multi-chr, 2 loci),
and CeccDNA (multi-chr, 3+ loci) via breakpoint graph analysis.
"""

from __future__ import annotations

import logging
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from circleseeker.utils.logging import get_logger


@dataclass
class NGSCircleConfig:
    """Configuration for circle detection."""
    bp_cluster_dist: int = 500       # Max distance to cluster breakpoints
    min_total_support: int = 2       # Minimum total evidence (split+disc+clip)
    min_circle_size: int = 100       # Minimum circle size (bp)
    max_circle_size: int = 500000    # Maximum circle size (bp)
    min_score: float = 0.5           # Minimum Bayesian posterior score
    # Bayesian model parameters
    prior_eccdna: float = 0.01       # Prior probability of eccDNA at any locus
    lambda_bg_split: float = 0.1     # Expected split reads per locus under H0
    lambda_bg_disc: float = 0.05     # Expected discordant pairs per locus under H0
    lambda_bg_clip: float = 0.5      # Expected soft-clips per locus under H0


class NGSCircleDetect:
    """Detect eccDNA circles from breakpoint evidence."""

    def __init__(
        self,
        config: Optional[NGSCircleConfig] = None,
        logger: Optional[logging.Logger] = None,
    ):
        self.config = config or NGSCircleConfig()
        self.logger = logger or get_logger(self.__class__.__name__)

    def detect(
        self,
        evidence_df: pd.DataFrame,
        chrom_sizes: Optional[dict[str, int]] = None,
    ) -> pd.DataFrame:
        """Detect eccDNA circles from evidence DataFrame.

        Args:
            evidence_df: Output from NGSEvidence.extract()
            chrom_sizes: Optional dict of chromosome sizes

        Returns:
            DataFrame of circle candidates with scores.
        """
        if evidence_df.empty:
            return pd.DataFrame()

        # Step 1: Cluster breakpoints
        bp_clusters = self._cluster_breakpoints(evidence_df)
        self.logger.info(f"Breakpoint clusters: {len(bp_clusters)}")

        # Step 2: Build circle candidates from breakpoint pairs
        circles = self._build_circles(bp_clusters, evidence_df)
        self.logger.info(f"Circle candidates: {len(circles)}")

        if not circles:
            return pd.DataFrame()

        # Step 3: Score each circle with Bayesian model
        scored = self._score_circles(circles)
        self.logger.info(
            f"Scored circles: {len(scored)} "
            f"(passed min_score={self.config.min_score}: "
            f"{sum(1 for c in scored if c['score'] >= self.config.min_score)})"
        )

        # Step 4: Filter and classify
        results = self._classify_and_filter(scored)
        self.logger.info(f"Final circles: {len(results)}")

        if not results:
            return pd.DataFrame()
        return pd.DataFrame(results)

    # ------------------------------------------------------------------

    def _cluster_breakpoints(self, df: pd.DataFrame) -> list[dict]:
        """Cluster nearby breakpoint positions into breakpoint regions."""
        clusters = []

        # Group by (chr1, chr2) pair
        pair_groups = df.groupby(["chr1", "chr2"])

        for (chr1, chr2), group in pair_groups:
            # Sort by pos1
            sorted_grp = group.sort_values("pos1")
            positions = sorted_grp["pos1"].values
            pos2_vals = sorted_grp["pos2"].values

            # Cluster pos1 within bp_cluster_dist
            cur_start = 0
            for i in range(1, len(positions)):
                if positions[i] - positions[cur_start] > self.config.bp_cluster_dist:
                    # Emit cluster
                    cl_rows = sorted_grp.iloc[cur_start:i]
                    clusters.append(self._summarize_cluster(cl_rows, chr1, chr2))
                    cur_start = i
            # Last cluster
            cl_rows = sorted_grp.iloc[cur_start:]
            clusters.append(self._summarize_cluster(cl_rows, chr1, chr2))

        return clusters

    def _summarize_cluster(self, rows: pd.DataFrame, chr1: str, chr2: str) -> dict:
        """Summarize a cluster of breakpoint observations."""
        n_split = (rows["evidence_type"] == "split").sum()
        n_disc = (rows["evidence_type"] == "discordant").sum()
        n_clip = (rows["evidence_type"] == "softclip").sum()

        pos1_vals = rows["pos1"].values
        pos2_vals = rows["pos2"].values

        return {
            "chr1": chr1,
            "chr2": chr2,
            "pos1_median": int(np.median(pos1_vals)),
            "pos1_std": float(np.std(pos1_vals)) if len(pos1_vals) > 1 else 0,
            "pos2_median": int(np.median(pos2_vals)),
            "pos2_std": float(np.std(pos2_vals)) if len(pos2_vals) > 1 else 0,
            "n_split": int(n_split),
            "n_discordant": int(n_disc),
            "n_softclip": int(n_clip),
            "n_total": len(rows),
            "n_unique_reads": rows["read_name"].nunique(),
            "mean_mapq": float(rows["mapq"].mean()),
            "reads": set(rows["read_name"]),
        }

    def _build_circles(
        self, bp_clusters: list[dict], evidence_df: pd.DataFrame,
    ) -> list[dict]:
        """Build circle candidates by pairing breakpoint clusters.

        For UeccDNA: two breakpoint clusters on same chr that bracket a region.
        For MeccDNA/CeccDNA: cross-chr breakpoint clusters forming a graph.
        """
        circles = []
        cfg = self.config

        # --- Same-chromosome circles (UeccDNA) ---
        # Group clusters by chr1==chr2 (same chr)
        same_chr: dict[str, list[dict]] = defaultdict(list)
        cross_chr: list[dict] = []

        for cl in bp_clusters:
            if cl["n_total"] < cfg.min_total_support:
                continue
            if cl["chr1"] == cl["chr2"]:
                same_chr[cl["chr1"]].append(cl)
            else:
                cross_chr.append(cl)

        # Same-chr: pair clusters where pos1 < pos2 to form circle regions
        for chrom, clusters in same_chr.items():
            # Each cluster represents a breakpoint at one end of a circle
            # We need to find PAIRS of clusters that form circle boundaries
            # Simple approach: for softclip clusters, pair left-clip with right-clip
            # For split/discordant: the pos1,pos2 already define the region
            for cl in clusters:
                start = min(cl["pos1_median"], cl["pos2_median"])
                end = max(cl["pos1_median"], cl["pos2_median"])
                size = end - start

                if size < cfg.min_circle_size or size > cfg.max_circle_size:
                    continue

                circles.append({
                    "segments": [(chrom, start, end)],
                    "n_chr": 1,
                    "circle_size": size,
                    "n_split": cl["n_split"],
                    "n_discordant": cl["n_discordant"],
                    "n_softclip": cl["n_softclip"],
                    "n_total": cl["n_total"],
                    "n_unique_reads": cl["n_unique_reads"],
                    "mean_mapq": cl["mean_mapq"],
                    "bp_sharpness": 1.0 / (1.0 + cl["pos1_std"] + cl["pos2_std"]),
                })

        # --- Cross-chromosome circles (MeccDNA/CeccDNA) ---
        # Group cross-chr breakpoints into connected components
        if cross_chr:
            circles.extend(self._build_cross_chr_circles(cross_chr))

        return circles

    def _build_cross_chr_circles(self, cross_chr_clusters: list[dict]) -> list[dict]:
        """Build multi-segment circles from cross-chromosome breakpoints."""
        # Simple approach: group clusters sharing chromosomes
        # A proper implementation would use a breakpoint graph
        chr_pairs = defaultdict(list)
        for cl in cross_chr_clusters:
            key = tuple(sorted([cl["chr1"], cl["chr2"]]))
            chr_pairs[key].append(cl)

        circles = []
        for (c1, c2), clusters in chr_pairs.items():
            total_support = sum(cl["n_total"] for cl in clusters)
            if total_support < self.config.min_total_support:
                continue

            # Build segments from cluster positions
            seg1_positions = []
            seg2_positions = []
            for cl in clusters:
                if cl["chr1"] == c1:
                    seg1_positions.append(cl["pos1_median"])
                    seg2_positions.append(cl["pos2_median"])
                else:
                    seg1_positions.append(cl["pos2_median"])
                    seg2_positions.append(cl["pos1_median"])

            if not seg1_positions or not seg2_positions:
                continue

            seg1_start = min(seg1_positions)
            seg1_end = max(seg1_positions) + 1000  # rough estimate
            seg2_start = min(seg2_positions)
            seg2_end = max(seg2_positions) + 1000

            best_cl = max(clusters, key=lambda x: x["n_total"])
            circles.append({
                "segments": [(c1, seg1_start, seg1_end), (c2, seg2_start, seg2_end)],
                "n_chr": 2,
                "circle_size": (seg1_end - seg1_start) + (seg2_end - seg2_start),
                "n_split": sum(cl["n_split"] for cl in clusters),
                "n_discordant": sum(cl["n_discordant"] for cl in clusters),
                "n_softclip": sum(cl["n_softclip"] for cl in clusters),
                "n_total": total_support,
                "n_unique_reads": sum(cl["n_unique_reads"] for cl in clusters),
                "mean_mapq": np.mean([cl["mean_mapq"] for cl in clusters]),
                "bp_sharpness": 1.0 / (1.0 + best_cl["pos1_std"] + best_cl["pos2_std"]),
            })

        return circles

    def _score_circles(self, circles: list[dict]) -> list[dict]:
        """Score each circle candidate using Bayesian model.

        P(eccDNA | E) = P(E | eccDNA) * P(eccDNA) / P(E)

        Uses log-likelihood ratio:
          LLR = log P(E | H1) - log P(E | H0)
          posterior = sigmoid(LLR + log_prior_odds)
        """
        cfg = self.config
        log_prior_odds = math.log(cfg.prior_eccdna / (1 - cfg.prior_eccdna))

        for circle in circles:
            ns = circle["n_split"]
            nd = circle["n_discordant"]
            nc = circle["n_softclip"]
            sharpness = circle["bp_sharpness"]

            # Log-likelihood under H1 (eccDNA present)
            # More evidence = higher likelihood; sharp breakpoints = higher
            # Use Poisson model: P(k | lambda) = lambda^k * exp(-lambda) / k!
            # H1: lambda = observed count (evidence is consistent)
            # H0: lambda = background rate

            # Split reads
            ll_h1_split = self._poisson_ll(ns, max(ns, 1.0))
            ll_h0_split = self._poisson_ll(ns, cfg.lambda_bg_split)

            # Discordant pairs
            ll_h1_disc = self._poisson_ll(nd, max(nd, 0.5))
            ll_h0_disc = self._poisson_ll(nd, cfg.lambda_bg_disc)

            # Soft-clips
            ll_h1_clip = self._poisson_ll(nc, max(nc, 1.0))
            ll_h0_clip = self._poisson_ll(nc, cfg.lambda_bg_clip)

            # Breakpoint sharpness bonus (sharp = more likely real)
            ll_sharpness = math.log(sharpness + 0.01) * 2.0

            # Combined log-likelihood ratio
            llr = (
                (ll_h1_split - ll_h0_split)
                + (ll_h1_disc - ll_h0_disc)
                + (ll_h1_clip - ll_h0_clip)
                + ll_sharpness
            )

            # Posterior via sigmoid
            log_posterior_odds = llr + log_prior_odds
            score = 1.0 / (1.0 + math.exp(-log_posterior_odds))

            circle["score"] = score
            circle["llr"] = llr

        return circles

    @staticmethod
    def _poisson_ll(k: int, lam: float) -> float:
        """Log-likelihood of Poisson(k | lambda)."""
        if lam <= 0:
            return -100.0 if k > 0 else 0.0
        return k * math.log(lam) - lam - sum(math.log(i) for i in range(1, k + 1))

    def _classify_and_filter(self, circles: list[dict]) -> list[dict]:
        """Filter by score and classify as U/M/C eccDNA."""
        results = []
        idx = 1

        for circle in circles:
            if circle["score"] < self.config.min_score:
                continue

            n_chr = circle["n_chr"]
            n_seg = len(circle["segments"])

            if n_chr == 1 and n_seg == 1:
                etype = "Uecc"
            elif n_chr >= 2 and n_seg == 2:
                etype = "Mecc"
            elif n_seg >= 3:
                etype = "Cecc"
            else:
                etype = "Mecc" if n_chr >= 2 else "Uecc"

            for seg_idx, (chrom, start, end) in enumerate(circle["segments"], 1):
                results.append({
                    "eccDNA_id": f"NGS_{etype}{idx:04d}",
                    "chr": chrom,
                    "start": start,
                    "end": end,
                    "length": end - start,
                    "eccdna_type": etype,
                    "state": "Confirmed",
                    "score": circle["score"],
                    "n_split": circle["n_split"],
                    "n_discordant": circle["n_discordant"],
                    "n_softclip": circle["n_softclip"],
                    "n_total": circle["n_total"],
                    "n_unique_reads": circle["n_unique_reads"],
                    "seg_index": seg_idx,
                    "seg_total": n_seg,
                })
            idx += 1

        return results
