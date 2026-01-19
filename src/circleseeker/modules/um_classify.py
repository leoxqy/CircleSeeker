"""
U/Mecc-Classify - Simple eccDNA Classification Module

This module performs initial classification of eccDNA from alignment results:
- Uecc: Unique/simple circular DNA
- Mecc: Multiple-copy repeat circular DNA
- Unclassified: All alignments from queries not classified as Uecc/Mecc

Migrated from step4_gatekeeper.py to the CircleSeeker architecture.
"""

from __future__ import annotations

from circleseeker.utils.logging import get_logger
from circleseeker.utils.column_standards import ColumnStandard
from pathlib import Path
from typing import Any, Optional, Iterable
import pandas as pd
import numpy as np
import logging
import json


class UMeccClassifier:
    """Classify eccDNA into Uecc/Mecc buckets."""

    # Import unified constants
    from circleseeker.constants import (
        MAPQ_LOW_THRESHOLD as _MAPQ_LOW,
        IDENTITY_LOW_THRESHOLD as _IDENTITY_LOW,
    )

    # Heuristic thresholds used for evidence flags and score scaling.
    # These do NOT affect classification unless mapq_u_min is explicitly enabled.
    MAPQ_MAX = 60  # minimap2 typically caps MAPQ at ~60
    MAPQ_LOW_THRESHOLD = _MAPQ_LOW
    IDENTITY_LOW_THRESHOLD = _IDENTITY_LOW

    # Alignment column names (BLAST outfmt 6 compatible + optional MAPQ)
    ALIGNMENT_COLUMNS = [
        "query_id",
        "subject_id",
        "identity",
        "alignment_length",
        "mismatches",
        "gap_opens",
        "q_start",
        "q_end",
        "s_start",
        "s_end",
        "evalue",
        "bit_score",
        "sstrand",
    ]
    ALIGNMENT_COLUMNS_WITH_MAPQ = ALIGNMENT_COLUMNS + ["mapq"]

    def __init__(
        self,
        gap_threshold: float = 10.0,
        min_full_length_coverage: float = 95.0,
        max_identity_gap_for_mecc: Optional[float] = 5.0,
        theta_full: Optional[float] = None,
        theta_u: Optional[float] = None,
        theta_m: Optional[float] = None,
        theta_u2_max: Optional[float] = None,
        mapq_u_min: int = 0,
        mapq_m_ambiguous_threshold: int = 0,
        mecc_identity_gap_threshold: float = 0,  # Disabled by default
        u_secondary_min_frac: float = 0.01,
        u_secondary_min_bp: int = 50,
        u_contig_gap_bp: int = 1000,
        u_secondary_max_ratio: float = 0.05,
        u_high_coverage_threshold: float = 0.98,
        u_high_mapq_threshold: int = 50,
        theta_locus: float = 0.95,
        pos_tol_bp: int = 50,
        logger: Optional[logging.Logger] = None,
    ) -> None:
        """
        Initialize U/Mecc classifier

        Args:
            gap_threshold: Maximum gap percentage for quality filtering (default: 10%)
            min_full_length_coverage: Legacy full-length threshold; accepts either percent
                (e.g. 95) or fraction (e.g. 0.95). Kept for backward compatibility.
            max_identity_gap_for_mecc: Legacy parameter kept for backward compatibility
                (not used by the coverage model).
            theta_full: Coverage fraction threshold for U/Mecc decisions (default: 0.95).
            theta_u: Coverage fraction threshold for Uecc decisions (default: theta_full).
            theta_m: Coverage fraction threshold for Mecc decisions (default: theta_full).
            theta_u2_max: Upper bound on the 2nd-best locus coverage for Uecc (default: 0.05).
                Set to 1.0 to disable the uniqueness constraint.
            mapq_u_min: If >0, require best-locus minimap2 MAPQ >= this threshold for Uecc.
                Intended to reduce false-positive Uecc calls in repetitive regions. Default 0 (off).
            mapq_m_ambiguous_threshold: When ALL loci in a potential Mecc have MAPQ below
                this threshold, the query is NOT classified as Mecc (left as unclassified).
                This can help avoid false-positive Mecc calls for Uecc from repetitive regions.
                WARNING: Most true MeccDNA also have low MAPQ (they come from repetitive regions).
                Enabling this filter (e.g., threshold=20) can remove ~85% misclassified UeccDNA
                but also loses ~67% of true MeccDNA. Use with caution.
                Default 0 (disabled). Set to 20 to enable a moderate filter.
            mecc_identity_gap_threshold: When the gap between max and second-max locus
                identity exceeds this threshold, reject Mecc classification.
                Originally designed to detect Uecc from repetitive regions that align
                to multiple loci but have one clearly superior match.
                However, true Mecc can also have identity gaps (~1%) due to mutation
                between repeat copies, causing false rejections. The coverage-based
                Mecc criteria (>= 2 loci with >= 95% coverage) is usually sufficient.
                Default 0 (disabled). Set to 1.5+ to enable a conservative filter.
            u_secondary_min_frac: Secondary coverage fraction threshold that triggers a
                contiguity check when attempting to call Uecc (default: 0.01).
            u_secondary_min_bp: Secondary coverage absolute bp threshold that triggers a
                contiguity check when attempting to call Uecc (default: 50).
            u_contig_gap_bp: Maximum allowed genomic gap (bp) between alignments considered
                part of the same locus when attempting to call Uecc (default: 1000).
            theta_locus: Locus clustering reciprocal-overlap threshold (default: 0.95).
            pos_tol_bp: Locus clustering boundary tolerance in bp (default: 50).
            logger: Optional logger instance
        """
        self.gap_threshold = gap_threshold
        self.min_full_length_coverage = min_full_length_coverage
        self.max_identity_gap_for_mecc = max_identity_gap_for_mecc

        theta_full_frac = (
            self._as_fraction(theta_full)
            if theta_full is not None
            else self._as_fraction(min_full_length_coverage)
        )
        self.theta_u = (
            self._as_fraction(theta_u) if theta_u is not None else float(theta_full_frac)
        )
        self.theta_m = (
            self._as_fraction(theta_m) if theta_m is not None else float(theta_full_frac)
        )
        # For backward compatibility: keep theta_full as an alias of the historical threshold.
        self.theta_full = float(theta_full_frac)

        # U-second-locus constraint: strict uniqueness by default; set to 1.0 to disable.
        self.theta_u2_max = (
            self._as_fraction(theta_u2_max) if theta_u2_max is not None else 0.05
        )
        self.mapq_u_min = int(mapq_u_min) if mapq_u_min is not None else 0
        self.mapq_m_ambiguous_threshold = (
            int(mapq_m_ambiguous_threshold) if mapq_m_ambiguous_threshold is not None else 0
        )
        self.mecc_identity_gap_threshold = (
            float(mecc_identity_gap_threshold)
            if mecc_identity_gap_threshold is not None
            else 1.0
        )
        self.u_secondary_min_frac = self._as_fraction(u_secondary_min_frac)
        self.u_secondary_min_bp = int(u_secondary_min_bp)
        self.u_contig_gap_bp = int(u_contig_gap_bp)
        # Adaptive thresholds for improved large-scale performance
        self.u_secondary_max_ratio = self._as_fraction(u_secondary_max_ratio)
        self.u_high_coverage_threshold = self._as_fraction(u_high_coverage_threshold)
        self.u_high_mapq_threshold = int(u_high_mapq_threshold)
        self.theta_locus = self._as_fraction(theta_locus)
        self.pos_tol_bp = int(pos_tol_bp)
        self.stats: dict[str, int] = {}

        # Setup logger
        self.logger = logger or get_logger(self.__class__.__name__)

    def _u_has_significant_secondary_mapping(
        self,
        all_alignments: pd.DataFrame,
        *,
        best_chr: str,
        best_start0: int,
        best_end0: int,
        cons_len: int,
        q_style: str,
        u_cov: float = 0.0,
        best_locus_mapq: Optional[int] = None,
    ) -> bool:
        """Return True if non-contiguous / multi-chr mappings are strong enough to veto Uecc.

        The veto logic uses adaptive thresholds:
        1. Relative threshold: secondary_cov / u_cov must not exceed u_secondary_max_ratio
        2. MAPQ weighting: high MAPQ increases tolerance (more confidence in primary mapping)
        3. High coverage tolerance: u_cov >= threshold gets relaxed thresholds
        """
        if all_alignments.empty or cons_len <= 0:
            return False

        gap_bp = int(self.u_contig_gap_bp)
        chr_col = ColumnStandard.CHR
        start_col = ColumnStandard.START0
        end_col = ColumnStandard.END0

        chrs = all_alignments.get(chr_col)
        if chrs is None:
            return False

        start0 = pd.to_numeric(all_alignments.get(start_col), errors="coerce")
        end0 = pd.to_numeric(all_alignments.get(end_col), errors="coerce")
        if start0 is None or end0 is None:
            return False

        best_chr_str = str(best_chr)
        other_chr_mask = chrs.astype(str) != best_chr_str
        same_chr_mask = ~other_chr_mask

        far_mask = same_chr_mask & (
            (end0 < (int(best_start0) - gap_bp)) | (start0 > (int(best_end0) + gap_bp))
        )
        secondary_df = all_alignments[other_chr_mask | far_mask]
        if secondary_df.empty:
            return False

        secondary_cov = self._coverage_fraction_for_alignments(secondary_df, cons_len, q_style)
        secondary_bp = secondary_cov * float(cons_len)

        # Compute adaptive thresholds
        effective_min_frac = float(self.u_secondary_min_frac)
        effective_min_bp = float(self.u_secondary_min_bp)

        # Improvement 1: High coverage tolerance
        # When primary coverage is very high, we trust it more and require stronger secondary
        # evidence to veto
        if u_cov >= float(self.u_high_coverage_threshold):
            effective_min_frac *= 2.0  # Increase from 1% to 2%
            effective_min_bp *= 2.0    # Increase from 50bp to 100bp

        # Improvement 2: MAPQ weighting
        # High MAPQ indicates unique mapping - increase tolerance for secondary
        if best_locus_mapq is not None and best_locus_mapq >= int(self.u_high_mapq_threshold):
            effective_min_frac *= 1.5
            effective_min_bp *= 1.5

        # Check absolute thresholds first
        if secondary_cov >= effective_min_frac or secondary_bp >= effective_min_bp:
            # Improvement 3: Relative threshold check
            # Even if absolute threshold is exceeded, allow if relative ratio is small
            if u_cov > 0:
                secondary_ratio = secondary_cov / u_cov
                if secondary_ratio < float(self.u_secondary_max_ratio):
                    # Secondary is small relative to primary - don't veto
                    return False
            return True

        return False

    @staticmethod
    def _as_fraction(value: float | None) -> float:
        """Accept either a fraction (0-1) or a percent (0-100) and return a fraction."""
        if value is None:
            return 0.0
        try:
            val = float(value)
        except (TypeError, ValueError):
            return 0.0
        if val > 1.0:
            return val / 100.0
        return val

    @staticmethod
    def _clamp01(value: float) -> float:
        try:
            x = float(value)
        except (TypeError, ValueError):
            return 0.0
        if x < 0.0:
            return 0.0
        if x > 1.0:
            return 1.0
        return x

    @classmethod
    def _norm_mapq(cls, mapq: float) -> float:
        """Normalize minimap2 MAPQ into [0,1]."""
        try:
            val = float(mapq)
        except (TypeError, ValueError):
            return 0.0
        cap = float(cls.MAPQ_MAX) if float(cls.MAPQ_MAX) > 0 else 60.0
        return cls._clamp01(val / cap)

    @classmethod
    def _norm_identity(cls, identity_pct: float) -> float:
        """Normalize identity (%) into [0,1] with a conservative floor."""
        try:
            val = float(identity_pct)
        except (TypeError, ValueError):
            return 0.0
        # Map [90,100] -> [0,1] (clamped). Below 90 is treated as very weak evidence.
        return cls._clamp01((val - 90.0) / 10.0)

    @staticmethod
    def _geom_mean(values: Iterable[float]) -> float:
        vals = [float(v) for v in values]
        if not vals:
            return 0.0
        prod = 1.0
        for v in vals:
            if v <= 0.0:
                return 0.0
            prod *= v
        return float(prod ** (1.0 / float(len(vals))))

    def _ensure_alignment_columns(self, df: pd.DataFrame, source: str) -> pd.DataFrame:
        """Ensure alignment DataFrame has the expected columns."""
        expected_cols = len(self.ALIGNMENT_COLUMNS)
        expected_cols_with_mapq = len(self.ALIGNMENT_COLUMNS_WITH_MAPQ)
        if df.columns.tolist() in (self.ALIGNMENT_COLUMNS, self.ALIGNMENT_COLUMNS_WITH_MAPQ):
            return df.copy()

        if len(df.columns) == 0:
            df = df.copy()
            df.columns = self.ALIGNMENT_COLUMNS
            return df

        if len(df.columns) not in (expected_cols, expected_cols_with_mapq):
            message = (
                f"Alignment results from {source} should have {expected_cols} columns "
                f"(or {expected_cols_with_mapq} with mapq), got {len(df.columns)}"
            )
            self.logger.error(message)
            raise ValueError(message)

        df = df.copy()
        df.columns = (
            self.ALIGNMENT_COLUMNS
            if len(df.columns) == expected_cols
            else self.ALIGNMENT_COLUMNS_WITH_MAPQ
        )
        return df

    def _preprocess_alignment_df(self, df: pd.DataFrame, source: str) -> pd.DataFrame:
        """Common preprocessing for alignment results."""
        df = self._ensure_alignment_columns(df, source)
        if df.empty:
            return df

        # Cast key numeric columns early (robust to minimap2/blast adapters)
        for col in ["identity", "alignment_length", "q_start", "q_end", "s_start", "s_end"]:
            df[col] = pd.to_numeric(df[col], errors="coerce")
        if "mapq" in df.columns:
            df["mapq"] = pd.to_numeric(df["mapq"], errors="coerce").fillna(0).astype(int)
        else:
            df["mapq"] = 0

        # Process strand information
        df["strand"] = df["sstrand"].apply(lambda x: "+" if x == "plus" else "-")
        neg_strand_mask = df["strand"] == "-"
        df.loc[neg_strand_mask, ["s_start", "s_end"]] = df.loc[
            neg_strand_mask, ["s_end", "s_start"]
        ].values

        # Parse query IDs with proper type conversion and standardize columns
        try:
            split_cols = df["query_id"].astype(str).str.split("|", expand=True)
            if split_cols.shape[1] < 4:
                message = (
                    f"Alignment query_id in {source} should have at least 4 '|' "
                    f"fields, got {split_cols.shape[1]}"
                )
                self.logger.error(message)
                raise ValueError(message)

            df[ColumnStandard.READS] = split_cols[0].astype(str)
            df[ColumnStandard.LENGTH] = (
                pd.to_numeric(split_cols[2], errors="coerce").fillna(0).astype(int)
            )
            df[ColumnStandard.COPY_NUMBER] = (
                pd.to_numeric(split_cols[3], errors="coerce").fillna(0).astype(float)
            )

            # Validate length to avoid division by zero later
            invalid_length = df[ColumnStandard.LENGTH] <= 0
            if invalid_length.any():
                self.logger.warning(
                    f"Found {invalid_length.sum()} entries with invalid length (<=0)"
                )
                # Filter out invalid entries
                df = df[~invalid_length].copy()

        except Exception as e:
            self.logger.error(f"Failed to parse query IDs: {e}")
            raise

        # Standardize coordinate columns
        df[ColumnStandard.CHR] = df["subject_id"]
        # Alignment inputs (BLAST/minimap2_align) use 1-based, strand-oriented coordinates.
        # Convert to 0-based, half-open coordinates for internal use.
        s_start_num = pd.to_numeric(df["s_start"], errors="coerce")
        s_end_num = pd.to_numeric(df["s_end"], errors="coerce")
        invalid_coords = s_start_num.isna() | s_end_num.isna()
        if invalid_coords.any():
            self.logger.warning(
                "Found %d entries with invalid subject coordinates; dropping",
                int(invalid_coords.sum()),
            )
            df = df[~invalid_coords].copy()
            s_start_num = s_start_num[~invalid_coords]
            s_end_num = s_end_num[~invalid_coords]

        df[ColumnStandard.START0] = (s_start_num - 1).astype(int)
        df[ColumnStandard.END0] = s_end_num.astype(int)
        df[ColumnStandard.STRAND] = df["strand"]

        # Calculate derived columns with safe division
        df["Rlength"] = df[ColumnStandard.END0] - df[ColumnStandard.START0]
        df["gap_Length"] = df[ColumnStandard.LENGTH] - df["Rlength"]

        # Safe division for gap_percentage
        df[ColumnStandard.GAP_PERCENTAGE] = 0.0
        valid_mask = df[ColumnStandard.LENGTH] > 0
        df.loc[valid_mask, ColumnStandard.GAP_PERCENTAGE] = (
            (df.loc[valid_mask, "gap_Length"].abs() / df.loc[valid_mask, ColumnStandard.LENGTH])
            * 100
        ).round(2)

        # Store statistics
        self.stats["total_alignments"] = len(df)
        self.stats["total_queries"] = df["query_id"].nunique()

        return df

    @staticmethod
    def _infer_query_coords_style(df: pd.DataFrame) -> str:
        """Infer whether q_start/q_end look 0-based half-open or 1-based inclusive (BLAST-like)."""
        if df.empty:
            return "blast"
        required = {"q_start", "q_end", "alignment_length"}
        if not required.issubset(df.columns):
            return "blast"

        sample = df[list(required)].dropna()
        if sample.empty:
            return "blast"

        qs = pd.to_numeric(sample["q_start"], errors="coerce")
        qe = pd.to_numeric(sample["q_end"], errors="coerce")
        alen = pd.to_numeric(sample["alignment_length"], errors="coerce")
        sample = pd.DataFrame({"qs": qs, "qe": qe, "alen": alen}).dropna()
        if sample.empty:
            return "blast"

        diff0 = (sample["qe"] - sample["qs"]).abs()
        diff1 = (sample["qe"] - sample["qs"] + 1).abs()
        err0 = (diff0 - sample["alen"]).abs().median()
        err1 = (diff1 - sample["alen"]).abs().median()
        return "0based" if err0 + 1e-6 < err1 else "blast"

    @staticmethod
    def _query_interval_0based_half_open(
        q_start: float, q_end: float, style: str
    ) -> tuple[int, int]:
        """Convert q_start/q_end to 0-based, half-open interval."""
        qs = int(q_start)
        qe = int(q_end)
        if style != "0based":
            qs -= 1
        if qe < qs:
            qs, qe = qe, qs
        return qs, qe

    @classmethod
    def _project_query_interval_to_ring(
        cls, q_start: float, q_end: float, cons_len: int, style: str
    ) -> list[tuple[int, int]]:
        """Project a query interval on doubled sequence onto ring coordinates [0, L).

        For a circular DNA of length L, queries are often mapped to a doubled
        reference (length 2L) to handle wrap-around. This function projects
        query intervals back onto the ring [0, L).

        Args:
            q_start: Query start position
            q_end: Query end position
            cons_len: Consensus length (ring circumference)
            style: Coordinate style ("0based" or "blast" for 1-based)

        Returns:
            List of (start, end) tuples representing ring intervals.
            May return multiple intervals if the projection wraps around.
        """
        L = int(cons_len)
        if L <= 0:
            return []

        qs, qe = cls._query_interval_0based_half_open(q_start, q_end, style)
        query_len = qe - qs

        # Validate query length
        if query_len <= 0:
            return []

        # If query covers full ring or more, return complete coverage
        if query_len >= L:
            return [(0, L)]

        # Project to ring coordinates using modulo
        u = qs % L
        v = qe % L

        # Handle edge case: if both project to same point and query is non-empty,
        # it means exactly one full wrap (shouldn't happen given query_len < L check above,
        # but handle defensively)
        if u == v:
            # This can happen if query_len == L exactly (handled above) or numerical edge case
            return [(0, L)] if query_len > 0 else []

        # Non-wrapping case: u < v means interval doesn't cross the origin
        if u < v:
            return [(u, v)]

        # Wrapping case: interval crosses the ring origin (position 0)
        # Split into two segments: [u, L) and [0, v)
        return [(u, L), (0, v)]

    @staticmethod
    def _union_len(segments: Iterable[tuple[int, int]]) -> int:
        """Compute union length of half-open segments."""
        segs = [(int(s), int(e)) for s, e in segments if e > s]
        if not segs:
            return 0
        segs.sort(key=lambda x: x[0])
        start, end = segs[0]
        total = 0
        for s, e in segs[1:]:
            if s <= end:
                end = max(end, e)
            else:
                total += end - start
                start, end = s, e
        total += end - start
        return int(total)

    @classmethod
    def _coverage_fraction_for_alignments(
        cls, df: pd.DataFrame, cons_len: int, style: str
    ) -> float:
        """Compute ring coverage fraction from q_start/q_end intervals."""
        L = int(cons_len)
        if df.empty or L <= 0:
            return 0.0
        segments: list[tuple[int, int]] = []
        q_start_series = df.get("q_start")
        q_end_series = df.get("q_end")
        if q_start_series is None or q_end_series is None:
            return 0.0
        # Avoid iterrows (slow); iterate over numpy arrays.
        q_start_vals = q_start_series.to_numpy()
        q_end_vals = q_end_series.to_numpy()
        for q_start, q_end in zip(q_start_vals, q_end_vals):
            if pd.isna(q_start) or pd.isna(q_end):
                continue
            segments.extend(cls._project_query_interval_to_ring(q_start, q_end, L, style))
        covered = cls._union_len(segments)
        return float(covered) / float(L) if L > 0 else 0.0

    def _cluster_loci(self, group: pd.DataFrame) -> dict[int, list[int]]:
        """Cluster alignments into loci for a single query."""
        if group.empty:
            return {}

        indices = list(group.index)
        parent: dict[int, int] = {idx: idx for idx in indices}

        def find(x: int) -> int:
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(a: int, b: int) -> None:
            ra = find(a)
            rb = find(b)
            if ra != rb:
                parent[rb] = ra

        for _, sub in group.groupby([ColumnStandard.CHR, ColumnStandard.STRAND], sort=False):
            sub_idx = list(sub.index)
            if len(sub_idx) <= 1:
                continue

            starts = sub[ColumnStandard.START0].to_numpy(dtype="int64", copy=False)
            ends = sub[ColumnStandard.END0].to_numpy(dtype="int64", copy=False)
            pos_tol = int(self.pos_tol_bp)
            theta = float(self.theta_locus)

            for i in range(len(sub_idx)):
                start_a = int(starts[i])
                end_a = int(ends[i])
                for j in range(i + 1, len(sub_idx)):
                    start_b = int(starts[j])
                    end_b = int(ends[j])

                    if abs(start_a - start_b) <= pos_tol and abs(end_a - end_b) <= pos_tol:
                        union(sub_idx[i], sub_idx[j])
                        continue

                    ov = max(0, min(end_a, end_b) - max(start_a, start_b))
                    min_len = min(max(0, end_a - start_a), max(0, end_b - start_b))
                    if min_len <= 0:
                        continue
                    if (ov / float(min_len)) >= theta:
                        union(sub_idx[i], sub_idx[j])

        clusters: dict[int, list[int]] = {}
        for idx in indices:
            root = find(idx)
            clusters.setdefault(root, []).append(idx)

        relabeled: dict[int, list[int]] = {}
        for locus_id, root in enumerate(sorted(clusters.keys(), key=lambda x: str(x))):
            relabeled[locus_id] = clusters[root]
        return relabeled

    def read_alignment_results(self, alignment_file: Path) -> pd.DataFrame:
        """
        Read and preprocess alignment results

        Args:
            alignment_file: Path to alignment results file

        Returns:
            Preprocessed DataFrame
        """
        self.logger.info(f"Reading alignment results from {alignment_file}")

        try:
            df = pd.read_csv(alignment_file, sep="\t", header=None)
        except pd.errors.EmptyDataError:
            self.logger.error(f"Alignment file is empty: {alignment_file}")
            raise
        except Exception as e:
            self.logger.error(f"Failed to read alignment file: {e}")
            raise

        df = self._preprocess_alignment_df(df, str(alignment_file))
        self.logger.info(
            f"Loaded {len(df):,} alignments from {df['query_id'].nunique():,} queries"
        )

        return df

    def classify_alignment_results(
        self, alignment_df: pd.DataFrame, output_prefix: Path
    ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        Classify eccDNA from pre-loaded alignment DataFrame

        Args:
            alignment_df: Raw alignment results DataFrame (without headers)
            output_prefix: Output file prefix path (used for saving files)

        Returns:
            Tuple of (uecc_df, mecc_df, unclassified_df)
        """
        self.logger.info("=" * 60)
        self.logger.info("Starting U/Mecc classification from DataFrame")
        self.logger.info("=" * 60)

        df = self._preprocess_alignment_df(alignment_df, "DataFrame")

        self.logger.info(
            f"Processed {len(df):,} alignments from {df['query_id'].nunique():,} queries"
        )

        # Check if we have valid data after preprocessing
        if df.empty:
            self.logger.warning("No valid data after preprocessing")
            empty_df = pd.DataFrame()
            # Save empty files
            self._save_classification_files(empty_df, empty_df, empty_df, output_prefix)
            return empty_df, empty_df, empty_df

        # Step 1: Classify Uecc and Mecc based on high-quality alignments
        uecc_df, mecc_df, classified_queries = self.classify_uecc_mecc(df)

        # Step 2: Extract ALL alignments from unclassified queries
        unclassified_df = self.extract_unclassified(df, classified_queries)

        # Format outputs
        uecc_df, mecc_df = self.format_outputs(uecc_df, mecc_df)

        # Save files
        self._save_classification_files(uecc_df, mecc_df, unclassified_df, output_prefix)

        # Final statistics
        uecc_count = uecc_df["query_id"].nunique() if not uecc_df.empty else 0
        mecc_count = mecc_df["query_id"].nunique() if not mecc_df.empty else 0
        unclassified_count = (
            unclassified_df["query_id"].nunique() if not unclassified_df.empty else 0
        )

        # Quality breakdown for unclassified
        unclass_hq_queries = 0
        unclass_lq_only_queries = 0
        if not unclassified_df.empty:
            for query_id, group in unclassified_df.groupby("query_id"):
                if (group["quality_category"] == "High_quality").any():
                    unclass_hq_queries += 1
                else:
                    unclass_lq_only_queries += 1

        self.logger.info("=" * 60)
        self.logger.info("Classification Summary:")
        self.logger.info(f"  Uecc: {uecc_count:,} queries")
        self.logger.info(f"  Mecc: {mecc_count:,} queries")
        self.logger.info(f"  Unclassified: {unclassified_count:,} queries")
        if unclassified_count > 0:
            self.logger.info(f"    - With high-quality alignments: {unclass_hq_queries:,}")
            self.logger.info(f"    - Only low-quality alignments: {unclass_lq_only_queries:,}")
        self.logger.info("=" * 60)

        return uecc_df, mecc_df, unclassified_df

    def _save_classification_files(
        self,
        uecc_df: pd.DataFrame,
        mecc_df: pd.DataFrame,
        unclassified_df: pd.DataFrame,
        output_prefix: Path,
    ) -> None:
        """Save classification results to CSV files."""
        uecc_out = output_prefix.parent / f"{output_prefix.name}.uecc.csv"
        mecc_out = output_prefix.parent / f"{output_prefix.name}.mecc.csv"
        unclass_out = output_prefix.parent / f"{output_prefix.name}.unclassified.csv"

        if uecc_df is not None and not uecc_df.empty:
            uecc_df.to_csv(uecc_out, index=False)
            self.logger.info(f"Saved Uecc results to {uecc_out}")
        else:
            # Write empty table for downstream processes
            pd.DataFrame().to_csv(uecc_out, index=False)
            self.logger.info(f"Saved empty Uecc results to {uecc_out}")

        if mecc_df is not None and not mecc_df.empty:
            mecc_df.to_csv(mecc_out, index=False)
            self.logger.info(f"Saved Mecc results to {mecc_out}")
        else:
            pd.DataFrame().to_csv(mecc_out, index=False)
            self.logger.info(f"Saved empty Mecc results to {mecc_out}")

        if unclassified_df is not None and not unclassified_df.empty:
            unclassified_df.to_csv(unclass_out, index=False)
            self.logger.info(f"Saved Unclassified results to {unclass_out}")
        else:
            pd.DataFrame().to_csv(unclass_out, index=False)
            self.logger.info(f"Saved empty Unclassified results to {unclass_out}")

    def classify_uecc_mecc(self, df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, set[str]]:
        """
        Classify Uecc and Mecc using a ring-coverage model on high-quality alignments.

        Args:
            df: Full DataFrame

        Returns:
            Tuple of (uecc_df, mecc_df, classified_query_ids)
        """
        self.logger.info("=" * 60)
        self.logger.info("Step 1: Classify using ring coverage on high-quality alignments")

        theta_u = float(self.theta_u)
        theta_m = float(self.theta_m)
        theta_u2_max = float(self.theta_u2_max)
        if theta_u <= 0.0 and theta_m <= 0.0:
            self.logger.warning("theta_u and theta_m are <= 0; nothing will be classified")
            return pd.DataFrame(), pd.DataFrame(), set()

        # Filter by quality (keep original df for downstream unclassified extraction)
        # Note: we still consult all alignments (including low-quality partial hits) when deciding
        # whether a U candidate has significant non-contiguous / multi-chr evidence.
        all_groups = df.groupby("query_id", sort=False).groups
        high_quality = df[df[ColumnStandard.GAP_PERCENTAGE] <= self.gap_threshold].copy()
        self.logger.info(f"High-quality alignments: {len(high_quality):,} / {len(df):,}")

        if high_quality.empty:
            self.logger.warning("No high-quality alignments found")
            return pd.DataFrame(), pd.DataFrame(), set()

        q_style = self._infer_query_coords_style(high_quality)

        uecc_rows: list[dict[str, Any]] = []
        mecc_rows: list[dict[str, Any]] = []
        classified_queries: set[str] = set()

        for query_id, group in high_quality.groupby("query_id", sort=False):
            cons_len = int(group[ColumnStandard.LENGTH].iloc[0])
            if cons_len <= 0:
                continue

            loci = self._cluster_loci(group)
            if not loci:
                continue

            locus_cov: dict[int, float] = {}
            for locus_id, idxs in loci.items():
                locus_df = group.loc[idxs]
                locus_cov[locus_id] = self._coverage_fraction_for_alignments(
                    locus_df, cons_len, q_style
                )

            cov_sorted = sorted(locus_cov.items(), key=lambda kv: kv[1], reverse=True)
            u_cov = cov_sorted[0][1] if cov_sorted else 0.0
            u_cov_2nd = cov_sorted[1][1] if len(cov_sorted) > 1 else 0.0
            best_lid = cov_sorted[0][0] if cov_sorted else None
            second_lid = cov_sorted[1][0] if len(cov_sorted) > 1 else None

            def locus_evidence(lid: Optional[int]) -> tuple[int, int, float, float]:
                if lid is None:
                    return 0, 0, 0.0, 0.0
                idxs = loci.get(lid)
                if not idxs:
                    return 0, 0, 0.0, 0.0
                locus_df_local = group.loc[idxs]
                if locus_df_local.empty:
                    return 0, 0, 0.0, 0.0
                # Columns were normalized to numeric in `_preprocess_alignment_df`.
                try:
                    mapq_best_local = int(locus_df_local["mapq"].max())
                    mapq_min_local = int(locus_df_local["mapq"].min())
                except (TypeError, ValueError):
                    mapq_best_local = 0
                    mapq_min_local = 0

                try:
                    id_best_local = float(locus_df_local["identity"].max())
                    id_min_local = float(locus_df_local["identity"].min())
                except (TypeError, ValueError):
                    id_best_local = 0.0
                    id_min_local = 0.0

                return int(mapq_best_local), int(mapq_min_local), float(id_best_local), float(id_min_local)

            full_loci = [lid for lid, cov in locus_cov.items() if cov >= theta_m]
            m_count = len(full_loci)

            if m_count >= 2:
                best_mapq_best, best_mapq_min, best_id_best, best_id_min = locus_evidence(best_lid)
                second_mapq_best, second_mapq_min, second_id_best, second_id_min = locus_evidence(
                    second_lid
                )
                mapq_best = max(best_mapq_best, second_mapq_best)
                mapq_min = (
                    min(best_mapq_min, second_mapq_min)
                    if second_lid is not None
                    else int(best_mapq_min)
                )
                identity_best = max(best_id_best, second_id_best)
                identity_min = (
                    min(best_id_min, second_id_min)
                    if second_lid is not None
                    else float(best_id_min)
                )
                low_mapq = mapq_best < int(self.MAPQ_LOW_THRESHOLD)
                low_identity = identity_best < float(self.IDENTITY_LOW_THRESHOLD)

                # MAPQ ambiguity check for Mecc: if ALL loci have low MAPQ,
                # this might be a Uecc from a repetitive region, not a true Mecc.
                # Skip classification and leave as unclassified.
                if self.mapq_m_ambiguous_threshold > 0:
                    all_loci_low_mapq = True
                    for lid_check in full_loci:
                        lid_mapq, _, _, _ = locus_evidence(lid_check)
                        if lid_mapq >= self.mapq_m_ambiguous_threshold:
                            all_loci_low_mapq = False
                            break
                    if all_loci_low_mapq:
                        self.stats["mecc_vetoed_low_mapq"] = self.stats.get(
                            "mecc_vetoed_low_mapq", 0
                        ) + 1
                        continue  # Skip Mecc classification, leave as unclassified

                # Identity gap check for Mecc: if one locus has significantly higher
                # identity than all others, this is likely a Uecc from a repetitive
                # region rather than a true Mecc (which should have similar identity
                # across all loci).
                if self.mecc_identity_gap_threshold > 0:
                    loci_max_identities = []
                    for lid_check in full_loci:
                        _, _, lid_id_best, _ = locus_evidence(lid_check)
                        loci_max_identities.append(float(lid_id_best))
                    if len(loci_max_identities) >= 2:
                        loci_max_identities.sort(reverse=True)
                        identity_gap = loci_max_identities[0] - loci_max_identities[1]
                        # Use > instead of >= to allow exact threshold matches
                        # (e.g., 1.0% gap with 1.0% threshold is allowed)
                        if identity_gap > self.mecc_identity_gap_threshold:
                            self.stats["mecc_vetoed_identity_gap"] = self.stats.get(
                                "mecc_vetoed_identity_gap", 0
                            ) + 1
                            continue  # Skip Mecc classification, leave as unclassified

                conf = self._geom_mean(
                    [
                        self._norm_mapq(mapq_best),
                        self._norm_identity(identity_best),
                        self._clamp01(u_cov),
                        self._clamp01(u_cov_2nd),
                    ]
                )
                for lid in full_loci:
                    idxs = loci[lid]
                    locus_df = group.loc[idxs]
                    if locus_df.empty:
                        continue
                    try:
                        best_idx = locus_df["alignment_length"].idxmax()
                        rep_row = locus_df.loc[best_idx]
                    except (ValueError, KeyError):
                        if locus_df.empty:
                            continue
                        rep_row = locus_df.iloc[0]

                    row = dict(rep_row.to_dict())
                    row["eccdna_type"] = "Mecc"
                    row["classification_reason"] = ">=2 loci with full ring coverage"
                    row["locus_id"] = int(lid)
                    row["locus_cov"] = round(float(locus_cov[lid]), 6)
                    row["U_cov"] = round(float(u_cov), 6)
                    row["U_cov_2nd"] = round(float(u_cov_2nd), 6)
                    row["M_count"] = int(m_count)
                    row[ColumnStandard.MAPQ_BEST] = int(mapq_best)
                    row[ColumnStandard.MAPQ_MIN] = int(mapq_min)
                    row[ColumnStandard.IDENTITY_BEST] = round(float(identity_best), 3)
                    row[ColumnStandard.IDENTITY_MIN] = round(float(identity_min), 3)
                    row[ColumnStandard.QUERY_COV_BEST] = round(float(u_cov), 6)
                    row[ColumnStandard.QUERY_COV_2ND] = round(float(u_cov_2nd), 6)
                    row[ColumnStandard.CONFIDENCE_SCORE] = round(float(conf), 6)
                    row[ColumnStandard.LOW_MAPQ] = bool(low_mapq)
                    row[ColumnStandard.LOW_IDENTITY] = bool(low_identity)
                    mecc_rows.append(row)
                classified_queries.add(query_id)
                continue

            if u_cov >= theta_u and u_cov_2nd <= theta_u2_max:
                # Uecc requires one strong locus explanation and no substantial 2nd locus.
                lid = best_lid if best_lid is not None else cov_sorted[0][0]
                locus_df = group.loc[loci[lid]]
                # Safety check: skip if locus_df is empty
                if locus_df.empty:
                    continue
                # Always compute best locus MAPQ for adaptive thresholds
                try:
                    best_locus_mapq = int(locus_df["mapq"].max())
                except (TypeError, ValueError):
                    best_locus_mapq = None
                # Check MAPQ minimum if configured
                if self.mapq_u_min > 0:
                    if best_locus_mapq is None or best_locus_mapq < int(self.mapq_u_min):
                        self.stats["uecc_vetoed_low_mapq"] = self.stats.get(
                            "uecc_vetoed_low_mapq", 0
                        ) + 1
                        continue
                best_chr = str(locus_df[ColumnStandard.CHR].iloc[0])
                try:
                    best_start0 = int(locus_df[ColumnStandard.START0].min())
                    best_end0 = int(locus_df[ColumnStandard.END0].max())
                except (TypeError, ValueError):
                    best_start0 = 0
                    best_end0 = 0
                all_idx = all_groups.get(query_id)
                all_alignments = df.loc[all_idx] if all_idx is not None else group
                if self._u_has_significant_secondary_mapping(
                    all_alignments,
                    best_chr=best_chr,
                    best_start0=best_start0,
                    best_end0=best_end0,
                    cons_len=cons_len,
                    q_style=q_style,
                    u_cov=u_cov,
                    best_locus_mapq=best_locus_mapq,
                ):
                    continue
                try:
                    best_idx = locus_df["alignment_length"].idxmax()
                    rep_row = locus_df.loc[best_idx]
                except (ValueError, KeyError):
                    if locus_df.empty:
                        continue
                    rep_row = locus_df.iloc[0]

                mapq_best, mapq_min, identity_best, identity_min = locus_evidence(lid)
                low_mapq = mapq_best < int(self.MAPQ_LOW_THRESHOLD)
                low_identity = identity_best < float(self.IDENTITY_LOW_THRESHOLD)
                denom = float(max(theta_u2_max, 0.05))
                uniq = self._clamp01(1.0 - (float(u_cov_2nd) / denom)) if denom > 0 else 1.0
                conf = self._geom_mean(
                    [
                        self._norm_mapq(mapq_best),
                        self._norm_identity(identity_best),
                        self._clamp01(u_cov),
                        uniq,
                    ]
                )
                row = dict(rep_row.to_dict())
                row["eccdna_type"] = "Uecc"
                row["classification_reason"] = "single locus with full ring coverage"
                row["locus_id"] = int(lid)
                row["locus_cov"] = round(float(locus_cov[lid]), 6)
                row["U_cov"] = round(float(u_cov), 6)
                row["U_cov_2nd"] = round(float(u_cov_2nd), 6)
                row["M_count"] = int(m_count)
                row[ColumnStandard.MAPQ_BEST] = int(mapq_best)
                row[ColumnStandard.MAPQ_MIN] = int(mapq_min)
                row[ColumnStandard.IDENTITY_BEST] = round(float(identity_best), 3)
                row[ColumnStandard.IDENTITY_MIN] = round(float(identity_min), 3)
                row[ColumnStandard.QUERY_COV_BEST] = round(float(u_cov), 6)
                row[ColumnStandard.QUERY_COV_2ND] = round(float(u_cov_2nd), 6)
                row[ColumnStandard.CONFIDENCE_SCORE] = round(float(conf), 6)
                row[ColumnStandard.LOW_MAPQ] = bool(low_mapq)
                row[ColumnStandard.LOW_IDENTITY] = bool(low_identity)
                uecc_rows.append(row)
                classified_queries.add(query_id)
                continue

        uecc_df = pd.DataFrame(uecc_rows) if uecc_rows else pd.DataFrame()
        mecc_df = pd.DataFrame(mecc_rows) if mecc_rows else pd.DataFrame()

        self.logger.info("Uecc: %d queries", uecc_df["query_id"].nunique() if not uecc_df.empty else 0)
        self.logger.info("Mecc: %d queries", mecc_df["query_id"].nunique() if not mecc_df.empty else 0)

        # Log filtering statistics
        if self.stats.get("mecc_vetoed_low_mapq", 0) > 0:
            self.logger.info(
                "Mecc vetoed due to low MAPQ (all loci < %d): %d queries",
                self.mapq_m_ambiguous_threshold,
                self.stats["mecc_vetoed_low_mapq"],
            )
        if self.stats.get("mecc_vetoed_identity_gap", 0) > 0:
            self.logger.info(
                "Mecc vetoed due to identity gap (>= %.1f): %d queries",
                self.mecc_identity_gap_threshold,
                self.stats["mecc_vetoed_identity_gap"],
            )
        if self.stats.get("uecc_vetoed_low_mapq", 0) > 0:
            self.logger.info(
                "Uecc vetoed due to low MAPQ (< %d): %d queries",
                self.mapq_u_min,
                self.stats["uecc_vetoed_low_mapq"],
            )

        return uecc_df, mecc_df, classified_queries

    def _process_overlaps_for_query(self, group: pd.DataFrame) -> pd.DataFrame:
        """
        Process overlaps for a single query's alignments.

        Args:
            group: DataFrame for single query

        Returns:
            Processed DataFrame with overlaps resolved
        """
        kept_alignments = []

        for _, chr_group in group.groupby(ColumnStandard.CHR):
            components = self._find_overlaps_sweepline(chr_group)
            for component in components:
                component_df = chr_group.loc[list(component)]
                if len(component_df) == 1:
                    kept_alignments.append(component_df)
                else:
                    # Select best from each overlap component
                    best_idx = component_df[ColumnStandard.GAP_PERCENTAGE].idxmin()
                    kept_alignments.append(chr_group.loc[[best_idx]])

        if kept_alignments:
            return pd.concat(kept_alignments, ignore_index=False)
        return pd.DataFrame()

    def _find_overlaps_sweepline(self, group: pd.DataFrame) -> list[set[int]]:
        """
        Find overlap components using sweep-line algorithm.

        Args:
            group: DataFrame with alignments on same chromosome

        Returns:
            List of sets, each containing indices in the same overlap component
        """
        if len(group) <= 1:
            return [set(group.index)] if len(group) else []

        parent = {idx: idx for idx in group.index}

        def find(x: int) -> int:
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(a: int, b: int) -> None:
            root_a = find(a)
            root_b = find(b)
            if root_a != root_b:
                parent[root_b] = root_a

        events = []
        for idx, row in group.iterrows():
            start = row[ColumnStandard.START0]
            end = row[ColumnStandard.END0]
            if end < start:
                start, end = end, start
            events.append((start, 0, idx))  # 0 for start event
            events.append((end, 1, idx))  # 1 for end event

        # Sort by position, then by event type (starts before ends)
        events.sort(key=lambda x: (x[0], x[1]))

        active_intervals: set[int] = set()
        for _, event_type, idx in events:
            if event_type == 0:  # start event
                for active_idx in active_intervals:
                    union(idx, active_idx)
                active_intervals.add(idx)
            else:  # end event
                active_intervals.discard(idx)

        components: dict[int, set[int]] = {}
        for idx in group.index:
            root = find(idx)
            components.setdefault(root, set()).add(idx)

        return list(components.values())

    def _is_full_length_repeat(self, group: pd.DataFrame) -> bool:
        """
        Check if alignments represent full-length repeats (Mecc)

        Args:
            group: DataFrame with multiple alignments for one query

        Returns:
            True if at least 2 full-length copies
        """
        # Filter out entries with zero or invalid length
        valid_group = group[group[ColumnStandard.LENGTH] > 0]

        if valid_group.empty or len(valid_group) < 2:
            return False

        # Calculate coverage for each alignment
        coverages = (valid_group["Rlength"] / valid_group[ColumnStandard.LENGTH]) * 100
        full_length_count = int((coverages >= self.min_full_length_coverage).sum())

        return full_length_count >= 2

    def _passes_identity_gap(self, group: pd.DataFrame) -> bool:
        """Check if identity gap between top two alignments is within Mecc threshold."""
        if self.max_identity_gap_for_mecc is None:
            return True

        if group.empty or "identity" not in group.columns:
            return False

        identities = pd.to_numeric(group["identity"], errors="coerce").dropna()
        if len(identities) < 2:
            return False

        identities = identities.sort_values(ascending=False)
        gap = float(identities.iloc[0] - identities.iloc[1])
        return gap <= float(self.max_identity_gap_for_mecc)

    def extract_unclassified(
        self, df_original: pd.DataFrame, classified_queries: set[str]
    ) -> pd.DataFrame:
        """
        Extract all alignments from unclassified queries

        Args:
            df_original: Original complete DataFrame
            classified_queries: Set of query IDs that were classified

        Returns:
            DataFrame with all alignments from unclassified queries
        """
        self.logger.info("=" * 60)
        self.logger.info("Step 2: Extract unclassified queries")

        all_queries = set(df_original["query_id"].unique())
        unclassified_queries = all_queries - classified_queries

        self.logger.info(f"Total queries: {len(all_queries):,}")
        self.logger.info(f"Classified queries: {len(classified_queries):,}")
        self.logger.info(f"Unclassified queries: {len(unclassified_queries):,}")

        # Extract ALL alignments for unclassified queries
        unclassified_df = df_original[df_original["query_id"].isin(unclassified_queries)].copy()

        if not unclassified_df.empty:
            # Add reason column for downstream analysis
            unclassified_df["unclass_reason"] = "Not_classified"

            # Add quality indicator
            unclassified_df["quality_category"] = unclassified_df[ColumnStandard.GAP_PERCENTAGE].apply(
                lambda x: "High_quality" if x <= self.gap_threshold else "Low_quality"
            )

            # Statistics
            hq_count = (unclassified_df["quality_category"] == "High_quality").sum()
            lq_count = (unclassified_df["quality_category"] == "Low_quality").sum()

            self.logger.info(f"Unclassified alignments: {len(unclassified_df):,}")
            self.logger.info(f"  - High quality (Gap <= {self.gap_threshold}%): {hq_count:,}")
            self.logger.info(f"  - Low quality (Gap > {self.gap_threshold}%): {lq_count:,}")

        return unclassified_df

    def format_outputs(
        self, uecc_df: pd.DataFrame, mecc_df: pd.DataFrame
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Format output DataFrames with appropriate columns

        Args:
            uecc_df: Uecc DataFrame
            mecc_df: Mecc DataFrame

        Returns:
            Tuple of formatted (uecc_df, mecc_df)
        """
        # Format Uecc (clean output)
        if not uecc_df.empty:
            uecc_formatted = uecc_df.copy()
            uecc_formatted[ColumnStandard.MATCH_DEGREE] = 100 - uecc_formatted[ColumnStandard.GAP_PERCENTAGE]
        else:
            uecc_formatted = uecc_df

        # Format Mecc (clean output)
        if not mecc_df.empty:
            mecc_formatted = mecc_df.copy()
            mecc_formatted[ColumnStandard.MATCH_DEGREE] = 100 - mecc_formatted[ColumnStandard.GAP_PERCENTAGE]
        else:
            mecc_formatted = mecc_df

        return uecc_formatted, mecc_formatted

    def run(self, alignment_file: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        Run complete U/Mecc classification

        Args:
            alignment_file: Path to alignment results

        Returns:
            Tuple of (uecc_df, mecc_df, unclassified_df)
        """
        self.logger.info("=" * 60)
        self.logger.info("Starting U/Mecc classification")
        self.logger.info("=" * 60)

        # Read and preprocess alignment results
        df_original = self.read_alignment_results(alignment_file)

        # Check if we have valid data after preprocessing
        if df_original.empty:
            self.logger.warning("No valid data after preprocessing")
            return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

        # Step 1: Classify Uecc and Mecc based on high-quality alignments
        uecc_df, mecc_df, classified_queries = self.classify_uecc_mecc(df_original)

        # Step 2: Extract ALL alignments from unclassified queries
        unclassified_df = self.extract_unclassified(df_original, classified_queries)

        # Format outputs
        uecc_df, mecc_df = self.format_outputs(uecc_df, mecc_df)

        # Final statistics
        uecc_count = uecc_df["query_id"].nunique() if not uecc_df.empty else 0
        mecc_count = mecc_df["query_id"].nunique() if not mecc_df.empty else 0
        unclassified_count = (
            unclassified_df["query_id"].nunique() if not unclassified_df.empty else 0
        )

        # Quality breakdown for unclassified
        unclass_hq_queries = 0
        unclass_lq_only_queries = 0
        if not unclassified_df.empty:
            for query_id, group in unclassified_df.groupby("query_id"):
                if (group["quality_category"] == "High_quality").any():
                    unclass_hq_queries += 1
                else:
                    unclass_lq_only_queries += 1

        self.logger.info("=" * 60)
        self.logger.info("Classification Summary:")
        self.logger.info(f"  Uecc: {uecc_count:,} queries")
        self.logger.info(f"  Mecc: {mecc_count:,} queries")
        self.logger.info(f"  Unclassified: {unclassified_count:,} queries")
        if unclassified_count > 0:
            self.logger.info(f"    - With high-quality alignments: {unclass_hq_queries:,}")
            self.logger.info(f"    - Only low-quality alignments: {unclass_lq_only_queries:,}")
        self.logger.info("=" * 60)

        return uecc_df, mecc_df, unclassified_df


def _parse_args():
    """Parse CLI arguments for direct script execution"""
    import argparse

    parser = argparse.ArgumentParser(
        description="U/Mecc-Classify - Classify eccDNA from alignment results"
    )
    parser.add_argument("-i", "--input", required=True, help="Input alignment result TSV file")
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output prefix (generates .uecc.csv/.mecc.csv/.unclassified.csv)",
    )
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
        help="Log level (default: INFO)",
    )
    return parser.parse_args()


def main():
    from pathlib import Path

    args = _parse_args()

    import logging as _logging

    _logging.getLogger().setLevel(getattr(_logging, args.log_level))

    input_path = Path(args.input)
    output_prefix = Path(args.output)
    # Remove suffix if user provided filename with extension
    if output_prefix.suffix:
        output_prefix = output_prefix.with_suffix("")

    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    clf = UMeccClassifier()
    uecc_df, mecc_df, unclassified_df = clf.run(input_path)

    uecc_out = output_prefix.with_suffix(".uecc.csv")
    mecc_out = output_prefix.with_suffix(".mecc.csv")
    un_out = output_prefix.with_suffix(".unclassified.csv")

    if uecc_df is not None and not uecc_df.empty:
        uecc_df.to_csv(uecc_out, index=False)
    else:
        # Write empty table for downstream
        pd.DataFrame().to_csv(uecc_out, index=False)

    if mecc_df is not None and not mecc_df.empty:
        mecc_df.to_csv(mecc_out, index=False)
    else:
        pd.DataFrame().to_csv(mecc_out, index=False)

    if unclassified_df is not None and not unclassified_df.empty:
        unclassified_df.to_csv(un_out, index=False)
    else:
        pd.DataFrame().to_csv(un_out, index=False)


if __name__ == "__main__":
    main()
