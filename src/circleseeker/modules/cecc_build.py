"""
Enhanced Cecc-Build - CECC DNA Analysis Pipeline Module with Overlap Filtering

Added post-processing step to filter out queries with overlapping genomic segments.
"""

from __future__ import annotations

import csv
import sys
import logging
from circleseeker.utils.logging import get_logger
from circleseeker.utils.column_standards import ColumnStandard
from pathlib import Path
from typing import Optional, Any, Iterable
import numpy as np
import pandas as pd

# Increase CSV field size limit to handle large sequence fields
csv.field_size_limit(sys.maxsize)

# Optional acceleration
try:
    from numba import jit

    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False


class CeccBuild:
    """Circular eccDNA analysis using segment rotation and path detection."""

    # Heuristic thresholds used for evidence flags and score scaling.
    # These do NOT affect Cecc detection logic (which remains thresholded by match_degree).
    MAPQ_MAX = 60
    MAPQ_LOW_THRESHOLD = 20
    IDENTITY_LOW_THRESHOLD = 95.0

    # Base required columns (independent of legacy vs standard names)
    BASE_REQUIRED_COLS = [
        "query_id",
        "subject_id",
        "q_start",
        "q_end",
        "s_start",
        "s_end",
        "strand",
        "alignment_length",
    ]
    STANDARD_EXTRA_COLS = ["reads", "length", "copy_number"]

    CANDIDATE_DELIMITERS = [",", "\t", ";", "|"]

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
        try:
            val = float(mapq)
        except (TypeError, ValueError):
            return 0.0
        cap = float(cls.MAPQ_MAX) if float(cls.MAPQ_MAX) > 0 else 60.0
        return cls._clamp01(val / cap)

    @classmethod
    def _norm_identity(cls, identity_pct: float) -> float:
        try:
            val = float(identity_pct)
        except (TypeError, ValueError):
            return 0.0
        return cls._clamp01((val - 90.0) / 10.0)

    @staticmethod
    def _geom_mean(values: list[float]) -> float:
        if not values:
            return 0.0
        prod = 1.0
        for v in values:
            if v <= 0.0:
                return 0.0
            prod *= float(v)
        return float(prod ** (1.0 / float(len(values))))

    def __init__(self, logger: Optional[logging.Logger] = None):
        """Initialize Cecc-Build analyzer."""
        self.logger = logger or get_logger(self.__class__.__name__)
        # Treat near-identical genomic intervals as the same locus in the final overlap filter.
        # This prevents false negatives caused by 1-2bp boundary jitter in alignments.
        self.locus_overlap_threshold = 0.95

        # Set up the optimized overlap function based on numba availability
        if NUMBA_AVAILABLE:

            @jit(nopython=True)
            def _overlap_numba_impl(qs, qe, thr=0.8):
                n = len(qs)
                max_r = 0.0
                for i in range(n):
                    for j in range(i + 1, n):
                        os = qs[i] if qs[i] > qs[j] else qs[j]
                        oe = qe[i] if qe[i] < qe[j] else qe[j]
                        if os < oe:
                            ov = oe - os
                            l1 = qe[i] - qs[i]
                            l2 = qe[j] - qs[j]
                            m = l1 if l1 < l2 else l2
                            if m > 0:
                                r = ov / m
                                if r > max_r:
                                    max_r = r
                                if r >= thr:
                                    return True, max_r
                return False, max_r

            self._overlap_numba = _overlap_numba_impl
        else:

            def _overlap_numba_impl(qs, qe, thr=0.8):
                n = len(qs)
                max_r = 0.0
                for i in range(n):
                    for j in range(i + 1, n):
                        os = max(qs[i], qs[j])
                        oe = min(qe[i], qe[j])
                        if os < oe:
                            ov = oe - os
                            l1 = qe[i] - qs[i]
                            l2 = qe[j] - qs[j]
                            m = min(l1, l2)
                            if m > 0:
                                r = ov / m
                                max_r = max(max_r, r)
                                if r >= thr:
                                    return True, max_r
                return False, max_r

            self._overlap_numba = _overlap_numba_impl

    def select_non_overlapping_alignments(
        self, group: pd.DataFrame, overlap_tolerance: int = 10
    ) -> pd.DataFrame:
        """
        Select non-overlapping alignments from a group using greedy algorithm.

        When BLAST produces multiple overlapping alignments for the same query region,
        this method selects the best non-overlapping subset by prioritizing longer
        alignments.

        Args:
            group: DataFrame containing alignments for a single query
            overlap_tolerance: Maximum allowed overlap in bp (default: 10)

        Returns:
            DataFrame with selected non-overlapping alignments
        """
        if len(group) <= 1:
            return group.copy()

        # Sort by alignment_length descending to prioritize longer alignments
        sorted_group = group.sort_values("alignment_length", ascending=False)

        selected_indices = []
        selected_intervals = []  # [(q_start, q_end), ...]

        for idx, row in sorted_group.iterrows():
            q_start = float(row["q_start"])
            q_end = float(row["q_end"])

            # Check overlap with already selected alignments
            is_overlapping = False
            for sel_start, sel_end in selected_intervals:
                overlap_start = max(q_start, sel_start)
                overlap_end = min(q_end, sel_end)
                overlap_len = max(0, overlap_end - overlap_start)

                if overlap_len > overlap_tolerance:
                    is_overlapping = True
                    break

            if not is_overlapping:
                selected_indices.append(idx)
                selected_intervals.append((q_start, q_end))

        return group.loc[selected_indices].copy()

    def preprocess_overlapping_alignments(
        self, df: pd.DataFrame, overlap_tolerance: int = 10
    ) -> pd.DataFrame:
        """
        Preprocess all queries to select non-overlapping alignment subsets.

        This fixes the bug where overlapping BLAST alignments cause the circular
        detection algorithm to fail due to negative gaps between adjacent records.

        Args:
            df: DataFrame with all alignment records
            overlap_tolerance: Maximum allowed overlap in bp (default: 10)

        Returns:
            Preprocessed DataFrame with non-overlapping alignments per query
        """
        original_count = len(df)
        original_queries = df["query_id"].nunique()

        processed_parts = []

        for query_id, group in df.groupby("query_id", sort=False):
            processed = self.select_non_overlapping_alignments(group, overlap_tolerance)
            processed_parts.append(processed)

        if not processed_parts:
            return df.iloc[0:0].copy()

        result = pd.concat(processed_parts, ignore_index=True)

        self.logger.info(
            f"Overlap preprocessing: {original_count} -> {len(result)} alignments "
            f"({original_queries} queries)"
        )

        return result

    def detect_genomic_overlaps_sweepline(self, segments: pd.DataFrame) -> bool:
        """
        Detect overlapping genomic segments on the reference.

        Notes
        -----
        We treat intervals as 0-based, half-open. Exact duplicate intervals are allowed
        (they often arise when a locus is visited twice on the doubled query).
        Near-duplicate intervals are also allowed when their reciprocal overlap
        (overlap / min(len_a, len_b)) is >= ``self.locus_overlap_threshold``.
        Other partial overlaps are considered conflicting and will be filtered.

        The reciprocal overlap is computed correctly as:
            overlap_start = max(cur_start, start)
            overlap_end = min(cur_end, end)
            overlap_len = max(0, overlap_end - overlap_start)

        Args:
            segments: DataFrame with chr, start, end columns

        Returns:
            True if any partial overlaps detected (excluding near-identical), False otherwise
        """
        if len(segments) <= 1:
            return False

        group_cols = [ColumnStandard.CHR]
        if ColumnStandard.STRAND in segments.columns:
            group_cols.append(ColumnStandard.STRAND)

        for _, chr_group in segments.groupby(group_cols, sort=False):
            if len(chr_group) <= 1:
                continue
            df_int = chr_group[[ColumnStandard.START0, ColumnStandard.END0]].copy()
            df_int[ColumnStandard.START0] = pd.to_numeric(
                df_int[ColumnStandard.START0], errors="coerce"
            )
            df_int[ColumnStandard.END0] = pd.to_numeric(df_int[ColumnStandard.END0], errors="coerce")
            df_int = df_int.dropna().astype(int)
            if len(df_int) <= 1:
                continue
            df_int = df_int.sort_values([ColumnStandard.START0, ColumnStandard.END0], kind="mergesort")

            cur_start = int(df_int.iloc[0][ColumnStandard.START0])
            cur_end = int(df_int.iloc[0][ColumnStandard.END0])

            for _, row in df_int.iloc[1:].iterrows():
                start = int(row[ColumnStandard.START0])
                end = int(row[ColumnStandard.END0])

                # No overlap with current interval
                if start >= cur_end:
                    cur_start, cur_end = start, end
                    continue

                # Calculate correct overlap using max of starts and min of ends
                overlap_start = max(cur_start, start)
                overlap_end = min(cur_end, end)
                overlap_len = max(0, overlap_end - overlap_start)

                # Calculate interval lengths (avoid division by zero)
                len_a = max(1, cur_end - cur_start)
                len_b = max(1, end - start)

                # Calculate reciprocal overlap for both intervals
                rec_ov_a = overlap_len / float(len_a)
                rec_ov_b = overlap_len / float(len_b)

                # Check if near-identical (both reciprocal overlaps exceed threshold)
                threshold = float(self.locus_overlap_threshold)
                if rec_ov_a >= threshold and rec_ov_b >= threshold:
                    # Merge into the current locus envelope
                    cur_start = min(cur_start, start)
                    cur_end = max(cur_end, end)
                    continue

                # Partial overlap detected - this is a conflict
                return True

        return False

    def filter_overlapping_queries(self, df_labeled: pd.DataFrame) -> pd.DataFrame:
        """
        Filter out queries that have overlapping genomic segments.

        Args:
            df_labeled: DataFrame with labeled circular patterns

        Returns:
            DataFrame with overlapping queries removed
        """
        if df_labeled.empty:
            return df_labeled

        self.logger.info("Filtering queries with overlapping genomic segments")

        valid_queries = set()
        overlapping_queries = set()

        # Check each query for genomic overlaps
        for query_id, query_group in df_labeled.groupby("query_id"):
            # Extract genomic coordinates for this query
            cols = [ColumnStandard.CHR, ColumnStandard.START0, ColumnStandard.END0]
            if ColumnStandard.STRAND in query_group.columns:
                cols.append(ColumnStandard.STRAND)
            segments = query_group[cols].copy()

            # Check for overlaps using sweep-line algorithm
            has_overlap = self.detect_genomic_overlaps_sweepline(segments)

            if has_overlap:
                overlapping_queries.add(query_id)
                # Log example of overlap for debugging
                if len(overlapping_queries) <= 5:  # Only log first few examples
                    self.logger.debug(f"Query {query_id} has overlapping segments:")
                    for _, row in segments.iterrows():
                        chrom = row[ColumnStandard.CHR]
                        start = row[ColumnStandard.START0]
                        end = row[ColumnStandard.END0]
                        self.logger.debug(f"  {chrom} {start} {end}")
            else:
                valid_queries.add(query_id)

        # Filter the DataFrame
        df_filtered = df_labeled[df_labeled["query_id"].isin(valid_queries)].copy()

        # Log statistics
        total_queries = len(df_labeled["query_id"].unique())
        overlapping_count = len(overlapping_queries)
        retained_count = len(valid_queries)

        self.logger.info("Overlap filtering results:")
        self.logger.info(f"  - Total queries before filtering: {total_queries}")
        self.logger.info(f"  - Queries with overlapping segments: {overlapping_count}")
        self.logger.info(f"  - Queries retained: {retained_count}")
        self.logger.info(f"  - Segments before filtering: {len(df_labeled)}")
        self.logger.info(f"  - Segments after filtering: {len(df_filtered)}")

        return df_filtered

    def ensure_required_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """Validate and prepare input DataFrame with column standardization."""
        # Strictly require standard columns
        required = self.BASE_REQUIRED_COLS + self.STANDARD_EXTRA_COLS
        missing = [c for c in required if c not in df.columns]
        if missing:
            raise ValueError(f"Input CSV missing required columns: {missing}")

        # Cast numeric columns
        numeric_cols = [
            "q_start",
            "q_end",
            "s_start",
            "s_end",
            "alignment_length",
            "length",
            "copy_number",
        ]
        for col in numeric_cols:
            df[col] = pd.to_numeric(df[col], errors="coerce")

        # Drop rows with invalid numeric data
        df = df.dropna(
            subset=["q_start", "q_end", "s_start", "s_end", "alignment_length", "length"]
        ).copy()

        # Standardize column names to internal standards
        df[ColumnStandard.CHR] = df["subject_id"]
        # Alignment inputs (BLAST/minimap2_align) use 1-based, strand-oriented coordinates.
        # Convert to 0-based, half-open coordinates for internal use.
        s_start = pd.to_numeric(df["s_start"], errors="coerce")
        s_end = pd.to_numeric(df["s_end"], errors="coerce")
        s_min = np.minimum(s_start, s_end)
        s_max = np.maximum(s_start, s_end)
        df[ColumnStandard.START0] = (s_min - 1).astype(int)
        df[ColumnStandard.END0] = s_max.astype(int)
        df[ColumnStandard.STRAND] = df["strand"]
        df[ColumnStandard.READS] = df["reads"]
        df[ColumnStandard.LENGTH] = df["length"]
        df[ColumnStandard.COPY_NUMBER] = df["copy_number"]

        return df

    def _read_input_dataframe(self, input_csv: Path) -> tuple[pd.DataFrame, str]:
        """Load the input CSV while detecting an appropriate delimiter."""
        sample = ""
        try:
            with input_csv.open("r", encoding="utf-8", errors="replace") as handle:
                sample = handle.read(65536)
        except OSError as exc:
            raise ValueError(f"Unable to read input file {input_csv}: {exc}") from exc

        candidate_seps: list[Optional[str]] = []
        if sample:
            try:
                dialect = csv.Sniffer().sniff(sample, delimiters="".join(self.CANDIDATE_DELIMITERS))
                candidate_seps.append(dialect.delimiter)
            except csv.Error:
                pass

        for candidate_sep in self.CANDIDATE_DELIMITERS:
            if candidate_sep not in candidate_seps:
                candidate_seps.append(candidate_sep)

        candidate_seps.append(None)  # auto-detect fallback

        errors: list[str] = []
        for sep in candidate_seps:
            try:
                if sep is None:
                    df = pd.read_csv(input_csv, sep=None, engine="python")
                    sep_label = "auto-detected"
                else:
                    df = pd.read_csv(input_csv, sep=sep)
                    sep_label = repr(sep)
            except Exception as exc:  # pragma: no cover - pandas-specific errors vary
                errors.append(f"sep={repr(sep)} failed: {exc}")
                continue

            # Validate required columns (base + standard extras)
            required = self.BASE_REQUIRED_COLS + self.STANDARD_EXTRA_COLS
            missing_required = [c for c in required if c not in df.columns]
            if missing_required:
                errors.append(
                    f"sep={repr(sep)} missing required columns: {', '.join(missing_required)}"
                )
                continue

            return df, sep_label

        raise ValueError(
            "Unable to parse input CSV with supported delimiters. " + "; ".join(errors)
        )

    def rotate_group(self, g: pd.DataFrame) -> pd.DataFrame:
        """Sort by q_start, move first row to the end, add 1-based segment_order."""
        s = g.sort_values("q_start").reset_index(drop=True)
        if len(s) > 1:
            # Move first segment to end (rotation)
            s = pd.concat([s.iloc[1:], s.iloc[:1]], ignore_index=True)
        s = s.copy()
        s["segment_order"] = np.arange(1, len(s) + 1)
        return s

    def rotate_all(self, df: pd.DataFrame) -> pd.DataFrame:
        """Rotate all groups efficiently in batch."""
        parts = [self.rotate_group(g) for _, g in df.groupby("query_id", sort=False)]
        return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()

    def vectorized_overlap_filter(
        self, df_rot: pd.DataFrame, overlap_threshold: float, min_segments: int
    ) -> pd.DataFrame:
        """
        Vectorized filtering by segment count and overlap criteria.
        This matches the exact algorithm from step5_trapeze.py.
        """
        keep_idx = []

        for qid, g in df_rot.groupby("query_id", sort=False):
            # Filter by minimum segments (Cecc must have at least 2 segments)
            if len(g) < min_segments:
                continue

            # Check for overlaps using optimized function
            qs = g["q_start"].to_numpy(np.float64)
            qe = g["q_end"].to_numpy(np.float64)
            has_overlap, _ = self._overlap_numba(qs, qe, overlap_threshold)

            if not has_overlap:
                keep_idx.append(g.index)

        if not keep_idx:
            return df_rot.iloc[0:0].copy()  # Return empty DataFrame with same structure

        filtered = df_rot.loc[np.concatenate(keep_idx)].reset_index(drop=True)
        self.logger.info(
            f"After filtering: {len(keep_idx)} queries retained, {len(filtered)} segments"
        )
        return filtered

    def position_match(self, a: pd.Series, b: pd.Series, tol: int) -> bool:
        """Check if two segments have matching chromosomal positions within tolerance."""
        if a[ColumnStandard.CHR] != b[ColumnStandard.CHR]:
            return False
        if str(a[ColumnStandard.STRAND]) != str(b[ColumnStandard.STRAND]):
            return False
        return (
            abs(float(a[ColumnStandard.START0]) - float(b[ColumnStandard.START0])) <= tol
            and abs(float(a[ColumnStandard.END0]) - float(b[ColumnStandard.END0])) <= tol
        )

    def analyze_gaps(self, path_df: pd.DataFrame) -> tuple[float, float, int, bool]:
        """Analyze gaps between segments in a path."""
        gaps = []
        for i in range(1, len(path_df)):
            prev = path_df.iloc[i - 1]
            cur = path_df.iloc[i]
            gaps.append(float(cur["q_start"]) - float(prev["q_end"]))

        if not gaps:
            return 0.0, 0.0, 0, True

        max_gap = max(abs(g) for g in gaps)
        avg_gap = sum(abs(g) for g in gaps) / len(gaps)
        return max_gap, avg_gap, len(gaps), max_gap <= 20

    def find_circular(
        self,
        group_df: pd.DataFrame,
        edge_tol: int,
        pos_tol: int,
        min_match_degree: Optional[float] = None,
    ) -> Optional[dict[str, Any]]:
        """Find circular path in a group of segments.

        Primary mode uses position-based closure (duplicate segment indicating wrap-around).
        Optional fallback mode accepts coverage-based paths when closure cannot be established.
        """
        g = group_df.sort_values("segment_order").reset_index(drop=True)
        if len(g) < 2:
            return None

        cons_len = float(g.loc[0, ColumnStandard.LENGTH])
        path = [0]
        cum_len = float(g.loc[0, "alignment_length"])

        for idx in range(1, len(g)):
            prev = g.iloc[path[-1]]
            cur = g.iloc[idx]
            gap = float(cur["q_start"]) - float(prev["q_end"])

            if abs(gap) > edge_tol:
                return None  # Break in continuity

            path.append(idx)
            cum_len += float(cur["alignment_length"])

            if cum_len >= cons_len:
                first = g.iloc[0]

                # Check closure at current segment
                if self.position_match(first, cur, pos_tol):
                    return {
                        # Keep the closing segment in the path. On doubled queries, the "duplicate"
                        # locus can still correspond to a distinct ring interval that is required
                        # to achieve full ring coverage after projection.
                        "path": path,
                        "closing_at": idx,
                        "cum_len": cum_len,
                        "mat_degree": round((cum_len / cons_len) * 100, 2),
                        "closure_found": True,
                        "closure_reason": "position_match",
                    }

                # Check closure at next segment (if still continuous)
                if idx < len(g) - 1:
                    nxt = g.iloc[idx + 1]
                    next_gap = float(nxt["q_start"]) - float(cur["q_end"])
                    if abs(next_gap) <= edge_tol and self.position_match(first, nxt, pos_tol):
                        cum_len_with_nxt = cum_len + float(nxt["alignment_length"])
                        return {
                            "path": path + [idx + 1],
                            "closing_at": idx + 1,
                            "cum_len": cum_len_with_nxt,
                            "mat_degree": round((cum_len_with_nxt / cons_len) * 100, 2),
                            "closure_found": True,
                            "closure_reason": "position_match",
                        }

        if min_match_degree is None or cons_len <= 0:
            return None

        mat_degree = (cum_len / cons_len) * 100
        if mat_degree < float(min_match_degree):
            return None

        return {
            "path": path,
            "closing_at": None,
            "cum_len": cum_len,
            "mat_degree": round(mat_degree, 2),
            "closure_found": False,
            "closure_reason": "coverage",
        }

    def detect_circles(
        self,
        df_filt: pd.DataFrame,
        edge_tol: int,
        pos_tol: int,
        min_match_degree: float = 90.0,
        max_rotations: int = 20,
    ) -> pd.DataFrame:
        """Detect circular patterns in filtered data."""
        rows: list[Dict] = []
        loci_cols = [
            ColumnStandard.CHR,
            ColumnStandard.START0,
            ColumnStandard.END0,
            ColumnStandard.STRAND,
        ]

        def loci_signature(unique_loci: pd.DataFrame) -> str:
            if unique_loci.empty:
                return ""
            df_sig = unique_loci.copy()
            df_sig = df_sig.sort_values(
                [ColumnStandard.CHR, ColumnStandard.START0, ColumnStandard.END0, ColumnStandard.STRAND],
                kind="mergesort",
            )
            parts: list[str] = []
            for _, row in df_sig.iterrows():
                parts.append(
                    f"{row[ColumnStandard.CHR]}:{int(row[ColumnStandard.START0])}-"
                    f"{int(row[ColumnStandard.END0])}:{row[ColumnStandard.STRAND]}"
                )
            return ";".join(parts)

        def infer_query_coords_style(df: pd.DataFrame) -> str:
            if df.empty:
                return "blast"
            sample = df[["q_start", "q_end", "alignment_length"]].dropna()
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

        def project_query_interval_to_ring(
            q_start: float, q_end: float, cons_len: int, style: str
        ) -> list[tuple[int, int]]:
            L = int(cons_len)
            if L <= 0:
                return []
            qs = int(q_start)
            qe = int(q_end)
            if style != "0based":
                qs -= 1
            if qe < qs:
                qs, qe = qe, qs
            if qe - qs >= L:
                return [(0, L)]
            u = qs % L
            v = qe % L
            if u < v:
                return [(u, v)]
            return [(u, L), (0, v)]

        def union_len(segments: Iterable[tuple[int, int]]) -> int:
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

        def ring_coverage_percent(path_df: pd.DataFrame, cons_len: int, style: str) -> float:
            L = int(cons_len)
            if path_df.empty or L <= 0:
                return 0.0
            segments: list[tuple[int, int]] = []
            for _, row in path_df.iterrows():
                q_start = row.get("q_start")
                q_end = row.get("q_end")
                if pd.isna(q_start) or pd.isna(q_end):
                    continue
                segments.extend(project_query_interval_to_ring(q_start, q_end, L, style))
            covered = union_len(segments)
            return (covered / float(L)) * 100.0 if L > 0 else 0.0

        for qid, g in df_filt.groupby("query_id", sort=False):
            group_sorted = g.sort_values("q_start", kind="mergesort").reset_index(drop=True)
            n_segments = len(group_sorted)
            if n_segments < 2:
                continue
            cons_len = int(group_sorted.loc[0, ColumnStandard.LENGTH])
            q_style = infer_query_coords_style(group_sorted)

            # Candidate rotations: try all for small groups; otherwise focus on likely breakpoints
            if n_segments <= max_rotations:
                starts = list(range(n_segments))
            else:
                starts = [0, 1] if n_segments > 1 else [0]
                gaps = []
                for i in range(n_segments - 1):
                    gap = abs(float(group_sorted.loc[i + 1, "q_start"]) - float(group_sorted.loc[i, "q_end"]))
                    gaps.append((gap, i + 1))
                gaps.sort(reverse=True, key=lambda x: x[0])
                for _, start in gaps:
                    if start not in starts:
                        starts.append(start)
                    if len(starts) >= max_rotations:
                        break

            best_rank: tuple[float, int, int] | None = None
            best_group: tuple[pd.DataFrame, dict[str, Any], pd.DataFrame] | None = None
            best_sig = ""
            best_degree: float | None = None

            second_rank: tuple[float, int, int] | None = None
            second_group: tuple[pd.DataFrame, dict[str, Any], pd.DataFrame] | None = None
            second_sig = ""
            second_degree: float | None = None

            for start in starts:
                rotated = (
                    group_sorted.copy()
                    if start == 0
                    else pd.concat(
                        [group_sorted.iloc[start:], group_sorted.iloc[:start]], ignore_index=True
                    )
                )
                rotated = rotated.copy()
                rotated["segment_order"] = np.arange(1, len(rotated) + 1)

                try:
                    res = self.find_circular(
                        rotated,
                        edge_tol=edge_tol,
                        pos_tol=pos_tol,
                        min_match_degree=min_match_degree,
                    )
                except TypeError:
                    # Backward-compatibility: allow patched/legacy implementations that
                    # don't accept the optional min_match_degree parameter.
                    res = self.find_circular(rotated, edge_tol=edge_tol, pos_tol=pos_tol)
                if not res:
                    continue

                path_df = rotated.iloc[res["path"]].reset_index(drop=True)
                segment_count = len(path_df)
                if segment_count < 2:
                    continue

                unique_loci = path_df[loci_cols].drop_duplicates()
                if unique_loci.shape[0] < 2:
                    continue

                cov_percent = ring_coverage_percent(path_df, cons_len, q_style)
                if min_match_degree is not None and cov_percent < float(min_match_degree):
                    continue

                sig = loci_signature(unique_loci)
                degree = float(cov_percent)
                rank = (
                    degree,
                    1 if res.get("closure_found") else 0,
                    -segment_count,
                )

                if best_rank is None or rank > best_rank:
                    # Old best may become a valid second-best candidate if distinct.
                    if best_rank is not None and best_group is not None and best_sig and best_sig != sig:
                        if second_rank is None or best_rank > second_rank:
                            second_rank = best_rank
                            second_group = best_group
                            second_sig = best_sig
                            second_degree = best_degree

                    best_rank = rank
                    best_group = (rotated, res, path_df)
                    best_sig = sig
                    best_degree = degree
                    continue

                if sig and sig != best_sig:
                    if second_rank is None or rank > second_rank:
                        second_rank = rank
                        second_group = (rotated, res, path_df)
                        second_sig = sig
                        second_degree = degree

            if not best_group:
                continue

            rotated, res, path_df = best_group
            chroms = path_df[ColumnStandard.CHR].unique()
            eclass = "Cecc-InterChr" if len(chroms) > 1 else "Cecc-IntraChr"
            max_gap, avg_gap, _, _ = self.analyze_gaps(path_df)

            mapq_series = (
                pd.to_numeric(path_df.get("mapq"), errors="coerce").dropna()
                if "mapq" in path_df.columns
                else pd.Series(dtype=float)
            )
            identity_series = (
                pd.to_numeric(path_df.get("identity"), errors="coerce").dropna()
                if "identity" in path_df.columns
                else pd.Series(dtype=float)
            )
            mapq_best = int(mapq_series.max()) if not mapq_series.empty else 0
            mapq_min = int(mapq_series.min()) if not mapq_series.empty else 0
            identity_best = float(identity_series.max()) if not identity_series.empty else 0.0
            identity_min = float(identity_series.min()) if not identity_series.empty else 0.0

            cov_best = (float(best_degree) / 100.0) if best_degree is not None else 0.0
            cov_2nd = (float(second_degree) / 100.0) if second_degree is not None else 0.0
            # If a distinct 2nd-best chain exists and is close, lower confidence.
            chain_unique = (
                1.0
                if second_degree is None
                else self._clamp01(max(0.0, cov_best - cov_2nd) / 0.05)
            )

            low_mapq = mapq_best < int(self.MAPQ_LOW_THRESHOLD)
            low_identity = identity_best < float(self.IDENTITY_LOW_THRESHOLD)
            conf = self._geom_mean(
                [
                    self._norm_mapq(mapq_best),
                    self._norm_identity(identity_best),
                    self._clamp01(cov_best),
                    self._clamp01(chain_unique),
                ]
            )

            base = {
                "query_id": qid,
                "reads": rotated.iloc[0][ColumnStandard.READS],
                "eccdna_type": "Cecc",
                "CeccClass": eclass,
                "length": rotated.iloc[0][ColumnStandard.LENGTH],
                "copy_number": rotated.iloc[0][ColumnStandard.COPY_NUMBER],
                "num_segments": len(path_df),
                "cumulative_length": res["cum_len"],
                "match_degree": best_degree,
                "match_degree_2nd": second_degree,
                "best_chain_signature": best_sig,
                "second_chain_signature": second_sig,
                ColumnStandard.MAPQ_BEST: int(mapq_best),
                ColumnStandard.MAPQ_MIN: int(mapq_min),
                ColumnStandard.IDENTITY_BEST: round(float(identity_best), 3),
                ColumnStandard.IDENTITY_MIN: round(float(identity_min), 3),
                ColumnStandard.QUERY_COV_BEST: round(float(cov_best), 6),
                ColumnStandard.QUERY_COV_2ND: round(float(cov_2nd), 6),
                ColumnStandard.CONFIDENCE_SCORE: round(float(conf), 6),
                ColumnStandard.LOW_MAPQ: bool(low_mapq),
                ColumnStandard.LOW_IDENTITY: bool(low_identity),
                "C_cov_best": round(float(cov_best), 6) if best_degree is not None else None,
                "C_cov_2nd": round(float(cov_2nd), 6) if second_degree is not None else None,
                "max_gap": max_gap,
                "avg_gap": round(avg_gap, 2),
                "chromosomes": ",".join(chroms.astype(str)),
            }

            for i, seg in path_df.iterrows():
                row = dict(base)
                row.update(
                    {
                        "segment_in_circle": int(i) + 1,
                        ColumnStandard.CHR: seg[ColumnStandard.CHR],
                        ColumnStandard.START0: seg[ColumnStandard.START0],
                        ColumnStandard.END0: seg[ColumnStandard.END0],
                        ColumnStandard.STRAND: seg[ColumnStandard.STRAND],
                        "q_start": seg["q_start"],
                        "q_end": seg["q_end"],
                        "alignment_length": seg["alignment_length"],
                    }
                )
                rows.append(row)

        result_df = pd.DataFrame(rows)
        if not result_df.empty:
            self.logger.info(
                f"Circular patterns detected: {result_df['query_id'].nunique()} queries"
            )

        return result_df

    def label_roles(self, df_circ: pd.DataFrame) -> pd.DataFrame:
        """Label segments with roles: head, middle, tail."""
        if df_circ.empty:
            return df_circ.copy()

        # Stable sort
        df = df_circ.sort_values(["query_id", "segment_in_circle"], kind="mergesort").copy()

        parts: list[pd.DataFrame] = []
        for _, g in df.groupby("query_id", sort=False):
            n = len(g)
            roles = ["middle"] * n
            if n >= 1:
                roles[0] = "head"
            if n >= 2:
                roles[-1] = "tail"

            g2 = g.copy()
            g2["segment_role"] = pd.Categorical(
                roles, categories=["head", "middle", "tail"], ordered=True
            )
            parts.append(g2)

        df = pd.concat(parts, ignore_index=True)

        # Reorder columns to place segment_role after segment_in_circle
        cols = df.columns.tolist()
        if "segment_role" in cols:
            cols.remove("segment_role")
            i = cols.index("segment_in_circle") + 1
            cols = cols[:i] + ["segment_role"] + cols[i:]
            df = df[cols]

        return df

    def run_pipeline(
        self,
        input_csv: Path,
        output_csv: Path,
        overlap_threshold: float = 0.95,
        min_segments: int = 2,
        edge_tolerance: int = 100,  # Increased from 20 to improve CeccDNA detection
        position_tolerance: int = 50,
        min_match_degree: float = 90.0,
        max_rotations: int = 20,
        locus_overlap_threshold: float = 0.95,
    ) -> pd.DataFrame:
        """
        Run the complete Trapeze circular analysis pipeline.

        Parameters
        ----------
        input_csv : Path
            Input CSV with per-segment alignments
        output_csv : Path
            Output path for labeled circular results
        overlap_threshold : float
            Overlap threshold to discard a query (relative to shorter segment)
        min_segments : int
            Minimum number of segments required (queries with fewer segments are discarded)
        edge_tolerance : int
            Max allowed |q-gap| between adjacent segments
        position_tolerance : int
            Genome position tolerance for closure (bp)

        Returns
        -------
        pd.DataFrame
            DataFrame containing labeled circular patterns
        """
        self.logger.info("=" * 60)
        self.logger.info("Cecc-Build - Circular eccDNA Analysis")
        self.logger.info("=" * 60)
        self.logger.info("Parameters:")
        self.logger.info(f"  - Overlap threshold: {overlap_threshold}")
        self.logger.info(f"  - Min segments: {min_segments}")
        self.logger.info(f"  - Edge tolerance: {edge_tolerance} bp")
        self.logger.info(f"  - Position tolerance: {position_tolerance} bp")
        self.logger.info(f"  - Min match degree (fallback): {min_match_degree}%")
        self.logger.info(f"  - Max rotations per query: {max_rotations}")
        self.logger.info("  - Final overlap filtering: Enabled")
        try:
            lot = float(locus_overlap_threshold)
        except (TypeError, ValueError):
            lot = 0.95
        self.locus_overlap_threshold = lot / 100.0 if lot > 1.0 else lot
        self.logger.info(
            "  - Final locus overlap threshold: %.3f", float(self.locus_overlap_threshold)
        )
        self.logger.info(f"  - Numba acceleration: {'Enabled' if NUMBA_AVAILABLE else 'Disabled'}")

        # Read and validate input
        self.logger.info(f"\nReading input from: {input_csv}")
        df, sep_label = self._read_input_dataframe(input_csv)
        self.logger.info(f"Detected delimiter: {sep_label}")
        self.logger.info(f"Loaded {len(df)} rows")

        # Ensure required columns and clean data
        df = self.ensure_required_columns(df)
        self.logger.info(f"Validated columns; {df['query_id'].nunique()} unique queries")

        # Preprocess: select non-overlapping alignments for each query
        # This fixes the bug where overlapping BLAST alignments cause circular
        # detection to fail due to negative gaps between adjacent records
        self.logger.info("Preprocessing: selecting non-overlapping alignments per query")
        df = self.preprocess_overlapping_alignments(df, overlap_tolerance=10)

        self.logger.info("Rotating segments by q_start")
        df_rot = self.rotate_all(df)

        if df_rot.empty:
            self.logger.warning("No data after rotation step")
            empty_result = pd.DataFrame()
            empty_result.to_csv(output_csv, index=False)
            return empty_result

        self.logger.info(f"  - Rotated {len(df_rot)} segments")

        self.logger.info("Filtering by segment count and overlaps")
        df_filt = self.vectorized_overlap_filter(df_rot, overlap_threshold, min_segments)

        if df_filt.empty:
            self.logger.warning("No data passed filtering criteria")
            empty_result = pd.DataFrame()
            empty_result.to_csv(output_csv, index=False)
            return empty_result

        self.logger.info("Detecting circular patterns")
        df_circles = self.detect_circles(
            df_filt,
            edge_tol=edge_tolerance,
            pos_tol=position_tolerance,
            min_match_degree=min_match_degree,
            max_rotations=max_rotations,
        )

        # Label roles
        if not df_circles.empty:
            self.logger.info("Labeling segment roles")
            df_labeled = self.label_roles(df_circles)
        else:
            self.logger.warning("No circular patterns detected")
            df_labeled = df_circles

        # Filter overlapping queries (mandatory)
        if not df_labeled.empty:
            df_labeled = self.filter_overlapping_queries(df_labeled)

        # Save results (clean output, no legacy columns)
        output_csv.parent.mkdir(parents=True, exist_ok=True)
        if not df_labeled.empty:
            # Ensure 'reads' is present
            if "reads" not in df_labeled.columns or df_labeled["reads"].isna().any():
                # Fallback from query_id prefix if necessary
                df_labeled["reads"] = df_labeled.get("reads", df_labeled.get("readName"))
                if "reads" not in df_labeled.columns or df_labeled["reads"].isna().any():
                    df_labeled["reads"] = (
                        df_labeled["query_id"].astype(str).str.split("|", n=1).str[0]
                    )

        df_labeled.to_csv(output_csv, index=False)

        # Generate summary statistics
        if not df_labeled.empty:
            unique_queries = df_labeled["query_id"].nunique()
            intra_chr = len(df_labeled[df_labeled["CeccClass"] == "Cecc-IntraChr"])
            inter_chr = len(df_labeled[df_labeled["CeccClass"] == "Cecc-InterChr"])

            # Calculate average metrics
            query_stats = df_labeled.groupby("query_id").first()
            avg_mat_degree = query_stats["match_degree"].mean()
            avg_segments = query_stats["num_segments"].mean()

            self.logger.info("\n" + "=" * 60)
            self.logger.info("RESULTS SUMMARY")
            self.logger.info("=" * 60)
            self.logger.info(f"Output saved to: {output_csv}")
            self.logger.info(f"Circular queries found: {unique_queries}")
            self.logger.info(f"Total segments in circles: {len(df_labeled)}")
            self.logger.info(f"  - Cecc-IntraChr: {intra_chr} segments")
            self.logger.info(f"  - Cecc-InterChr: {inter_chr} segments")
            self.logger.info(f"Average match_degree: {avg_mat_degree:.2f}%")
            self.logger.info(f"Average segments per circle: {avg_segments:.1f}")
        else:
            self.logger.info("\n" + "=" * 60)
            self.logger.info(f"No circular patterns found. Empty file saved to: {output_csv}")

        self.logger.info("=" * 60)

        return df_labeled


def _parse_args():
    """Parse CLI arguments for direct script execution"""
    import argparse

    parser = argparse.ArgumentParser(
        description="Cecc-Build - Identify circular eccDNA via fragment rotation"
    )
    parser.add_argument("-i", "--input", required=True, help="Input CSV with alignment info")
    parser.add_argument("-o", "--output", required=True, help="Output CSV with circular results")
    parser.add_argument(
        "--overlap-threshold", type=float, default=0.95, help="Overlap threshold (default: 0.95)"
    )
    parser.add_argument(
        "--locus-overlap-threshold",
        type=float,
        default=0.95,
        help="Reciprocal overlap threshold to treat two genomic intervals as the same locus "
        "in the final overlap filter (default: 0.95)",
    )
    parser.add_argument("--min-segments", type=int, default=2, help="Min segments (default: 2)")
    parser.add_argument(
        "--edge-tolerance", type=int, default=20, help="Edge tolerance in bp (default: 20)"
    )
    parser.add_argument(
        "--position-tolerance", type=int, default=50, help="Position tolerance in bp (default: 50)"
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
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    runner = CeccBuild()
    runner.run_pipeline(
        input_csv=input_path,
        output_csv=output_path,
        overlap_threshold=args.overlap_threshold,
        min_segments=args.min_segments,
        edge_tolerance=args.edge_tolerance,
        position_tolerance=args.position_tolerance,
        locus_overlap_threshold=args.locus_overlap_threshold,
    )


if __name__ == "__main__":
    main()
