"""
Enhanced Cecc-Build - CECC DNA Analysis Pipeline Module with Overlap Filtering

Added post-processing step to filter out queries with overlapping genomic segments.
"""

from __future__ import annotations

import csv
import logging
from circleseeker.utils.logging import get_logger
from circleseeker.utils.column_standards import ColumnStandard
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any, Set
import numpy as np
import pandas as pd

# Optional acceleration
try:
    from numba import jit
    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False


class CeccBuild:
    """Circular eccDNA analysis using segment rotation and path detection."""
    
    # Base required columns (independent of legacy vs standard names)
    BASE_REQUIRED_COLS = [
        "query_id", "subject_id", "q_start", "q_end", "s_start", "s_end",
        "strand", "alignment_length",
    ]
    STANDARD_EXTRA_COLS = ["reads", "length", "copy_number"]

    CANDIDATE_DELIMITERS = [",", "\t", ";", "|"]

    def __init__(self, logger: Optional[logging.Logger] = None):
        """Initialize Cecc-Build analyzer."""
        self.logger = logger or get_logger(self.__class__.__name__)
        
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
    
    def detect_genomic_overlaps_sweepline(self, segments: pd.DataFrame) -> bool:
        """
        Use sweep-line algorithm to detect overlapping genomic segments.
        
        Args:
            segments: DataFrame with chr, start, end columns
            
        Returns:
            True if any overlaps detected, False otherwise
        """
        if len(segments) <= 1:
            return False
        
        # Group by chromosome to check overlaps within each chromosome
        for chr_name, chr_group in segments.groupby(ColumnStandard.CHR):
            if len(chr_group) <= 1:
                continue
                
            # Create events for sweep line algorithm
            events = []
            for idx, row in chr_group.iterrows():
                start = int(row[ColumnStandard.START0])
                end = int(row[ColumnStandard.END0])
                events.append((start, 0, idx))  # 0 for start event
                events.append((end, 1, idx))    # 1 for end event
            
            # Sort by position, then by event type (starts before ends)
            events.sort(key=lambda x: (x[0], x[1]))
            
            # Sweep line to detect overlaps
            active_intervals = set()
            
            for pos, event_type, idx in events:
                if event_type == 0:  # start event
                    if active_intervals:
                        # Found overlap - there's already an active interval
                        return True
                    active_intervals.add(idx)
                else:  # end event
                    active_intervals.discard(idx)
        
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
        for query_id, query_group in df_labeled.groupby('query_id'):
            # Extract genomic coordinates for this query
            segments = query_group[[ColumnStandard.CHR, ColumnStandard.START0, ColumnStandard.END0]].copy()
            
            # Check for overlaps using sweep-line algorithm
            has_overlap = self.detect_genomic_overlaps_sweepline(segments)
            
            if has_overlap:
                overlapping_queries.add(query_id)
                # Log example of overlap for debugging
                if len(overlapping_queries) <= 5:  # Only log first few examples
                    self.logger.debug(f"Query {query_id} has overlapping segments:")
                    for _, row in segments.iterrows():
                        self.logger.debug(f"  {row[ColumnStandard.CHR]} {row[ColumnStandard.START0]} {row[ColumnStandard.END0]}")
            else:
                valid_queries.add(query_id)
        
        # Filter the DataFrame
        df_filtered = df_labeled[df_labeled['query_id'].isin(valid_queries)].copy()
        
        # Log statistics
        total_queries = len(df_labeled['query_id'].unique())
        overlapping_count = len(overlapping_queries)
        retained_count = len(valid_queries)
        
        self.logger.info(f"Overlap filtering results:")
        self.logger.info(f"  - Total queries before filtering: {total_queries}")
        self.logger.info(f"  - Queries with overlapping segments: {overlapping_count}")
        self.logger.info(f"  - Queries retained: {retained_count}")
        self.logger.info(f"  - Segments before filtering: {len(df_labeled)}")
        self.logger.info(f"  - Segments after filtering: {len(df_filtered)}")
        
        return df_filtered
    
    def ensure_required_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """Validate and prepare input DataFrame with column standardization (standard names only)."""
        # Strictly require standard columns
        required = self.BASE_REQUIRED_COLS + self.STANDARD_EXTRA_COLS
        missing = [c for c in required if c not in df.columns]
        if missing:
            raise ValueError(f"Input CSV missing required columns: {missing}")

        # Cast numeric columns
        numeric_cols = ["q_start", "q_end", "s_start", "s_end", "alignment_length", "length", "copy_number"]
        for col in numeric_cols:
            df[col] = pd.to_numeric(df[col], errors="coerce")

        # Drop rows with invalid numeric data
        df = df.dropna(subset=["q_start", "q_end", "s_start", "s_end", "alignment_length", "length"]).copy()

        # Standardize column names to internal standards
        df[ColumnStandard.CHR] = df['subject_id']
        df[ColumnStandard.START0] = df['s_start']
        df[ColumnStandard.END0] = df['s_end']
        df[ColumnStandard.STRAND] = df['strand']
        df[ColumnStandard.READS] = df['reads']
        df[ColumnStandard.LENGTH] = df['length']
        df[ColumnStandard.COPY_NUMBER] = df['copy_number']

        return df

    def _read_input_dataframe(self, input_csv: Path) -> Tuple[pd.DataFrame, str]:
        """Load the input CSV while detecting an appropriate delimiter."""
        sample = ""
        try:
            with input_csv.open("r", encoding="utf-8", errors="replace") as handle:
                sample = handle.read(65536)
        except OSError as exc:
            raise ValueError(f"Unable to read input file {input_csv}: {exc}") from exc

        candidate_seps: List[Optional[str]] = []
        if sample:
            try:
                dialect = csv.Sniffer().sniff(sample, delimiters=self.CANDIDATE_DELIMITERS)
                candidate_seps.append(dialect.delimiter)
            except csv.Error:
                pass

        for sep in self.CANDIDATE_DELIMITERS:
            if sep not in candidate_seps:
                candidate_seps.append(sep)

        candidate_seps.append(None)  # auto-detect fallback

        errors: List[str] = []
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
            "Unable to parse input CSV with supported delimiters. "
            + "; ".join(errors)
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
        self,
        df_rot: pd.DataFrame,
        overlap_threshold: float,
        min_segments: int
    ) -> pd.DataFrame:
        """
        Vectorized filtering by segment count and overlap criteria.
        This matches the exact algorithm from step5_trapeze.py.
        """
        keep_idx = []
        
        for qid, g in df_rot.groupby("query_id", sort=False):
            # Filter by minimum segments (using <= as in original)
            if len(g) <= min_segments:
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
        self.logger.info(f"After filtering: {len(keep_idx)} queries retained, {len(filtered)} segments")
        return filtered
    
    def position_match(self, a: pd.Series, b: pd.Series, tol: int) -> bool:
        """Check if two segments have matching chromosomal positions within tolerance."""
        if a[ColumnStandard.CHR] != b[ColumnStandard.CHR]:
            return False
        if str(a[ColumnStandard.STRAND]) != str(b[ColumnStandard.STRAND]):
            return False
        return (abs(float(a[ColumnStandard.START0]) - float(b[ColumnStandard.START0])) <= tol and
                abs(float(a[ColumnStandard.END0]) - float(b[ColumnStandard.END0])) <= tol)
    
    def analyze_gaps(self, path_df: pd.DataFrame) -> Tuple[float, float, int, bool]:
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
        pos_tol: int
    ) -> Optional[Dict[str, Any]]:
        """Find circular path in a group of segments."""
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
                        "path": path[:-1],
                        "closing_at": idx,
                        "cum_len": cum_len - float(cur["alignment_length"]),
                        "mat_degree": round(((cum_len - float(cur["alignment_length"])) / cons_len) * 100, 2),
                    }
                
                # Check closure at next segment (if still continuous)
                if idx < len(g) - 1:
                    nxt = g.iloc[idx + 1]
                    next_gap = float(nxt["q_start"]) - float(cur["q_end"])
                    if abs(next_gap) <= edge_tol and self.position_match(first, nxt, pos_tol):
                        return {
                            "path": path,
                            "closing_at": idx + 1,
                            "cum_len": cum_len,
                            "mat_degree": round((cum_len / cons_len) * 100, 2),
                        }
        
        return None
    
    def detect_circles(
        self,
        df_filt: pd.DataFrame,
        edge_tol: int,
        pos_tol: int
    ) -> pd.DataFrame:
        """Detect circular patterns in filtered data."""
        rows: List[Dict] = []
        
        for qid, g in df_filt.groupby("query_id", sort=False):
            res = self.find_circular(g, edge_tol, pos_tol)
            if not res:
                continue
            
            path_df = g.iloc[res["path"]].reset_index(drop=True)
            chroms = path_df[ColumnStandard.CHR].unique()
            eclass = "Cecc-InterChr" if len(chroms) > 1 else "Cecc-IntraChr"
            max_gap, avg_gap, _, _ = self.analyze_gaps(path_df)
            
            base = {
                "query_id": qid,
                "reads": g.iloc[0][ColumnStandard.READS],
                "eccdna_type": "Cecc",
                "CeccClass": eclass,
                "length": g.iloc[0][ColumnStandard.LENGTH],
                "copy_number": g.iloc[0][ColumnStandard.COPY_NUMBER],
                "num_segments": len(path_df),
                "cumulative_length": res["cum_len"],
                "match_degree": res["mat_degree"],
                "max_gap": max_gap,
                "avg_gap": round(avg_gap, 2),
                "chromosomes": ",".join(chroms.astype(str)),
            }
            
            for i in range(len(path_df)):
                seg = path_df.iloc[i]
                row = dict(base)
                row.update({
                    "segment_in_circle": i + 1,
                    ColumnStandard.CHR: seg[ColumnStandard.CHR],
                    ColumnStandard.START0: seg[ColumnStandard.START0],
                    ColumnStandard.END0: seg[ColumnStandard.END0],
                    ColumnStandard.STRAND: seg[ColumnStandard.STRAND],
                    "q_start": seg["q_start"],
                    "q_end": seg["q_end"],
                    "alignment_length": seg["alignment_length"],
                })
                rows.append(row)
        
        result_df = pd.DataFrame(rows)
        if not result_df.empty:
            self.logger.info(f"Circular patterns detected: {result_df['query_id'].nunique()} queries")
        
        return result_df
    
    def label_roles(self, df_circ: pd.DataFrame) -> pd.DataFrame:
        """Label segments with roles: head, middle, tail."""
        if df_circ.empty:
            return df_circ.copy()
        
        # Stable sort
        df = df_circ.sort_values(["query_id", "segment_in_circle"], kind="mergesort").copy()
        
        parts: List[pd.DataFrame] = []
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
        edge_tolerance: int = 20,
        position_tolerance: int = 50,
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
            Discard queries with segments <= this value
        edge_tolerance : int
            Max allowed |q-gap| between adjacent segments
        position_tolerance : int
            Genome position tolerance for closure (bp)
        (Final overlap filtering is always enabled)
        
        Returns
        -------
        pd.DataFrame
            DataFrame containing labeled circular patterns
        """
        self.logger.info("=" * 60)
        self.logger.info("Cecc-Build - Circular eccDNA Analysis")
        self.logger.info("=" * 60)
        self.logger.info(f"Parameters:")
        self.logger.info(f"  - Overlap threshold: {overlap_threshold}")
        self.logger.info(f"  - Min segments: {min_segments}")
        self.logger.info(f"  - Edge tolerance: {edge_tolerance} bp")
        self.logger.info(f"  - Position tolerance: {position_tolerance} bp")
        self.logger.info(f"  - Final overlap filtering: Enabled")
        self.logger.info(f"  - Numba acceleration: {'Enabled' if NUMBA_AVAILABLE else 'Disabled'}")
        
        # Read and validate input
        self.logger.info(f"\nReading input from: {input_csv}")
        df, sep_label = self._read_input_dataframe(input_csv)
        self.logger.info(f"Detected delimiter: {sep_label}")
        self.logger.info(f"Loaded {len(df)} rows")

        # Ensure required columns and clean data
        df = self.ensure_required_columns(df)
        self.logger.info(f"Validated columns; {df['query_id'].nunique()} unique queries")
        
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
        df_circles = self.detect_circles(df_filt, edge_tolerance, position_tolerance)
        
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
                    df_labeled["reads"] = df_labeled["query_id"].astype(str).str.split("|", n=1).str[0]
        
        df_labeled.to_csv(output_csv, index=False)
        
        # Generate summary statistics
        if not df_labeled.empty:
            unique_queries = df_labeled['query_id'].nunique()
            intra_chr = len(df_labeled[df_labeled["CeccClass"] == "Cecc-IntraChr"])
            inter_chr = len(df_labeled[df_labeled["CeccClass"] == "Cecc-InterChr"])
            
            # Calculate average metrics
            query_stats = df_labeled.groupby('query_id').first()
            avg_mat_degree = query_stats['match_degree'].mean()
            avg_segments = query_stats['num_segments'].mean()
            
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
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input CSV with alignment info"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output CSV with circular results"
    )
    parser.add_argument("--overlap-threshold", type=float, default=0.95, help="Overlap threshold (default: 0.95)")
    parser.add_argument("--min-segments", type=int, default=2, help="Min segments (default: 2)")
    parser.add_argument("--edge-tolerance", type=int, default=20, help="Edge tolerance in bp (default: 20)")
    parser.add_argument("--position-tolerance", type=int, default=50, help="Position tolerance in bp (default: 50)")
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
        help="Log level (default: INFO)"
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
    )

if __name__ == "__main__":
    main()
