"""
Trapeze - CeccDNA Analysis Pipeline Module

This module processes unclassified eccDNA candidates to identify circular patterns.
Based on the correct implementation from step5_trapeze.py, adapted for integration
into the CircleSeeker2 architecture.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
import numpy as np
import pandas as pd

# Optional acceleration
try:
    from numba import jit
    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False


class Trapeze:
    """Circular eccDNA analysis using segment rotation and path detection."""
    
    # Required columns for the algorithm
    REQUIRED_COLS = [
        "query_id", "subject_id", "q_start", "q_end", "s_start", "s_end",
        "strand", "alignment_length", "consLen", "readName", "copyNum",
    ]
    
    def __init__(self, logger: Optional[logging.Logger] = None):
        """Initialize Trapeze analyzer."""
        self.logger = logger or logging.getLogger(self.__class__.__name__)
        
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
    
    def ensure_required_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """Validate and prepare input DataFrame."""
        missing = [c for c in self.REQUIRED_COLS if c not in df.columns]
        if missing:
            raise ValueError(f"Input CSV missing required columns: {missing}")
        
        # Cast numeric columns
        numeric_cols = ["q_start", "q_end", "s_start", "s_end", "alignment_length", "consLen", "copyNum"]
        for col in numeric_cols:
            df[col] = pd.to_numeric(df[col], errors="coerce")
        
        # Drop rows with invalid numeric data
        df = df.dropna(subset=["q_start", "q_end", "s_start", "s_end", "alignment_length", "consLen"]).copy()
        
        return df
    
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
        if a["subject_id"] != b["subject_id"]:
            return False
        if str(a["strand"]) != str(b["strand"]):
            return False
        return (abs(float(a["s_start"]) - float(b["s_start"])) <= tol and
                abs(float(a["s_end"]) - float(b["s_end"])) <= tol)
    
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
        
        cons_len = float(g.loc[0, "consLen"])
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
            chroms = path_df["subject_id"].unique()
            eclass = "Cecc-InterChr" if len(chroms) > 1 else "Cecc-IntraChr"
            max_gap, avg_gap, _, _ = self.analyze_gaps(path_df)
            
            base = {
                "query_id": qid,
                "readName": g.iloc[0]["readName"],
                "eClass": eclass,
                "consLen": g.iloc[0]["consLen"],
                "copyNum": g.iloc[0]["copyNum"],
                "num_segments": len(path_df),
                "cumulative_length": res["cum_len"],
                "MatDegree": res["mat_degree"],
                "max_gap": max_gap,
                "avg_gap": round(avg_gap, 2),
                "chromosomes": ",".join(chroms.astype(str)),
            }
            
            for i in range(len(path_df)):
                seg = path_df.iloc[i]
                row = dict(base)
                row.update({
                    "segment_in_circle": i + 1,
                    "chr": seg["subject_id"],
                    "start": seg["s_start"],
                    "end": seg["s_end"],
                    "strand": seg["strand"],
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
        """Label segments with roles: head, body, tail."""
        if df_circ.empty:
            return df_circ.copy()
        
        # Stable sort
        df = df_circ.sort_values(["query_id", "segment_in_circle"], kind="mergesort").copy()
        
        parts: List[pd.DataFrame] = []
        for _, g in df.groupby("query_id", sort=False):
            n = len(g)
            roles = ["body"] * n
            if n >= 1:
                roles[0] = "head"
            if n >= 2:
                roles[-1] = "tail"
            
            g2 = g.copy()
            g2["segment_role"] = pd.Categorical(
                roles, categories=["head", "body", "tail"], ordered=True
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
        overlap_threshold: float = 0.8,
        min_segments: int = 3,
        edge_tolerance: int = 20,
        position_tolerance: int = 50
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
        
        Returns
        -------
        pd.DataFrame
            DataFrame containing labeled circular patterns
        """
        self.logger.info("=" * 60)
        self.logger.info("Trapeze - Circular eccDNA Analysis")
        self.logger.info("=" * 60)
        self.logger.info(f"Parameters:")
        self.logger.info(f"  - Overlap threshold: {overlap_threshold}")
        self.logger.info(f"  - Min segments: {min_segments}")
        self.logger.info(f"  - Edge tolerance: {edge_tolerance} bp")
        self.logger.info(f"  - Position tolerance: {position_tolerance} bp")
        self.logger.info(f"  - Numba acceleration: {'Enabled' if NUMBA_AVAILABLE else 'Disabled'}")
        
        # Read and validate input
        self.logger.info(f"\nReading input from: {input_csv}")
        df = pd.read_csv(input_csv, sep=None, engine="python")
        self.logger.info(f"Loaded {len(df)} rows, {df['query_id'].nunique()} unique queries")
        
        # Ensure required columns and clean data
        df = self.ensure_required_columns(df)
        
        # Step 1: Rotate segments within each query
        self.logger.info("\nStep 1: Rotating segments by q_start")
        df_rot = self.rotate_all(df)
        
        if df_rot.empty:
            self.logger.warning("No data after rotation step")
            empty_result = pd.DataFrame()
            empty_result.to_csv(output_csv, index=False)
            return empty_result
        
        self.logger.info(f"  - Rotated {len(df_rot)} segments")
        
        # Step 2: Filter by segments and overlap (vectorized)
        self.logger.info("\nStep 2: Filtering by segment count and overlaps")
        df_filt = self.vectorized_overlap_filter(df_rot, overlap_threshold, min_segments)
        
        if df_filt.empty:
            self.logger.warning("No data passed filtering criteria")
            empty_result = pd.DataFrame()
            empty_result.to_csv(output_csv, index=False)
            return empty_result
        
        # Step 3: Detect circular patterns
        self.logger.info("\nStep 3: Detecting circular patterns")
        df_circles = self.detect_circles(df_filt, edge_tolerance, position_tolerance)
        
        # Step 4: Label roles
        if not df_circles.empty:
            self.logger.info("\nStep 4: Labeling segment roles")
            df_labeled = self.label_roles(df_circles)
        else:
            self.logger.warning("No circular patterns detected")
            df_labeled = df_circles
        
        # Step 5: Save results
        output_csv.parent.mkdir(parents=True, exist_ok=True)
        df_labeled.to_csv(output_csv, index=False)
        
        # Generate summary statistics
        if not df_labeled.empty:
            unique_queries = df_labeled['query_id'].nunique()
            intra_chr = len(df_labeled[df_labeled["eClass"] == "Cecc-IntraChr"])
            inter_chr = len(df_labeled[df_labeled["eClass"] == "Cecc-InterChr"])
            
            # Calculate average metrics
            query_stats = df_labeled.groupby('query_id').first()
            avg_mat_degree = query_stats['MatDegree'].mean()
            avg_segments = query_stats['num_segments'].mean()
            
            self.logger.info("\n" + "=" * 60)
            self.logger.info("RESULTS SUMMARY")
            self.logger.info("=" * 60)
            self.logger.info(f"Output saved to: {output_csv}")
            self.logger.info(f"Circular queries found: {unique_queries}")
            self.logger.info(f"Total segments in circles: {len(df_labeled)}")
            self.logger.info(f"  - Cecc-IntraChr: {intra_chr} segments")
            self.logger.info(f"  - Cecc-InterChr: {inter_chr} segments")
            self.logger.info(f"Average MatDegree: {avg_mat_degree:.2f}%")
            self.logger.info(f"Average segments per circle: {avg_segments:.1f}")
        else:
            self.logger.info("\n" + "=" * 60)
            self.logger.info(f"No circular patterns found. Empty file saved to: {output_csv}")
        
        self.logger.info("=" * 60)
        
        return df_labeled