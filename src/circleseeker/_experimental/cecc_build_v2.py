"""
Cecc-Build v2 - Simplified CeccDNA Detection

Based on the key insight that TideHunter produces doubled circular sequences:
- Original ring: B1 -> B2 -> B3 -> (back to B1)
- Doubled query: [B1'][B2][B3][B1][B2][B3']
                  A1   A2  A3  A4  A5  A6

Detection logic:
1. Sort alignments by query position (A)
2. Check continuity of A (gaps <= tolerance)
3. Extract genomic regions (B)
4. Find first repeat starting from B2 (skip B1 as it may be incomplete)
5. If B2 repeats -> cycle found -> Cecc

Output format is fully compatible with cecc_build.py (V1).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Any

import pandas as pd
import numpy as np

from circleseeker.utils.logging import get_logger
from circleseeker.utils.column_standards import ColumnStandard


@dataclass
class AlignmentSegment:
    """Represents an alignment segment with all metadata."""
    # Genomic coordinates (0-based, half-open)
    chr: str
    start0: int
    end0: int
    strand: str
    # Query coordinates
    q_start: int
    q_end: int
    alignment_length: int
    # Original row data for metadata access
    row_data: dict = field(default_factory=dict)

    @property
    def length(self) -> int:
        return self.end0 - self.start0

    def matches_position(self, other: "AlignmentSegment", tolerance: int = 50) -> bool:
        """Check if two segments represent the same genomic location."""
        if self.chr != other.chr:
            return False
        if self.strand != other.strand:
            return False
        if abs(self.start0 - other.start0) > tolerance:
            return False
        if abs(self.end0 - other.end0) > tolerance:
            return False
        return True

    def __repr__(self) -> str:
        return f"{self.chr}:{self.start0}-{self.end0}:{self.strand}"


class CeccBuildV2:
    """Simplified CeccDNA detection based on doubled sequence structure.

    Output format is fully compatible with CeccBuild (V1).
    """

    # Heuristic thresholds (same as V1)
    MAPQ_MAX = 60
    MAPQ_LOW_THRESHOLD = 20
    IDENTITY_LOW_THRESHOLD = 95.0

    # Required columns from alignment data
    REQUIRED_COLS = [
        "query_id", "subject_id", "q_start", "q_end",
        "s_start", "s_end", "alignment_length"
    ]

    # Additional columns we try to preserve
    METADATA_COLS = ["reads", "length", "copy_number", "mapq", "identity"]

    def __init__(
        self,
        gap_tolerance: int = 20,
        position_tolerance: int = 50,
        min_segments: int = 2,
        min_match_degree: float = 95.0,
        logger: Optional[logging.Logger] = None,
    ):
        """
        Initialize CeccBuild v2.

        Args:
            gap_tolerance: Max gap between adjacent alignments on query (bp)
            position_tolerance: Tolerance for matching genomic positions (bp)
            min_segments: Minimum number of distinct segments for Cecc
            min_match_degree: Minimum coverage percentage to accept as Cecc
            logger: Optional logger
        """
        self.gap_tolerance = gap_tolerance
        self.position_tolerance = position_tolerance
        self.min_segments = min_segments
        self.min_match_degree = min_match_degree
        self.logger = logger or get_logger(self.__class__.__name__)

    @staticmethod
    def _clamp01(value: float) -> float:
        """Clamp value to [0, 1] range."""
        try:
            x = float(value)
        except (TypeError, ValueError):
            return 0.0
        return max(0.0, min(1.0, x))

    @classmethod
    def _norm_mapq(cls, mapq: float) -> float:
        """Normalize MAPQ to [0, 1]."""
        try:
            val = float(mapq)
        except (TypeError, ValueError):
            return 0.0
        cap = float(cls.MAPQ_MAX) if cls.MAPQ_MAX > 0 else 60.0
        return cls._clamp01(val / cap)

    @classmethod
    def _norm_identity(cls, identity_pct: float) -> float:
        """Normalize identity percentage to [0, 1]."""
        try:
            val = float(identity_pct)
        except (TypeError, ValueError):
            return 0.0
        return cls._clamp01((val - 90.0) / 10.0)

    @staticmethod
    def _geom_mean(values: list[float]) -> float:
        """Calculate geometric mean of values."""
        if not values:
            return 0.0
        prod = 1.0
        for v in values:
            if v <= 0.0:
                return 0.0
            prod *= float(v)
        return float(prod ** (1.0 / float(len(values))))

    def _parse_segment(self, row: pd.Series) -> AlignmentSegment:
        """Extract alignment segment from a row."""
        s_start = int(row["s_start"])
        s_end = int(row["s_end"])
        # Normalize to 0-based, half-open coordinates
        start0 = min(s_start, s_end) - 1  # Convert from 1-based
        end0 = max(s_start, s_end)

        # Handle strand column (may be "strand" or "sstrand")
        strand = str(row.get("strand", row.get("sstrand", "+")))
        if strand in ("plus", "1"):
            strand = "+"
        elif strand in ("minus", "-1"):
            strand = "-"
        elif strand not in ("+", "-"):
            strand = "+"

        # Preserve all row data for metadata
        row_data = row.to_dict()

        return AlignmentSegment(
            chr=str(row["subject_id"]),
            start0=start0,
            end0=end0,
            strand=strand,
            q_start=int(row["q_start"]),
            q_end=int(row["q_end"]),
            alignment_length=int(row["alignment_length"]),
            row_data=row_data,
        )

    def _check_continuity(self, segments: list[AlignmentSegment]) -> tuple[bool, float, float]:
        """Check if segments are continuous on query sequence.

        Returns:
            Tuple of (is_continuous, max_gap, avg_gap)
        """
        if len(segments) < 2:
            return True, 0.0, 0.0

        gaps = []
        for i in range(1, len(segments)):
            gap = segments[i].q_start - segments[i - 1].q_end
            gaps.append(gap)
            if gap > self.gap_tolerance:
                return False, max(gaps), sum(gaps) / len(gaps)

        return True, max(gaps) if gaps else 0.0, sum(gaps) / len(gaps) if gaps else 0.0

    def _find_cycle(
        self,
        segments: list[AlignmentSegment],
        cons_len: int,
    ) -> Optional[tuple[list[int], bool, str]]:
        """
        Find circular pattern in segments.

        For doubled sequences from TideHunter, there MUST be at least one
        repeated genomic region to confirm circularity. Without a repeat,
        we cannot confirm this is a circular DNA.

        Returns:
            Tuple of (path_indices, closure_found, closure_reason) or None
        """
        if len(segments) < 2:
            return None

        # Step 1: Find repeated genomic regions (required for Cecc)
        # For doubled sequence [A'][B][C][A][B][C'], we need to find
        # where a region repeats to identify the cycle boundary.
        repeat_found = False

        for i in range(len(segments) - 1):
            for j in range(i + 1, len(segments)):
                if segments[j].matches_position(segments[i], self.position_tolerance):
                    repeat_found = True
                    break
            if repeat_found:
                break

        if not repeat_found:
            # No repeated region found - cannot be a valid Cecc
            return None

        # Step 2: Build path (include all segments) and check coverage
        # For doubled sequences, we include all segments as they represent
        # the complete doubled structure
        path = list(range(len(segments)))
        cum_len = sum(seg.alignment_length for seg in segments)

        # Step 3: Verify coverage is sufficient
        if cons_len > 0:
            match_degree = (cum_len / cons_len) * 100
            if match_degree >= self.min_match_degree:
                return (path, True, "position_match")

        return None

    def _detect_single_query(
        self,
        query_id: str,
        group: pd.DataFrame,
    ) -> Optional[dict[str, Any]]:
        """Detect Cecc for a single query.

        Returns dict with detection results or None.
        """
        # Sort by query position
        sorted_df = group.sort_values("q_start").reset_index(drop=True)

        if len(sorted_df) < 2:
            return None

        # Parse segments
        segments = [self._parse_segment(row) for _, row in sorted_df.iterrows()]

        # Check continuity
        is_continuous, max_gap, avg_gap = self._check_continuity(segments)
        if not is_continuous:
            return None

        # Get consensus length from metadata
        cons_len = 0
        if "length" in sorted_df.columns:
            cons_len = int(sorted_df.iloc[0]["length"])
        if cons_len <= 0:
            # Estimate from query span
            cons_len = int(sorted_df["q_end"].max() - sorted_df["q_start"].min())

        # Find circular pattern
        cycle_result = self._find_cycle(segments, cons_len)
        if cycle_result is None:
            return None

        path_indices, closure_found, closure_reason = cycle_result
        path_segments = [segments[i] for i in path_indices]

        # Need at least min_segments distinct loci
        unique_loci: list[AlignmentSegment] = []
        for seg in path_segments:
            is_dup = any(seg.matches_position(u, self.position_tolerance) for u in unique_loci)
            if not is_dup:
                unique_loci.append(seg)

        if len(unique_loci) < self.min_segments:
            return None

        # Calculate match degree
        cum_len = sum(seg.alignment_length for seg in path_segments)
        match_degree = (cum_len / cons_len * 100) if cons_len > 0 else 0.0

        # Determine class
        chroms = list(set(seg.chr for seg in unique_loci))
        cecc_class = "Cecc-InterChr" if len(chroms) > 1 else "Cecc-IntraChr"

        # Extract quality metrics
        mapq_values = []
        identity_values = []
        for seg in path_segments:
            if "mapq" in seg.row_data and pd.notna(seg.row_data.get("mapq")):
                mapq_values.append(float(seg.row_data["mapq"]))
            if "identity" in seg.row_data and pd.notna(seg.row_data.get("identity")):
                identity_values.append(float(seg.row_data["identity"]))

        mapq_best = int(max(mapq_values)) if mapq_values else 0
        mapq_min = int(min(mapq_values)) if mapq_values else 0
        identity_best = max(identity_values) if identity_values else 0.0
        identity_min = min(identity_values) if identity_values else 0.0

        # Calculate coverage metrics
        query_span = sorted_df["q_end"].max() - sorted_df["q_start"].min()
        cov_best = cum_len / query_span if query_span > 0 else 0.0
        cov_2nd = 0.0  # V2 doesn't track second-best chain

        # Calculate confidence score (same as V1)
        chain_unique = len(unique_loci) / len(path_segments) if path_segments else 0.0
        conf = self._geom_mean([
            self._norm_mapq(mapq_best),
            self._norm_identity(identity_best),
            self._clamp01(cov_best),
            self._clamp01(chain_unique),
        ])

        # Quality flags
        low_mapq = mapq_best < self.MAPQ_LOW_THRESHOLD
        low_identity = identity_best < self.IDENTITY_LOW_THRESHOLD

        # Get metadata from first row
        first_row = sorted_df.iloc[0]
        reads = first_row.get("reads", query_id)
        copy_number = first_row.get("copy_number", 1.0)
        length = int(first_row.get("length", cons_len))

        return {
            "query_id": query_id,
            "path_segments": path_segments,
            "unique_loci": unique_loci,
            "closure_found": closure_found,
            "closure_reason": closure_reason,
            "cum_len": cum_len,
            "match_degree": match_degree,
            "cecc_class": cecc_class,
            "chromosomes": chroms,
            "max_gap": max_gap,
            "avg_gap": avg_gap,
            "mapq_best": mapq_best,
            "mapq_min": mapq_min,
            "identity_best": identity_best,
            "identity_min": identity_min,
            "cov_best": cov_best,
            "cov_2nd": cov_2nd,
            "confidence_score": conf,
            "low_mapq": low_mapq,
            "low_identity": low_identity,
            "reads": reads,
            "length": length,
            "copy_number": copy_number,
        }

    def detect_circles(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Detect CeccDNA from alignment results.

        Args:
            df: DataFrame with alignment results (must have required columns)

        Returns:
            DataFrame with Cecc results in V1-compatible format
        """
        # Validate columns
        missing = [c for c in self.REQUIRED_COLS if c not in df.columns]
        if missing:
            raise ValueError(f"Missing required columns: {missing}")

        rows: list[dict] = []

        for query_id, group in df.groupby("query_id", sort=False):
            result = self._detect_single_query(str(query_id), group)
            if result is None:
                continue

            # Build output rows (one per segment, same as V1)
            base = {
                "query_id": result["query_id"],
                "reads": result["reads"],
                "eccdna_type": "Cecc",
                "CeccClass": result["cecc_class"],
                "length": result["length"],
                "copy_number": result["copy_number"],
                "num_segments": len(result["path_segments"]),
                "cumulative_length": result["cum_len"],
                "match_degree": round(result["match_degree"], 2),
                "match_degree_2nd": None,  # V2 doesn't track 2nd best
                "best_chain_signature": None,
                "second_chain_signature": None,
                ColumnStandard.MAPQ_BEST: result["mapq_best"],
                ColumnStandard.MAPQ_MIN: result["mapq_min"],
                ColumnStandard.IDENTITY_BEST: round(result["identity_best"], 3),
                ColumnStandard.IDENTITY_MIN: round(result["identity_min"], 3),
                ColumnStandard.QUERY_COV_BEST: round(result["cov_best"], 6),
                ColumnStandard.QUERY_COV_2ND: round(result["cov_2nd"], 6),
                ColumnStandard.CONFIDENCE_SCORE: round(result["confidence_score"], 6),
                ColumnStandard.LOW_MAPQ: result["low_mapq"],
                ColumnStandard.LOW_IDENTITY: result["low_identity"],
                "C_cov_best": round(result["cov_best"], 6),
                "C_cov_2nd": round(result["cov_2nd"], 6),
                "max_gap": result["max_gap"],
                "avg_gap": round(result["avg_gap"], 2),
                "chromosomes": ",".join(str(c) for c in result["chromosomes"]),
            }

            for i, seg in enumerate(result["path_segments"]):
                row = dict(base)
                row.update({
                    "segment_in_circle": i + 1,
                    ColumnStandard.CHR: seg.chr,
                    ColumnStandard.START0: seg.start0,
                    ColumnStandard.END0: seg.end0,
                    ColumnStandard.STRAND: seg.strand,
                    "q_start": seg.q_start,
                    "q_end": seg.q_end,
                    "alignment_length": seg.alignment_length,
                })
                rows.append(row)

        result_df = pd.DataFrame(rows)
        if not result_df.empty:
            self.logger.info(
                f"Circular patterns detected: {result_df['query_id'].nunique()} queries"
            )

        return result_df

    def run(
        self,
        input_csv: Path,
        output_csv: Path,
        edge_tolerance: Optional[int] = None,
        position_tolerance: Optional[int] = None,
        min_match_degree: Optional[float] = None,
    ) -> pd.DataFrame:
        """
        Run the complete Cecc detection pipeline.

        This method signature is compatible with CeccBuild.run() from V1.

        Args:
            input_csv: Input alignment CSV/TSV file
            output_csv: Output CSV file
            edge_tolerance: Override gap tolerance (optional)
            position_tolerance: Override position tolerance (optional)
            min_match_degree: Override minimum match degree (optional)

        Returns:
            DataFrame with Cecc results
        """
        # Apply overrides if provided
        if edge_tolerance is not None:
            self.gap_tolerance = edge_tolerance
        if position_tolerance is not None:
            self.position_tolerance = position_tolerance
        if min_match_degree is not None:
            self.min_match_degree = min_match_degree

        self.logger.info("=" * 60)
        self.logger.info("CeccBuild v2 - Simplified CeccDNA Detection")
        self.logger.info("=" * 60)
        self.logger.info(f"Parameters:")
        self.logger.info(f"  - Gap tolerance: {self.gap_tolerance} bp")
        self.logger.info(f"  - Position tolerance: {self.position_tolerance} bp")
        self.logger.info(f"  - Min segments: {self.min_segments}")
        self.logger.info(f"  - Min match degree: {self.min_match_degree}%")

        # Read input (auto-detect separator)
        self.logger.info(f"Reading: {input_csv}")
        try:
            df = pd.read_csv(input_csv)
        except Exception:
            df = pd.read_csv(input_csv, sep="\t")

        if df.empty:
            self.logger.warning("Input file is empty")
            result_df = pd.DataFrame()
            output_csv.parent.mkdir(parents=True, exist_ok=True)
            result_df.to_csv(output_csv, index=False)
            return result_df

        self.logger.info(f"Loaded {len(df)} alignments, {df['query_id'].nunique()} queries")

        # Detect circles
        result_df = self.detect_circles(df)

        # Save output
        output_csv.parent.mkdir(parents=True, exist_ok=True)
        result_df.to_csv(output_csv, index=False)

        # Summary
        self.logger.info("=" * 60)
        self.logger.info("Results:")
        if not result_df.empty:
            n_cecc = result_df["query_id"].nunique()
            self.logger.info(f"  - Cecc detected: {n_cecc}")
            if "CeccClass" in result_df.columns:
                class_counts = result_df.groupby("query_id")["CeccClass"].first().value_counts()
                for cls, cnt in class_counts.items():
                    self.logger.info(f"  - {cls}: {cnt}")
        else:
            self.logger.info("  - Cecc detected: 0")
        self.logger.info(f"Output: {output_csv}")
        self.logger.info("=" * 60)

        return result_df


def main():
    """CLI entry point."""
    import argparse

    parser = argparse.ArgumentParser(description="CeccBuild v2 - Simplified CeccDNA Detection")
    parser.add_argument("-i", "--input", required=True, help="Input alignment TSV/CSV")
    parser.add_argument("-o", "--output", required=True, help="Output CSV")
    parser.add_argument("--gap-tolerance", type=int, default=20, help="Gap tolerance (bp)")
    parser.add_argument("--position-tolerance", type=int, default=50, help="Position tolerance (bp)")
    parser.add_argument("--min-segments", type=int, default=2, help="Min distinct segments")
    parser.add_argument("--min-match-degree", type=float, default=95.0, help="Min match degree (%%)")

    args = parser.parse_args()

    builder = CeccBuildV2(
        gap_tolerance=args.gap_tolerance,
        position_tolerance=args.position_tolerance,
        min_segments=args.min_segments,
        min_match_degree=args.min_match_degree,
    )

    builder.run(
        input_csv=Path(args.input),
        output_csv=Path(args.output),
    )


if __name__ == "__main__":
    main()
