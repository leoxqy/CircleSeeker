"""
Cecc-Build v2 - Simplified CeccDNA Detection

Based on the key insight that TideHunter produces doubled circular sequences:
- Original ring: B1 → B2 → B3 → (back to B1)
- Doubled query: [B1'][B2][B3][B1][B2][B3']
                  A1   A2  A3  A4  A5  A6

Detection logic:
1. Sort alignments by query position (A)
2. Check continuity of A (gaps <= tolerance)
3. Extract genomic regions (B)
4. Find first repeat starting from B2 (skip B1 as it may be incomplete)
5. If B2 repeats → cycle found → Cecc
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import pandas as pd
import numpy as np

from circleseeker.utils.logging import get_logger


@dataclass
class GenomicRegion:
    """Represents a genomic region."""
    chr: str
    start: int
    end: int
    strand: str = "+"

    @property
    def length(self) -> int:
        return self.end - self.start

    def matches(self, other: "GenomicRegion", tolerance: int = 50) -> bool:
        """Check if two regions represent the same genomic location."""
        if self.chr != other.chr:
            return False
        if self.strand != other.strand:
            return False
        if abs(self.start - other.start) > tolerance:
            return False
        if abs(self.end - other.end) > tolerance:
            return False
        return True

    def __repr__(self) -> str:
        return f"{self.chr}:{self.start}-{self.end}:{self.strand}"


@dataclass
class CeccResult:
    """Result of Cecc detection for a single query."""
    query_id: str
    regions: list[GenomicRegion]
    cycle_length: int
    is_inter_chr: bool
    coverage: float
    # Quality metrics from alignment data
    mapq_best: int = 0
    identity_best: float = 0.0
    query_cov_best: float = 0.0
    query_cov_2nd: float = 0.0
    # Sequence extraction info (for umc_process)
    q_start: int = 0
    length: int = 0

    @property
    def num_regions(self) -> int:
        return len(self.regions)

    @property
    def cecc_class(self) -> str:
        return "Cecc-InterChr" if self.is_inter_chr else "Cecc-IntraChr"


class CeccBuildV2:
    """Simplified CeccDNA detection based on doubled sequence structure."""

    # Required columns from alignment data
    # Note: sstrand may be called "strand" in unclassified.csv
    REQUIRED_COLS = [
        "query_id", "subject_id", "q_start", "q_end",
        "s_start", "s_end", "alignment_length"
    ]

    def __init__(
        self,
        gap_tolerance: int = 50,
        region_tolerance: int = 50,
        min_regions: int = 2,
        logger: Optional[logging.Logger] = None,
    ):
        """
        Initialize CeccBuild v2.

        Args:
            gap_tolerance: Max gap between adjacent alignments on query (bp)
            region_tolerance: Tolerance for matching genomic regions (bp)
            min_regions: Minimum number of distinct regions for Cecc
            logger: Optional logger
        """
        self.gap_tolerance = gap_tolerance
        self.region_tolerance = region_tolerance
        self.min_regions = min_regions
        self.logger = logger or get_logger(self.__class__.__name__)

    def _parse_region(self, row: pd.Series) -> GenomicRegion:
        """Extract genomic region from alignment row."""
        s_start = int(row["s_start"])
        s_end = int(row["s_end"])
        # Normalize coordinates (ensure start < end)
        start = min(s_start, s_end)
        end = max(s_start, s_end)

        # Handle both "sstrand" and "strand" column names
        strand = str(row.get("sstrand", row.get("strand", "+")))
        if strand not in ("+", "-", "plus", "minus", "1", "-1"):
            strand = "+"
        if strand in ("plus", "1"):
            strand = "+"
        elif strand in ("minus", "-1"):
            strand = "-"

        return GenomicRegion(
            chr=str(row["subject_id"]),
            start=start,
            end=end,
            strand=strand,
        )

    def _check_continuity(self, sorted_alns: pd.DataFrame) -> bool:
        """Check if alignments are continuous on query sequence."""
        if len(sorted_alns) < 2:
            return True

        q_starts = sorted_alns["q_start"].values
        q_ends = sorted_alns["q_end"].values

        for i in range(1, len(sorted_alns)):
            gap = q_starts[i] - q_ends[i - 1]
            if gap > self.gap_tolerance:
                return False

        return True

    def _merge_adjacent_regions(
        self,
        regions: list[GenomicRegion],
    ) -> list[GenomicRegion]:
        """
        Merge genomically adjacent regions.

        If two consecutive regions are on the same chromosome, same strand,
        and their genomic coordinates are adjacent (within tolerance),
        merge them into one region.

        This handles cases where a long genomic region is split into
        multiple alignment segments.
        """
        if len(regions) <= 1:
            return regions

        merged = [regions[0]]
        for r in regions[1:]:
            last = merged[-1]
            # Check if same chr, same strand, and genomically adjacent
            if (last.chr == r.chr and
                last.strand == r.strand and
                abs(r.start - last.end) <= self.region_tolerance):
                # Merge: extend the last region
                merged[-1] = GenomicRegion(
                    chr=last.chr,
                    start=min(last.start, r.start),
                    end=max(last.end, r.end),
                    strand=last.strand,
                )
            else:
                merged.append(r)

        return merged

    def _find_region_repeat(
        self,
        regions: list[GenomicRegion],
    ) -> Optional[tuple[int, int]]:
        """
        Find the first repeat by checking all regions starting from region[1].

        We skip region[0] as it may be incomplete (the ring's start point
        could be anywhere). We then check region[1], region[2], etc.
        until we find one that repeats later in the sequence.

        Returns:
            Tuple of (start_index, end_index) if repeat found, None otherwise.
            The cycle is regions[start_index:end_index]
        """
        if len(regions) < 3:
            return None

        # Check regions[1], regions[2], ... until we find a repeat
        for start_idx in range(1, len(regions) - 1):
            target = regions[start_idx]
            for end_idx in range(start_idx + 1, len(regions)):
                if regions[end_idx].matches(target, self.region_tolerance):
                    # Found repeat! Cycle is from start_idx to end_idx
                    return (start_idx, end_idx)

        return None

    def _detect_single_query(
        self,
        query_id: str,
        alignments: pd.DataFrame,
    ) -> Optional[CeccResult]:
        """Detect Cecc for a single query."""
        # Sort by query position
        sorted_alns = alignments.sort_values("q_start").reset_index(drop=True)

        # Check continuity
        if not self._check_continuity(sorted_alns):
            return None

        # Extract genomic regions
        regions = [self._parse_region(row) for _, row in sorted_alns.iterrows()]

        # Merge adjacent regions (handles split alignments)
        regions = self._merge_adjacent_regions(regions)

        if len(regions) < 3:
            return None

        # Find repeat (check all regions starting from region[1])
        repeat_info = self._find_region_repeat(regions)

        if repeat_info is None:
            return None

        start_idx, end_idx = repeat_info
        cycle_regions = regions[start_idx:end_idx]

        # Check minimum distinct regions
        unique_regions: list[GenomicRegion] = []
        for r in cycle_regions:
            is_dup = False
            for u in unique_regions:
                if r.matches(u, self.region_tolerance):
                    is_dup = True
                    break
            if not is_dup:
                unique_regions.append(r)

        if len(unique_regions) < self.min_regions:
            return None

        # Calculate coverage
        total_aln_len = sorted_alns["alignment_length"].sum()
        query_len = sorted_alns["q_end"].max() - sorted_alns["q_start"].min()
        coverage = total_aln_len / query_len if query_len > 0 else 0.0

        # Check if inter-chromosomal
        chrs = set(r.chr for r in unique_regions)
        is_inter_chr = len(chrs) > 1

        # Extract quality metrics from alignment data
        mapq_best = 0
        identity_best = 0.0
        if "mapq" in sorted_alns.columns:
            mapq_best = int(sorted_alns["mapq"].max())
        if "identity" in sorted_alns.columns:
            identity_best = float(sorted_alns["identity"].max())

        # Calculate query coverage (best and 2nd best by subject)
        query_cov_best = coverage  # Use overall coverage as best
        query_cov_2nd = 0.0
        if len(unique_regions) > 1:
            # Approximate 2nd coverage from region lengths
            region_lengths = sorted([r.length for r in unique_regions], reverse=True)
            if len(region_lengths) > 1 and query_len > 0:
                query_cov_2nd = region_lengths[1] / query_len

        # Extract q_start and length for sequence extraction
        # q_start: first alignment's query start position
        q_start = int(sorted_alns["q_start"].min())
        # length: try to get from input data, otherwise use total region length
        if "length" in sorted_alns.columns:
            length = int(sorted_alns["length"].iloc[0])
        else:
            # Use sum of unique region lengths as approximation
            length = sum(r.length for r in unique_regions)

        return CeccResult(
            query_id=query_id,
            regions=unique_regions,
            cycle_length=len(unique_regions),
            is_inter_chr=is_inter_chr,
            coverage=coverage,
            mapq_best=mapq_best,
            identity_best=identity_best,
            query_cov_best=query_cov_best,
            query_cov_2nd=query_cov_2nd,
            q_start=q_start,
            length=length,
        )

    def detect(self, alignments_df: pd.DataFrame) -> list[CeccResult]:
        """
        Detect CeccDNA from alignment results.

        Args:
            alignments_df: DataFrame with alignment results

        Returns:
            List of CeccResult for detected Cecc
        """
        # Validate columns
        missing = [c for c in self.REQUIRED_COLS if c not in alignments_df.columns]
        if missing:
            raise ValueError(f"Missing required columns: {missing}")

        results = []

        for query_id, group in alignments_df.groupby("query_id", sort=False):
            result = self._detect_single_query(str(query_id), group)
            if result is not None:
                results.append(result)

        self.logger.info(f"Detected {len(results)} CeccDNA from {alignments_df['query_id'].nunique()} queries")

        return results

    # Thresholds for low quality flags
    MAPQ_LOW_THRESHOLD = 20
    IDENTITY_LOW_THRESHOLD = 90.0

    def results_to_dataframe(self, results: list[CeccResult]) -> pd.DataFrame:
        """Convert CeccResult list to DataFrame with all required columns."""
        rows = []

        for res in results:
            # Calculate confidence score (geometric mean of normalized metrics)
            norm_mapq = min(res.mapq_best / 60.0, 1.0) if res.mapq_best > 0 else 0.0
            norm_identity = res.identity_best / 100.0 if res.identity_best > 0 else 0.0
            norm_cov = min(res.query_cov_best, 1.0)

            if norm_mapq > 0 and norm_identity > 0 and norm_cov > 0:
                confidence_score = (norm_mapq * norm_identity * norm_cov) ** (1/3)
            else:
                confidence_score = 0.0

            # Quality flags
            low_mapq = res.mapq_best < self.MAPQ_LOW_THRESHOLD
            low_identity = res.identity_best < self.IDENTITY_LOW_THRESHOLD

            for i, region in enumerate(res.regions):
                rows.append({
                    # Required columns per contract
                    "query_id": res.query_id,
                    "eccdna_type": "Cecc",
                    "chr": region.chr,
                    "start0": region.start,  # 0-based
                    "end0": region.end,      # 0-based, half-open
                    "strand": region.strand,
                    "segment_in_circle": i + 1,
                    "confidence_score": round(confidence_score, 4),
                    "mapq_best": res.mapq_best,
                    "identity_best": round(res.identity_best, 3),
                    "query_cov_best": round(res.query_cov_best, 4),
                    "query_cov_2nd": round(res.query_cov_2nd, 4),
                    "low_mapq": low_mapq,
                    "low_identity": low_identity,
                    # Sequence extraction columns (for umc_process)
                    "q_start": res.q_start,
                    "length": res.length,
                    # Extra informational columns
                    "cecc_class": res.cecc_class,
                    "num_regions": res.num_regions,
                    "coverage": round(res.coverage, 4),
                    "region_length": region.length,
                })

        return pd.DataFrame(rows)

    def run_pipeline(
        self,
        input_csv: Path,
        output_csv: Path,
    ) -> pd.DataFrame:
        """
        Run the complete Cecc detection pipeline.

        Args:
            input_csv: Input alignment CSV/TSV file
            output_csv: Output CSV file

        Returns:
            DataFrame with Cecc results
        """
        self.logger.info("=" * 60)
        self.logger.info("CeccBuild v2 - Simplified CeccDNA Detection")
        self.logger.info("=" * 60)
        self.logger.info(f"Parameters:")
        self.logger.info(f"  - Gap tolerance: {self.gap_tolerance} bp")
        self.logger.info(f"  - Region tolerance: {self.region_tolerance} bp")
        self.logger.info(f"  - Min regions: {self.min_regions}")

        # Read input (auto-detect separator)
        self.logger.info(f"Reading: {input_csv}")
        try:
            df = pd.read_csv(input_csv)  # Try CSV first
        except Exception:
            df = pd.read_csv(input_csv, sep="\t")  # Fallback to TSV

        # Define required columns for empty DataFrame
        required_columns = [
            "query_id", "eccdna_type", "chr", "start0", "end0", "strand",
            "segment_in_circle", "confidence_score", "mapq_best", "identity_best",
            "query_cov_best", "query_cov_2nd", "low_mapq", "low_identity",
            "q_start", "length",  # For sequence extraction
            "cecc_class", "num_regions", "coverage", "region_length",
        ]

        if df.empty:
            self.logger.warning("Input file is empty")
            result_df = pd.DataFrame(columns=required_columns)
            output_csv.parent.mkdir(parents=True, exist_ok=True)
            result_df.to_csv(output_csv, index=False)
            return result_df

        self.logger.info(f"Loaded {len(df)} alignments, {df['query_id'].nunique()} queries")

        # Detect
        results = self.detect(df)

        # Convert to DataFrame
        if results:
            result_df = self.results_to_dataframe(results)
        else:
            result_df = pd.DataFrame(columns=required_columns)

        # Save
        output_csv.parent.mkdir(parents=True, exist_ok=True)
        result_df.to_csv(output_csv, index=False)

        # Summary
        self.logger.info("=" * 60)
        self.logger.info("Results:")
        self.logger.info(f"  - Cecc detected: {len(results)}")
        if results:
            inter_chr = sum(1 for r in results if r.is_inter_chr)
            intra_chr = len(results) - inter_chr
            self.logger.info(f"  - Cecc-IntraChr: {intra_chr}")
            self.logger.info(f"  - Cecc-InterChr: {inter_chr}")
            avg_regions = sum(r.num_regions for r in results) / len(results)
            self.logger.info(f"  - Avg regions per Cecc: {avg_regions:.1f}")
        self.logger.info(f"Output: {output_csv}")
        self.logger.info("=" * 60)

        return result_df


def main():
    """CLI entry point."""
    import argparse

    parser = argparse.ArgumentParser(description="CeccBuild v2 - Simplified CeccDNA Detection")
    parser.add_argument("-i", "--input", required=True, help="Input alignment TSV")
    parser.add_argument("-o", "--output", required=True, help="Output CSV")
    parser.add_argument("--gap-tolerance", type=int, default=50, help="Gap tolerance (bp)")
    parser.add_argument("--region-tolerance", type=int, default=50, help="Region tolerance (bp)")
    parser.add_argument("--min-regions", type=int, default=2, help="Min distinct regions")

    args = parser.parse_args()

    builder = CeccBuildV2(
        gap_tolerance=args.gap_tolerance,
        region_tolerance=args.region_tolerance,
        min_regions=args.min_regions,
    )

    builder.run_pipeline(
        input_csv=Path(args.input),
        output_csv=Path(args.output),
    )


if __name__ == "__main__":
    main()
