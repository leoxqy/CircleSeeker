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
from typing import Tuple, Set, Optional
import pandas as pd
import logging


class UMeccClassifier:
    """Classify eccDNA into Uecc/Mecc buckets."""

    # Alignment column names (BLAST outfmt 6 compatible)
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

    def __init__(
        self,
        gap_threshold: float = 10.0,
        min_full_length_coverage: float = 95.0,
        logger: Optional[logging.Logger] = None,
    ):
        """
        Initialize U/Mecc classifier

        Args:
            gap_threshold: Maximum gap percentage for quality filtering (default: 10%)
            min_full_length_coverage: Minimum coverage for full-length classification (default: 95%)
            logger: Optional logger instance
        """
        self.gap_threshold = gap_threshold
        self.min_full_length_coverage = min_full_length_coverage
        self.stats = {}

        # Setup logger
        self.logger = logger or get_logger(self.__class__.__name__)

    def _ensure_alignment_columns(self, df: pd.DataFrame, source: str) -> pd.DataFrame:
        """Ensure alignment DataFrame has the expected columns."""
        expected_cols = len(self.ALIGNMENT_COLUMNS)
        if df.columns.tolist() == self.ALIGNMENT_COLUMNS:
            return df.copy()

        if len(df.columns) == 0:
            df = df.copy()
            df.columns = self.ALIGNMENT_COLUMNS
            return df

        if len(df.columns) != expected_cols:
            message = (
                f"Alignment results from {source} should have {expected_cols} columns, "
                f"got {len(df.columns)}"
            )
            self.logger.error(message)
            raise ValueError(message)

        df = df.copy()
        df.columns = self.ALIGNMENT_COLUMNS
        return df

    def _preprocess_alignment_df(self, df: pd.DataFrame, source: str) -> pd.DataFrame:
        """Common preprocessing for alignment results."""
        df = self._ensure_alignment_columns(df, source)
        if df.empty:
            return df

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
        df[ColumnStandard.START0] = df["s_start"]
        df[ColumnStandard.END0] = df["s_end"]
        df[ColumnStandard.STRAND] = df["strand"]

        # Calculate derived columns with safe division
        df["Rlength"] = df[ColumnStandard.END0] - df[ColumnStandard.START0] + 1
        df["gap_Length"] = df[ColumnStandard.LENGTH] - df["Rlength"]

        # Safe division for Gap_Percentage
        df["Gap_Percentage"] = 0.0
        valid_mask = df[ColumnStandard.LENGTH] > 0
        df.loc[valid_mask, "Gap_Percentage"] = (
            (df.loc[valid_mask, "gap_Length"].abs() / df.loc[valid_mask, ColumnStandard.LENGTH])
            * 100
        ).round(2)

        # Store statistics
        self.stats["total_alignments"] = len(df)
        self.stats["total_queries"] = df["query_id"].nunique()

        return df

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
    ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
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

    def classify_uecc_mecc(self, df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, Set[str]]:
        """
        Classify Uecc and Mecc based on high-quality alignments

        Args:
            df: Full DataFrame

        Returns:
            Tuple of (uecc_df, mecc_df, classified_query_ids)
        """
        self.logger.info("=" * 60)
        self.logger.info("Step 1: Classify using high-quality alignments")

        # Filter by quality
        high_quality = df[df["Gap_Percentage"] <= self.gap_threshold].copy()
        self.logger.info(f"High-quality alignments: {len(high_quality):,} / {len(df):,}")

        if high_quality.empty:
            self.logger.warning("No high-quality alignments found")
            return pd.DataFrame(), pd.DataFrame(), set()

        # Process each query based on high-quality alignments
        uecc_list = []
        mecc_list = []
        classified_queries = set()

        for query_id, group in high_quality.groupby("query_id"):
            if len(group) == 1:
                # Single alignment -> Uecc
                group_copy = group.copy()
                group_copy["eccdna_type"] = "Uecc"
                group_copy["classification_reason"] = "Single alignment"
                uecc_list.append(group_copy)
                classified_queries.add(query_id)
            else:
                # Multiple alignments - check for overlaps
                processed_group = self._process_overlaps_for_query(group)

                if len(processed_group) == 1:
                    # After overlap removal -> Uecc
                    processed_copy = processed_group.copy()
                    processed_copy["eccdna_type"] = "Uecc"
                    processed_copy["classification_reason"] = "Single after overlap removal"
                    uecc_list.append(processed_copy)
                    classified_queries.add(query_id)
                else:
                    # Check if it's Mecc (full-length repeats)
                    if self._is_full_length_repeat(processed_group):
                        processed_copy = processed_group.copy()
                        processed_copy["eccdna_type"] = "Mecc"
                        mecc_list.append(processed_copy)
                        classified_queries.add(query_id)
                    # If not Mecc, it remains unclassified

        # Combine results
        uecc_df = pd.concat(uecc_list, ignore_index=False) if uecc_list else pd.DataFrame()
        mecc_df = pd.concat(mecc_list, ignore_index=False) if mecc_list else pd.DataFrame()

        uecc_queries = set(uecc_df["query_id"].unique()) if not uecc_df.empty else set()
        mecc_queries = set(mecc_df["query_id"].unique()) if not mecc_df.empty else set()
        self.logger.info(f"Uecc: {len(classified_queries & uecc_queries)} queries")
        self.logger.info(f"Mecc: {len(classified_queries & mecc_queries)} queries")

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
                    best_idx = component_df["Gap_Percentage"].idxmin()
                    kept_alignments.append(chr_group.loc[[best_idx]])

        if kept_alignments:
            return pd.concat(kept_alignments, ignore_index=False)
        return pd.DataFrame()

    def _find_overlaps_sweepline(self, group: pd.DataFrame) -> list[Set[int]]:
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

        active_intervals: Set[int] = set()
        for _, event_type, idx in events:
            if event_type == 0:  # start event
                for active_idx in active_intervals:
                    union(idx, active_idx)
                active_intervals.add(idx)
            else:  # end event
                active_intervals.discard(idx)

        components: dict[int, Set[int]] = {}
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
        full_length_count = (coverages >= self.min_full_length_coverage).sum()

        return full_length_count >= 2

    def extract_unclassified(
        self, df_original: pd.DataFrame, classified_queries: Set[str]
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
            unclassified_df["quality_category"] = unclassified_df["Gap_Percentage"].apply(
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
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
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
            uecc_formatted["match_degree"] = 100 - uecc_formatted["Gap_Percentage"]
        else:
            uecc_formatted = uecc_df

        # Format Mecc (clean output)
        if not mecc_df.empty:
            mecc_formatted = mecc_df.copy()
            mecc_formatted["match_degree"] = 100 - mecc_formatted["Gap_Percentage"]
        else:
            mecc_formatted = mecc_df

        return uecc_formatted, mecc_formatted

    def run(self, alignment_file: Path) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
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
