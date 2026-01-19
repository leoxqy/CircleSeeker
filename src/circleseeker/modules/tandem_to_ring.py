"""
Tandem-to-Ring - Enhanced Circular DNA Analysis Pipeline (CircleSeeker)

This module processes tandem repeat sequences to identify and classify circular DNA.
Rewritten to exactly match the original step2_carousel.py algorithm.

Terminology (user-facing):
- Ctc = Concatemeric tandem copies
- CtcReads = reads carrying Ctc signals (tracked as CtcR-* classes in `tandem_to_ring.csv`)
- This module is a core component of the CtcReads-Caller workflow.
"""

from __future__ import annotations

import logging
from circleseeker.utils.logging import get_logger
from circleseeker.utils.column_standards import ColumnStandard
import pandas as pd
import numpy as np
import networkx as nx
from pathlib import Path
from typing import Any, Optional
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

try:
    from intervaltree import IntervalTree

    INTERVALTREE_AVAILABLE = True
except ImportError:
    INTERVALTREE_AVAILABLE = False


class TandemToRing:

    def __init__(
        self,
        input_file: str | Path,
        output_file: str | Path,
        circular_fasta: str | Path,
        min_ave_match: float = 99.0,
        logger: Optional[logging.Logger] = None,
    ):
        self.input_file = Path(input_file)
        self.output_file = Path(output_file)
        self.circular_fasta = Path(circular_fasta)
        self.min_ave_match = float(min_ave_match)

        # Setup logger
        self.logger = logger or get_logger(self.__class__.__name__)

        if not INTERVALTREE_AVAILABLE:
            self.logger.warning(
                "intervaltree not installed, falling back to O(n^2) overlap detection"
            )

    def read_and_preprocess_data(self) -> pd.DataFrame:
        """Read and preprocess tandem repeat data"""
        self.logger.info("Reading and preprocessing data...")

        # Define column names (input format)
        columns = [
            "readName",
            "repN",
            "copyNum",
            "readLen",
            "start",
            "end",
            "consLen",
            "aveMatch",
            "fullLen",
            "subPos",
            "consSeq",
        ]

        # Read TSV file
        df = pd.read_csv(self.input_file, sep="\t", header=None, names=columns)

        # Convert subPos to string and replace commas with pipes
        df["subPos"] = df["subPos"].astype(str).str.replace(",", "|")

        # Calculate Effective_Length
        df["Effective_Length"] = ((df["end"] - df["start"] + 1) / df["readLen"]) * 100

        # Drop fullLen column
        df = df.drop(columns=["fullLen"])

        # Filter low quality (aveMatch < threshold)
        df = df[df["aveMatch"] >= self.min_ave_match]

        # Reset index
        df = df.reset_index(drop=True)

        # Standardize internal columns while preserving input compatibility
        standardize_mapping = {
            "readName": ColumnStandard.READS,
            "start": ColumnStandard.START0,
            "end": ColumnStandard.END0,
        }

        # Apply standardization to internal columns
        for old_col, new_col in standardize_mapping.items():
            if old_col in df.columns and old_col != new_col:
                df[new_col] = df[old_col]

        self.logger.info(
            f"Loaded {len(df)} regions from {df[ColumnStandard.READS].nunique()} reads"
        )

        return df

    def separate_simple_complex(self, df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Separate reads into simple and complex categories"""
        self.logger.info("Separating simple and complex reads...")

        # Count occurrences of each read name
        read_counts = df[ColumnStandard.READS].value_counts()

        # Simple reads: appear only once
        simple_reads = read_counts[read_counts == 1].index
        # Complex reads: appear multiple times
        complex_reads = read_counts[read_counts > 1].index

        # Separate data
        df_simple = df[df[ColumnStandard.READS].isin(simple_reads)].copy()
        df_complex = df[df[ColumnStandard.READS].isin(complex_reads)].copy()

        complex_count = df_complex[ColumnStandard.READS].nunique()
        self.logger.info(f"Simple reads: {len(df_simple)}, Complex reads: {complex_count}")

        return df_simple, df_complex

    def classify_simple_reads(self, df_simple: pd.DataFrame) -> pd.DataFrame:
        """Classify simple reads based on effective length"""
        self.logger.info("Classifying simple reads...")

        # Create classification column
        df_simple["classification"] = "Other"  # Default to Other

        # Classify based on effective length
        # CtcR-perfect: >=99%
        df_simple.loc[df_simple["Effective_Length"] >= 99, "classification"] = "CtcR-perfect"

        # CtcR-hybrid: 70-99%
        df_simple.loc[
            (df_simple["Effective_Length"] >= 70) & (df_simple["Effective_Length"] < 99),
            "classification",
        ] = "CtcR-hybrid"

        return df_simple

    def process_complex_group_with_graph_optimized(self, group_df: pd.DataFrame) -> pd.DataFrame:
        """Process single complex read group using optimized overlap detection"""
        # Reset index for operation
        group_df = group_df.reset_index(drop=True)

        # Build overlap graph
        G = nx.Graph()

        # Add nodes
        for idx in group_df.index:
            G.add_node(idx)

        if INTERVALTREE_AVAILABLE and len(group_df) > 10:
            # Use interval tree for better performance on larger groups
            tree = IntervalTree()
            for idx, row in group_df.iterrows():
                tree[row[ColumnStandard.START0] : row[ColumnStandard.END0] + 1] = idx

            # Find overlaps efficiently
            # Note: G.add_edge is idempotent, so we only check idx < other_idx to avoid
            # processing each pair twice (saves memory by not tracking processed pairs)
            for idx, row in group_df.iterrows():
                overlaps = tree[row[ColumnStandard.START0] : row[ColumnStandard.END0] + 1]
                for interval in overlaps:
                    other_idx = interval.data
                    if idx < other_idx:  # Process each pair only once
                        other_row = group_df.loc[other_idx]
                        overlap = (
                            min(row[ColumnStandard.END0], other_row[ColumnStandard.END0])
                            - max(row[ColumnStandard.START0], other_row[ColumnStandard.START0])
                            + 1
                        )

                        if overlap > 50:
                            G.add_edge(idx, other_idx, overlap=overlap)
        else:
            # Fallback to O(n^2) method for small groups or when intervaltree not available
            for i in group_df.index:
                for j in group_df.index:
                    if i < j:
                        # Check for overlap (fixed calculation with +1)
                        overlap = (
                            min(
                                group_df.loc[i, ColumnStandard.END0],
                                group_df.loc[j, ColumnStandard.END0],
                            )
                            - max(
                                group_df.loc[i, ColumnStandard.START0],
                                group_df.loc[j, ColumnStandard.START0],
                            )
                            + 1
                        )

                        # Significant overlap if > 50bp
                        if overlap > 50:
                            G.add_edge(i, j, overlap=overlap)

        # Find all connected components
        selected_indices = []

        for component in nx.connected_components(G):
            if len(component) == 1:
                # No overlap, keep directly
                selected_indices.extend(component)
            else:
                # Has overlap, keep the longest coverage
                component_df = group_df.loc[list(component)]
                component_df = component_df.copy()
                component_df["coverage_length"] = (
                    component_df[ColumnStandard.END0] - component_df[ColumnStandard.START0] + 1
                )
                best_idx = component_df["coverage_length"].idxmax()
                selected_indices.append(best_idx)

        # Return processed data
        processed_df = group_df.loc[selected_indices].copy()

        return processed_df

    def process_all_complex_reads(self, df_complex: pd.DataFrame) -> pd.DataFrame:
        """Process all complex reads"""
        self.logger.info("Processing complex reads with optimized algorithm...")

        processed_groups = []

        # Process by read name groups
        for read_name, group in df_complex.groupby(ColumnStandard.READS):
            # Process this group
            processed_group = self.process_complex_group_with_graph_optimized(group)

            # Calculate total effective length for the group
            group_total_effective_length = processed_group["Effective_Length"].sum()

            # Add group total effective length to each row
            processed_group["group_total_effective_length"] = group_total_effective_length

            # Save processed group
            processed_groups.append(processed_group)

        # Merge all processed groups
        if processed_groups:
            df_complex_processed = pd.concat(processed_groups, ignore_index=True)
            self.logger.info(f"Processed: {len(df_complex)} -> {len(df_complex_processed)} regions")
        else:
            df_complex_processed = pd.DataFrame()

        return df_complex_processed

    def analyze_consLen_consistency(self, df_complex_processed: pd.DataFrame) -> pd.DataFrame:
        """Analyze consLen consistency within each group"""
        self.logger.info("Analyzing consLen consistency...")

        consistency_results = []

        # Analyze by read name groups
        for read_name, group in df_complex_processed.groupby(ColumnStandard.READS):
            consLen_values = group["consLen"].values

            # Calculate statistics
            mean_consLen = np.mean(consLen_values)
            std_consLen = np.std(consLen_values)
            min_consLen = np.min(consLen_values)
            max_consLen = np.max(consLen_values)

            # Coefficient of variation (CV)
            cv = (std_consLen / mean_consLen * 100) if mean_consLen > 0 else 0

            # Max/min ratio (with zero protection)
            ratio_max_min = max_consLen / min_consLen if min_consLen > 0 else np.inf

            # Check for integer multiple relationship (with zero protection)
            if min_consLen > 0:
                base_unit = min_consLen
                is_integer_multiple = True
                multiples = []
                for consLen in consLen_values:
                    multiple = consLen / base_unit
                    multiples.append(multiple)
                    if abs(multiple - round(multiple)) > 0.05:  # 5% tolerance
                        is_integer_multiple = False
            else:
                is_integer_multiple = False
                multiples = None

            consistency_results.append(
                {
                    "readName": read_name,
                    "num_regions": len(group),
                    "consLen_list": list(consLen_values),
                    "mean_consLen": mean_consLen,
                    "std_consLen": std_consLen,
                    "cv_percent": cv,
                    "min_consLen": min_consLen,
                    "max_consLen": max_consLen,
                    "ratio_max_min": ratio_max_min,
                    "is_integer_multiple": is_integer_multiple,
                    "multiples": multiples if is_integer_multiple else None,
                    "group_total_effective_length": group["group_total_effective_length"].iloc[0],
                }
            )

        # Convert to DataFrame
        consistency_df = pd.DataFrame(consistency_results)

        # Add consistency classification
        def classify_consistency(row):
            if row["cv_percent"] < 1:  # CV < 1%
                return "highly_consistent"
            elif row["cv_percent"] < 5:  # CV < 5%
                return "consistent"
            elif row["is_integer_multiple"]:
                return "integer_multiple"
            else:
                return "variable"

        consistency_df["consistency_type"] = consistency_df.apply(classify_consistency, axis=1)

        # Statistics
        highly_consistent = len(
            consistency_df[consistency_df["consistency_type"] == "highly_consistent"]
        )
        self.logger.info(f"Highly consistent groups: {highly_consistent}")

        return consistency_df

    def process_and_classify_complex_reads(
        self, df_complex_processed: pd.DataFrame, consistency_analysis: pd.DataFrame
    ) -> pd.DataFrame:
        """Process and classify complex reads"""
        self.logger.info("Classifying complex reads...")

        # Get highly consistent groups with multiple records
        highly_consistent_multi = consistency_analysis[
            (consistency_analysis["consistency_type"] == "highly_consistent")
            & (consistency_analysis["num_regions"] > 1)
        ]["readName"].tolist()

        processed_results = []

        # Process by read name groups
        for read_name, group in df_complex_processed.groupby(ColumnStandard.READS):

            if read_name in highly_consistent_multi:
                # Highly consistent multi-record groups: merge to single record
                first_row = group.iloc[0].copy()

                # Sum Effective_Length and copyNum (with proper rounding)
                first_row["Effective_Length"] = group["Effective_Length"].sum()
                first_row["copyNum"] = np.round(group["copyNum"].sum())

                # Classify
                if first_row["Effective_Length"] >= 70:
                    first_row["classification"] = "CtcR-inversion"
                else:
                    first_row["classification"] = "Other"

                processed_results.append(first_row)

            else:
                # Other groups
                if len(group) == 1:
                    # Single record, use simple read rules
                    row = group.iloc[0].copy()
                    eff_length = row["Effective_Length"]

                    if eff_length >= 99:
                        row["classification"] = "CtcR-perfect"
                    elif eff_length >= 70:
                        row["classification"] = "CtcR-hybrid"
                    else:
                        row["classification"] = "Other"

                    processed_results.append(row)

                else:
                    # Multiple records but not highly consistent, mark all as hybrid
                    for _, row in group.iterrows():
                        row_copy = row.copy()
                        row_copy["classification"] = "CtcR-hybrid"
                        processed_results.append(row_copy)

        # Merge results
        df_complex_final = pd.DataFrame(processed_results)

        return df_complex_final

    def create_main_results_table(
        self, df_simple_classified: pd.DataFrame, df_complex_final: pd.DataFrame
    ) -> pd.DataFrame:
        """Create main results table"""
        # Merge simple and complex read results
        df_main = pd.concat([df_simple_classified, df_complex_final], ignore_index=True)

        # Convert copyNum to integer with proper rounding
        df_main["copyNum"] = np.round(df_main["copyNum"]).astype(int)

        # Create unique ID column (with proper type conversion)
        df_main["unique_id"] = (
            df_main[ColumnStandard.READS].astype(str)
            + "|"
            + df_main["repN"].astype(str)
            + "|"
            + df_main["consLen"].astype(str)
            + "|"
            + df_main["copyNum"].astype(str)
        )

        # Move unique_id to first column
        cols = df_main.columns.tolist()
        cols = ["unique_id"] + [col for col in cols if col != "unique_id"]
        df_main = df_main[cols]

        return df_main

    def circularize_sequences(self, df_main: pd.DataFrame) -> list[SeqRecord]:
        """Generate circularized sequences"""
        circular_sequences = []

        for _, row in df_main.iterrows():
            # Get sequence
            seq = row["consSeq"]

            # Circularize: concatenate sequence with itself
            circular_seq = seq + seq

            # Create SeqRecord object
            record = SeqRecord(Seq(circular_seq), id=row["unique_id"] + "|circular", description="")

            circular_sequences.append(record)

        return circular_sequences

    def write_fasta(self, sequences: list[SeqRecord], output_file: Path) -> None:
        """Write sequences to FASTA file"""
        SeqIO.write(sequences, output_file, "fasta")

    def create_readname_classification(self, df_main: pd.DataFrame) -> pd.DataFrame:
        """
        Create readName-based classification table
        Each readName appears only once with its classification
        """
        # Get read name and classification
        df_readname = df_main[[ColumnStandard.READS, "classification"]].copy()

        # For reads with multiple records, determine final classification
        # Priority: CtcR-perfect > CtcR-inversion > CtcR-hybrid > Other
        priority_map = {"CtcR-perfect": 1, "CtcR-inversion": 2, "CtcR-hybrid": 3, "Other": 4}

        # Add priority column
        df_readname["priority"] = df_readname["classification"].map(priority_map)

        # Sort by read name and priority
        df_readname = df_readname.sort_values([ColumnStandard.READS, "priority"])

        # Keep the highest priority classification for each read name
        df_readname = df_readname.drop_duplicates(subset=[ColumnStandard.READS], keep="first")

        # Drop priority column
        df_readname = df_readname.drop(columns=["priority"])

        # Create legacy output format for backward compatibility
        legacy_output = df_readname.copy()
        legacy_output.columns = ["readName", "readClass"]

        return legacy_output

    def process(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Main processing pipeline"""
        self.logger.info("=" * 60)
        self.logger.info("TandemToRing")
        self.logger.info("=" * 60)

        # 1. Read and preprocess
        df = self.read_and_preprocess_data()

        # 2. Separate simple and complex reads
        df_simple, df_complex = self.separate_simple_complex(df)

        # 3. Process simple reads
        df_simple = self.classify_simple_reads(df_simple)

        # 4. Process complex reads
        if not df_complex.empty:
            # Process overlaps with optimized algorithm
            df_complex_processed = self.process_all_complex_reads(df_complex)

            # Analyze consistency
            consistency_analysis = self.analyze_consLen_consistency(df_complex_processed)

            # Classify complex reads
            df_complex_final = self.process_and_classify_complex_reads(
                df_complex_processed, consistency_analysis
            )
        else:
            df_complex_final = pd.DataFrame()

        # 5. Create main results table
        df_main = self.create_main_results_table(df_simple, df_complex_final)

        # 6. Create readName-based classification (one entry per readName)
        df_classification = self.create_readname_classification(df_main)

        # 7. Save classification results
        self.output_file.parent.mkdir(parents=True, exist_ok=True)
        df_classification.to_csv(self.output_file, index=False)
        self.logger.info(f"Classification results saved to '{self.output_file}'")

        # 8. Generate circular FASTA
        circular_sequences = self.circularize_sequences(df_main)
        self.write_fasta(circular_sequences, self.circular_fasta)
        self.logger.info(f"Circular sequences saved to '{self.circular_fasta}'")

        # 9. Display statistics
        self.logger.info("=" * 60)
        self.logger.info("Final Classification Statistics:")
        self.logger.info("-" * 60)

        # Statistics from df_classification (unique readNames)
        for class_name, count in df_classification["readClass"].value_counts().items():
            self.logger.info(f"{class_name}: {count}")

        self.logger.info("-" * 60)
        self.logger.info(f"Total unique reads: {len(df_classification)}")
        self.logger.info(f"Total records in main table: {len(df_main)}")
        self.logger.info("=" * 60)

        # Additional detailed statistics
        if len(df_main) > len(df_classification):
            self.logger.info("Note: Some reads have multiple records in the main table")
            self.logger.info("The classification file contains one entry per unique readName")

        return df_main, df_classification

    def run(self) -> None:
        """Run the tandem-to-ring processing pipeline."""
        self.process()


def _parse_args() -> Any:
    """Parse CLI arguments for direct script execution"""
    import argparse

    parser = argparse.ArgumentParser(
        description="Tandem-to-Ring - Process tandem repeats and classify circular DNA"
    )
    parser.add_argument("-i", "--input", required=True, help="Input TideHunter result TSV file")
    parser.add_argument("-o", "--output", required=True, help="Output CSV file for classification")
    parser.add_argument(
        "-c",
        "--circular-fasta",
        required=False,
        default=None,
        help="Optional: Output FASTA file (default: same as output with .fasta)",
    )
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
        help="Log level (default: INFO)",
    )
    return parser.parse_args()


def main() -> None:
    from pathlib import Path

    args = _parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)

    if args.circular_fasta is None:
        # Same directory and name as output, with .fasta extension
        circular_fasta_path = output_path.with_suffix(".fasta")
    else:
        circular_fasta_path = Path(args.circular_fasta)

    # Set root log level
    import logging as _logging

    _logging.getLogger().setLevel(getattr(_logging, args.log_level))

    runner = TandemToRing(
        input_file=input_path,
        output_file=output_path,
        circular_fasta=circular_fasta_path,
    )
    runner.process()


if __name__ == "__main__":
    main()
