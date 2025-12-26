#!/usr/bin/env python3
"""
Unified UMC Process Module with Correct Logic
- U type: Single location per query_id, cluster by location signature
- M type: Multiple locations per query_id, cluster by multi-location signature
- C type: Multiple segments per query_id, cluster by ordered segment signature
All using new column names: chr, start0, end0, length, etc.
"""

import logging
from circleseeker.utils.logging import get_logger, setup_logging
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Union

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


@dataclass
class UMCProcessConfig:
    """Configuration for UMC processing pipeline."""

    process_xecc: bool = True
    cluster_uecc: bool = True
    cluster_mecc: bool = True
    cluster_cecc: bool = True
    debug: bool = False


def _sanitize_fasta_id(id_str: str) -> str:
    """Sanitize a string to be a valid FASTA identifier."""
    id_str = re.sub(r"\s+", "_", str(id_str))
    id_str = re.sub(r"[^A-Za-z0-9._-]", "_", id_str)
    return id_str


def _coerce_to_paths(value: Optional[Union[Path, Sequence[Path]]]) -> List[Path]:
    """Normalize path inputs to a list of Path objects."""
    if value is None:
        return []
    if isinstance(value, (str, Path)):
        return [Path(value)]
    paths: List[Path] = []
    for item in value:
        if not item:
            continue
        paths.append(Path(item))
    return paths


class SequenceLibrary:
    """Manages FASTA sequences with flexible lookup strategies."""

    def __init__(self, logger: Optional[logging.Logger] = None):
        """Initialize sequence library."""
        self.logger = logger if logger else get_logger(self.__class__.__name__)
        self.fasta_sequences: Dict[str, str] = {}
        self.primary_ids: Set[str] = set()

    def load_fasta(self, fasta_file: Path) -> None:
        """Load FASTA sequences into memory."""
        if not fasta_file.exists():
            raise FileNotFoundError(f"FASTA file not found: {fasta_file}")

        self.fasta_sequences.clear()
        self.primary_ids.clear()

        record_count = 0
        for record in SeqIO.parse(str(fasta_file), "fasta"):
            self.fasta_sequences[record.id] = str(record.seq)
            self.primary_ids.add(record.id)

            if "|" in record.id:
                base_id = record.id.split("|")[0]
                if base_id not in self.fasta_sequences:
                    self.fasta_sequences[base_id] = str(record.seq)

            record_count += 1

        self.logger.info(f"Loaded {record_count:,} sequences from {fasta_file.name}")

    def find_sequence(self, query_id: str) -> Optional[str]:
        """Find sequence by query_id with multiple fallback strategies."""
        query_id = str(query_id)

        if query_id in self.fasta_sequences:
            return self.fasta_sequences[query_id]

        if "|" in query_id:
            base_id = query_id.split("|")[0]
            if base_id in self.fasta_sequences:
                return self.fasta_sequences[base_id]

        return None


def extract_ring_sequence(seq_str: str, q_start: int, cons_len: int) -> str:
    """Extract circular sequence with proper wrapping."""
    seq_len = len(seq_str)
    start_index = q_start - 1

    if q_start > cons_len:
        start_index = (q_start - cons_len) - 1

    end_index = start_index + cons_len

    if end_index <= seq_len:
        return seq_str[start_index:end_index]
    else:
        part1 = seq_str[start_index:]
        part2 = seq_str[: (end_index - seq_len)]
        return part1 + part2


class UeccProcessor:
    """Process UeccDNA sequences (single location per query_id)."""

    def __init__(
        self,
        seq_library: SequenceLibrary,
        config: UMCProcessConfig,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize UeccDNA processor."""
        self.seq_library = seq_library
        self.config = config
        self.logger = logger if logger else get_logger(self.__class__.__name__)
        self.fasta_records: List[SeqRecord] = []
        self.counter = 0

    def cluster_by_location(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Cluster UeccDNA by single location signature.
        Each query_id should have only one row with location chr:start0-end0.
        """
        if not all(col in df.columns for col in ["chr", "start0", "end0"]):
            df["location_signature"] = ""
            df["cluster_id"] = 0
            df["cluster_size"] = 1
            return df

        # Generate signature for each row (UeccDNA has one row per query_id)
        df["location_signature"] = df.apply(
            lambda row: f"{row['chr']}:{row['start0']}-{row['end0']}", axis=1
        )

        # Group by signature
        signature_groups = df.groupby("location_signature")

        result_rows = []
        cluster_id = 0

        for signature, group in signature_groups:
            if signature and signature != ":0-0":
                cluster_id += 1
                cluster_size = len(group)

                # Use first row as representative
                representative = group.iloc[0].copy()

                # Aggregate copy_number
                if "copy_number" in group.columns:
                    total_copynum = group["copy_number"].sum()
                    representative["copy_number"] = total_copynum

                # Average Gap_Percentage and match_degree
                if "Gap_Percentage" in group.columns:
                    gaps = group["Gap_Percentage"].dropna()
                    if len(gaps) > 0:
                        avg_gap = gaps.mean()
                        representative["Gap_Percentage"] = round(avg_gap, 2)
                        representative["match_degree"] = round(100 - avg_gap, 2)

                # Aggregate reads if present
                if "reads" in group.columns:
                    all_reads = []
                    for reads_str in group["reads"].dropna():
                        if str(reads_str) not in ["", "NA"]:
                            reads_list = str(reads_str).split(";")
                            all_reads.extend(r.strip() for r in reads_list if r.strip())
                    unique_reads = list(dict.fromkeys(all_reads))
                    representative["reads"] = ";".join(unique_reads) if unique_reads else ""

                # Add cluster metadata
                representative["cluster_id"] = cluster_id
                representative["cluster_size"] = cluster_size
                representative["cluster_members"] = ";".join(group["query_id"].astype(str))

                result_rows.append(representative)
            else:
                # Unclustered
                for _, row in group.iterrows():
                    row_copy = row.copy()
                    row_copy["cluster_id"] = 0
                    row_copy["cluster_size"] = 1
                    row_copy["cluster_members"] = row["query_id"]
                    result_rows.append(row_copy)

        result_df = pd.DataFrame(result_rows)

        self.logger.debug(f"UeccDNA clustering: {len(df)} → {len(result_df)} rows")
        if cluster_id > 0:
            self.logger.debug(f"Formed {cluster_id} clusters")

        return result_df

    def compute_sequences(self, df: pd.DataFrame) -> pd.DataFrame:
        """Extract sequences for UeccDNA (one per row)."""
        self.logger.debug("Extracting UeccDNA sequences...")

        df["eSeq"] = ""

        if "length" not in df.columns or "q_start" not in df.columns:
            self.logger.warning("Missing required columns: length and/or q_start")
            return df

        sequences_found = 0
        sequences_missing = 0

        # For UeccDNA, process each row directly
        for idx, row in df.iterrows():
            query_id = row["query_id"]
            q_start = row["q_start"]
            cons_len = row["length"]

            seq_str = self.seq_library.find_sequence(query_id)

            if seq_str and q_start > 0 and cons_len > 0:
                try:
                    extracted_seq = extract_ring_sequence(seq_str, int(q_start), int(cons_len))
                    df.loc[idx, "eSeq"] = extracted_seq
                    sequences_found += 1
                except Exception as e:
                    self.logger.debug(f"Failed to extract sequence for {query_id}: {e}")
                    sequences_missing += 1
            else:
                sequences_missing += 1

        self.logger.debug(f"Sequences extracted: {sequences_found}, missing: {sequences_missing}")
        return df

    def add_numbering_and_export(self, df: pd.DataFrame) -> pd.DataFrame:
        """Add UeccDNA IDs and prepare FASTA records."""
        df = df.copy()
        df["eccDNA_id"] = "NA"

        if "eSeq" not in df.columns:
            return df

        for idx, row in df.iterrows():
            seq = row["eSeq"]
            if seq and len(seq) > 0:
                self.counter += 1
                eccDNA_id = f"U{self.counter}"
                df.loc[idx, "eccDNA_id"] = eccDNA_id

                query_id = row["query_id"]
                if "cluster_id" in df.columns and row["cluster_id"] > 0:
                    description = f"UeccDNA_cluster{row['cluster_id']}_{query_id}"
                else:
                    description = f"UeccDNA_{query_id}"

                record = SeqRecord(Seq(seq), id=eccDNA_id, description=description)
                self.fasta_records.append(record)

        return df

    def process(
        self,
        csv_files: Sequence[Path],
        output_dir: Path,
        prefix: str = "sample",
        cluster: bool = True,
    ) -> Optional[pd.DataFrame]:
        """Main processing pipeline for UeccDNA."""
        self.fasta_records = []
        self.counter = 0

        # Load all CSV files
        dataframes: List[pd.DataFrame] = []
        for csv_file in csv_files:
            if not csv_file:
                continue
            csv_path = Path(csv_file)
            if not csv_path.exists() or csv_path.stat().st_size == 0:
                continue
            try:
                df = pd.read_csv(csv_path)
                if not df.empty:
                    dataframes.append(df)
            except Exception as exc:
                self.logger.warning(f"Unable to read {csv_path}: {exc}")

        if not dataframes:
            self.logger.info("No UeccDNA inputs provided")
            return None

        df = pd.concat(dataframes, ignore_index=True)
        initial_count = len(df)

        # Apply clustering if requested
        if cluster and self.config.cluster_uecc:
            df = self.cluster_by_location(df)

        # Extract sequences and add IDs
        df = self.compute_sequences(df)
        df = self.add_numbering_and_export(df)

        # Save outputs
        output_dir.mkdir(parents=True, exist_ok=True)
        output_csv = output_dir / f"{prefix}_UeccDNA_processed.csv"
        df.to_csv(output_csv, index=False)

        if self.fasta_records:
            output_fasta = output_dir / f"{prefix}_UeccDNA_pre.fasta"
            with open(output_fasta, "w") as handle:
                SeqIO.write(self.fasta_records, handle, "fasta")
            self.logger.debug(f"Saved {len(self.fasta_records)} UeccDNA sequences")

        # Summary
        sequences_count = len(self.fasta_records)
        self.logger.info(
            f"UeccDNA: original sequences={initial_count}, after clustering={sequences_count}"
        )

        return df


class MeccProcessor:
    """Process MeccDNA sequences (multiple locations per query_id)."""

    def __init__(
        self,
        seq_library: SequenceLibrary,
        config: UMCProcessConfig,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize MeccDNA processor."""
        self.seq_library = seq_library
        self.config = config
        self.logger = logger if logger else get_logger(self.__class__.__name__)
        self.fasta_records: List[SeqRecord] = []
        self.counter = 0

    def generate_mecc_signature(self, df_group: pd.DataFrame) -> str:
        """
        Generate location signature for MeccDNA.
        Signature format: chr1:start1-end1;chr2:start2-end2...
        Sorted by (chr, start, end)
        """
        if not all(col in df_group.columns for col in ["chr", "start0", "end0"]):
            return ""

        locations = []
        for _, row in df_group.iterrows():
            chrom = row["chr"]
            start = row["start0"]
            end = row["end0"]

            try:
                start = int(start)
                end = int(end)
                locations.append((chrom, start, end))
            except (ValueError, TypeError):
                continue

        if not locations:
            return ""

        # Sort to ensure consistency
        locations.sort(key=lambda x: (str(x[0]), x[1], x[2]))

        # Generate signature
        signature = ";".join([f"{chr}:{start}-{end}" for chr, start, end in locations])
        return signature

    def cluster_by_signature(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Cluster MeccDNA by multi-location signature.
        Each query_id has multiple rows (multiple locations).
        """
        self.logger.debug("Clustering MeccDNA by location signatures...")

        if not all(col in df.columns for col in ["chr", "start0", "end0"]):
            self.logger.warning("Location columns not found, skipping clustering")
            df["mecc_signature"] = ""
            df["cluster_id"] = 0
            df["cluster_size"] = 1
            return df

        # Generate signatures for each query_id
        signature_map = {}
        query_id_to_group = {}

        for query_id, group in df.groupby("query_id"):
            signature = self.generate_mecc_signature(group)
            signature_map[query_id] = signature
            query_id_to_group[query_id] = group
            df.loc[group.index, "mecc_signature"] = signature

        # Group query_ids by signature
        signature_to_queries = defaultdict(list)
        for query_id, signature in signature_map.items():
            if signature:
                signature_to_queries[signature].append(query_id)

        rows_to_keep = []
        cluster_id = 0

        # Process each cluster
        for signature, query_ids in signature_to_queries.items():
            cluster_id += 1
            cluster_size = len(query_ids)

            # Use first query_id as representative
            representative_qid = query_ids[0]
            representative_group = query_id_to_group[representative_qid].copy()

            # For each location in the representative group
            for idx, rep_row in representative_group.iterrows():
                chr_match = rep_row["chr"]
                start_match = rep_row["start0"]
                end_match = rep_row["end0"]

                # Aggregate copy_number
                if "copy_number" in representative_group.columns:
                    total_copynum = 0
                    for qid in query_ids:
                        group = query_id_to_group[qid]
                        matching_rows = group[
                            (group["chr"] == chr_match)
                            & (group["start0"] == start_match)
                            & (group["end0"] == end_match)
                        ]
                        if not matching_rows.empty:
                            copynum_val = matching_rows.iloc[0].get("copy_number", 0)
                            if pd.notna(copynum_val):
                                total_copynum += copynum_val
                    representative_group.loc[idx, "copy_number"] = total_copynum

                # Aggregate Gap_Percentage
                if "Gap_Percentage" in representative_group.columns:
                    gaps = []
                    for qid in query_ids:
                        group = query_id_to_group[qid]
                        matching_rows = group[
                            (group["chr"] == chr_match)
                            & (group["start0"] == start_match)
                            & (group["end0"] == end_match)
                        ]
                        if not matching_rows.empty:
                            gap = matching_rows.iloc[0].get("Gap_Percentage")
                            if pd.notna(gap):
                                gaps.append(float(gap))

                    if gaps:
                        avg_gap = np.mean(gaps)
                        representative_group.loc[idx, "Gap_Percentage"] = round(avg_gap, 2)
                        representative_group.loc[idx, "match_degree"] = round(100 - avg_gap, 2)

            # Add cluster information
            representative_group["cluster_id"] = cluster_id
            representative_group["cluster_size"] = cluster_size
            representative_group["cluster_members"] = ";".join(query_ids)

            rows_to_keep.append(representative_group)

        # Handle unclustered sequences
        for query_id, signature in signature_map.items():
            if not signature:
                unclustered_group = query_id_to_group[query_id].copy()
                unclustered_group["cluster_id"] = 0
                unclustered_group["cluster_size"] = 1
                unclustered_group["cluster_members"] = query_id
                rows_to_keep.append(unclustered_group)

        if rows_to_keep:
            result_df = pd.concat(rows_to_keep, ignore_index=True)
        else:
            result_df = df.copy()
            result_df["cluster_id"] = 0
            result_df["cluster_size"] = 1
            result_df["cluster_members"] = result_df["query_id"]

        self.logger.debug(f"MeccDNA clustering: {len(df)} → {len(result_df)} rows")
        return result_df

    def compute_sequences(self, df: pd.DataFrame) -> pd.DataFrame:
        """Extract sequences for all unique query_ids."""
        self.logger.debug("Extracting MeccDNA sequences...")

        df["eSeq"] = ""

        if "length" not in df.columns or "q_start" not in df.columns:
            self.logger.warning("Missing required columns: length and/or q_start")
            return df

        sequences_found = 0
        sequences_missing = 0

        # Process each unique query_id (use first row for each)
        for query_id in df["query_id"].unique():
            query_mask = df["query_id"] == query_id
            first_idx = df[query_mask].index[0]

            q_start = df.loc[first_idx, "q_start"]
            cons_len = df.loc[first_idx, "length"]

            seq_str = self.seq_library.find_sequence(query_id)

            if seq_str and q_start > 0 and cons_len > 0:
                try:
                    extracted_seq = extract_ring_sequence(seq_str, int(q_start), int(cons_len))
                    df.loc[query_mask, "eSeq"] = extracted_seq
                    sequences_found += 1
                except Exception:
                    sequences_missing += 1
            else:
                sequences_missing += 1

        self.logger.debug(f"Sequences extracted: {sequences_found}, missing: {sequences_missing}")
        return df

    def add_numbering_and_export(self, df: pd.DataFrame) -> pd.DataFrame:
        """Add MeccDNA IDs and prepare FASTA records."""
        df = df.copy()
        df["eccDNA_id"] = "NA"

        if "eSeq" not in df.columns:
            return df

        # If we have clusters, assign one ID per cluster
        if "cluster_id" in df.columns and df["cluster_id"].max() > 0:
            processed_clusters = set()

            for cluster_id in df[df["cluster_id"] > 0]["cluster_id"].unique():
                if cluster_id in processed_clusters:
                    continue

                cluster_mask = df["cluster_id"] == cluster_id
                cluster_df = df[cluster_mask]

                first_idx = cluster_df.index[0]
                seq = cluster_df.loc[first_idx, "eSeq"]

                if seq and len(seq) > 0:
                    self.counter += 1
                    eccDNA_id = f"M{self.counter}"
                    df.loc[cluster_mask, "eccDNA_id"] = eccDNA_id

                    query_ids = cluster_df["query_id"].unique()
                    description = f"MeccDNA_cluster{cluster_id}_{len(query_ids)}queries"

                    record = SeqRecord(Seq(seq), id=eccDNA_id, description=description)
                    self.fasta_records.append(record)

                processed_clusters.add(cluster_id)

            # Handle unclustered
            unclustered_mask = df["cluster_id"] == 0
            if unclustered_mask.any():
                for query_id in df[unclustered_mask]["query_id"].unique():
                    query_mask = (df["query_id"] == query_id) & unclustered_mask
                    first_idx = df[query_mask].index[0]
                    seq = df.loc[first_idx, "eSeq"]

                    if seq and len(seq) > 0:
                        self.counter += 1
                        eccDNA_id = f"M{self.counter}"
                        df.loc[query_mask, "eccDNA_id"] = eccDNA_id

                        record = SeqRecord(
                            Seq(seq), id=eccDNA_id, description=f"MeccDNA_{query_id}"
                        )
                        self.fasta_records.append(record)
        else:
            # No clustering
            for query_id, group in df.groupby("query_id"):
                first_idx = group.index[0]
                seq = group.loc[first_idx, "eSeq"]

                if seq and len(seq) > 0:
                    self.counter += 1
                    eccDNA_id = f"M{self.counter}"
                    df.loc[group.index, "eccDNA_id"] = eccDNA_id

                    record = SeqRecord(Seq(seq), id=eccDNA_id, description=f"MeccDNA_{query_id}")
                    self.fasta_records.append(record)

        return df

    def process(
        self,
        csv_files: Sequence[Path],
        output_dir: Path,
        prefix: str = "sample",
        cluster: bool = True,
    ) -> Optional[pd.DataFrame]:
        """Main processing pipeline for MeccDNA."""
        self.fasta_records = []
        self.counter = 0

        # Load all CSV files
        dataframes: List[pd.DataFrame] = []
        for csv_file in csv_files:
            if not csv_file:
                continue
            csv_path = Path(csv_file)
            if not csv_path.exists() or csv_path.stat().st_size == 0:
                continue
            try:
                df = pd.read_csv(csv_path)
                if not df.empty:
                    dataframes.append(df)
            except Exception as exc:
                self.logger.warning(f"Unable to read {csv_path}: {exc}")

        if not dataframes:
            self.logger.info("No MeccDNA inputs provided")
            return None

        df = pd.concat(dataframes, ignore_index=True)
        initial_query_ids = df["query_id"].nunique()

        # Apply clustering if requested
        if cluster and self.config.cluster_mecc:
            df = self.cluster_by_signature(df)

        # Calculate match_degree if Gap_Percentage exists
        if "Gap_Percentage" in df.columns:
            df["match_degree"] = (100 - pd.to_numeric(df["Gap_Percentage"], errors="coerce")).round(
                2
            )

        # Extract sequences and add IDs
        df = self.compute_sequences(df)
        df = self.add_numbering_and_export(df)

        # Save outputs
        output_dir.mkdir(parents=True, exist_ok=True)
        output_csv = output_dir / f"{prefix}_MeccDNA_processed.csv"
        df.to_csv(output_csv, index=False)

        if self.fasta_records:
            output_fasta = output_dir / f"{prefix}_MeccDNA_pre.fasta"
            with open(output_fasta, "w") as handle:
                SeqIO.write(self.fasta_records, handle, "fasta")
            self.logger.debug(f"Saved {len(self.fasta_records)} MeccDNA sequences")

        # Summary
        sequences_count = len(self.fasta_records)
        self.logger.info(
            f"MeccDNA: original sequences={initial_query_ids}, after clustering={sequences_count}"
        )

        return df


class CeccProcessor:
    """Process CeccDNA sequences (complex multi-segment chimeras)."""

    def __init__(
        self,
        seq_library: SequenceLibrary,
        config: UMCProcessConfig,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize CeccDNA processor."""
        self.seq_library = seq_library
        self.config = config
        self.logger = logger if logger else get_logger(self.__class__.__name__)
        self.fasta_records: List[SeqRecord] = []
        self.counter = 0

    def generate_cecc_signature(self, df_group: pd.DataFrame) -> str:
        """
        Generate location signature for CeccDNA.
        Signature format: chr1:start1-end1;chr2:start2-end2...
        Sorted by segment_in_circle (the order of segments in the circular DNA)
        """
        if not all(col in df_group.columns for col in ["chr", "start0", "end0"]):
            return ""

        # Sort by segment_in_circle if available
        if "segment_in_circle" in df_group.columns:
            df_sorted = df_group.sort_values("segment_in_circle")
        else:
            df_sorted = df_group

        segments = []
        for _, row in df_sorted.iterrows():
            chrom = row["chr"]
            start = row["start0"]
            end = row["end0"]

            try:
                start = int(start)
                end = int(end)
                segments.append(f"{chrom}:{start}-{end}")
            except (ValueError, TypeError):
                continue

        return ";".join(segments) if segments else ""

    def cluster_by_signature(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Cluster CeccDNA by multi-segment signature.
        Each query_id has multiple rows (segments).
        """
        self.logger.debug("Clustering CeccDNA by segment signatures...")

        if not all(col in df.columns for col in ["chr", "start0", "end0"]):
            self.logger.warning("Location columns not found, skipping clustering")
            df["cecc_signature"] = ""
            df["cluster_id"] = 0
            df["cluster_size"] = 1
            df["num_segments"] = df.groupby("query_id")["query_id"].transform("count")
            return df

        # Generate signatures for each query_id
        signature_map = {}
        segment_counts = {}
        query_id_to_group = {}

        for query_id, group in df.groupby("query_id"):
            signature = self.generate_cecc_signature(group)
            signature_map[query_id] = signature
            segment_counts[query_id] = len(group)
            query_id_to_group[query_id] = group
            df.loc[group.index, "cecc_signature"] = signature
            df.loc[group.index, "num_segments"] = len(group)

        # Group query_ids by signature
        signature_to_queries = defaultdict(list)
        for query_id, signature in signature_map.items():
            if signature:
                signature_to_queries[signature].append(query_id)

        rows_to_keep = []
        cluster_id = 0

        # Process each cluster
        for signature, query_ids in signature_to_queries.items():
            cluster_id += 1
            cluster_size = len(query_ids)

            representative_qid = query_ids[0]
            representative_group = query_id_to_group[representative_qid].copy()

            # Aggregate values for each segment
            for idx, rep_row in representative_group.iterrows():
                chr_match = rep_row["chr"]
                start_match = rep_row["start0"]
                end_match = rep_row["end0"]

                # Aggregate copy_number
                if "copy_number" in representative_group.columns:
                    total_copynum = 0
                    for qid in query_ids:
                        group = query_id_to_group[qid]
                        matching_rows = group[
                            (group["chr"] == chr_match)
                            & (group["start0"] == start_match)
                            & (group["end0"] == end_match)
                        ]
                        if not matching_rows.empty:
                            copynum_val = matching_rows.iloc[0].get("copy_number", 0)
                            if pd.notna(copynum_val):
                                total_copynum += copynum_val

                    if total_copynum > 0:
                        representative_group.loc[idx, "copy_number"] = total_copynum

            representative_group["cluster_id"] = cluster_id
            representative_group["cluster_size"] = cluster_size
            representative_group["cluster_members"] = ";".join(query_ids)

            rows_to_keep.append(representative_group)

        # Handle unclustered sequences
        for query_id, signature in signature_map.items():
            if not signature:
                unclustered_group = query_id_to_group[query_id].copy()
                unclustered_group["cluster_id"] = 0
                unclustered_group["cluster_size"] = 1
                unclustered_group["cluster_members"] = query_id
                rows_to_keep.append(unclustered_group)

        if rows_to_keep:
            result_df = pd.concat(rows_to_keep, ignore_index=True)
        else:
            result_df = df.copy()
            result_df["cluster_id"] = 0
            result_df["cluster_size"] = 1
            result_df["cluster_members"] = result_df["query_id"]

        self.logger.debug(f"CeccDNA clustering: {len(df)} → {len(result_df)} rows")
        return result_df

    def compute_sequences(self, df: pd.DataFrame) -> pd.DataFrame:
        """Extract sequences for all unique query_ids."""
        self.logger.debug("Extracting CeccDNA sequences...")

        df["eSeq"] = ""

        if "length" not in df.columns or "q_start" not in df.columns:
            self.logger.warning("Missing required columns: length and/or q_start")
            return df

        sequences_found = 0
        sequences_missing = 0

        for query_id in df["query_id"].unique():
            query_mask = df["query_id"] == query_id
            first_idx = df[query_mask].index[0]

            q_start = df.loc[first_idx, "q_start"]
            cons_len = df.loc[first_idx, "length"]

            seq_str = self.seq_library.find_sequence(query_id)

            if seq_str and q_start > 0 and cons_len > 0:
                try:
                    extracted_seq = extract_ring_sequence(seq_str, int(q_start), int(cons_len))
                    df.loc[query_mask, "eSeq"] = extracted_seq
                    sequences_found += 1
                except Exception as e:
                    self.logger.debug(f"Failed to extract sequence for {query_id}: {e}")
                    sequences_missing += 1
            else:
                sequences_missing += 1

        self.logger.debug(f"Sequences extracted: {sequences_found}, missing: {sequences_missing}")
        return df

    def add_numbering_and_export(self, df: pd.DataFrame) -> pd.DataFrame:
        """Add CeccDNA IDs and prepare FASTA records."""
        df = df.copy()
        df["eccDNA_id"] = "NA"

        if "eSeq" not in df.columns:
            return df

        # If we have clusters, assign one ID per cluster
        if "cluster_id" in df.columns and df["cluster_id"].max() > 0:
            processed_clusters = set()

            for cluster_id in df[df["cluster_id"] > 0]["cluster_id"].unique():
                if cluster_id in processed_clusters:
                    continue

                cluster_mask = df["cluster_id"] == cluster_id
                cluster_df = df[cluster_mask]

                first_idx = cluster_df.index[0]
                seq = cluster_df.loc[first_idx, "eSeq"]

                if seq and len(seq) > 0:
                    self.counter += 1
                    eccDNA_id = f"C{self.counter}"
                    df.loc[cluster_mask, "eccDNA_id"] = eccDNA_id

                    query_ids = cluster_df["query_id"].unique()
                    num_segments = cluster_df.iloc[0].get("num_segments", 0)

                    if len(query_ids) == 1:
                        description = f"CeccDNA_{query_ids[0]}_segments{num_segments}"
                    else:
                        n_queries = len(query_ids)
                        description = (
                            f"CeccDNA_cluster{cluster_id}_{n_queries}queries_"
                            f"segments{num_segments}"
                        )

                    record = SeqRecord(Seq(seq), id=eccDNA_id, description=description)
                    self.fasta_records.append(record)

                processed_clusters.add(cluster_id)

            # Handle unclustered
            unclustered_mask = df["cluster_id"] == 0
            if unclustered_mask.any():
                for query_id in df[unclustered_mask]["query_id"].unique():
                    query_mask = (df["query_id"] == query_id) & unclustered_mask
                    first_idx = df[query_mask].index[0]
                    seq = df.loc[first_idx, "eSeq"]
                    num_segments = df.loc[first_idx].get("num_segments", 0)

                    if seq and len(seq) > 0:
                        self.counter += 1
                        eccDNA_id = f"C{self.counter}"
                        df.loc[query_mask, "eccDNA_id"] = eccDNA_id

                        record = SeqRecord(
                            Seq(seq),
                            id=eccDNA_id,
                            description=f"CeccDNA_{query_id}_segments{num_segments}",
                        )
                        self.fasta_records.append(record)
        else:
            # No clustering
            for query_id, group in df.groupby("query_id"):
                first_idx = group.index[0]
                seq = group.loc[first_idx, "eSeq"]
                num_segments = len(group)

                if seq and len(seq) > 0:
                    self.counter += 1
                    eccDNA_id = f"C{self.counter}"
                    df.loc[group.index, "eccDNA_id"] = eccDNA_id

                    record = SeqRecord(
                        Seq(seq),
                        id=eccDNA_id,
                        description=f"CeccDNA_{query_id}_segments{num_segments}",
                    )
                    self.fasta_records.append(record)

        return df

    def process(
        self,
        csv_files: Sequence[Path],
        output_dir: Path,
        prefix: str = "sample",
        cluster: bool = True,
    ) -> Optional[pd.DataFrame]:
        """Main processing pipeline for CeccDNA."""
        self.fasta_records = []
        self.counter = 0

        # Load all CSV files
        dataframes: List[pd.DataFrame] = []
        for csv_file in csv_files:
            if not csv_file:
                continue
            csv_path = Path(csv_file)
            if not csv_path.exists() or csv_path.stat().st_size == 0:
                continue
            try:
                df = pd.read_csv(csv_path)
                if not df.empty:
                    dataframes.append(df)
            except Exception as exc:
                self.logger.warning(f"Unable to read {csv_path}: {exc}")

        if not dataframes:
            self.logger.info("No CeccDNA inputs provided")
            return None

        df = pd.concat(dataframes, ignore_index=True)
        initial_query_ids = df["query_id"].nunique()

        # Apply clustering if requested
        if cluster and self.config.cluster_cecc:
            df = self.cluster_by_signature(df)

        # Extract sequences and add IDs
        df = self.compute_sequences(df)
        df = self.add_numbering_and_export(df)

        # Save outputs
        output_dir.mkdir(parents=True, exist_ok=True)
        output_csv = output_dir / f"{prefix}_CeccDNA_processed.csv"
        df.to_csv(output_csv, index=False)

        if self.fasta_records:
            output_fasta = output_dir / f"{prefix}_CeccDNA_pre.fasta"
            with open(output_fasta, "w") as handle:
                SeqIO.write(self.fasta_records, handle, "fasta")
            self.logger.debug(f"Saved {len(self.fasta_records)} CeccDNA sequences")

        # Summary
        sequences_count = len(self.fasta_records)
        self.logger.info(
            f"CeccDNA: original sequences={initial_query_ids}, after clustering={sequences_count}"
        )

        return df


class XeccExporter:
    """Identifies and extracts unclassified sequences (XeccDNA)."""

    def __init__(
        self,
        seq_library: SequenceLibrary,
        config: UMCProcessConfig,
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize XeccDNA exporter."""
        self.seq_library = seq_library
        self.config = config
        self.logger = logger if logger else get_logger(self.__class__.__name__)

    def _get_classified_ids(self, classified_csvs: List[Path]) -> Set[str]:
        """Reads all classified query_ids from the U/M/C CSV files."""
        classified_ids: Set[str] = set()

        for csv_file in classified_csvs:
            if not csv_file.exists():
                continue

            try:
                df = pd.read_csv(csv_file)
                if "query_id" in df.columns:
                    unique_ids = df["query_id"].dropna().unique()
                    classified_ids.update(unique_ids)
            except Exception as e:
                self.logger.warning(f"Error reading {csv_file}: {e}")

        return classified_ids

    def generate(
        self, classified_csvs: Sequence[Path], output_dir: Path, prefix: str
    ) -> Optional[Path]:
        """Generate XeccDNA FASTA from unclassified sequences."""
        all_fasta_ids = self.seq_library.primary_ids
        classified_ids = self._get_classified_ids(classified_csvs)
        xecc_ids = all_fasta_ids - classified_ids

        if not xecc_ids:
            self.logger.info("No unclassified sequences found")
            return None

        self.logger.info(f"Found {len(xecc_ids):,} unclassified sequences (XeccDNA)")

        xecc_records = []
        for i, seq_id in enumerate(sorted(list(xecc_ids)), 1):
            seq_str = self.seq_library.fasta_sequences.get(seq_id)
            if not seq_str:
                continue

            # Take first half of sequence
            half_length = len(seq_str) // 2
            trimmed_seq = seq_str[:half_length] if half_length else seq_str

            new_id = f"{seq_id}__X{i}"

            record = SeqRecord(
                Seq(trimmed_seq),
                id=_sanitize_fasta_id(new_id),
                description=f"XeccDNA_original_id={seq_id}",
            )
            xecc_records.append(record)

        if xecc_records:
            output_dir.mkdir(parents=True, exist_ok=True)
            output_fasta_file = output_dir / f"{prefix}_XeccDNA.fasta"

            with open(output_fasta_file, "w") as handle:
                SeqIO.write(xecc_records, handle, "fasta")

            self.logger.info(
                f"Saved XeccDNA: {output_fasta_file.name} ({len(xecc_records)} sequences)"
            )
            return output_fasta_file

        return None


class UMCProcess:
    """Integrated U/M/C/X processing pipeline."""

    def __init__(
        self, config: Optional[UMCProcessConfig] = None, logger: Optional[logging.Logger] = None
    ):
        """Initialize the UMC process."""
        self.config = config or UMCProcessConfig()
        self.logger = logger if logger else get_logger(self.__class__.__name__)
        self.seq_library = SequenceLibrary(self.logger)

    def process_all(
        self,
        fasta_file: Path,
        uecc_csv: Optional[Union[Path, Sequence[Path]]],
        mecc_csv: Optional[Union[Path, Sequence[Path]]],
        cecc_csv: Optional[Union[Path, Sequence[Path]]],
        output_dir: Path,
        prefix: str = "sample",
    ) -> Dict[str, Union[pd.DataFrame, Path]]:
        """Process all eccDNA types with correct logic."""
        self.logger.info("Starting UMC processing module")

        # Load FASTA sequences
        self.seq_library.load_fasta(fasta_file)

        # Normalize paths
        uecc_files = _coerce_to_paths(uecc_csv)
        mecc_files = _coerce_to_paths(mecc_csv)
        cecc_files = _coerce_to_paths(cecc_csv)

        results: Dict[str, Union[pd.DataFrame, Path]] = {}

        # Process X (unclassified) first
        if self.config.process_xecc:
            classified_csvs = [*uecc_files, *mecc_files, *cecc_files]
            xecc_exporter = XeccExporter(self.seq_library, self.config, self.logger)
            xecc_path = xecc_exporter.generate(classified_csvs, output_dir, prefix)
            if xecc_path:
                results["xecc"] = xecc_path

        # Process U (single location)
        if uecc_files:
            uecc_processor = UeccProcessor(self.seq_library, self.config, self.logger)
            uecc_df = uecc_processor.process(
                uecc_files, output_dir, prefix, cluster=self.config.cluster_uecc
            )
            if uecc_df is not None:
                results["uecc"] = uecc_df

        # Process M (multiple locations)
        if mecc_files:
            mecc_processor = MeccProcessor(self.seq_library, self.config, self.logger)
            mecc_df = mecc_processor.process(
                mecc_files, output_dir, prefix, cluster=self.config.cluster_mecc
            )
            if mecc_df is not None:
                results["mecc"] = mecc_df

        # Process C (complex segments)
        if cecc_files:
            cecc_processor = CeccProcessor(self.seq_library, self.config, self.logger)
            cecc_df = cecc_processor.process(
                cecc_files, output_dir, prefix, cluster=self.config.cluster_cecc
            )
            if cecc_df is not None:
                results["cecc"] = cecc_df

        # Summary
        processed_keys = ", ".join(results.keys()) or "none"
        self.logger.info(f"Module complete. Processed types: {processed_keys}")

        return results

    def run(
        self,
        fasta_file: Optional[Path] = None,
        uecc_csv: Optional[Path] = None,
        mecc_csv: Optional[Path] = None,
        cecc_csv: Optional[Path] = None,
        output_dir: Optional[Path] = None,
        prefix: str = "sample",
    ) -> Dict[str, Union[pd.DataFrame, Path]]:
        """Run the UMC processing pipeline with auto-discovery of input files."""

        # Auto-discover file paths if not provided
        if output_dir is None:
            output_dir = Path.cwd()

        if fasta_file is None:
            # Look for common FASTA file names
            possible_fasta = [
                output_dir / "tandem_to_ring.fasta",
                output_dir / f"{prefix}_circular.fasta",
                output_dir / f"{prefix}_processed.fasta",
            ]
            for fasta_path in possible_fasta:
                if fasta_path.exists():
                    fasta_file = fasta_path
                    break

            if fasta_file is None:
                raise FileNotFoundError("No FASTA file found. Please specify fasta_file parameter.")

        if uecc_csv is None:
            uecc_path = output_dir / "um_classify.uecc.csv"
            if uecc_path.exists():
                uecc_csv = uecc_path

        if mecc_csv is None:
            mecc_path = output_dir / "um_classify.mecc.csv"
            if mecc_path.exists():
                mecc_csv = mecc_path

        if cecc_csv is None:
            cecc_path = output_dir / "cecc_build.csv"
            if cecc_path.exists():
                cecc_csv = cecc_path

        self.logger.info("Running UMC process with:")
        self.logger.info(f"  FASTA file: {fasta_file}")
        self.logger.info(f"  UECC CSV: {uecc_csv or 'None'}")
        self.logger.info(f"  MECC CSV: {mecc_csv or 'None'}")
        self.logger.info(f"  CECC CSV: {cecc_csv or 'None'}")
        self.logger.info(f"  Output dir: {output_dir}")

        # Call the main processing method
        return self.process_all(
            fasta_file=fasta_file,
            uecc_csv=uecc_csv,
            mecc_csv=mecc_csv,
            cecc_csv=cecc_csv,
            output_dir=output_dir,
            prefix=prefix,
        )


def main():
    """Command-line interface."""
    import argparse

    parser = argparse.ArgumentParser(
        description="UMC Process: Process eccDNA files with correct clustering logic"
    )

    parser.add_argument(
        "-f", "--fasta", type=Path, required=True, help="Master FASTA file with all sequences"
    )
    parser.add_argument("-u", "--uecc", type=Path, help="UeccDNA CSV file")
    parser.add_argument("-m", "--mecc", type=Path, help="MeccDNA CSV file")
    parser.add_argument("-c", "--cecc", type=Path, help="CeccDNA CSV file")
    parser.add_argument(
        "-o", "--output", type=Path, required=True, help="Output directory for results"
    )
    parser.add_argument(
        "-p", "--prefix", type=str, default="sample", help="Prefix for output files"
    )
    parser.add_argument("--no-cluster", action="store_true", help="Skip clustering for all types")
    parser.add_argument("--debug", action="store_true", help="Enable debug logging")

    args = parser.parse_args()

    # Setup logging
    log_level = logging.DEBUG if args.debug else logging.INFO
    setup_logging(level=log_level)
    logger = get_logger("umc_process")

    # Validate inputs
    if not args.fasta.exists():
        logger.error(f"FASTA file not found: {args.fasta}")
        sys.exit(1)

    # Create config
    config = UMCProcessConfig(
        process_xecc=True,
        cluster_uecc=not args.no_cluster,
        cluster_cecc=not args.no_cluster,
        cluster_mecc=not args.no_cluster,
        debug=args.debug,
    )

    try:
        # Run pipeline
        processor = UMCProcess(config, logger)
        processor.process_all(
            fasta_file=args.fasta,
            uecc_csv=args.uecc,
            mecc_csv=args.mecc,
            cecc_csv=args.cecc,
            output_dir=args.output,
            prefix=args.prefix,
        )

        logger.info(f"Processing completed. Results: {args.output}")

    except Exception as e:
        logger.error(f"Processing failed: {e}")
        if args.debug:
            import traceback

            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
