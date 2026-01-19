"""
Read Filter - FASTA Sequence Filtering Module

This module filters FASTA sequences based on TandemToRing (formerly Carousel)
classification results, removing sequences classified as CtcReads variants
(legacy label: CtcR-* in `tandem_to_ring.csv`).

Key features:
- Reads TandemToRing output CSV with readName and readClass columns
- Filters out sequences with specific CtcR-* classifications (CtcReads variants)
- Memory-efficient processing for large files
- Generates samtools faidx index for output

Migrated from step9_sieve.py to the CircleSeeker architecture.
"""

from __future__ import annotations

import csv
import sys
import logging
from circleseeker.utils.logging import get_logger
import subprocess
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Optional
from circleseeker.exceptions import PipelineError

# Increase CSV field size limit to handle large sequence fields
csv.field_size_limit(sys.maxsize)


@dataclass
class FilterStats:
    """Statistics for filtering operation."""

    total_reads: int = 0
    filtered_reads: int = 0
    retained_reads: int = 0
    csv_total_reads: int = 0
    csv_ctcr_reads: int = 0

    @property
    def filtered_percentage(self) -> float:
        if self.total_reads == 0:
            return 0.0
        return (self.filtered_reads / self.total_reads) * 100

    @property
    def retained_percentage(self) -> float:
        if self.total_reads == 0:
            return 0.0
        return (self.retained_reads / self.total_reads) * 100


class Sieve:
    """Filter FASTA sequences based on TandemToRing classification."""

    # Default CtcR classes to filter out
    DEFAULT_CTCR_CLASSES = {"CtcR-perfect", "CtcR-inversion", "CtcR-hybrid"}

    def __init__(self, logger: Optional[logging.Logger] = None) -> None:
        """Initialize Sieve filter."""
        self.logger = logger or get_logger(self.__class__.__name__)

        # Data structures
        self.reads_to_filter: set[str] = set()

        # Statistics
        self.stats = FilterStats()

        # Check for samtools availability
        self.has_samtools = shutil.which("samtools") is not None
        if not self.has_samtools:
            self.logger.warning("samtools not found - faidx index will not be generated")

    def load_tandem_to_ring_classification(
        self, csv_file: Path, ctcr_classes: Optional[set[str]] = None
    ) -> None:
        """
        Load read classification from TandemToRing output CSV.

        Args:
            csv_file: Path to TandemToRing output CSV (tandem_to_ring.csv)
            ctcr_classes: Set of classes to filter out (default: CtcR variants)
        """
        if not csv_file.exists():
            raise FileNotFoundError(f"Classification CSV not found: {csv_file}")

        ctcr_classes = ctcr_classes or self.DEFAULT_CTCR_CLASSES
        self.logger.info(f"Loading classifications from: {csv_file}")
        self.logger.info(f"Filtering classes: {ctcr_classes}")

        # Track statistics
        class_counts: dict[str, int] = {}

        try:
            with open(csv_file, "r", encoding="utf-8") as f:
                # Try to detect delimiter
                first_line = f.readline()
                f.seek(0)

                if "\t" in first_line:
                    delimiter = "\t"
                else:
                    delimiter = ","

                reader = csv.DictReader(f, delimiter=delimiter)

                # Find the correct columns (case-insensitive)
                if not reader.fieldnames:
                    raise ValueError("Classification CSV is missing a header row")
                fieldnames_lower = {fn.lower(): fn for fn in reader.fieldnames}

                if "readname" not in fieldnames_lower or "readclass" not in fieldnames_lower:
                    raise ValueError(
                        f"CSV must contain 'readName' and 'readClass' columns. "
                        f"Found: {reader.fieldnames}"
                    )

                readname_col = fieldnames_lower["readname"]
                readclass_col = fieldnames_lower["readclass"]

                # Process each row
                for row in reader:
                    self.stats.csv_total_reads += 1
                    read_name = row[readname_col].strip()
                    read_class = row[readclass_col].strip()

                    # Count by class
                    class_counts[read_class] = class_counts.get(read_class, 0) + 1

                    # Check if this read should be filtered
                    if read_class in ctcr_classes:
                        self.reads_to_filter.add(read_name)
                        self.stats.csv_ctcr_reads += 1

        except Exception as e:
            raise PipelineError(f"Error loading classification CSV: {e}")

        # Log statistics
        self.logger.info(f"Loaded {self.stats.csv_total_reads} classifications")
        self.logger.info("Class distribution:")
        for cls, count in sorted(class_counts.items()):
            if cls in ctcr_classes:
                self.logger.info(f"  {cls}: {count} (will filter)")
            else:
                self.logger.info(f"  {cls}: {count}")

        self.logger.info(f"Total CtcR reads to filter: {self.stats.csv_ctcr_reads}")

    def filter_fasta_files(
        self, input_fastas: list[Path], output_fasta: Path, generate_index: bool = True
    ) -> FilterStats:
        """
        Filter and combine multiple FASTA files.

        Args:
            input_fastas: List of input FASTA files to filter
            output_fasta: Path for combined filtered output
            generate_index: Whether to generate samtools faidx index

        Returns:
            FilterStats object with filtering statistics
        """
        output_fasta = Path(output_fasta)
        output_fasta.parent.mkdir(parents=True, exist_ok=True)

        self.logger.info(f"Filtering {len(input_fastas)} FASTA file(s)...")

        # Open output file
        with open(output_fasta, "w") as outfile:
            for input_fasta in input_fastas:
                if not input_fasta.exists():
                    self.logger.warning(f"Input file not found, skipping: {input_fasta}")
                    continue

                self.logger.info(f"Processing: {input_fasta}")

                # Process each input file
                with open(input_fasta, "r") as infile:
                    current_header: Optional[str] = None
                    current_sequence: list[str] = []
                    write_current = True

                    for line in infile:
                        line = line.rstrip("\n")

                        if line.startswith(">"):
                            # Write previous sequence if needed
                            if current_header is not None:
                                self.stats.total_reads += 1
                                if write_current:
                                    outfile.write(f"{current_header}\n")
                                    outfile.write("\n".join(current_sequence) + "\n")
                                    self.stats.retained_reads += 1
                                else:
                                    self.stats.filtered_reads += 1

                            # Start new sequence
                            current_header = line
                            current_sequence = []

                            # Extract read ID (first part after >)
                            read_id = line[1:].split()[0]

                            # Check if should be filtered
                            write_current = read_id not in self.reads_to_filter

                        elif line:  # Sequence line
                            current_sequence.append(line)

                    # Write last sequence
                    if current_header is not None:
                        self.stats.total_reads += 1
                        if write_current:
                            outfile.write(f"{current_header}\n")
                            outfile.write("\n".join(current_sequence) + "\n")
                            self.stats.retained_reads += 1
                        else:
                            self.stats.filtered_reads += 1

        # Generate samtools index if requested
        if generate_index and self.has_samtools:
            self.generate_faidx_index(output_fasta)

        # Log final statistics
        self.logger.info("Filtering complete:")
        self.logger.info(f"  Total sequences: {self.stats.total_reads}")
        self.logger.info(
            f"  Filtered out: {self.stats.filtered_reads} ({self.stats.filtered_percentage:.2f}%)"
        )
        self.logger.info(
            f"  Retained: {self.stats.retained_reads} ({self.stats.retained_percentage:.2f}%)"
        )
        self.logger.info(f"  CtcR classifications: {self.stats.csv_ctcr_reads}")
        self.logger.info(f"  Output: {output_fasta}")

        return self.stats

    def generate_faidx_index(self, fasta_file: Path) -> bool:
        """
        Generate samtools faidx index for FASTA file.

        Args:
            fasta_file: Path to FASTA file

        Returns:
            True if successful, False otherwise
        """
        try:
            cmd = ["samtools", "faidx", str(fasta_file)]
            subprocess.run(cmd, capture_output=True, text=True, check=True)

            fai_file = Path(f"{fasta_file}.fai")
            if fai_file.exists():
                self.logger.info(f"Generated index: {fai_file}")
                return True
            else:
                self.logger.warning("Index file was not created")
                return False

        except subprocess.CalledProcessError as e:
            self.logger.error(f"Failed to generate index: {e}")
            return False
        except FileNotFoundError:
            self.logger.warning("samtools not found, skipping index generation")
            return False

    def run_sieve(
        self,
        tandem_to_ring_csv: Path,
        input_fastas: list[Path],
        output_fasta: Path,
        ctcr_classes: Optional[set[str]] = None,
    ) -> FilterStats:
        """
        Run complete sieve filtering pipeline.

        Args:
            tandem_to_ring_csv: Path to TandemToRing classification CSV (tandem_to_ring.csv)
            input_fastas: List of FASTA files to filter (from previous steps)
            output_fasta: Path for output filtered FASTA
            ctcr_classes: Optional set of classes to filter (default: CtcR variants)

        Returns:
            Filtering statistics
        """
        self.logger.info("=" * 60)
        self.logger.info("Sieve - FASTA Filtering Pipeline")
        self.logger.info("=" * 60)

        try:
            # Step 1: Load classifications
            self.load_tandem_to_ring_classification(tandem_to_ring_csv, ctcr_classes)

            # Step 2: Filter FASTA files
            stats = self.filter_fasta_files(input_fastas, output_fasta)

            # Step 3: Summary
            self.logger.info("=" * 60)
            self.logger.info("Filtering Complete")
            self.logger.info(f"Output file: {output_fasta}")
            if Path(f"{output_fasta}.fai").exists():
                self.logger.info(f"Index file: {output_fasta}.fai")
            self.logger.info("=" * 60)

            return stats

        except Exception as e:
            self.logger.error(f"Sieve pipeline failed: {e}")
            raise PipelineError(f"Sieve filtering failed: {e}")
