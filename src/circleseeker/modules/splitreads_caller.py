"""SplitReads-Caller: Built-in split-reads based eccDNA detection module.

This module provides a fully internal implementation for detecting eccDNA from
long-read sequencing data using split-read alignment patterns. It replaces the
external Cresil tool dependency.

Algorithm overview:
1. Trim phase: Align reads to reference using mappy, detect CTC (Circular Tandem
   Copy) patterns where the same genomic region appears multiple times.
2. Identify phase: Build a graph from aligned segments, detect circular subgraphs
   using networkx, and validate eccDNA candidates.

The output format is compatible with the existing Cresil adapter, allowing seamless
integration with the CircleSeeker pipeline.
"""

from __future__ import annotations

from bisect import bisect_right
import logging
import os
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterator, Optional

import mappy as mp
import networkx as nx
import numpy as np
import pandas as pd
import pybedtools as bt
import tempfile

from circleseeker.config import SplitReadsConfig
from circleseeker.utils.logging import get_logger


# ============================================================================
# Data Classes
# ============================================================================


@dataclass
class AlignmentSegment:
    """Represents a single alignment segment from mappy."""

    read_name: str
    read_start: int  # 0-based start on query
    read_end: int  # 0-based end on query
    chrom: str
    ref_start: int  # 0-based start on reference
    ref_end: int  # 0-based end on reference
    strand: str  # '+' or '-'
    mapq: int
    cigar: str
    nm: int  # edit distance

    @property
    def read_length(self) -> int:
        return self.read_end - self.read_start

    @property
    def ref_length(self) -> int:
        return self.ref_end - self.ref_start


@dataclass
class TrimResult:
    """Result from the trim phase for a single read."""

    read_name: str
    segments: list[AlignmentSegment]
    is_ctc: bool  # Circular Tandem Copy pattern detected
    ctc_regions: list[tuple[str, int, int]]  # Repeated genomic regions


@dataclass
class MergedRegion:
    """A merged genomic region from multiple overlapping alignments."""

    merge_id: int
    chrom: str
    start: int
    end: int
    strand: str
    supporting_reads: set[str] = field(default_factory=set)
    coverage_depth: float = 0.0
    has_5prime_break: bool = False
    has_3prime_break: bool = False

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass
class EccDNACandidate:
    """An eccDNA candidate detected from split-read patterns."""

    ecc_id: str
    regions: list[MergedRegion]
    is_circular: bool
    is_ctc: bool
    num_reads: int
    total_base: int
    coverage: float

    @property
    def merge_region_str(self) -> str:
        """Format regions as Cresil-compatible string."""
        parts = []
        for r in self.regions:
            parts.append(f"{r.chrom}:{r.start}-{r.end}_{r.strand}")
        return ";".join(parts)

    @property
    def merge_len(self) -> int:
        return sum(r.length for r in self.regions)


# ============================================================================
# Interval Coverage (Pure Python Scanline Algorithm)
# ============================================================================


class IntervalCoverage:
    """Calculate coverage depth using a scanline algorithm.

    This is a pure Python implementation that replaces bedtools genomecov.
    It uses an event-based approach to efficiently compute coverage.
    """

    def __init__(self) -> None:
        self._events: dict[str, list[tuple[int, int]]] = defaultdict(list)

    def add_interval(self, chrom: str, start: int, end: int) -> None:
        """Add an interval to track."""
        self._events[chrom].append((start, 1))  # +1 at start
        self._events[chrom].append((end, -1))  # -1 at end

    def calculate(self) -> dict[str, list[tuple[int, int, int]]]:
        """Calculate coverage for all intervals.

        Returns:
            Dict mapping chrom -> list of (start, end, depth) tuples
        """
        result: dict[str, list[tuple[int, int, int]]] = {}

        for chrom, events in self._events.items():
            if not events:
                continue

            # Sort by position, then by event type (-1 before +1 for same pos)
            events.sort(key=lambda x: (x[0], -x[1]))

            coverage_spans: list[tuple[int, int, int]] = []
            current_depth = 0
            prev_pos = events[0][0]

            for pos, delta in events:
                if pos > prev_pos and current_depth > 0:
                    coverage_spans.append((prev_pos, pos, current_depth))
                current_depth += delta
                prev_pos = pos

            result[chrom] = coverage_spans

        return result

    def get_regions_above_threshold(
        self, min_depth: float
    ) -> list[tuple[str, int, int, float]]:
        """Get regions with coverage >= min_depth.

        Returns:
            List of (chrom, start, end, avg_depth) tuples
        """
        coverage = self.calculate()
        regions: list[tuple[str, int, int, float]] = []

        for chrom, spans in coverage.items():
            # Merge consecutive spans above threshold
            current_start: Optional[int] = None
            current_end: Optional[int] = None
            total_depth = 0.0
            total_len = 0

            for start, end, depth in spans:
                if depth >= min_depth:
                    if current_start is None:
                        current_start = start
                        current_end = end
                        total_depth = depth * (end - start)
                        total_len = end - start
                    elif start <= current_end:  # Overlapping/adjacent
                        current_end = max(current_end, end)
                        total_depth += depth * (end - start)
                        total_len += end - start
                    else:
                        # Gap found, emit current region
                        avg = total_depth / total_len if total_len > 0 else 0
                        regions.append((chrom, current_start, current_end, avg))
                        current_start = start
                        current_end = end
                        total_depth = depth * (end - start)
                        total_len = end - start
                else:
                    if current_start is not None:
                        avg = total_depth / total_len if total_len > 0 else 0
                        regions.append((chrom, current_start, current_end, avg))
                        current_start = None
                        current_end = None
                        total_depth = 0.0
                        total_len = 0

            # Emit last region
            if current_start is not None and current_end is not None:
                avg = total_depth / total_len if total_len > 0 else 0
                regions.append((chrom, current_start, current_end, avg))

        return regions


# ============================================================================
# Trim Processor
# ============================================================================


class TrimProcessor:
    """Process reads using mappy alignment to detect CTC patterns.

    The trim phase aligns reads to the reference genome and identifies reads
    that show Circular Tandem Copy (CTC) patterns - where the same genomic
    region appears multiple times in a single read.
    """

    # mappy flag to disable long-join (similar to Cresil behavior)
    # Correct value is 0x400 (matches minimap2's MM_F_NO_LJOIN definition)
    MM_F_NO_LJOIN = 0x400

    def __init__(
        self,
        reference_path: Path,
        config: SplitReadsConfig,
        logger: Optional[logging.Logger] = None,
        threads: int = 1,
    ) -> None:
        """Initialize the TrimProcessor.

        Args:
            reference_path: Path to reference FASTA or MMI index
            config: SplitReads configuration
            logger: Logger instance
            threads: Number of threads for alignment
        """
        self.config = config
        self.logger = logger or get_logger(self.__class__.__name__)
        self.threads = threads

        # Parse excluded chromosomes
        self.exclude_chrs: set[str] = set()
        if config.exclude_chrs:
            self.exclude_chrs = {
                c.strip() for c in config.exclude_chrs.split(",") if c.strip()
            }
            self.logger.info(f"Excluding chromosomes: {self.exclude_chrs}")

        # Initialize mappy aligner with extra_flags to control behavior
        self.logger.info(f"Loading reference: {reference_path}")
        self.aligner = mp.Aligner(
            str(reference_path),
            preset="map-hifi",
            n_threads=threads,
            extra_flags=self.MM_F_NO_LJOIN,  # Disable long-join
        )
        if not self.aligner:
            raise RuntimeError(f"Failed to load reference: {reference_path}")
        self.logger.info("Reference loaded successfully")

        # Extract chromosome sizes from aligner
        self.chrom_sizes: dict[str, int] = {}
        if hasattr(self.aligner, "seq_names") and hasattr(self.aligner, "seq"):
            for seq_name in self.aligner.seq_names:
                seq = self.aligner.seq(seq_name)
                if seq:
                    self.chrom_sizes[seq_name] = len(seq)
            self.logger.info(f"Extracted {len(self.chrom_sizes)} chromosome sizes from reference")

    def process_read(
        self, read_name: str, sequence: str
    ) -> Optional[TrimResult]:
        """Process a single read and detect CTC patterns.

        Args:
            read_name: Name of the read
            sequence: Read sequence

        Returns:
            TrimResult if valid alignments found, None otherwise
        """
        if len(sequence) < 100:
            return None

        # Align read to reference
        alignments = list(self.aligner.map(sequence))
        if not alignments:
            return None

        # Filter by MAPQ and excluded chromosomes
        segments: list[AlignmentSegment] = []
        for hit in alignments:
            # Match Cresil behavior: keep only primary alignments.
            if getattr(hit, "is_primary", True) is False:
                continue
            if hit.mapq < self.config.mapq_threshold:
                continue
            if hit.ctg in self.exclude_chrs:
                continue

            seg = AlignmentSegment(
                read_name=read_name,
                read_start=hit.q_st,
                read_end=hit.q_en,
                chrom=hit.ctg,
                ref_start=hit.r_st,
                ref_end=hit.r_en,
                strand="+" if hit.strand == 1 else "-",
                mapq=hit.mapq,
                cigar=hit.cigar_str if hit.cigar_str else "",
                nm=hit.NM if hasattr(hit, "NM") else 0,
            )
            segments.append(seg)

        if not segments:
            return None

        # Sort segments by read position
        segments.sort(key=lambda s: s.read_start)

        # Merge adjacent alignments
        segments = self._merge_adjacent_alignments(segments)

        # Check for CTC pattern
        is_ctc, ctc_regions = self._check_ctc_pattern(segments)

        return TrimResult(
            read_name=read_name,
            segments=segments,
            is_ctc=is_ctc,
            ctc_regions=ctc_regions,
        )

    def _merge_adjacent_alignments(
        self, segments: list[AlignmentSegment]
    ) -> list[AlignmentSegment]:
        """Merge adjacent alignment segments with small gaps/overlaps.

        Args:
            segments: List of alignment segments sorted by read position

        Returns:
            Merged list of segments
        """
        if len(segments) <= 1:
            return segments

        merged: list[AlignmentSegment] = []
        current = segments[0]

        for next_seg in segments[1:]:
            # Check if segments can be merged
            gap = next_seg.read_start - current.read_end
            # Also ensure reference coordinates are adjacent on the reference.
            # Otherwise, two distant loci that appear adjacent on the read (e.g., split alignments)
            # could be incorrectly merged into a huge genomic span.
            current_ref_end = current.ref_end if current.strand == "+" else current.ref_start
            next_ref_start = next_seg.ref_start if next_seg.strand == "+" else next_seg.ref_end
            ref_gap = next_ref_start - current_ref_end

            # Same chromosome and strand, small gap/overlap
            if (
                current.chrom == next_seg.chrom
                and current.strand == next_seg.strand
                and -self.config.overlap_tolerance <= gap <= self.config.gap_tolerance
                and -self.config.overlap_tolerance <= ref_gap <= self.config.gap_tolerance
            ):
                # Merge: extend current segment
                current = AlignmentSegment(
                    read_name=current.read_name,
                    read_start=current.read_start,
                    read_end=max(current.read_end, next_seg.read_end),
                    chrom=current.chrom,
                    ref_start=min(current.ref_start, next_seg.ref_start),
                    ref_end=max(current.ref_end, next_seg.ref_end),
                    strand=current.strand,
                    mapq=max(current.mapq, next_seg.mapq),
                    cigar="",  # Merged cigar is complex, skip
                    nm=current.nm + next_seg.nm,
                )
            else:
                merged.append(current)
                current = next_seg

        merged.append(current)
        return merged

    def _check_ctc_pattern(
        self, segments: list[AlignmentSegment]
    ) -> tuple[bool, list[tuple[str, int, int]]]:
        """Check if segments show a CTC (Circular Tandem Copy) pattern.

        CTC pattern: The same genomic region appears multiple times in the read,
        indicating the read traversed the circular DNA multiple times.

        Args:
            segments: List of alignment segments

        Returns:
            Tuple of (is_ctc, list of repeated regions)
        """
        if len(segments) < 2:
            return False, []

        # Group segments by genomic region (chrom:start-end with tolerance)
        region_groups: dict[str, list[AlignmentSegment]] = defaultdict(list)
        tolerance = 500  # bp tolerance for grouping

        for seg in segments:
            # Create region key with rounded coordinates
            start_round = (seg.ref_start // tolerance) * tolerance
            end_round = ((seg.ref_end + tolerance - 1) // tolerance) * tolerance
            key = f"{seg.chrom}:{start_round}-{end_round}"
            region_groups[key].append(seg)

        # Find repeated regions
        ctc_regions: list[tuple[str, int, int]] = []
        for key, group in region_groups.items():
            if len(group) >= 2:
                # This region appears multiple times - CTC pattern
                chrom = group[0].chrom
                start = min(s.ref_start for s in group)
                end = max(s.ref_end for s in group)
                ctc_regions.append((chrom, start, end))

        is_ctc = len(ctc_regions) > 0
        return is_ctc, ctc_regions

    def process_fasta(
        self, fasta_path: Path
    ) -> Iterator[TrimResult]:
        """Process all reads in a FASTA file.

        Args:
            fasta_path: Path to input FASTA file

        Yields:
            TrimResult for each read with valid alignments
        """
        self.logger.info(f"Processing FASTA: {fasta_path}")
        processed = 0
        valid = 0

        for name, seq, _ in mp.fastx_read(str(fasta_path)):
            processed += 1
            result = self.process_read(name, seq)
            if result:
                valid += 1
                yield result

            if processed % 10000 == 0:
                self.logger.debug(f"Processed {processed} reads, {valid} with valid alignments")

        self.logger.info(f"Processed {processed} reads, {valid} with valid alignments")

    def write_trim_tsv(
        self, results: list[TrimResult], output_path: Path
    ) -> None:
        """Write trim results to TSV file.

        Args:
            results: List of trim results
            output_path: Output TSV path
        """
        rows = []
        for result in results:
            for seg in result.segments:
                rows.append({
                    "read_name": seg.read_name,
                    "read_start": seg.read_start,
                    "read_end": seg.read_end,
                    "chrom": seg.chrom,
                    "ref_start": seg.ref_start,
                    "ref_end": seg.ref_end,
                    "strand": seg.strand,
                    "mapq": seg.mapq,
                    "is_ctc": result.is_ctc,
                })

        df = pd.DataFrame(rows)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_path, sep="\t", index=False)
        self.logger.info(f"Wrote trim results to {output_path}")


# ============================================================================
# Identify Processor
# ============================================================================


class IdentifyProcessor:
    """Identify eccDNA candidates from trim results using graph-based detection.

    This processor:
    1. Calculates coverage depth across all aligned regions
    2. Filters low-coverage regions
    3. Merges overlapping regions
    4. Detects breakpoints (5'/3' ends)
    5. Builds an adjacency graph
    6. Finds circular subgraphs using networkx
    """

    def __init__(
        self,
        config: SplitReadsConfig,
        logger: Optional[logging.Logger] = None,
        chrom_sizes: Optional[dict[str, int]] = None,
    ) -> None:
        """Initialize the IdentifyProcessor.

        Args:
            config: SplitReads configuration
            logger: Logger instance
            chrom_sizes: Dictionary mapping chromosome names to sizes (required for bedtools)
        """
        self.config = config
        self.logger = logger or get_logger(self.__class__.__name__)
        self.chrom_sizes = chrom_sizes or {}

    def identify(
        self, trim_results: list[TrimResult]
    ) -> list[EccDNACandidate]:
        """Identify eccDNA candidates from trim results.

        Args:
            trim_results: List of trim results from TrimProcessor

        Returns:
            List of eccDNA candidates
        """
        if not trim_results:
            return []

        # Separate all reads (for coverage) from split-reads (for graph inference)
        all_results = trim_results
        split_results = [r for r in trim_results if len(r.segments) >= 2]

        # Determine which reads to use for coverage calculation based on config
        coverage_read_mode = getattr(self.config, "coverage_read_mode", "all_reads")
        if coverage_read_mode == "all_reads":
            coverage_results = all_results
            self.logger.info(
                f"Coverage calculation using all {len(all_results)} reads "
                f"({len(split_results)} split-reads, {len(all_results) - len(split_results)} single-segment)"
            )
        else:
            # Legacy behavior: only use split-reads
            coverage_results = split_results
            skipped = len(all_results) - len(split_results)
            if skipped:
                self.logger.info(
                    "Coverage calculation: %d/%d single-segment reads ignored; %d split-read reads used",
                    skipped,
                    len(all_results),
                    len(split_results),
                )

        if not split_results:
            self.logger.info("No split-read alignments (>=2 segments) available; skipping inference candidates")
            return []

        self.logger.info(f"Identifying eccDNA from {len(split_results)} split-read results")

        # Step 1: Calculate coverage using configured strategy
        merge_strategy = getattr(self.config, "merge_strategy", "cresil")
        if merge_strategy == "cresil":
            regions = self._calculate_coverage_cresil_style(coverage_results)
            self.logger.info(f"Found {len(regions)} merged regions using Cresil-style merge")
        else:
            regions = self._calculate_coverage_bedtools(coverage_results)
            self.logger.info(f"Found {len(regions)} merged regions from bedtools")

        if not regions:
            return []

        # Step 2: Filter by minimum depth (already filtered in cresil_style, but double-check)
        regions = [r for r in regions if r[3] >= self.config.min_avg_depth]
        self.logger.info(f"After depth filter: {len(regions)} regions with depth >= {self.config.min_avg_depth}")

        if not regions:
            return []

        # Step 3: Filter by minimum size
        regions = [
            r for r in regions
            if (r[2] - r[1]) >= self.config.min_region_size
        ]
        self.logger.info(f"After size filter: {len(regions)} regions >= {self.config.min_region_size}bp")

        if not regions:
            return []

        # Step 4: Create merged regions with metadata (use all reads for supporting_reads)
        merged_regions = self._create_merged_regions(regions, all_results)

        # Step 5: Annotate breakpoints (use split-reads only for breakpoint evidence)
        self._annotate_breakpoints(merged_regions, split_results)

        # Step 6: Build adjacency graph and find circles (use split-reads only)
        candidates = self._find_circular_candidates(merged_regions, split_results)

        self.logger.info(f"Identified {len(candidates)} eccDNA candidates")
        return candidates

    def _calculate_coverage(self, trim_results: list[TrimResult]) -> IntervalCoverage:
        """Calculate raw per-base coverage using a scanline approach.

        This is primarily used as a robust fallback when bedtools/genome sizes
        are unavailable, and also makes unit testing easier.
        """
        coverage = IntervalCoverage()
        for result in trim_results:
            for seg in result.segments:
                coverage.add_interval(seg.chrom, seg.ref_start, seg.ref_end)
        return coverage

    def _calculate_coverage_scanline_regions(
        self, trim_results: list[TrimResult]
    ) -> list[tuple[str, int, int, float]]:
        """Fallback: merge scanline coverage spans into regions with avg depth."""
        coverage = self._calculate_coverage(trim_results)
        spans_by_chrom = coverage.calculate()

        regions: list[tuple[str, int, int, float]] = []
        for chrom, spans in spans_by_chrom.items():
            if not spans:
                continue

            current_start, current_end, depth = spans[0]
            total_depth = float(depth) * (current_end - current_start)
            total_len = current_end - current_start

            for start, end, d in spans[1:]:
                if start <= current_end:  # overlapping/adjacent
                    current_end = max(current_end, end)
                    total_depth += float(d) * (end - start)
                    total_len += end - start
                else:
                    avg = total_depth / total_len if total_len > 0 else 0.0
                    regions.append((chrom, current_start, current_end, avg))
                    current_start, current_end = start, end
                    total_depth = float(d) * (end - start)
                    total_len = end - start

            avg = total_depth / total_len if total_len > 0 else 0.0
            regions.append((chrom, current_start, current_end, avg))

        return regions

    def _calculate_coverage_bedtools(
        self, trim_results: list[TrimResult]
    ) -> list[tuple[str, int, int, float]]:
        """Calculate coverage using bedtools merge and genomecov.

        This method uses pybedtools to:
        1. Merge only adjacent/overlapping regions (unlike the scanline approach)
        2. Calculate accurate coverage depth using bedtools genomecov

        Args:
            trim_results: List of trim results

        Returns:
            List of (chrom, start, end, avg_depth) tuples
        """
        # Collect all alignment segments into BED format
        bed_rows = []
        for result in trim_results:
            for seg in result.segments:
                bed_rows.append({
                    "chrom": seg.chrom,
                    "start": seg.ref_start,
                    "end": seg.ref_end,
                    "name": seg.read_name,
                    "score": seg.mapq,
                    "strand": seg.strand,
                })

        if not bed_rows:
            return []

        # Create DataFrame and BedTool
        df = pd.DataFrame(bed_rows)
        df = df.sort_values(["chrom", "start", "end"])

        try:
            aln_bed = bt.BedTool.from_dataframe(df).sort()

            # Step 1: Merge adjacent/overlapping regions
            # This is the KEY difference from the scanline approach:
            # bedtools merge only combines regions that actually overlap or are adjacent
            merged_bed = aln_bed.merge()

            # Step 2: Calculate coverage depth using genomecov
            # We need chrom sizes for genomecov
            genome_file = None
            chrom_sizes = self.chrom_sizes
            if not chrom_sizes:
                # Infer minimal chrom sizes from alignment ends for testability/fallback.
                inferred = df.groupby("chrom")["end"].max()
                chrom_sizes = {str(chrom): int(end) for chrom, end in inferred.items() if pd.notna(end)}

            if chrom_sizes:
                with tempfile.NamedTemporaryFile(mode="w", suffix=".genome", delete=False) as f:
                    for chrom, size in chrom_sizes.items():
                        f.write(f"{chrom}\t{int(size)}\n")
                    genome_file = f.name
                genome_cov = aln_bed.genome_coverage(bg=True, g=genome_file)
            else:
                # Extremely defensive fallback; should not happen for non-empty inputs.
                genome_cov = aln_bed.genome_coverage(bg=True)

            # Convert to DataFrames
            merged_df = merged_bed.to_dataframe(names=["chrom", "start", "end"])
            cov_df = genome_cov.to_dataframe(names=["chrom", "start", "end", "depth"])

            # Step 3: For each merged region, calculate average depth
            regions: list[tuple[str, int, int, float]] = []

            for _, row in merged_df.iterrows():
                chrom = row["chrom"]
                start = row["start"]
                end = row["end"]

                # Find coverage entries that overlap this merged region
                mask = (
                    (cov_df["chrom"] == chrom)
                    & (cov_df["start"] < end)
                    & (cov_df["end"] > start)
                )
                overlapping = cov_df[mask]

                if overlapping.empty:
                    regions.append((chrom, int(start), int(end), 0.0))
                    continue

                # Cresil-like trimming: shrink merged regions to the min/max span
                # where coverage depth is >= min_avg_depth.
                high = overlapping[overlapping["depth"] >= float(self.config.min_avg_depth)]
                if high.empty:
                    regions.append((chrom, int(start), int(end), 0.0))
                    continue

                trimmed_start = int(max(start, high["start"].min()))
                trimmed_end = int(min(end, high["end"].max()))
                if trimmed_end <= trimmed_start:
                    continue

                # Weighted average depth within trimmed region.
                total_depth = 0.0
                total_len = 0
                for _, cov_row in high.iterrows():
                    cov_start = max(int(cov_row["start"]), trimmed_start)
                    cov_end = min(int(cov_row["end"]), trimmed_end)
                    length = cov_end - cov_start
                    if length > 0:
                        total_depth += float(cov_row["depth"]) * length
                        total_len += length

                avg_depth = total_depth / total_len if total_len > 0 else 0.0
                regions.append((chrom, trimmed_start, trimmed_end, avg_depth))

            # Cleanup bedtools temp files
            bt.cleanup()
            if genome_file:
                try:
                    os.unlink(genome_file)
                except OSError:
                    pass

            return regions

        except Exception as e:
            self.logger.error(f"bedtools error: {e}")
            # Fallback: scanline-based regions with avg depth.
            bt.cleanup()
            if genome_file:
                try:
                    os.unlink(genome_file)
                except OSError:
                    pass
            return self._calculate_coverage_scanline_regions(trim_results)

    def _calculate_coverage_cresil_style(
        self, trim_results: list[TrimResult]
    ) -> list[tuple[str, int, int, float]]:
        """Calculate coverage using Cresil-style groupby(mergeid).agg(min/max) strategy.

        This method preserves large regions even when coverage gaps exist within them.
        Unlike the bedtools merge approach which splits regions at coverage gaps,
        this method:
        1. First merges overlapping segments to assign merge IDs
        2. Calculates per-base coverage
        3. Filters by depth threshold
        4. Uses groupby(mergeid).agg({'start': 'min', 'end': 'max'}) to get the
           full extent of high-coverage areas within each original merged region

        This is critical for detecting large eccDNA (>50kb) where internal coverage
        gaps are common but the overall circular structure should be preserved.

        Args:
            trim_results: List of trim results

        Returns:
            List of (chrom, start, end, avg_depth) tuples
        """
        # Collect all alignment segments into BED format
        bed_rows = []
        for result in trim_results:
            for seg in result.segments:
                bed_rows.append({
                    "chrom": seg.chrom,
                    "start": seg.ref_start,
                    "end": seg.ref_end,
                    "name": seg.read_name,
                    "score": seg.mapq,
                    "strand": seg.strand,
                })

        if not bed_rows:
            return []

        df = pd.DataFrame(bed_rows)
        df = df.sort_values(["chrom", "start", "end"])

        genome_file = None
        try:
            aln_bed = bt.BedTool.from_dataframe(df).sort()

            # Step 1: Merge adjacent/overlapping regions to assign merge IDs
            merged_bed = aln_bed.merge()
            merged_df = merged_bed.to_dataframe(names=["chrom", "start", "end"])
            merged_df["mergeid"] = range(len(merged_df))

            # Build chrom sizes for genomecov
            chrom_sizes = self.chrom_sizes
            if not chrom_sizes:
                inferred = df.groupby("chrom")["end"].max()
                chrom_sizes = {str(chrom): int(end) for chrom, end in inferred.items() if pd.notna(end)}

            if chrom_sizes:
                with tempfile.NamedTemporaryFile(mode="w", suffix=".genome", delete=False) as f:
                    for chrom, size in chrom_sizes.items():
                        f.write(f"{chrom}\t{int(size)}\n")
                    genome_file = f.name
                genome_cov = aln_bed.genome_coverage(bg=True, g=genome_file)
            else:
                genome_cov = aln_bed.genome_coverage(bg=True)

            cov_df = genome_cov.to_dataframe(names=["chrom", "start", "end", "depth"])

            # Step 2: Assign mergeid to each coverage interval based on which merged region it overlaps
            # For efficiency, build an interval index per chromosome
            cov_df["mergeid"] = -1

            for chrom in cov_df["chrom"].unique():
                chrom_merged = merged_df[merged_df["chrom"] == chrom].copy()
                chrom_cov_mask = cov_df["chrom"] == chrom

                if chrom_merged.empty:
                    continue

                # Simple interval assignment: for each coverage interval, find overlapping merge region
                chrom_cov_indices = cov_df[chrom_cov_mask].index
                for idx in chrom_cov_indices:
                    cov_start = cov_df.loc[idx, "start"]
                    cov_end = cov_df.loc[idx, "end"]

                    # Find merged region that contains this coverage interval
                    overlaps = chrom_merged[
                        (chrom_merged["start"] <= cov_start) & (chrom_merged["end"] >= cov_end)
                    ]
                    if not overlaps.empty:
                        cov_df.loc[idx, "mergeid"] = overlaps.iloc[0]["mergeid"]
                    else:
                        # Partial overlap: find best match
                        partial = chrom_merged[
                            (chrom_merged["start"] < cov_end) & (chrom_merged["end"] > cov_start)
                        ]
                        if not partial.empty:
                            cov_df.loc[idx, "mergeid"] = partial.iloc[0]["mergeid"]

            # Step 3: Filter by minimum depth threshold
            high_cov = cov_df[cov_df["depth"] >= self.config.min_avg_depth].copy()

            if high_cov.empty:
                bt.cleanup()
                if genome_file:
                    try:
                        os.unlink(genome_file)
                    except OSError:
                        pass
                return []

            # Step 4: Key Cresil strategy - groupby mergeid and take min(start), max(end)
            # This preserves the full extent of the region even with internal gaps
            high_cov = high_cov[high_cov["mergeid"] >= 0]

            if high_cov.empty:
                bt.cleanup()
                if genome_file:
                    try:
                        os.unlink(genome_file)
                    except OSError:
                        pass
                return []

            # Group by (chrom, mergeid) and aggregate
            grouped = high_cov.groupby(["chrom", "mergeid"]).agg({
                "start": "min",
                "end": "max",
                "depth": "mean"  # Average depth across the region
            }).reset_index()

            # Build result list
            regions: list[tuple[str, int, int, float]] = []
            for _, row in grouped.iterrows():
                regions.append((
                    str(row["chrom"]),
                    int(row["start"]),
                    int(row["end"]),
                    float(row["depth"])
                ))

            # Cleanup
            bt.cleanup()
            if genome_file:
                try:
                    os.unlink(genome_file)
                except OSError:
                    pass

            return regions

        except Exception as e:
            self.logger.error(f"Cresil-style coverage error: {e}")
            bt.cleanup()
            if genome_file:
                try:
                    os.unlink(genome_file)
                except OSError:
                    pass
            # Fallback to scanline method
            return self._calculate_coverage_scanline_regions(trim_results)

    def _create_merged_regions(
        self,
        regions: list[tuple[str, int, int, float]],
        trim_results: list[TrimResult],
    ) -> list[MergedRegion]:
        """Create MergedRegion objects with supporting reads.

        Args:
            regions: List of (chrom, start, end, avg_depth) tuples
            trim_results: Original trim results

        Returns:
            List of MergedRegion objects
        """
        merged = []

        for i, (chrom, start, end, avg_depth) in enumerate(regions):
            region = MergedRegion(
                merge_id=i,
                chrom=chrom,
                start=start,
                end=end,
                strand="+",  # Will be determined by majority
                coverage_depth=avg_depth,
            )

            # Find supporting reads
            strand_votes = defaultdict(int)
            for result in trim_results:
                for seg in result.segments:
                    if (
                        seg.chrom == chrom
                        and seg.ref_start < end
                        and seg.ref_end > start
                    ):
                        region.supporting_reads.add(seg.read_name)
                        strand_votes[seg.strand] += 1

            # Set strand by majority vote
            if strand_votes:
                region.strand = max(strand_votes, key=strand_votes.get)

            merged.append(region)

        return merged

    def _annotate_breakpoints(
        self,
        regions: list[MergedRegion],
        trim_results: list[TrimResult],
    ) -> None:
        """Annotate breakpoint evidence on each merged region.

        Cresil-style: define a 5' and 3' end window for every merged region and
        mark a region as having a breakpoint when enough *distinct reads* overlap
        that end window.

        This is intentionally junction-oriented (end-window overlap), rather
        than clustering raw start/end coordinates, to reduce false cycles.
        """
        min_support = int(self.config.min_breakpoint_depth)
        end_size = int(self.config.overlap_check_size)
        if self.config.min_region_size < 200:
            end_size = max(15, int(round(self.config.min_region_size * 0.3)))
        end_size = max(1, end_size)

        # Build a per-chrom index for fast seg->region mapping.
        regions_by_chrom: dict[str, list[MergedRegion]] = defaultdict(list)
        for region in regions:
            regions_by_chrom[region.chrom].append(region)
        starts_by_chrom: dict[str, list[int]] = {}
        for chrom, lst in regions_by_chrom.items():
            lst.sort(key=lambda r: r.start)
            starts_by_chrom[chrom] = [r.start for r in lst]

        def _find_region(seg: AlignmentSegment) -> Optional[MergedRegion]:
            lst = regions_by_chrom.get(seg.chrom)
            if not lst:
                return None
            starts = starts_by_chrom[seg.chrom]
            idx = bisect_right(starts, seg.ref_start) - 1
            for j in (idx, idx + 1):
                if 0 <= j < len(lst):
                    candidate = lst[j]
                    if seg.ref_start < candidate.end and seg.ref_end > candidate.start:
                        return candidate
            return None

        ovl_5_reads: dict[int, set[str]] = defaultdict(set)
        ovl_3_reads: dict[int, set[str]] = defaultdict(set)

        for result in trim_results:
            for seg in result.segments:
                region = _find_region(seg)
                if region is None:
                    continue

                this_end = min(end_size, max(1, region.length))
                end5_start = region.start
                end5_end = min(region.end, region.start + this_end)
                end3_start = max(region.start, region.end - this_end)
                end3_end = region.end

                if seg.ref_start < end5_end and seg.ref_end > end5_start:
                    ovl_5_reads[region.merge_id].add(result.read_name)
                if seg.ref_start < end3_end and seg.ref_end > end3_start:
                    ovl_3_reads[region.merge_id].add(result.read_name)

        for region in regions:
            region.has_5prime_break = len(ovl_5_reads.get(region.merge_id, set())) >= min_support
            region.has_3prime_break = len(ovl_3_reads.get(region.merge_id, set())) >= min_support

    def _find_circular_candidates(
        self,
        regions: list[MergedRegion],
        trim_results: list[TrimResult],
    ) -> list[EccDNACandidate]:
        """Find circular eccDNA candidates using graph-based detection.

        Args:
            regions: List of merged regions
            trim_results: Original trim results

        Returns:
            List of eccDNA candidates
        """
        # Build a Cresil-like breakpoint graph:
        # - Nodes: merged regions
        # - Edges: only when a read traverses region end windows with a
        #   strand-consistent breakpoint direction (junction-level evidence).
        regions_by_id: dict[int, MergedRegion] = {r.merge_id: r for r in regions}

        regions_by_chrom: dict[str, list[MergedRegion]] = defaultdict(list)
        for region in regions:
            regions_by_chrom[region.chrom].append(region)
        starts_by_chrom: dict[str, list[int]] = {}
        for chrom, lst in regions_by_chrom.items():
            lst.sort(key=lambda r: r.start)
            starts_by_chrom[chrom] = [r.start for r in lst]

        def _find_region_id(seg: AlignmentSegment) -> Optional[int]:
            lst = regions_by_chrom.get(seg.chrom)
            if not lst:
                return None
            starts = starts_by_chrom[seg.chrom]
            idx = bisect_right(starts, seg.ref_start) - 1
            for j in (idx, idx + 1):
                if 0 <= j < len(lst):
                    candidate = lst[j]
                    if seg.ref_start < candidate.end and seg.ref_end > candidate.start:
                        return candidate.merge_id
            return None

        min_support = int(self.config.min_breakpoint_depth)
        end_size = int(self.config.overlap_check_size)
        if self.config.min_region_size < 200:
            end_size = max(15, int(round(self.config.min_region_size * 0.3)))
        end_size = max(1, end_size)

        # Cresil uses a fixed +50bp query adjacency window.
        query_link_window = 50

        edge_weight: dict[tuple[int, int], int] = defaultdict(int)
        edge_strand_pairs: dict[tuple[int, int], set[str]] = defaultdict(set)

        for result in trim_results:
            segs = result.segments
            if len(segs) < 2:
                continue

            mapped: list[Optional[tuple[int, bool, bool, str]]] = []
            for seg in segs:
                rid = _find_region_id(seg)
                if rid is None:
                    mapped.append(None)
                    continue

                region = regions_by_id[rid]
                this_end = min(end_size, max(1, region.length))
                end5_end = min(region.end, region.start + this_end)
                end3_start = max(region.start, region.end - this_end)

                ovl5 = seg.ref_start < end5_end and seg.ref_end > region.start
                ovl3 = seg.ref_start < region.end and seg.ref_end > end3_start
                mapped.append((rid, ovl5, ovl3, seg.strand))

            for i in range(len(segs) - 1):
                left = mapped[i]
                right = mapped[i + 1]
                if left is None or right is None:
                    continue

                # Require adjacency on the read (no large gaps between consecutive segments).
                if segs[i + 1].read_start > segs[i].read_end + query_link_window:
                    break

                l_rid, l_ovl5, l_ovl3, l_strand = left
                r_rid, r_ovl5, r_ovl3, r_strand = right

                src: Optional[int] = None
                dst: Optional[int] = None
                strand_pair: Optional[str] = None

                if l_strand == "+" and r_strand == "+":
                    if l_ovl3 and r_ovl5:
                        src, dst, strand_pair = l_rid, r_rid, "+_+"
                elif l_strand == "-" and r_strand == "-":
                    # Reverse direction on reference.
                    if l_ovl5 and r_ovl3:
                        src, dst, strand_pair = r_rid, l_rid, "-_-"
                elif l_strand == "-" and r_strand == "+":
                    if l_ovl5 and r_ovl5:
                        src, dst, strand_pair = l_rid, r_rid, "-_+"
                elif l_strand == "+" and r_strand == "-":
                    if l_ovl3 and r_ovl3:
                        src, dst, strand_pair = l_rid, r_rid, "+_-"

                if src is None or dst is None or strand_pair is None:
                    continue

                edge_weight[(src, dst)] += 1
                edge_strand_pairs[(src, dst)].add(strand_pair)

        # Materialize graph with edge support filtering (breakpoint depth).
        G = nx.DiGraph()
        for region in regions:
            G.add_node(
                region.merge_id,
                chrom=region.chrom,
                start=region.start,
                end=region.end,
                strand=region.strand,
            )

        for (src, dst), weight in edge_weight.items():
            if weight < min_support:
                continue
            G.add_edge(src, dst, weight=weight, strand_pairs=edge_strand_pairs[(src, dst)])

        candidates: list[EccDNACandidate] = []
        ecc_counter = 0

        def _node_key(node_id: int) -> tuple[str, int, int]:
            r = regions_by_id[node_id]
            return (r.chrom, r.start, r.end)

        UG = G.to_undirected()
        for comp_nodes in nx.connected_components(UG):
            nodes = sorted(comp_nodes, key=_node_key)
            if not nodes:
                continue

            contain_selfloop = any(G.has_edge(n, n) for n in nodes)

            comp_graph = UG.subgraph(nodes).copy()
            comp_graph.remove_edges_from(nx.selfloop_edges(comp_graph))

            can_be_solved = all(deg <= 2 for _, deg in comp_graph.degree())

            is_cyclic = False
            traversal: list[int] = []

            if len(nodes) == 1:
                is_cyclic = contain_selfloop
                traversal = nodes
            elif len(nodes) == 2:
                a, b = nodes
                if G.has_edge(a, b) and G.has_edge(b, a):
                    pairs_ab = G[a][b].get("strand_pairs", set())
                    pairs_ba = G[b][a].get("strand_pairs", set())
                    compat = False
                    for p in pairs_ab:
                        a_s, b_s = p.split("_")
                        for q in pairs_ba:
                            b2_s, a2_s = q.split("_")
                            if a_s == a2_s and b_s == b2_s:
                                compat = True
                                break
                        if compat:
                            break
                    is_cyclic = compat
                traversal = nodes
            else:
                basis = nx.cycle_basis(comp_graph)
                if basis:
                    traversal = max(basis, key=len)
                    # Canonicalize traversal start for stability.
                    min_node = min(traversal, key=_node_key)
                    shift = traversal.index(min_node)
                    traversal = traversal[shift:] + traversal[:shift]
                    is_cyclic = True

            if not is_cyclic:
                continue
            if len(nodes) > 1 and not can_be_solved:
                continue

            cycle_regions = [regions_by_id[n] for n in traversal]

            # Breakpoint validation based on config
            strict_validation = getattr(self.config, "strict_breakpoint_validation", False)

            if strict_validation:
                # Legacy behavior: require each region to have both ends supported
                if not all(r.has_5prime_break and r.has_3prime_break for r in cycle_regions):
                    continue
            else:
                # Cresil-style relaxed validation:
                # - Single-region circles: still require both ends (strong evidence needed)
                # - Multi-region circles: edges already filtered by min_support,
                #   only require at least one region to have breakpoint evidence
                if len(cycle_regions) == 1:
                    if not (cycle_regions[0].has_5prime_break and cycle_regions[0].has_3prime_break):
                        continue
                else:
                    # For multi-region circles, at least one region must have breakpoint evidence
                    if not any(r.has_5prime_break or r.has_3prime_break for r in cycle_regions):
                        continue

            all_reads: set[str] = set()
            for r in cycle_regions:
                all_reads.update(r.supporting_reads)

            total_base = sum(r.length * len(r.supporting_reads) for r in cycle_regions)
            merge_len = sum(r.length for r in cycle_regions)
            coverage = total_base / merge_len if merge_len > 0 else 0.0

            is_ctc = any(
                result.is_ctc
                for result in trim_results
                if result.read_name in all_reads
            )

            ecc_counter += 1
            candidates.append(
                EccDNACandidate(
                    ecc_id=f"ec{ecc_counter}",
                    regions=cycle_regions,
                    is_circular=True,
                    is_ctc=is_ctc,
                    num_reads=len(all_reads),
                    total_base=total_base,
                    coverage=coverage,
                )
            )

        # Method 2: Single-region circles (self-loops via CTC)
        for result in trim_results:
            if not result.is_ctc:
                continue

            for chrom, start, end in result.ctc_regions:
                # Find corresponding merged region
                for region in regions:
                    if (
                        region.chrom == chrom
                        and abs(region.start - start) < 500
                        and abs(region.end - end) < 500
                    ):
                        # Check if this region is already in a cycle
                        already_in_cycle = any(
                            region in c.regions for c in candidates
                        )
                        if already_in_cycle:
                            continue

                        all_reads = region.supporting_reads
                        total_base = region.length * len(all_reads)
                        coverage = total_base / region.length if region.length > 0 else 0

                        ecc_counter += 1
                        candidate = EccDNACandidate(
                            ecc_id=f"ec{ecc_counter}",
                            regions=[region],
                            is_circular=True,
                            is_ctc=True,
                            num_reads=len(all_reads),
                            total_base=total_base,
                            coverage=coverage,
                        )
                        candidates.append(candidate)
                        break

        return candidates


# ============================================================================
# Main SplitReadsCaller Class
# ============================================================================


class SplitReadsCaller:
    """Main entry point for the SplitReads-Caller module.

    This class orchestrates the full pipeline:
    1. Initialize mappy aligner
    2. Run trim phase on input FASTA
    3. Run identify phase to detect eccDNA
    4. Write output in Cresil-compatible format
    """

    def __init__(
        self,
        reference_path: Path,
        config: SplitReadsConfig,
        logger: Optional[logging.Logger] = None,
        threads: int = 1,
    ) -> None:
        """Initialize SplitReadsCaller.

        Args:
            reference_path: Path to reference FASTA or MMI index
            config: SplitReads configuration
            logger: Logger instance
            threads: Number of threads
        """
        self.reference_path = Path(reference_path)
        self.config = config
        self.logger = logger or get_logger(self.__class__.__name__)
        self.threads = threads

    def run(
        self,
        input_fasta: Path,
        output_dir: Path,
    ) -> Path:
        """Run the full SplitReads-Caller pipeline.

        Args:
            input_fasta: Input FASTA file with reads
            output_dir: Output directory

        Returns:
            Path to eccDNA_final.txt output file
        """
        input_fasta = Path(input_fasta)
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        self.logger.info("Starting SplitReads-Caller pipeline")
        self.logger.info(f"  Input: {input_fasta}")
        self.logger.info(f"  Reference: {self.reference_path}")
        self.logger.info(f"  Output: {output_dir}")

        # Phase 1: Trim
        self.logger.info("Phase 1/2: Trim (alignment and CTC detection)")
        trim_processor = TrimProcessor(
            reference_path=self.reference_path,
            config=self.config,
            logger=self.logger.getChild("trim"),
            threads=self.threads,
        )

        trim_results = list(trim_processor.process_fasta(input_fasta))
        self.logger.info(f"Trim phase complete: {len(trim_results)} reads with alignments")

        # Write trim results
        trim_tsv = output_dir / "trim.tsv"
        trim_processor.write_trim_tsv(trim_results, trim_tsv)

        # Phase 2: Identify
        self.logger.info("Phase 2/2: Identify (eccDNA detection)")
        identify_processor = IdentifyProcessor(
            config=self.config,
            logger=self.logger.getChild("identify"),
            chrom_sizes=trim_processor.chrom_sizes,
        )

        candidates = identify_processor.identify(trim_results)
        self.logger.info(f"Identify phase complete: {len(candidates)} eccDNA candidates")

        # Write output in Cresil-compatible format
        output_file = output_dir / "eccDNA_final.txt"
        self._write_output(candidates, output_file)

        self.logger.info(f"SplitReads-Caller complete. Output: {output_file}")
        return output_file

    def _write_output(
        self, candidates: list[EccDNACandidate], output_path: Path
    ) -> None:
        """Write eccDNA candidates in Cresil-compatible format.

        Args:
            candidates: List of eccDNA candidates
            output_path: Output file path
        """
        rows = []
        for candidate in candidates:
            rows.append({
                "id": candidate.ecc_id,
                "merge_region": candidate.merge_region_str,
                "merge_len": candidate.merge_len,
                "num_region": len(candidate.regions),
                "ctc": candidate.is_ctc,
                "numreads": candidate.num_reads,
                "totalbase": candidate.total_base,
                "coverage": round(candidate.coverage, 2),
            })

        df = pd.DataFrame(rows)

        # Ensure column order matches Cresil output
        columns = [
            "id", "merge_region", "merge_len", "num_region",
            "ctc", "numreads", "totalbase", "coverage"
        ]
        df = df.reindex(columns=columns)

        output_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_path, sep="\t", index=False)
        self.logger.info(f"Wrote {len(df)} eccDNA candidates to {output_path}")
