"""
NGS Evidence Extraction — Extract eccDNA-supporting evidence from short-read BAM.

Extracts three types of evidence:
  1. Split reads (SA tag): reads split across a circle junction
  2. Discordant pairs: PE reads with outward orientation or large insert size
  3. Soft-clipped reads: reads partially aligned at junction boundaries

Output: a DataFrame of breakpoint observations with evidence type and quality.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pysam
import pandas as pd

from circleseeker.utils.logging import get_logger


@dataclass
class Breakpoint:
    """A single breakpoint observation from one read."""
    chr1: str
    pos1: int
    chr2: str
    pos2: int
    evidence_type: str  # "split", "discordant", "softclip"
    read_name: str
    mapq: int = 0
    strand1: str = "+"
    strand2: str = "+"


@dataclass
class NGSEvidenceConfig:
    """Configuration for NGS evidence extraction."""
    min_mapq: int = 10
    min_clip_len: int = 20
    max_normal_insert: int = 1000
    min_split_align_len: int = 30


class NGSEvidence:
    """Extract eccDNA evidence from a sorted, indexed BAM file."""

    def __init__(
        self,
        config: Optional[NGSEvidenceConfig] = None,
        logger: Optional[logging.Logger] = None,
    ):
        self.config = config or NGSEvidenceConfig()
        self.logger = logger or get_logger(self.__class__.__name__)

    def extract(self, bam_path: Path) -> pd.DataFrame:
        """Extract all eccDNA evidence from BAM.

        Returns DataFrame with columns:
            chr1, pos1, strand1, chr2, pos2, strand2,
            evidence_type, read_name, mapq
        """
        bam_path = Path(bam_path)
        breakpoints: list[dict] = []

        bam = pysam.AlignmentFile(str(bam_path), "rb")

        n_reads = 0
        n_split = 0
        n_discordant = 0
        n_softclip = 0

        for read in bam.fetch():
            n_reads += 1

            if read.is_unmapped or read.is_secondary or read.is_duplicate:
                continue
            if read.mapping_quality < self.config.min_mapq:
                continue

            chrom = read.reference_name
            start = read.reference_start
            end = read.reference_end

            # --- Split reads (SA tag) ---
            if read.has_tag("SA"):
                for bp in self._parse_sa(read, chrom, start, end):
                    breakpoints.append(bp)
                    n_split += 1

            # --- Discordant pairs ---
            if read.is_paired and not read.mate_is_unmapped and read.is_read1:
                bp = self._check_discordant(read, chrom, start, end)
                if bp:
                    breakpoints.append(bp)
                    n_discordant += 1

            # --- Soft-clipped reads ---
            for bp in self._check_softclip(read, chrom, start, end):
                breakpoints.append(bp)
                n_softclip += 1

        bam.close()

        self.logger.info(
            f"NGS evidence: {n_reads} reads → "
            f"{n_split} split, {n_discordant} discordant, {n_softclip} softclip"
        )

        if not breakpoints:
            return pd.DataFrame()

        df = pd.DataFrame(breakpoints)
        return df

    def _parse_sa(self, read, chrom: str, start: int, end: int) -> list[dict]:
        """Parse supplementary alignments from SA tag."""
        results = []
        sa_str = read.get_tag("SA")
        for entry in sa_str.strip(";").split(";"):
            parts = entry.split(",")
            if len(parts) < 6:
                continue
            sa_chr = parts[0]
            sa_pos = int(parts[1]) - 1  # 1-based to 0-based
            sa_strand = parts[2]
            sa_mapq = int(parts[4])

            if sa_mapq < self.config.min_mapq:
                continue

            # Parse SA CIGAR to get alignment length
            sa_cigar = parts[3]
            sa_alen = sum(int(n) for n, op in
                         __import__("re").findall(r"(\d+)([MIDNSHP=X])", sa_cigar)
                         if op in "MDN=X")
            if sa_alen < self.config.min_split_align_len:
                continue

            # Determine breakpoint: which end of primary connects to SA?
            # Use read query positions to determine junction order
            pri_strand = "-" if read.is_reverse else "+"

            results.append({
                "chr1": chrom,
                "pos1": end if not read.is_reverse else start,
                "strand1": pri_strand,
                "chr2": sa_chr,
                "pos2": sa_pos,
                "strand2": sa_strand,
                "evidence_type": "split",
                "read_name": read.query_name,
                "mapq": min(read.mapping_quality, sa_mapq),
            })
        return results

    def _check_discordant(self, read, chrom: str, start: int, end: int) -> Optional[dict]:
        """Check if PE pair is discordant (outward-facing or cross-chr)."""
        mate_chr = read.next_reference_name
        mate_pos = read.next_reference_start

        # Cross-chromosome = definitely interesting
        if chrom != mate_chr:
            return {
                "chr1": chrom,
                "pos1": start if not read.is_reverse else end,
                "strand1": "-" if read.is_reverse else "+",
                "chr2": mate_chr,
                "pos2": mate_pos,
                "strand2": "-" if read.mate_is_reverse else "+",
                "evidence_type": "discordant",
                "read_name": read.query_name,
                "mapq": read.mapping_quality,
            }

        # Same chromosome: check for outward-facing or large insert
        insert = abs(read.template_length)

        # Outward-facing pairs (both same strand) = eccDNA signature
        if read.is_reverse == read.mate_is_reverse:
            return {
                "chr1": chrom, "pos1": min(start, mate_pos),
                "strand1": "-" if read.is_reverse else "+",
                "chr2": chrom, "pos2": max(end, mate_pos),
                "strand2": "-" if read.mate_is_reverse else "+",
                "evidence_type": "discordant",
                "read_name": read.query_name,
                "mapq": read.mapping_quality,
            }

        # Very large insert size
        if insert > self.config.max_normal_insert:
            return {
                "chr1": chrom, "pos1": min(start, mate_pos),
                "strand1": "-" if read.is_reverse else "+",
                "chr2": chrom, "pos2": max(end, mate_pos + 150),
                "strand2": "-" if read.mate_is_reverse else "+",
                "evidence_type": "discordant",
                "read_name": read.query_name,
                "mapq": read.mapping_quality,
            }

        return None

    def _check_softclip(self, read, chrom: str, start: int, end: int) -> list[dict]:
        """Extract breakpoint positions from soft-clipped reads."""
        results = []
        if read.cigartuples is None:
            return results

        # Left soft-clip
        if read.cigartuples[0][0] == 4 and read.cigartuples[0][1] >= self.config.min_clip_len:
            results.append({
                "chr1": chrom, "pos1": start,
                "strand1": "-" if read.is_reverse else "+",
                "chr2": chrom, "pos2": start,
                "strand2": "-" if read.is_reverse else "+",
                "evidence_type": "softclip",
                "read_name": read.query_name,
                "mapq": read.mapping_quality,
            })

        # Right soft-clip
        if read.cigartuples[-1][0] == 4 and read.cigartuples[-1][1] >= self.config.min_clip_len:
            results.append({
                "chr1": chrom, "pos1": end,
                "strand1": "-" if read.is_reverse else "+",
                "chr2": chrom, "pos2": end,
                "strand2": "-" if read.is_reverse else "+",
                "evidence_type": "softclip",
                "read_name": read.query_name,
                "mapq": read.mapping_quality,
            })

        return results
