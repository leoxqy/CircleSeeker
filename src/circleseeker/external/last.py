"""LAST aligner wrapper for CircleSeeker.

LAST is a genome-scale sequence comparison tool that can:
- Handle repeat-rich genomes
- Find many-to-many alignments
- Use different scoring schemes for different sequence types

This module provides:
1. LastAligner - wrapper for running LAST alignment
2. last_tab_to_blast_tsv - convert LAST TAB format to BLAST outfmt 6 + MAPQ
"""

from __future__ import annotations

import math
import subprocess
import tempfile
from pathlib import Path
from typing import Optional, Tuple, List
import logging

from circleseeker.external.base import ExternalTool
from circleseeker.exceptions import ExternalToolError


def parse_alignment_blocks(blocks_str: str) -> Tuple[int, int, int, int]:
    """Parse LAST alignment blocks string to extract alignment statistics.

    LAST TAB alignment blocks format:
    - Simple: "778" means 778bp continuous match
    - Complex: "size,ins:del,size,ins:del,..." where:
      - size: gapless match length
      - ins: insertions in query after this block
      - del: deletions in query (insertions in ref) after this block

    Example: "29,6:0,15,1:0,5,0:4,178" means:
      - 29bp match, 6bp insertion in query, 0bp deletion
      - 15bp match, 1bp insertion in query, 0bp deletion
      - 5bp match, 0bp insertion, 4bp deletion
      - 178bp match (final block, no following gaps)

    Args:
        blocks_str: LAST alignment blocks string

    Returns:
        Tuple of (total_matches, total_mismatches, gap_opens, gap_bases)
    """
    if not blocks_str or blocks_str == "0":
        return 0, 0, 0, 0

    total_matches = 0
    gap_opens = 0
    gap_bases = 0

    # Replace ':' with ',' to get flat list of numbers
    # Format: size,ins,del,size,ins,del,...
    parts = blocks_str.replace(":", ",").split(",")

    for i, p in enumerate(parts):
        try:
            val = int(p)
        except ValueError:
            continue

        if i % 3 == 0:
            # Position 0, 3, 6, ... = match size
            total_matches += val
        else:
            # Position 1, 2, 4, 5, ... = gap sizes
            if val > 0:
                gap_opens += 1
                gap_bases += val

    # LAST doesn't directly report mismatches in TAB format
    # Mismatches are embedded in the score but not reported separately
    mismatches = 0

    return total_matches, mismatches, gap_opens, gap_bases


def mismap_to_mapq(mismap: float) -> int:
    """Convert LAST mismap probability to MAPQ score.

    MAPQ = -10 * log10(mismap)

    Args:
        mismap: Mismap probability (0 to 1)

    Returns:
        MAPQ score (0 to 60, capped)
    """
    if mismap <= 0:
        return 60  # Perfect mapping
    if mismap >= 1:
        return 0  # Completely ambiguous

    mapq = int(-10 * math.log10(mismap))
    return min(60, max(0, mapq))


def last_tab_to_blast_tsv(
    tab_path: Path,
    output_path: Path,
    logger: Optional[logging.Logger] = None,
    min_identity: float = 0.0,
    min_alignment_length: int = 0,
) -> int:
    """Convert LAST TAB format to BLAST outfmt 6-like TSV + MAPQ.

    LAST TAB columns (0-indexed):
    0: score
    1: subject_id (reference)
    2: s_start (0-based)
    3: s_aln_size
    4: s_strand
    5: s_seq_size
    6: query_id
    7: q_start (0-based)
    8: q_aln_size
    9: q_strand
    10: q_seq_size
    11: alignment_blocks
    12+: mismap=value (optional)

    Output columns (TSV, no header):
    query_id, subject_id, identity, alignment_length, mismatches, gap_opens,
    q_start, q_end, s_start, s_end, evalue, bit_score, sstrand, mapq

    Args:
        tab_path: Path to LAST TAB file
        output_path: Path to output TSV file
        logger: Optional logger
        min_identity: Minimum identity threshold (%)
        min_alignment_length: Minimum alignment length filter

    Returns:
        Number of alignments written
    """
    from circleseeker.utils.logging import get_logger

    log = logger or get_logger("last")
    tab_path = Path(tab_path)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    written = 0
    skipped = 0
    filtered_identity = 0
    filtered_length = 0

    with open(tab_path, "r") as infile, open(output_path, "w") as outfile:
        for line_num, line in enumerate(infile, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 12:
                skipped += 1
                if skipped <= 5:
                    log.warning(f"Line {line_num}: insufficient fields ({len(parts)} < 12)")
                continue

            try:
                score = int(parts[0])
                subject_id = parts[1]
                s_start = int(parts[2])  # 0-based
                s_aln_size = int(parts[3])
                s_strand = parts[4]
                s_seq_size = int(parts[5])
                query_id = parts[6]
                q_start = int(parts[7])  # 0-based
                q_aln_size = int(parts[8])
                q_strand = parts[9]
                q_seq_size = int(parts[10])
                alignment_blocks = parts[11]
            except (ValueError, IndexError) as e:
                skipped += 1
                if skipped <= 5:
                    log.warning(f"Line {line_num}: parse error: {e}")
                continue

            # Parse mismap probability
            mismap = 1.0  # Default: unknown
            for part in parts[12:]:
                if part.startswith("mismap="):
                    try:
                        mismap = float(part.split("=")[1])
                    except (ValueError, IndexError):
                        pass
                    break

            # Parse alignment blocks
            total_matches, mismatches, gap_opens, gap_bases = parse_alignment_blocks(alignment_blocks)

            # Calculate alignment length (use q_aln_size as reference)
            alignment_length = q_aln_size
            if alignment_length <= 0:
                alignment_length = total_matches + gap_bases

            # Filter by alignment length
            if min_alignment_length > 0 and alignment_length < min_alignment_length:
                filtered_length += 1
                continue

            # Calculate identity
            # For LAST, identity = matches / alignment_length
            if alignment_length > 0:
                identity = (total_matches / alignment_length) * 100.0
            else:
                identity = 0.0

            # Filter by identity
            if min_identity > 0 and identity < min_identity:
                filtered_identity += 1
                continue

            # Convert to MAPQ
            mapq = mismap_to_mapq(mismap)

            # Convert coordinates to 1-based BLAST-like format
            # Query coordinates (always forward strand in output)
            q_start_1based = q_start + 1
            q_end_1based = q_start + q_aln_size

            # Subject coordinates (strand-oriented like BLAST)
            # Determine effective strand (combination of q_strand and s_strand)
            if q_strand == s_strand:
                # Same strand - forward alignment
                s_start_1based = s_start + 1
                s_end_1based = s_start + s_aln_size
                sstrand = "plus"
            else:
                # Opposite strand - reverse alignment
                s_start_1based = s_start + s_aln_size
                s_end_1based = s_start + 1
                sstrand = "minus"

            # Output fields
            out_fields = [
                query_id,
                subject_id,
                f"{identity:.2f}",
                str(alignment_length),
                str(mismatches),
                str(gap_opens),
                str(q_start_1based),
                str(q_end_1based),
                str(s_start_1based),
                str(s_end_1based),
                f"{mismap:.2e}",  # Use mismap as evalue
                str(score),       # Use score as bit_score
                sstrand,
                str(mapq),
            ]
            outfile.write("\t".join(out_fields) + "\n")
            written += 1

    if skipped > 0:
        log.warning(f"Skipped {skipped} malformed lines in TAB file")
    if filtered_identity > 0:
        log.info(f"Filtered {filtered_identity} alignments below {min_identity}% identity")
    if filtered_length > 0:
        log.info(f"Filtered {filtered_length} alignments below {min_alignment_length}bp length")

    log.info(f"Converted {written} alignments from LAST TAB to BLAST TSV")
    return written


class LastAligner(ExternalTool):
    """LAST genome aligner wrapper.

    LAST workflow:
    1. lastdb - build database from reference
    2. lastal - align queries to database
    3. last-split - split alignments (optional)
    4. Convert TAB to BLAST TSV format
    """

    def __init__(
        self,
        threads: int = 1,
        scoring_matrix: str = "HOXD70",  # Good for ~95% identity
        logger: Optional[logging.Logger] = None,
    ):
        """Initialize LAST aligner.

        Args:
            threads: Number of parallel threads
            scoring_matrix: Scoring matrix preset (HOXD70, BLOSUM62, etc.)
            logger: Optional logger
        """
        super().__init__(logger)
        self.threads = threads
        self.scoring_matrix = scoring_matrix
        self._check_tools()

    def _check_tools(self) -> None:
        """Check if LAST tools are available."""
        required_tools = ["lastdb", "lastal"]
        missing = []
        for tool in required_tools:
            try:
                subprocess.run(
                    [tool, "--version"],
                    capture_output=True,
                    check=False,
                )
            except FileNotFoundError:
                missing.append(tool)

        if missing:
            raise ExternalToolError(
                f"Required LAST tools not found: {', '.join(missing)}. "
                "Install LAST: conda install -c bioconda last"
            )

    def build_database(
        self,
        reference_fasta: Path,
        db_prefix: Path,
    ) -> Path:
        """Build LAST database from reference.

        Args:
            reference_fasta: Path to reference FASTA
            db_prefix: Output database prefix

        Returns:
            Path to database prefix
        """
        reference_fasta = Path(reference_fasta)
        db_prefix = Path(db_prefix)
        db_prefix.parent.mkdir(parents=True, exist_ok=True)

        cmd = ["lastdb"]
        cmd.extend(["-P", str(self.threads)])
        cmd.extend([str(db_prefix), str(reference_fasta)])

        self.logger.info(f"Building LAST database: {' '.join(cmd)}")

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
        )

        if result.returncode != 0:
            raise ExternalToolError(f"lastdb failed: {result.stderr}")

        self.logger.info(f"LAST database built: {db_prefix}")
        return db_prefix

    def align(
        self,
        query_fasta: Path,
        db_prefix: Path,
        output_tab: Path,
        max_multiplicity: int = 10,
    ) -> Path:
        """Run LAST alignment.

        Args:
            query_fasta: Path to query FASTA
            db_prefix: LAST database prefix
            output_tab: Output TAB file path
            max_multiplicity: Maximum alignments per query position (-m)

        Returns:
            Path to output TAB file
        """
        query_fasta = Path(query_fasta)
        db_prefix = Path(db_prefix)
        output_tab = Path(output_tab)
        output_tab.parent.mkdir(parents=True, exist_ok=True)

        # Build lastal command
        cmd = ["lastal"]
        cmd.extend(["-P", str(self.threads)])
        cmd.extend(["-Q", "0"])  # FASTA input format
        cmd.extend(["-m", str(max_multiplicity)])  # Max alignments per query position
        # Output format: TAB with mismap probabilities
        cmd.extend(["-f", "TAB"])
        cmd.extend([str(db_prefix), str(query_fasta)])

        self.logger.info(f"Running LAST alignment: {' '.join(cmd)}")

        with open(output_tab, "w") as outfile:
            result = subprocess.run(
                cmd,
                stdout=outfile,
                stderr=subprocess.PIPE,
                text=True,
            )

        if result.returncode != 0:
            raise ExternalToolError(f"lastal failed: {result.stderr}")

        # Count alignments
        with open(output_tab, "r") as f:
            n_alignments = sum(1 for line in f if line.strip() and not line.startswith("#"))

        self.logger.info(f"LAST alignment complete: {n_alignments} alignments")
        return output_tab

    def run_alignment(
        self,
        query_fasta: Path,
        reference_fasta: Path,
        output_tsv: Path,
        db_prefix: Optional[Path] = None,
        min_identity: float = 90.0,
        min_alignment_length: int = 50,
    ) -> Path:
        """Run complete LAST alignment workflow.

        1. Build database (if needed)
        2. Run alignment
        3. Convert to BLAST TSV format

        Args:
            query_fasta: Query FASTA file
            reference_fasta: Reference FASTA file
            output_tsv: Output TSV file (BLAST format)
            db_prefix: Optional pre-built database prefix
            min_identity: Minimum identity filter (%)
            min_alignment_length: Minimum alignment length filter

        Returns:
            Path to output TSV file
        """
        query_fasta = Path(query_fasta)
        reference_fasta = Path(reference_fasta)
        output_tsv = Path(output_tsv)

        # Build database if not provided
        if db_prefix is None:
            db_prefix = output_tsv.parent / f"{reference_fasta.stem}_lastdb"

        db_prefix = Path(db_prefix)

        # Check if database exists
        if not db_prefix.with_suffix(".suf").exists():
            self.build_database(reference_fasta, db_prefix)
        else:
            self.logger.info(f"Using existing LAST database: {db_prefix}")

        # Run alignment
        tab_file = output_tsv.with_suffix(".last.tab")
        self.align(query_fasta, db_prefix, tab_file)

        # Convert to BLAST TSV format
        last_tab_to_blast_tsv(
            tab_file,
            output_tsv,
            logger=self.logger,
            min_identity=min_identity,
            min_alignment_length=min_alignment_length,
        )

        return output_tsv
