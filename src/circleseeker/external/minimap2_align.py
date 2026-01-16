"""Minimap2 PAF -> alignment TSV adapter for CircleSeeker."""

from __future__ import annotations

import os
import re
import shlex
import subprocess
import tempfile
from pathlib import Path
from typing import Iterable, List, Optional

from circleseeker.external.base import ExternalTool
from circleseeker.exceptions import ExternalToolError

_CS_OPS = {":", "*", "+", "-", "~"}


def _parse_cs_tag(cs: str) -> tuple[int, int, int]:
    """Parse minimap2 cs tag into mismatches, gap opens, gap bases."""
    mismatches = 0
    gap_opens = 0
    gap_bases = 0
    i = 0
    while i < len(cs):
        op = cs[i]
        if op == ":":
            i += 1
            while i < len(cs) and cs[i].isdigit():
                i += 1
        elif op == "*":
            mismatches += 1
            # cs mismatch is two bases (ref/query)
            i += 3 if i + 2 < len(cs) else len(cs)
        elif op in "+-":
            gap_opens += 1
            i += 1
            start = i
            while i < len(cs) and cs[i] not in _CS_OPS:
                i += 1
            gap_bases += i - start
        elif op == "~":
            gap_opens += 1
            i += 1
            start = i
            while i < len(cs) and cs[i] not in _CS_OPS:
                i += 1
            segment = cs[start:i]
            match = re.search(r"(\d+)", segment)
            if match:
                gap_bases += int(match.group(1))
        else:
            i += 1
    return mismatches, gap_opens, gap_bases


def _extract_cs_tag(tags: Iterable[str]) -> Optional[str]:
    for tag in tags:
        if tag.startswith("cs:Z:"):
            return tag[5:]
    return None


def paf_to_alignment_tsv(
    paf_path: Path,
    output_path: Path,
    logger=None,
    min_identity: float = 0.0,
    identity_decay_per_10kb: float = 0.0,
    min_identity_floor: float = 97.0,
) -> int:
    """Convert minimap2 PAF to alignment TSV (BLAST outfmt 6-like + MAPQ).

    Output columns (TSV, no header):
      query_id, subject_id, identity, alignment_length, mismatches, gap_opens,
      q_start, q_end, s_start, s_end, evalue, bit_score, sstrand, mapq

    Args:
        paf_path: Path to input PAF file
        output_path: Path to output TSV file
        logger: Optional logger for warnings
        min_identity: Base identity threshold (%). For HiFi data, 99.0 is recommended.
        identity_decay_per_10kb: Identity threshold decay per 10kb of query length (%).
            This implements a length-based compensation: longer sequences are allowed
            slightly lower identity because sequencing errors accumulate with length.
            Default 0.0 (no decay). Recommended: 0.5 for HiFi data.
        min_identity_floor: Minimum identity threshold floor (%). The threshold will
            never go below this value regardless of sequence length.
            Default 97.0.

    Returns:
        Number of alignments written
    """
    from circleseeker.utils.logging import get_logger

    log = logger or get_logger("minimap2_align")
    paf_path = Path(paf_path)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    REQUIRED_FIELDS = 12
    written = 0
    skipped = 0
    identity_filtered = 0

    with open(paf_path, "r") as infile, open(output_path, "w") as outfile:
        for line_num, line in enumerate(infile, 1):
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < REQUIRED_FIELDS:
                skipped += 1
                if skipped <= 5:
                    log.warning(
                        f"Line {line_num}: insufficient fields ({len(parts)} < {REQUIRED_FIELDS})"
                    )
                continue

            try:
                qname = parts[0]
                qlen = int(parts[1])  # query length for compensation calc
                qstart = int(parts[2])
                qend = int(parts[3])
                strand = parts[4]
                tname = parts[5]
                tstart = int(parts[7])
                tend = int(parts[8])
                nmatch = int(parts[9])
                alen = int(parts[10])
                mapq = int(parts[11])
            except (TypeError, ValueError, IndexError) as e:
                skipped += 1
                if skipped <= 5:
                    log.warning(f"Line {line_num}: parse error: {e}")
                continue

            cs_tag = _extract_cs_tag(parts[12:])
            if cs_tag:
                mismatches, gap_opens, gap_bases = _parse_cs_tag(cs_tag)
            else:
                gap_bases = 0
                mismatches = max(alen - nmatch, 0)
                gap_opens = 0

            identity = (nmatch / alen * 100.0) if alen else 0.0

            # Filter by identity threshold with length-based compensation
            # Longer sequences allow slightly lower identity (sequencing errors accumulate)
            if min_identity > 0.0:
                if identity_decay_per_10kb > 0.0:
                    # Calculate length-compensated threshold
                    length_kb = qlen / 1000.0
                    decay = (length_kb / 10.0) * identity_decay_per_10kb
                    identity_threshold = max(min_identity - decay, min_identity_floor)
                else:
                    identity_threshold = min_identity

                if identity < identity_threshold:
                    identity_filtered += 1
                    continue

            q_start = qstart + 1
            q_end = qend

            # Emit BLAST-like coordinates for downstream compatibility:
            # - q_start/q_end and s_start/s_end are 1-based and strand-oriented
            # - Downstream modules convert to 0-based half-open (start0/end0)
            if strand == "+":
                s_start = tstart + 1
                s_end = tend
                sstrand = "plus"
            else:
                s_start = tend
                s_end = tstart + 1
                sstrand = "minus"

            out_fields = [
                qname,
                tname,
                f"{identity:.2f}",
                str(alen),
                str(mismatches),
                str(gap_opens),
                str(q_start),
                str(q_end),
                str(s_start),
                str(s_end),
                "0",  # evalue unavailable in minimap2
                "0",  # bit_score unavailable in minimap2
                sstrand,
                str(mapq),
            ]
            outfile.write("\t".join(out_fields) + "\n")
            written += 1

    if skipped > 0:
        log.warning(f"Skipped {skipped} malformed lines in PAF file")

    if identity_filtered > 0:
        if identity_decay_per_10kb > 0.0:
            log.info(
                f"Filtered {identity_filtered} alignments with identity below threshold "
                f"(base={min_identity:.1f}%, decay={identity_decay_per_10kb:.1f}%/10kb, floor={min_identity_floor:.1f}%)"
            )
        else:
            log.info(
                f"Filtered {identity_filtered} alignments with identity < {min_identity:.1f}%"
            )

    return written


class Minimap2Aligner(ExternalTool):
    """Run minimap2 and emit alignment TSV for downstream steps."""

    tool_name = "minimap2"

    def run_alignment(
        self,
        reference: Path,
        query: Path,
        output_tsv: Path,
        preset: str = "sr",
        max_target_seqs: int = 200,
        additional_args: str = "",
        output_paf: Optional[Path] = None,
        min_identity: float = 0.0,
        identity_decay_per_10kb: float = 0.0,
        min_identity_floor: float = 97.0,
        split_by_length: bool = False,
        split_length: int = 5000,
        preset_short: str = "sr",
        preset_long: str = "asm5",
    ) -> Path:
        """Run minimap2 alignment with optional length-based preset splitting.

        Args:
            reference: Reference genome FASTA file
            query: Query FASTA file
            output_tsv: Output alignment TSV file
            preset: Minimap2 preset (used when split_by_length=False)
            max_target_seqs: Maximum secondary alignments (-N)
            additional_args: Additional minimap2 arguments
            output_paf: Output PAF file (default: output_tsv with .paf suffix)
            min_identity: Base identity threshold (%) for filtering
            identity_decay_per_10kb: Identity decay per 10kb (%). Longer sequences
                get a lower threshold. Set to 0.0 to disable (default).
                Recommended: 0.5 for HiFi data.
            min_identity_floor: Minimum identity floor (%). Threshold never goes
                below this value. Default: 97.0.
            split_by_length: Split queries by length and use different presets
            split_length: Length threshold for splitting (default: 5000bp)
            preset_short: Preset for short sequences (default: "sr")
            preset_long: Preset for long sequences (default: "asm5")

        Returns:
            Path to output TSV file
        """
        reference = Path(reference)
        query = Path(query)
        output_tsv = Path(output_tsv)

        if output_paf is None:
            output_paf = output_tsv.with_suffix(".paf")
        output_paf = Path(output_paf)
        output_paf.parent.mkdir(parents=True, exist_ok=True)

        if split_by_length:
            self._run_alignment_split(
                reference=reference,
                query=query,
                output_paf=output_paf,
                max_target_seqs=max_target_seqs,
                additional_args=additional_args,
                split_length=split_length,
                preset_short=preset_short,
                preset_long=preset_long,
            )
        else:
            self._run_minimap2(
                reference=reference,
                query=query,
                output_paf=output_paf,
                preset=preset,
                max_target_seqs=max_target_seqs,
                additional_args=additional_args,
            )

        written = paf_to_alignment_tsv(
            output_paf,
            output_tsv,
            logger=self.logger,
            min_identity=min_identity,
            identity_decay_per_10kb=identity_decay_per_10kb,
            min_identity_floor=min_identity_floor,
        )
        self.logger.info(f"Minimap2 alignments converted: {written} records")

        return output_tsv

    def _run_minimap2(
        self,
        reference: Path,
        query: Path,
        output_paf: Path,
        preset: str,
        max_target_seqs: int,
        additional_args: str,
    ) -> None:
        """Run minimap2 with a single preset."""
        cmd = [
            self.tool_name,
            "-x",
            preset,
            "-t",
            str(self.threads),
            "-N",
            str(max_target_seqs),
            "--secondary=yes",
            "-c",
            "--cs",
        ]

        if additional_args:
            cmd.extend(shlex.split(additional_args))

        cmd.extend([str(reference), str(query)])

        self.logger.info(f"Running minimap2 (preset={preset})")
        self.logger.debug(f"Command: {' '.join(cmd)}")

        try:
            log_file = output_paf.parent / "minimap2_align.log"
            with open(output_paf, "w") as out_handle, open(log_file, "w") as log_handle:
                subprocess.run(cmd, stdout=out_handle, stderr=log_handle, text=True, check=True)
        except subprocess.CalledProcessError as exc:
            stderr_tail = self._read_log_tail(output_paf.parent / "minimap2_align.log")
            raise ExternalToolError(
                "minimap2 failed",
                command=cmd,
                returncode=exc.returncode,
                stderr=stderr_tail or None,
            ) from exc

    def _run_alignment_split(
        self,
        reference: Path,
        query: Path,
        output_paf: Path,
        max_target_seqs: int,
        additional_args: str,
        split_length: int,
        preset_short: str,
        preset_long: str,
    ) -> None:
        """Run minimap2 with length-based preset splitting."""
        work_dir = output_paf.parent / "minimap2_split_tmp"
        work_dir.mkdir(parents=True, exist_ok=True)

        query_short = work_dir / "query.short.fa"
        query_long = work_dir / "query.long.fa"
        paf_short = work_dir / "align.short.paf"
        paf_long = work_dir / "align.long.paf"

        # Split query sequences by length
        n_short, n_long = self._split_fasta_by_length(
            query, query_short, query_long, split_length
        )
        self.logger.info(
            f"Split queries by length ({split_length}bp): "
            f"{n_short} short, {n_long} long"
        )

        paf_files: List[Path] = []

        # Align short sequences
        if n_short > 0:
            self._run_minimap2(
                reference=reference,
                query=query_short,
                output_paf=paf_short,
                preset=preset_short,
                max_target_seqs=max_target_seqs,
                additional_args=additional_args,
            )
            paf_files.append(paf_short)

        # Align long sequences
        if n_long > 0:
            self._run_minimap2(
                reference=reference,
                query=query_long,
                output_paf=paf_long,
                preset=preset_long,
                max_target_seqs=max_target_seqs,
                additional_args=additional_args,
            )
            paf_files.append(paf_long)

        # Merge PAF files
        self._merge_paf_files(paf_files, output_paf)

        # Cleanup temporary files
        try:
            import shutil
            shutil.rmtree(work_dir)
        except Exception:
            pass

    def _split_fasta_by_length(
        self,
        input_fasta: Path,
        output_short: Path,
        output_long: Path,
        split_length: int,
    ) -> tuple[int, int]:
        """Split FASTA file by sequence length.

        Returns:
            Tuple of (n_short, n_long) counts
        """
        n_short = 0
        n_long = 0

        with (
            open(input_fasta, "r") as fin,
            open(output_short, "w") as f_short,
            open(output_long, "w") as f_long,
        ):
            header = ""
            seq_parts: List[str] = []

            for line in fin:
                line = line.rstrip("\n")
                if line.startswith(">"):
                    # Write previous sequence
                    if header and seq_parts:
                        seq = "".join(seq_parts)
                        if len(seq) >= split_length:
                            f_long.write(f"{header}\n{seq}\n")
                            n_long += 1
                        else:
                            f_short.write(f"{header}\n{seq}\n")
                            n_short += 1
                    header = line
                    seq_parts = []
                else:
                    seq_parts.append(line)

            # Write last sequence
            if header and seq_parts:
                seq = "".join(seq_parts)
                if len(seq) >= split_length:
                    f_long.write(f"{header}\n{seq}\n")
                    n_long += 1
                else:
                    f_short.write(f"{header}\n{seq}\n")
                    n_short += 1

        return n_short, n_long

    def _merge_paf_files(self, paf_files: List[Path], output_paf: Path) -> None:
        """Merge multiple PAF files into one."""
        with open(output_paf, "w") as fout:
            for paf_file in paf_files:
                if paf_file.exists() and paf_file.stat().st_size > 0:
                    with open(paf_file, "r") as fin:
                        for line in fin:
                            fout.write(line)

    def _read_log_tail(self, log_path: Path, max_bytes: int = 32000) -> str:
        """Read the tail of a log file."""
        try:
            with open(log_path, "rb") as handle:
                try:
                    handle.seek(0, 2)
                    size = handle.tell()
                    handle.seek(max(0, size - max_bytes))
                except OSError:
                    pass
                return handle.read().decode(errors="replace")
        except OSError:
            return ""
