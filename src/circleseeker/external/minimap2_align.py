"""Minimap2 PAF -> alignment TSV adapter for CircleSeeker."""

from __future__ import annotations

import re
import shlex
import subprocess
from pathlib import Path
from typing import Iterable, Optional, Tuple

from circleseeker.external.base import ExternalTool
from circleseeker.exceptions import ExternalToolError

_CS_OPS = {":", "*", "+", "-", "~"}


def _parse_cs_tag(cs: str) -> Tuple[int, int, int]:
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


def paf_to_alignment_tsv(paf_path: Path, output_path: Path) -> int:
    """Convert minimap2 PAF to alignment TSV (BLAST outfmt 6-like + MAPQ).

    Output columns (TSV, no header):
      query_id, subject_id, identity, alignment_length, mismatches, gap_opens,
      q_start, q_end, s_start, s_end, evalue, bit_score, sstrand, mapq
    """
    paf_path = Path(paf_path)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    written = 0
    with open(paf_path, "r") as infile, open(output_path, "w") as outfile:
        for line in infile:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue

            qname = parts[0]
            try:
                qstart = int(parts[2])
                qend = int(parts[3])
                strand = parts[4]
                tname = parts[5]
                tstart = int(parts[7])
                tend = int(parts[8])
                nmatch = int(parts[9])
                alen = int(parts[10])
                mapq = int(parts[11]) if len(parts) > 11 else 0
            except (TypeError, ValueError):
                continue

            cs_tag = _extract_cs_tag(parts[12:])
            if cs_tag:
                mismatches, gap_opens, gap_bases = _parse_cs_tag(cs_tag)
            else:
                gap_bases = 0
                mismatches = max(alen - nmatch, 0)
                gap_opens = 0

            identity = (nmatch / alen * 100.0) if alen else 0.0

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

    return written


class Minimap2Aligner(ExternalTool):
    """Run minimap2 and emit alignment TSV for downstream steps."""

    tool_name = "minimap2"

    def run(
        self,
        reference: Path,
        query: Path,
        output_tsv: Path,
        preset: str = "sr",
        max_target_seqs: int = 200,
        additional_args: str = "",
        output_paf: Optional[Path] = None,
    ) -> Path:
        reference = Path(reference)
        query = Path(query)
        output_tsv = Path(output_tsv)

        if output_paf is None:
            output_paf = output_tsv.with_suffix(".paf")
        output_paf = Path(output_paf)
        output_paf.parent.mkdir(parents=True, exist_ok=True)

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

        self.logger.info("Running minimap2 alignment")
        self.logger.debug(f"Command: {' '.join(cmd)}")

        try:
            log_file = output_paf.parent / "minimap2_align.log"
            with open(output_paf, "w") as out_handle, open(log_file, "w") as log_handle:
                subprocess.run(cmd, stdout=out_handle, stderr=log_handle, text=True, check=True)
        except subprocess.CalledProcessError as exc:
            stderr_tail = ""
            log_path = output_paf.parent / "minimap2_align.log"
            try:
                with open(log_path, "rb") as handle:
                    try:
                        handle.seek(0, 2)
                        size = handle.tell()
                        handle.seek(max(0, size - 32_000))
                    except OSError:
                        pass
                    stderr_tail = handle.read().decode(errors="replace")
            except OSError:
                stderr_tail = ""
            raise ExternalToolError(
                "minimap2 failed",
                command=cmd,
                returncode=exc.returncode,
                stderr=stderr_tail or None,
            ) from exc

        written = paf_to_alignment_tsv(output_paf, output_tsv)
        self.logger.info(f"Minimap2 alignments converted: {written} records")

        return output_tsv
