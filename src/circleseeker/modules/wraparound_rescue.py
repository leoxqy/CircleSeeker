"""
Wrap-around Rescue — Detect eccDNA from reads with >1 but <2 tandem copies.

For ONT RCA data, many reads are too short to contain ≥2 copies of the
eccDNA circle.  TideHunter requires ≥2 copies for tandem repeat detection
and misses these reads entirely.

However, reads with >1 copy exhibit a detectable "wrap-around" pattern:
the read's tail maps to the same reference region as its head, because
the read extends past one full circle into the beginning of the next copy.

    Circle: [--A--B--C--D--]
    Read (1.2 copies): A B C D A'B'
                        ↑             ↑
                        Both map to same reference region

This module detects wrap-around reads, extracts the consensus circle unit,
and doubles it to create input compatible with the existing pipeline
(alignment → um_classify / cecc_build).
"""

from __future__ import annotations

import logging
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Optional

from circleseeker.utils.logging import get_logger


def _ofrac(s1: int, e1: int, s2: int, e2: int) -> float:
    ov = max(0, min(e1, e2) - max(s1, s2))
    return min(ov / max(e1 - s1, 1), ov / max(e2 - s2, 1))


class WraparoundRescue:
    """Detect and rescue eccDNA reads that TideHunter missed."""

    def __init__(
        self,
        reference: Path,
        logger: Optional[logging.Logger] = None,
        threads: int = 4,
        min_ref_overlap: float = 0.3,
        min_query_gap: int = 100,
        min_unit_len: int = 200,
        minimap2_preset: str = "map-ont",
    ):
        self.reference = Path(reference)
        self.logger = logger or get_logger(self.__class__.__name__)
        self.threads = threads
        self.min_ref_overlap = min_ref_overlap
        self.min_query_gap = min_query_gap
        self.min_unit_len = min_unit_len
        self.minimap2_preset = minimap2_preset

    def run(
        self,
        input_fastx: Path,
        tidehunter_read_ids: set[str],
        output_dir: Path,
        prefix: str = "rescue",
    ) -> Optional[Path]:
        """Run wrap-around rescue pipeline.

        Args:
            input_fastx: Original ONT reads (FASTQ/FASTA).
            tidehunter_read_ids: Read IDs already detected by TideHunter.
            output_dir: Output directory.
            prefix: File prefix.

        Returns:
            Path to doubled FASTA of rescued reads, or None.
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Step 1: Align all non-TH reads to reference
        missed_fa = output_dir / f"{prefix}_missed.fasta"
        n_missed = self._extract_missed(input_fastx, tidehunter_read_ids, missed_fa)
        if n_missed == 0:
            self.logger.info("No TideHunter-missed reads to rescue")
            return None
        self.logger.info(f"Wrap-around rescue: {n_missed} reads to scan")

        paf = output_dir / f"{prefix}_rescue.paf"
        self._align_to_ref(missed_fa, paf)

        # Step 2: Detect wrap-around pattern
        read_aligns = self._parse_paf(paf)
        wr_reads = self._detect_wraparound(read_aligns)
        self.logger.info(f"Wrap-around rescue: {len(wr_reads)} wrap-around reads detected")

        if not wr_reads:
            self._cleanup(missed_fa, paf)
            return None

        # Step 3: Extract consensus units and double them
        read_seqs = self._load_seqs(input_fastx, set(wr_reads.keys()))
        doubled_fa = output_dir / f"{prefix}_wraparound_doubled.fasta"
        n_doubled = self._write_doubled(wr_reads, read_seqs, doubled_fa)
        self.logger.info(f"Wrap-around rescue: wrote {n_doubled} doubled sequences")

        self._cleanup(missed_fa, paf)
        return doubled_fa if n_doubled > 0 else None

    # ---- Private ----

    def _extract_missed(self, input_fastx: Path, th_ids: set[str], output: Path) -> int:
        n = 0
        is_fq = str(input_fastx).endswith((".fastq", ".fq"))
        with open(input_fastx) as fin, open(output, "w") as fout:
            if is_fq:
                while True:
                    h = fin.readline()
                    if not h:
                        break
                    s = fin.readline()
                    fin.readline()
                    fin.readline()
                    rid = h[1:].split()[0]
                    if rid not in th_ids:
                        fout.write(f">{rid}\n{s}")
                        n += 1
            else:
                write = False
                for line in fin:
                    if line.startswith(">"):
                        rid = line[1:].split()[0]
                        write = rid not in th_ids
                        if write:
                            fout.write(line)
                            n += 1
                    elif write:
                        fout.write(line)
        return n

    def _align_to_ref(self, query: Path, output_paf: Path) -> None:
        cmd = [
            "minimap2", "-cx", self.minimap2_preset,
            "--secondary=yes", "-N", "5",
            "-t", str(self.threads),
            str(self.reference), str(query),
        ]
        with open(output_paf, "w") as fout:
            subprocess.run(cmd, stdout=fout, stderr=subprocess.PIPE, timeout=3600)

    def _parse_paf(self, paf: Path) -> dict[str, list[dict]]:
        read_aligns: dict[str, list[dict]] = defaultdict(list)
        with open(paf) as f:
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 12:
                    continue
                read_aligns[cols[0]].append({
                    "rlen": int(cols[1]),
                    "qs": int(cols[2]),
                    "qe": int(cols[3]),
                    "chr": cols[5],
                    "cs": int(cols[7]),
                    "ce": int(cols[8]),
                    "strand": cols[4],
                    "mapq": int(cols[11]),
                })
        return dict(read_aligns)

    def _detect_wraparound(self, read_aligns: dict[str, list[dict]]) -> dict[str, dict]:
        """Detect reads whose tail wraps around to overlap with the head."""
        results = {}
        for rid, aligns in read_aligns.items():
            if len(aligns) < 2:
                continue
            sa = sorted(aligns, key=lambda a: a["qs"])
            first, last = sa[0], sa[-1]

            # Wrap-around: last alignment overlaps first on reference
            if first["chr"] != last["chr"]:
                continue
            if _ofrac(first["cs"], first["ce"], last["cs"], last["ce"]) < self.min_ref_overlap:
                continue
            # Must have query gap (= the circle body between head and tail)
            if last["qs"] - first["qe"] < self.min_query_gap:
                continue

            circle_size = last["qs"] - first["qs"]
            if circle_size < self.min_unit_len:
                continue

            results[rid] = {
                "circle_size": circle_size,
                "first_qs": first["qs"],
                "last_qs": last["qs"],
                "aligns": sa,
            }
        return results

    def _load_seqs(self, fastx: Path, rids: set[str]) -> dict[str, str]:
        seqs = {}
        is_fq = str(fastx).endswith((".fastq", ".fq"))
        with open(fastx) as f:
            if is_fq:
                while True:
                    h = f.readline()
                    if not h:
                        break
                    s = f.readline().strip()
                    f.readline()
                    f.readline()
                    rid = h[1:].split()[0]
                    if rid in rids:
                        seqs[rid] = s
            else:
                cur_rid, cur_seq = None, []
                for line in f:
                    if line.startswith(">"):
                        if cur_rid and cur_rid in rids:
                            seqs[cur_rid] = "".join(cur_seq)
                        cur_rid = line[1:].split()[0]
                        cur_seq = []
                    else:
                        cur_seq.append(line.strip())
                if cur_rid and cur_rid in rids:
                    seqs[cur_rid] = "".join(cur_seq)
        return seqs

    def _write_doubled(
        self, wr_reads: dict[str, dict], read_seqs: dict[str, str], output: Path,
    ) -> int:
        n = 0
        with open(output, "w") as f:
            for rid, info in wr_reads.items():
                seq = read_seqs.get(rid)
                if not seq:
                    continue
                # Extract one consensus unit from the read
                unit = seq[info["first_qs"]:info["last_qs"]]
                if len(unit) < self.min_unit_len:
                    continue
                # Double it (like tandem_to_ring does)
                doubled = unit + unit
                f.write(f">{rid}|rep0|{len(unit)}|2|circular\n{doubled}\n")
                n += 1
        return n

    def _cleanup(self, *files: Path) -> None:
        for f in files:
            f.unlink(missing_ok=True)
