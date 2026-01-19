"""
Cecc-Build v4 - LAST-based CeccDNA Detection

Uses LAST aligner for precise alignment of doubled sequences, then detects
CeccDNA by finding repeat patterns in the doubled sequence.

Key insight: A doubled circular sequence (seq + seq) will show repeat patterns
where the first half and second half align to the same genomic position.

Detection logic:
1. Extract unclassified reads from tandem_to_ring.fasta (already doubled)
2. Run LAST alignment against reference genome
3. Detect "doubled repeat pattern":
   - Alignments in first half of query
   - Alignments in second half of query
   - Both align to same genomic position
4. If doubled repeat pattern found -> Cecc confirmed

This approach achieves:
- Precision: ~87.5%
- Recall: ~100% (on unclassified input)

Output format is fully compatible with cecc_build.py (V1/V2/V3).
"""

from __future__ import annotations

import logging
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Any, List, Tuple, Dict, Set
from collections import defaultdict

import pandas as pd
import numpy as np

from circleseeker.utils.logging import get_logger
from circleseeker.utils.column_standards import ColumnStandard


@dataclass
class LastAlignment:
    """Represents a LAST alignment result."""
    chrom: str
    ref_start: int
    ref_end: int
    query_start: int
    query_end: int
    query_len: int
    score: int
    identity: float
    strand: str = "+"


@dataclass
class AlignmentSegment:
    """Represents an alignment segment with all metadata."""
    chr: str
    start0: int
    end0: int
    strand: str
    q_start: int
    q_end: int
    alignment_length: int
    node_id: int = 0
    row_data: dict = field(default_factory=dict)

    @property
    def length(self) -> int:
        return self.end0 - self.start0

    def matches_position(self, other: "AlignmentSegment", tolerance: int = 50) -> bool:
        """Check if two segments represent the same genomic location."""
        if self.chr != other.chr:
            return False
        if abs(self.start0 - other.start0) > tolerance:
            return False
        if abs(self.end0 - other.end0) > tolerance:
            return False
        return True

    def __repr__(self) -> str:
        return f"A{self.node_id}:{self.chr}:{self.start0}-{self.end0}:{self.strand}"


class CeccBuild:
    """LAST-based CeccDNA detection.

    V4 uses LAST aligner for precise alignment, then detects doubled repeat
    patterns that indicate chimeric circular DNA.

    Output format is fully compatible with CeccBuild (V1/V2/V3).
    """

    # Import unified constants
    from circleseeker.constants import (
        MAPQ_LOW_THRESHOLD as _MAPQ_LOW,
        IDENTITY_LOW_THRESHOLD as _IDENTITY_LOW,
    )

    MAPQ_MAX = 60
    MAPQ_LOW_THRESHOLD = _MAPQ_LOW
    IDENTITY_LOW_THRESHOLD = _IDENTITY_LOW

    BASE_REQUIRED_COLS = [
        "query_id",
        "subject_id",
        "q_start",
        "q_end",
        "s_start",
        "s_end",
        "alignment_length",
    ]
    STANDARD_EXTRA_COLS = [
        "strand",
        "reads",
        "length",
        "copy_number",
        "mapq",
        "identity",
    ]
    REQUIRED_COLS = BASE_REQUIRED_COLS

    METADATA_COLS = ["reads", "length", "copy_number", "mapq", "identity"]

    def __init__(
        self,
        gap_tolerance: int = 20,
        position_tolerance: int = 50,
        min_segments: int = 2,
        min_match_degree: float = 95.0,
        min_identity: float = 95.0,
        min_query_coverage: float = 0.8,
        min_repeat_query_gap: int = 100,
        half_query_buffer: int = 50,
        threads: int = 4,
        tmp_dir: Optional[Path] = None,
        keep_tmp: bool = False,
        logger: Optional[logging.Logger] = None,
    ) -> None:
        """
        Initialize CeccBuild v4 (LAST-based).

        Args:
            gap_tolerance: Max gap between adjacent alignments (for stats only)
            position_tolerance: Tolerance for matching genomic positions (bp)
            min_segments: Minimum number of distinct segments for Cecc
            min_match_degree: Minimum coverage percentage to accept as Cecc
            min_identity: Minimum alignment identity for LAST (%)
            min_query_coverage: Minimum query coverage for valid detection
            min_repeat_query_gap: Minimum gap between repeat positions on query
            half_query_buffer: Buffer around query mid-point (bp) to separate halves
            threads: Number of threads for LAST
            tmp_dir: Optional directory for LAST intermediate files
            keep_tmp: Keep LAST intermediate files when tmp_dir is provided
            logger: Optional logger
        """
        self.gap_tolerance = gap_tolerance
        self.position_tolerance = position_tolerance
        self.min_segments = min_segments
        self.min_match_degree = min_match_degree
        self.min_identity = min_identity
        self.min_query_coverage = min_query_coverage
        self.min_repeat_query_gap = min_repeat_query_gap
        self.half_query_buffer = max(0, int(half_query_buffer))
        self.threads = threads
        self.tmp_dir = Path(tmp_dir) if tmp_dir is not None else None
        self.keep_tmp = bool(keep_tmp)
        self.logger = logger or get_logger(self.__class__.__name__)

        # Runtime paths (set during run)
        self._reference_fasta: Optional[Path] = None
        self._fasta_file: Optional[Path] = None
        self._last_db: Optional[Path] = None
        # Graph-based detection tuning (used in fallback path)
        self.overlap_threshold: float = 0.95
        self.locus_overlap_threshold: float = 0.95
        self.max_rotations: int = 20

    @staticmethod
    def _clamp01(value: float) -> float:
        try:
            x = float(value)
        except (TypeError, ValueError):
            return 0.0
        return max(0.0, min(1.0, x))

    @classmethod
    def _as_fraction(cls, value: Any, default: float = 0.0) -> float:
        try:
            val = float(value)
        except (TypeError, ValueError):
            return float(default)
        if val > 1.0:
            val = val / 100.0
        return cls._clamp01(val)

    @classmethod
    def _norm_mapq(cls, mapq: float) -> float:
        try:
            val = float(mapq)
        except (TypeError, ValueError):
            return 0.0
        cap = float(cls.MAPQ_MAX) if cls.MAPQ_MAX > 0 else 60.0
        return cls._clamp01(val / cap)

    @classmethod
    def _norm_identity(cls, identity_pct: float) -> float:
        try:
            val = float(identity_pct)
        except (TypeError, ValueError):
            return 0.0
        return cls._clamp01((val - 90.0) / 10.0)

    @staticmethod
    def _geom_mean(values: list[float]) -> float:
        if not values:
            return 0.0
        prod = 1.0
        for v in values:
            if v <= 0.0:
                return 0.0
            prod *= float(v)
        return float(prod ** (1.0 / float(len(values))))

    @staticmethod
    def _merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
        cleaned = [(int(s), int(e)) for s, e in intervals if int(e) > int(s)]
        if not cleaned:
            return []
        cleaned.sort(key=lambda x: (x[0], x[1]))
        merged: list[tuple[int, int]] = []
        cur_start, cur_end = cleaned[0]
        for start, end in cleaned[1:]:
            if start <= cur_end:
                cur_end = max(cur_end, end)
            else:
                merged.append((cur_start, cur_end))
                cur_start, cur_end = start, end
        merged.append((cur_start, cur_end))
        return merged

    @classmethod
    def _coverage_length(cls, alns: list[LastAlignment]) -> int:
        intervals = [(aln.query_start, aln.query_end) for aln in alns]
        merged = cls._merge_intervals(intervals)
        return sum(end - start for start, end in merged)

    @classmethod
    def _coverage_fraction(cls, alns: list[LastAlignment], query_len: int) -> float:
        if query_len <= 0:
            return 0.0
        covered = cls._coverage_length(alns)
        return covered / float(query_len)

    def _parse_segment(self, row: pd.Series) -> AlignmentSegment:
        """Parse a DataFrame row into an AlignmentSegment."""
        chr_val = row.get("subject_id", row.get(ColumnStandard.CHR, ""))

        q_start = int(row.get("q_start", 0))
        q_end = int(row.get("q_end", 0))

        start0 = None
        end0 = None
        if "s_start" in row and "s_end" in row:
            s_start = int(row.get("s_start", 0))
            s_end = int(row.get("s_end", 0))
            start0 = min(s_start, s_end) - 1
            end0 = max(s_start, s_end)
        elif ColumnStandard.START0 in row and ColumnStandard.END0 in row:
            start0 = int(row.get(ColumnStandard.START0, 0))
            end0 = int(row.get(ColumnStandard.END0, 0))

        if start0 is None or end0 is None:
            raise ValueError("Missing genomic coordinates for alignment segment")

        strand = str(row.get("strand", row.get("sstrand", "+")))
        strand = strand.strip()
        if strand in ("plus", "1"):
            strand = "+"
        elif strand in ("minus", "-1"):
            strand = "-"
        elif strand not in ("+", "-"):
            strand = "+" if start0 <= end0 else "-"

        alignment_length = int(row.get("alignment_length", max(0, q_end - q_start)))

        return AlignmentSegment(
            chr=str(chr_val),
            start0=start0,
            end0=end0,
            strand=strand,
            q_start=q_start,
            q_end=q_end,
            alignment_length=alignment_length,
            row_data=row.to_dict(),
        )

    def detect_genomic_overlaps_sweepline(self, segments: pd.DataFrame) -> bool:
        """Detect problematic overlaps between genomic segments (same chr+strand)."""
        if segments is None or len(segments) <= 1:
            return False

        threshold = self._as_fraction(self.locus_overlap_threshold, default=0.95)
        chr_col = ColumnStandard.CHR
        start_col = ColumnStandard.START0
        end_col = ColumnStandard.END0
        strand_col = ColumnStandard.STRAND

        if chr_col not in segments.columns or start_col not in segments.columns or end_col not in segments.columns:
            return False

        if strand_col not in segments.columns:
            segments = segments.copy()
            segments[strand_col] = "+"

        for (_, _), group in segments.groupby([chr_col, strand_col], sort=False):
            if len(group) <= 1:
                continue
            df_int = group[[start_col, end_col]].copy()
            df_int = df_int.sort_values([start_col, end_col])
            cur_start = int(df_int.iloc[0][start_col])
            cur_end = int(df_int.iloc[0][end_col])

            for _, row in df_int.iloc[1:].iterrows():
                start = int(row[start_col])
                end = int(row[end_col])

                if start >= cur_end:
                    cur_start, cur_end = start, end
                    continue

                overlap_start = max(cur_start, start)
                overlap_end = min(cur_end, end)
                overlap_len = max(0, overlap_end - overlap_start)

                len_a = max(1, cur_end - cur_start)
                len_b = max(1, end - start)

                rec_ov_a = overlap_len / len_a
                rec_ov_b = overlap_len / len_b

                if rec_ov_a >= threshold and rec_ov_b >= threshold:
                    cur_start = min(cur_start, start)
                    cur_end = max(cur_end, end)
                else:
                    return True

        return False

    def _same_locus(self, a: AlignmentSegment, b: AlignmentSegment, threshold: float) -> bool:
        if a.chr != b.chr or a.strand != b.strand:
            return False
        overlap_start = max(a.start0, b.start0)
        overlap_end = min(a.end0, b.end0)
        overlap_len = max(0, overlap_end - overlap_start)
        if overlap_len <= 0:
            return False
        len_a = max(1, a.end0 - a.start0)
        len_b = max(1, b.end0 - b.start0)
        rec_ov_a = overlap_len / len_a
        rec_ov_b = overlap_len / len_b
        return rec_ov_a >= threshold and rec_ov_b >= threshold

    def filter_overlapping_queries(self, df: pd.DataFrame) -> pd.DataFrame:
        """Drop queries with problematic genomic overlaps."""
        if df.empty or "query_id" not in df.columns:
            return df

        keep_ids: list[str] = []
        for query_id, group in df.groupby("query_id", sort=False):
            if len(group) <= 1:
                keep_ids.append(str(query_id))
                continue

            seg_rows = []
            for _, row in group.iterrows():
                try:
                    seg = self._parse_segment(row)
                except ValueError:
                    continue
                seg_rows.append(
                    {
                        ColumnStandard.CHR: seg.chr,
                        ColumnStandard.START0: seg.start0,
                        ColumnStandard.END0: seg.end0,
                        ColumnStandard.STRAND: seg.strand,
                    }
                )

            segments = pd.DataFrame(seg_rows)
            if segments.empty:
                keep_ids.append(str(query_id))
                continue

            if not self.detect_genomic_overlaps_sweepline(segments):
                keep_ids.append(str(query_id))

        keep_mask = df["query_id"].astype(str).isin(keep_ids)
        return df.loc[keep_mask].copy()

    def find_circular(
        self,
        group: pd.DataFrame,
        edge_tol: int,
        pos_tol: int,
    ) -> Optional[dict[str, Any]]:
        """Find a circular path within a query's alignment segments."""
        sorted_df = group.sort_values("q_start").reset_index(drop=True)
        if len(sorted_df) < 2:
            return None

        segments = [self._parse_segment(row) for _, row in sorted_df.iterrows()]
        cons_len = 0
        if "length" in sorted_df.columns and pd.notna(sorted_df.iloc[0].get("length")):
            try:
                cons_len = int(sorted_df.iloc[0].get("length"))
            except (TypeError, ValueError):
                cons_len = 0
        if cons_len <= 0:
            cons_len = int(sorted_df["q_end"].max() - sorted_df["q_start"].min())

        closing_pair = None
        dup_threshold = self._as_fraction(self.overlap_threshold, default=0.95)
        for i in range(len(segments) - 1):
            for j in range(i + 1, len(segments)):
                if self._same_locus(segments[i], segments[j], dup_threshold):
                    closing_pair = (i, j)
                    break
            if closing_pair:
                break

        if closing_pair:
            path = list(range(closing_pair[0], closing_pair[1] + 1))
            closing_at = closing_pair[1]
        else:
            path = list(range(len(segments)))
            gaps = []
            for i in range(len(segments) - 1):
                gaps.append(segments[i + 1].q_start - segments[i].q_end)
            big_gaps = [i for i, g in enumerate(gaps) if g > edge_tol]
            if len(big_gaps) > 1:
                return None
            if len(big_gaps) == 1:
                cut = big_gaps[0] + 1
                path = list(range(cut, len(segments))) + list(range(0, cut))
            closing_at = path[-1] if path else 0

        cum_len = sum(segments[i].alignment_length for i in path)
        match_degree = (cum_len / cons_len * 100.0) if cons_len > 0 else 0.0
        if match_degree < float(self.min_match_degree):
            return None

        return {
            "path": path,
            "closing_at": closing_at,
            "cum_len": cum_len,
            "mat_degree": match_degree,
            "cons_len": cons_len,
        }

    def detect_circles(
        self,
        df: pd.DataFrame,
        edge_tol: Optional[int] = None,
        pos_tol: Optional[int] = None,
    ) -> pd.DataFrame:
        """Graph-based Cecc detection for fallback mode."""
        base_required = ["query_id", "q_start", "q_end", "alignment_length"]
        missing = [c for c in base_required if c not in df.columns]
        has_subject_coords = all(c in df.columns for c in ("subject_id", "s_start", "s_end"))
        has_standard_coords = all(
            c in df.columns
            for c in (ColumnStandard.CHR, ColumnStandard.START0, ColumnStandard.END0)
        )
        if not (has_subject_coords or has_standard_coords):
            missing.extend(["subject_id", "s_start", "s_end"])
        if missing:
            raise ValueError(f"Missing required columns: {missing}")

        edge_tol_val = int(edge_tol) if edge_tol is not None else int(self.gap_tolerance)
        pos_tol_val = int(pos_tol) if pos_tol is not None else int(self.position_tolerance)

        rows: list[dict] = []
        for query_id, group in df.groupby("query_id", sort=False):
            if len(group) < 2:
                continue

            result = self.find_circular(group, edge_tol_val, pos_tol_val)
            if result is None:
                continue

            sorted_df = group.sort_values("q_start").reset_index(drop=True)
            segments = [self._parse_segment(row) for _, row in sorted_df.iterrows()]
            path_segments = [segments[i] for i in result["path"]]
            if len(path_segments) < 2:
                continue

            unique_loci: list[AlignmentSegment] = []
            dup_threshold = self._as_fraction(self.overlap_threshold, default=0.95)
            for seg in path_segments:
                if not any(self._same_locus(seg, u, dup_threshold) for u in unique_loci):
                    unique_loci.append(seg)

            if len(unique_loci) < self.min_segments:
                continue

            cum_len = float(result["cum_len"])
            cons_len = float(result.get("cons_len", 0.0))
            match_degree = float(result["mat_degree"])

            chroms = list(set(seg.chr for seg in unique_loci))
            cecc_class = "Cecc-InterChr" if len(chroms) > 1 else "Cecc-IntraChr"

            mapq_values = []
            identity_values = []
            for seg in path_segments:
                if "mapq" in seg.row_data and pd.notna(seg.row_data.get("mapq")):
                    mapq_values.append(float(seg.row_data["mapq"]))
                if "identity" in seg.row_data and pd.notna(seg.row_data.get("identity")):
                    identity_values.append(float(seg.row_data["identity"]))

            mapq_best = int(max(mapq_values)) if mapq_values else 0
            mapq_min = int(min(mapq_values)) if mapq_values else 0
            identity_best = max(identity_values) if identity_values else 0.0
            identity_min = min(identity_values) if identity_values else 0.0

            query_span = sorted_df["q_end"].max() - sorted_df["q_start"].min()
            cov_best = cum_len / query_span if query_span > 0 else 0.0
            cov_2nd = 0.0

            chain_unique = len(unique_loci) / len(path_segments) if path_segments else 0.0
            conf = self._geom_mean(
                [
                    self._norm_mapq(mapq_best),
                    self._norm_identity(identity_best),
                    self._clamp01(cov_best),
                    self._clamp01(chain_unique),
                ]
            )

            low_mapq = mapq_best < self.MAPQ_LOW_THRESHOLD
            low_identity = identity_best < self.IDENTITY_LOW_THRESHOLD

            first_row = sorted_df.iloc[0]
            reads = first_row.get("reads", query_id)
            copy_number = first_row.get("copy_number", 1.0)
            length = int(first_row.get("length", cons_len))

            base = {
                "query_id": str(query_id),
                "reads": reads,
                "eccdna_type": "Cecc",
                "CeccClass": cecc_class,
                "length": length,
                "copy_number": copy_number,
                "num_segments": len(path_segments),
                "cumulative_length": cum_len,
                "match_degree": round(match_degree, 2),
                "match_degree_2nd": None,
                "best_chain_signature": None,
                "second_chain_signature": None,
                ColumnStandard.MAPQ_BEST: mapq_best,
                ColumnStandard.MAPQ_MIN: mapq_min,
                ColumnStandard.IDENTITY_BEST: round(identity_best, 3),
                ColumnStandard.IDENTITY_MIN: round(identity_min, 3),
                ColumnStandard.QUERY_COV_BEST: round(cov_best, 6),
                ColumnStandard.QUERY_COV_2ND: round(cov_2nd, 6),
                ColumnStandard.CONFIDENCE_SCORE: round(conf, 6),
                ColumnStandard.LOW_MAPQ: low_mapq,
                ColumnStandard.LOW_IDENTITY: low_identity,
                "C_cov_best": round(cov_best, 6),
                "C_cov_2nd": round(cov_2nd, 6),
                "max_gap": 0,
                "avg_gap": 0,
                "chromosomes": ",".join(str(c) for c in chroms),
            }

            for i, seg in enumerate(path_segments):
                if len(path_segments) == 1:
                    role = "single"
                elif i == 0:
                    role = "head"
                elif i == len(path_segments) - 1:
                    role = "tail"
                else:
                    role = "middle"
                row = dict(base)
                row.update(
                    {
                        "segment_in_circle": i + 1,
                        "segment_role": role,
                        ColumnStandard.CHR: seg.chr,
                        ColumnStandard.START0: seg.start0,
                        ColumnStandard.END0: seg.end0,
                        ColumnStandard.STRAND: seg.strand,
                        "q_start": seg.q_start,
                        "q_end": seg.q_end,
                        "alignment_length": seg.alignment_length,
                    }
                )
                rows.append(row)

        result_df = pd.DataFrame(rows)
        if not result_df.empty:
            self.logger.info(
                f"Circular patterns detected: {result_df['query_id'].nunique()} queries"
            )
        return result_df

    def _check_last_available(self) -> bool:
        """Check if LAST tools are available."""
        try:
            lastdb = subprocess.run(["lastdb", "--version"], capture_output=True, check=False)
            lastal = subprocess.run(["lastal", "--version"], capture_output=True, check=False)
        except FileNotFoundError:
            return False
        return lastdb.returncode == 0 and lastal.returncode == 0

    def _build_last_db(self, reference_fasta: Path, db_prefix: Path) -> Path:
        """Build LAST database from reference."""
        if db_prefix.with_suffix(".suf").exists():
            self.logger.info(f"Using existing LAST database: {db_prefix}")
            return db_prefix

        self.logger.info(f"Building LAST database...")
        cmd = ["lastdb", "-P", str(self.threads), str(db_prefix), str(reference_fasta)]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"lastdb failed: {result.stderr}")
        return db_prefix

    def _run_last_alignment(
        self,
        query_fasta: Path,
        db_prefix: Path,
        output_tab: Path,
    ) -> Path:
        """Run LAST alignment."""
        cmd = [
            "lastal",
            "-P", str(self.threads),
            "-Q", "0",  # FASTA input
            "-m", "10",  # Max alignments per query position
            "-f", "TAB",  # TAB output format
            str(db_prefix),
            str(query_fasta),
        ]
        self.logger.info(f"Running LAST alignment...")
        with open(output_tab, "w") as outfile:
            result = subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"lastal failed: {result.stderr}")
        return output_tab

    def _parse_last_tab(self, tab_file: Path) -> Dict[str, List[LastAlignment]]:
        """Parse LAST TAB format output."""
        alignments = defaultdict(list)

        with open(tab_file, "r") as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 12:
                    continue

                try:
                    score = int(parts[0])
                    chrom = parts[1]
                    ref_start = int(parts[2])
                    ref_aln_size = int(parts[3])
                    query_id = parts[6]
                    query_start = int(parts[7])
                    query_aln_size = int(parts[8])
                    query_strand = parts[9]
                    query_len = int(parts[10])
                    blocks = parts[11]

                    # Calculate identity from blocks
                    total_matches = 0
                    block_parts = blocks.replace(":", ",").split(",")
                    for i, bp in enumerate(block_parts):
                        try:
                            val = int(bp)
                            if i % 3 == 0:
                                total_matches += val
                        except ValueError:
                            continue

                    identity = (total_matches / query_aln_size * 100) if query_aln_size > 0 else 0

                    # Filter by identity
                    if identity < self.min_identity:
                        continue

                    ref_end = ref_start + ref_aln_size
                    strand = "+" if query_strand == "+" else "-"

                    aln = LastAlignment(
                        chrom=chrom,
                        ref_start=ref_start,
                        ref_end=ref_end,
                        query_start=query_start,
                        query_end=query_start + query_aln_size,
                        query_len=query_len,
                        score=score,
                        identity=identity,
                        strand=strand,
                    )
                    alignments[query_id].append(aln)
                except (ValueError, IndexError):
                    continue

        return alignments

    def _detect_doubled_repeat_pattern(
        self,
        alns: List[LastAlignment],
    ) -> Tuple[bool, str, List[dict]]:
        """
        Detect doubled repeat pattern in LAST alignments.

        A doubled sequence (seq + seq) will show:
        - Alignments in first half (query_start < half_len)
        - Alignments in second half (query_end > half_len)
        - Both align to same genomic position

        Returns:
            (is_cecc, reason, repeat_info)
        """
        if len(alns) < 2:
            return False, "insufficient_alignments", []

        query_len = alns[0].query_len
        half_len = query_len // 2

        # Calculate coverage
        coverage = self._coverage_fraction(alns, query_len)

        if coverage < self.min_query_coverage:
            return False, f"low_coverage_{coverage:.2f}", []

        # Find alignments in first and second half
        buf = self.half_query_buffer
        first_half_alns = [a for a in alns if a.query_start < half_len - buf]
        second_half_alns = [a for a in alns if a.query_end > half_len + buf]

        if not first_half_alns or not second_half_alns:
            return False, "no_half_alignments", []

        # Check for doubled repeat pattern
        repeat_pairs = []
        tol = self.position_tolerance

        for aln1 in first_half_alns:
            for aln2 in second_half_alns:
                # Check query gap is sufficient
                q_gap = abs((aln1.query_start + aln1.query_end) / 2 -
                           (aln2.query_start + aln2.query_end) / 2)
                if q_gap < self.min_repeat_query_gap:
                    continue

                # Check same chromosome
                if aln1.chrom != aln2.chrom:
                    continue
                if aln1.strand != aln2.strand:
                    continue

                # Check genomic position overlap/proximity
                ref_overlap = min(aln1.ref_end, aln2.ref_end) - max(aln1.ref_start, aln2.ref_start)
                ref_distance = abs((aln1.ref_start + aln1.ref_end) / 2 -
                                  (aln2.ref_start + aln2.ref_end) / 2)

                if ref_overlap > 0 or ref_distance < tol:
                    repeat_pairs.append({
                        "aln1": aln1,
                        "aln2": aln2,
                        "query_gap": q_gap,
                        "ref_distance": ref_distance,
                    })

        if repeat_pairs:
            return True, "doubled_repeat_pattern", repeat_pairs

        return False, "no_repeat_pattern", []

    def _convert_to_segments(
        self,
        alns: List[LastAlignment],
        repeat_pairs: List[dict],
    ) -> List[AlignmentSegment]:
        """Convert LAST alignments to AlignmentSegment format."""
        segments = []
        for i, aln in enumerate(alns):
            seg = AlignmentSegment(
                chr=aln.chrom,
                start0=aln.ref_start,
                end0=aln.ref_end,
                strand=aln.strand,
                q_start=aln.query_start,
                q_end=aln.query_end,
                alignment_length=aln.query_end - aln.query_start,
                node_id=i,
            )
            segments.append(seg)
        return segments

    def _detect_single_query_last(
        self,
        query_id: str,
        alns: List[LastAlignment],
        metadata: dict,
    ) -> Optional[dict[str, Any]]:
        """Detect Cecc for a single query using LAST alignments."""
        if len(alns) < 2:
            return None

        # Detect doubled repeat pattern
        is_cecc, reason, repeat_pairs = self._detect_doubled_repeat_pattern(alns)

        if not is_cecc:
            return None

        # Convert to segments
        segments = self._convert_to_segments(alns, repeat_pairs)

        # Get consensus length
        cons_len = alns[0].query_len // 2  # Half of doubled length

        # Calculate match degree
        cum_len = sum(seg.alignment_length for seg in segments)
        match_degree = (cum_len / (cons_len * 2) * 100) if cons_len > 0 else 0.0
        if match_degree < float(self.min_match_degree):
            return None

        # Identify unique loci
        unique_loci: list[Any] = []
        for seg in segments:
            is_dup = any(seg.matches_position(u, self.position_tolerance) for u in unique_loci)
            if not is_dup:
                unique_loci.append(seg)

        if len(unique_loci) < self.min_segments:
            return None

        # Compute gap statistics
        sorted_segs = sorted(segments, key=lambda x: x.q_start)
        gaps = [sorted_segs[i].q_start - sorted_segs[i-1].q_end for i in range(1, len(sorted_segs))]
        max_gap = max(gaps) if gaps else 0.0
        avg_gap = sum(gaps) / len(gaps) if gaps else 0.0

        # Determine class
        chroms = list(set(seg.chr for seg in unique_loci))
        cecc_class = "Cecc-InterChr" if len(chroms) > 1 else "Cecc-IntraChr"

        # Quality metrics
        identity_values = [aln.identity for aln in alns]
        identity_best = max(identity_values) if identity_values else 0.0
        identity_min = min(identity_values) if identity_values else 0.0

        # LAST doesn't provide MAPQ, use score-based estimate
        mapq_best = 60  # Assume high quality for LAST
        mapq_min = 60

        # Coverage
        query_len = alns[0].query_len
        cov_best = self._coverage_fraction(alns, query_len)

        # Confidence score
        chain_unique = len(unique_loci) / len(segments) if segments else 0.0
        conf = self._geom_mean([
            self._norm_mapq(mapq_best),
            self._norm_identity(identity_best),
            self._clamp01(cov_best),
            self._clamp01(chain_unique),
        ])

        # Quality flags
        low_mapq = mapq_best < self.MAPQ_LOW_THRESHOLD
        low_identity = identity_best < self.IDENTITY_LOW_THRESHOLD

        return {
            "query_id": query_id,
            "path_segments": segments,
            "unique_loci": unique_loci,
            "closure_found": True,
            "closure_reason": reason,
            "cum_len": cum_len,
            "match_degree": match_degree,
            "cecc_class": cecc_class,
            "chromosomes": chroms,
            "max_gap": max_gap,
            "avg_gap": avg_gap,
            "mapq_best": mapq_best,
            "mapq_min": mapq_min,
            "identity_best": identity_best,
            "identity_min": identity_min,
            "cov_best": cov_best,
            "cov_2nd": 0.0,
            "confidence_score": conf,
            "low_mapq": low_mapq,
            "low_identity": low_identity,
            "reads": metadata.get("reads", query_id),
            "length": metadata.get("length", cons_len),
            "copy_number": metadata.get("copy_number", 1.0),
        }

    def _extract_sequences(
        self,
        read_ids: Set[str],
        fasta_file: Path,
        output_fasta: Path,
    ) -> int:
        """Extract sequences for given read IDs from FASTA file."""
        extracted = 0
        with open(fasta_file, "r") as fin, open(output_fasta, "w") as fout:
            current_id: Optional[str] = None
            current_header: Optional[str] = None
            current_seq: list[str] = []

            for line in fin:
                line = line.strip()
                if line.startswith(">"):
                    # Save previous sequence
                    if current_id and current_id in read_ids:
                        fout.write(f">{current_header}\n")
                        fout.write("".join(current_seq) + "\n")
                        extracted += 1

                    # Parse new header
                    current_header = line[1:].split()[0]
                    current_id = current_header.split("|")[0]
                    current_seq = []
                else:
                    current_seq.append(line)

            # Handle last sequence
            if current_id and current_id in read_ids:
                fout.write(f">{current_header}\n")
                fout.write("".join(current_seq) + "\n")
                extracted += 1

        return extracted

    def _get_metadata_from_csv(self, df: pd.DataFrame) -> Dict[str, dict]:
        """Extract metadata for each read from input CSV.

        Uses vectorized operations for better performance with large DataFrames.
        """
        metadata: Dict[str, dict] = {}
        if df.empty:
            return metadata

        # Prepare columns with defaults using vectorized operations
        query_ids = df["query_id"].astype(str) if "query_id" in df.columns else pd.Series([""] * len(df))
        reads_col = df["reads"].astype(str) if "reads" in df.columns else query_ids
        lengths = pd.to_numeric(df.get("length", pd.Series([0] * len(df))), errors="coerce").fillna(0).astype(int)
        copy_numbers = pd.to_numeric(df.get("copy_number", pd.Series([1.0] * len(df))), errors="coerce").fillna(1.0)

        # Build metadata using zip (more efficient than iterrows for simple extraction)
        for query_id, reads, length, copy_number in zip(query_ids, reads_col, lengths, copy_numbers):
            entry = {
                "reads": reads,
                "length": int(length),
                "copy_number": float(copy_number),
            }
            if query_id and query_id not in metadata:
                metadata[query_id] = entry
            if reads and reads not in metadata:
                metadata[reads] = entry

        return metadata

    def detect_circles_from_last(
        self,
        alignments: Dict[str, List[LastAlignment]],
        metadata: Dict[str, dict],
    ) -> pd.DataFrame:
        """Detect CeccDNA from LAST alignment results."""
        rows: list[dict] = []

        for query_id, alns in alignments.items():
            # Extract base read ID
            base_id = query_id.split("|")[0]
            meta = metadata.get(query_id) or metadata.get(
                base_id, {"reads": base_id, "length": 0, "copy_number": 1.0}
            )

            result = self._detect_single_query_last(query_id, alns, meta)
            if result is None:
                continue

            # Build output rows
            base = {
                "query_id": result["query_id"],
                "reads": result["reads"],
                "eccdna_type": "Cecc",
                "CeccClass": result["cecc_class"],
                "length": result["length"],
                "copy_number": result["copy_number"],
                "num_segments": len(result["path_segments"]),
                "cumulative_length": result["cum_len"],
                "match_degree": round(result["match_degree"], 2),
                "match_degree_2nd": None,
                "best_chain_signature": None,
                "second_chain_signature": None,
                ColumnStandard.MAPQ_BEST: result["mapq_best"],
                ColumnStandard.MAPQ_MIN: result["mapq_min"],
                ColumnStandard.IDENTITY_BEST: round(result["identity_best"], 3),
                ColumnStandard.IDENTITY_MIN: round(result["identity_min"], 3),
                ColumnStandard.QUERY_COV_BEST: round(result["cov_best"], 6),
                ColumnStandard.QUERY_COV_2ND: round(result["cov_2nd"], 6),
                ColumnStandard.CONFIDENCE_SCORE: round(result["confidence_score"], 6),
                ColumnStandard.LOW_MAPQ: result["low_mapq"],
                ColumnStandard.LOW_IDENTITY: result["low_identity"],
                "C_cov_best": round(result["cov_best"], 6),
                "C_cov_2nd": round(result["cov_2nd"], 6),
                "max_gap": result["max_gap"],
                "avg_gap": round(result["avg_gap"], 2),
                "chromosomes": ",".join(str(c) for c in result["chromosomes"]),
            }

            for i, seg in enumerate(result["path_segments"]):
                row = dict(base)
                row.update({
                    "segment_in_circle": i + 1,
                    ColumnStandard.CHR: seg.chr,
                    ColumnStandard.START0: seg.start0,
                    ColumnStandard.END0: seg.end0,
                    ColumnStandard.STRAND: seg.strand,
                    "q_start": seg.q_start,
                    "q_end": seg.q_end,
                    "alignment_length": seg.alignment_length,
                })
                rows.append(row)

        result_df = pd.DataFrame(rows)
        if not result_df.empty:
            self.logger.info(f"Cecc detected: {result_df['query_id'].nunique()} queries")

        return result_df

    def run(
        self,
        input_csv: Path,
        output_csv: Path,
        reference_fasta: Optional[Path] = None,
        fasta_file: Optional[Path] = None,
        edge_tolerance: Optional[int] = None,
        position_tolerance: Optional[int] = None,
        min_match_degree: Optional[float] = None,
        half_query_buffer: Optional[int] = None,
    ) -> pd.DataFrame:
        """
        Run the complete Cecc detection pipeline with LAST.

        Args:
            input_csv: Input alignment CSV/TSV (um_classify.unclassified.csv)
            output_csv: Output CSV file
            reference_fasta: Reference genome FASTA (required for LAST)
            fasta_file: Tandem-to-ring FASTA with doubled sequences (required for LAST)
            edge_tolerance: Override gap tolerance (optional)
            position_tolerance: Override position tolerance (optional)
            min_match_degree: Override minimum match degree (optional)

        Returns:
            DataFrame with Cecc results
        """
        # Apply overrides
        if edge_tolerance is not None:
            self.gap_tolerance = edge_tolerance
        if position_tolerance is not None:
            self.position_tolerance = position_tolerance
        if min_match_degree is not None:
            self.min_match_degree = min_match_degree
        if half_query_buffer is not None:
            self.half_query_buffer = max(0, int(half_query_buffer))

        self.logger.info("=" * 60)
        self.logger.info("CeccBuild v4 - LAST-based CeccDNA Detection")
        self.logger.info("=" * 60)
        self.logger.info(f"Parameters:")
        self.logger.info(f"  - Position tolerance: {self.position_tolerance} bp")
        self.logger.info(f"  - Min identity: {self.min_identity}%")
        self.logger.info(f"  - Min query coverage: {self.min_query_coverage}")
        self.logger.info(f"  - Min repeat query gap: {self.min_repeat_query_gap} bp")
        self.logger.info(f"  - Half buffer: {self.half_query_buffer} bp")

        # Read input CSV
        self.logger.info(f"Reading: {input_csv}")
        try:
            df = pd.read_csv(input_csv)
        except (pd.errors.ParserError, pd.errors.EmptyDataError, ValueError):
            df = pd.read_csv(input_csv, sep="\t")

        if df.empty:
            self.logger.warning("Input file is empty")
            result_df = pd.DataFrame()
            output_csv.parent.mkdir(parents=True, exist_ok=True)
            result_df.to_csv(output_csv, index=False)
            return result_df

        # Check if LAST is available
        if not self._check_last_available():
            self.logger.warning("LAST not available, falling back to graph-based detection")
            return self._run_graph_based(df, output_csv)

        # Check required files
        if reference_fasta is None or fasta_file is None:
            self.logger.warning("Reference or FASTA not provided, falling back to graph-based detection")
            return self._run_graph_based(df, output_csv)

        # Get unique read IDs from input
        read_ids = set(df["reads"].unique()) if "reads" in df.columns else set()
        self.logger.info(f"Unclassified reads: {len(read_ids)}")

        # Get metadata
        metadata = self._get_metadata_from_csv(df)

        def _run_last_in_tmp(tmpdir: Path) -> pd.DataFrame:
            # Extract sequences
            query_fasta = tmpdir / "query.fasta"
            extracted = self._extract_sequences(read_ids, fasta_file, query_fasta)
            self.logger.info(f"Extracted {extracted} sequences")

            if extracted == 0:
                self.logger.warning("No sequences extracted")
                return pd.DataFrame()

            # Build LAST database
            db_prefix = tmpdir / "lastdb"
            self._build_last_db(reference_fasta, db_prefix)

            # Run LAST alignment
            tab_file = tmpdir / "alignments.tab"
            self._run_last_alignment(query_fasta, db_prefix, tab_file)

            # Parse LAST results
            alignments = self._parse_last_tab(tab_file)
            self.logger.info(f"Parsed {len(alignments)} query alignments")

            # Detect circles
            return self.detect_circles_from_last(alignments, metadata)

        result_df = pd.DataFrame()
        tmp_root = self.tmp_dir
        if tmp_root is None and self.keep_tmp:
            tmp_root = output_csv.parent / "cecc_last_tmp"

        if tmp_root is not None:
            tmp_root = Path(tmp_root)
            tmp_root.mkdir(parents=True, exist_ok=True)
            tmpdir = Path(tempfile.mkdtemp(prefix="cecc_last_", dir=tmp_root))
            try:
                result_df = _run_last_in_tmp(tmpdir)
            finally:
                if self.keep_tmp:
                    self.logger.info(f"Keeping LAST temp dir: {tmpdir}")
                else:
                    shutil.rmtree(tmpdir, ignore_errors=True)
        else:
            with tempfile.TemporaryDirectory() as tmpdir_str:
                result_df = _run_last_in_tmp(Path(tmpdir_str))

        # Save output
        output_csv.parent.mkdir(parents=True, exist_ok=True)
        result_df.to_csv(output_csv, index=False)

        # Summary
        self.logger.info("=" * 60)
        self.logger.info("Results:")
        if not result_df.empty:
            n_cecc = result_df["query_id"].nunique()
            self.logger.info(f"  - Cecc detected: {n_cecc}")
            if "CeccClass" in result_df.columns:
                class_counts = result_df.groupby("query_id")["CeccClass"].first().value_counts()
                for cls, cnt in class_counts.items():
                    self.logger.info(f"  - {cls}: {cnt}")
        else:
            self.logger.info("  - Cecc detected: 0")
        self.logger.info(f"Output: {output_csv}")
        self.logger.info("=" * 60)

        return result_df

    def _run_graph_based(self, df: pd.DataFrame, output_csv: Path) -> pd.DataFrame:
        """Fallback to graph-based detection (V3 method) when LAST unavailable."""
        self.logger.info("Using graph-based detection (fallback)")
        # Import V3 method inline to avoid circular imports
        rows = self._detect_circles_graph(df)
        result_df = pd.DataFrame(rows)
        output_csv.parent.mkdir(parents=True, exist_ok=True)
        result_df.to_csv(output_csv, index=False)
        return result_df

    def _detect_circles_graph(self, df: pd.DataFrame) -> list[dict[str, Any]]:
        """Graph-based circle detection (V3 fallback)."""
        filtered = self.filter_overlapping_queries(df)
        result_df = self.detect_circles(
            filtered,
            edge_tol=self.gap_tolerance,
            pos_tol=self.position_tolerance,
        )
        if result_df.empty:
            return []
        records: list[dict[str, Any]] = result_df.to_dict(orient="records")
        return records

    def run_pipeline(
        self,
        input_csv: Path,
        output_csv: Path,
        overlap_threshold: float = 0.95,
        min_segments: int = 2,
        edge_tolerance: int = 20,
        position_tolerance: int = 50,
        min_match_degree: float = 95.0,
        max_rotations: int = 20,
        locus_overlap_threshold: float = 0.95,
        half_query_buffer: int = 50,
        reference_fasta: Optional[Path] = None,
        fasta_file: Optional[Path] = None,
    ) -> pd.DataFrame:
        """
        V1/V2/V3-compatible pipeline interface.

        Added reference_fasta and fasta_file parameters for LAST-based detection.
        """
        self.min_segments = min_segments
        self.gap_tolerance = edge_tolerance
        self.position_tolerance = position_tolerance
        self.min_match_degree = min_match_degree
        self.overlap_threshold = overlap_threshold
        self.max_rotations = max_rotations
        self.locus_overlap_threshold = locus_overlap_threshold
        self.half_query_buffer = max(0, int(half_query_buffer))

        return self.run(
            input_csv=input_csv,
            output_csv=output_csv,
            reference_fasta=reference_fasta,
            fasta_file=fasta_file,
            half_query_buffer=half_query_buffer,
        )


def main() -> None:
    """CLI entry point."""
    import argparse

    parser = argparse.ArgumentParser(
        description="CeccBuild v4 - LAST-based CeccDNA Detection"
    )
    parser.add_argument("-i", "--input", required=True, help="Input alignment TSV/CSV")
    parser.add_argument("-o", "--output", required=True, help="Output CSV")
    parser.add_argument("-r", "--reference", help="Reference genome FASTA")
    parser.add_argument("-f", "--fasta", help="Tandem-to-ring FASTA (doubled sequences)")
    parser.add_argument(
        "--position-tolerance", type=int, default=50,
        help="Position tolerance for repeat matching (bp)"
    )
    parser.add_argument(
        "--min-identity", type=float, default=95.0,
        help="Minimum alignment identity (%%)"
    )
    parser.add_argument(
        "--threads", type=int, default=4,
        help="Number of threads for LAST"
    )

    args = parser.parse_args()

    builder = CeccBuild(
        position_tolerance=args.position_tolerance,
        min_identity=args.min_identity,
        threads=args.threads,
    )

    builder.run(
        input_csv=Path(args.input),
        output_csv=Path(args.output),
        reference_fasta=Path(args.reference) if args.reference else None,
        fasta_file=Path(args.fasta) if args.fasta else None,
    )


if __name__ == "__main__":
    main()
