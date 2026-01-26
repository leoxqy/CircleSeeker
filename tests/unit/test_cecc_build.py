from pathlib import Path
import shutil
import sys

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.modules.cecc_build import CeccBuild
from circleseeker.utils.column_standards import ColumnStandard

import numpy as np
from circleseeker.modules.cecc_build import LastAlignment, AlignmentSegment

# Check if all LAST tools are available for integration tests
HAS_LAST = all(
    shutil.which(t) is not None for t in ("lastal", "lastdb", "last-split")
)


def _random_dna(n: int, seed: int = 42) -> str:
    """Generate reproducible random DNA sequence."""
    import random
    rng = random.Random(seed)
    return "".join(rng.choices("ACGT", k=n))


def _write_fasta(path, records: dict[str, str]) -> None:
    """Write dict of {header: sequence} to FASTA file."""
    with open(path, "w") as f:
        for header, seq in records.items():
            f.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")


def _sorted(df: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    if df.empty:
        return df.copy()
    return df.sort_values(keys).reset_index(drop=True)


def _normalise_roles(df: pd.DataFrame) -> pd.DataFrame:
    if "segment_role" in df.columns:
        df = df.copy()
        df["segment_role"] = df["segment_role"].astype(str)
    return df


@pytest.mark.skipif(not HAS_LAST, reason="LAST aligner tools (lastal/lastdb/last-split) not installed")
def test_cecc_build_intra_chr(tmp_path):
    """Doubled chimeric read from two chr1 loci → should detect Cecc-IntraChr."""
    # Reference genome: chr1 with 5kb random DNA
    chr1_seq = _random_dna(5000, seed=42)
    ref_fasta = tmp_path / "ref.fasta"
    _write_fasta(ref_fasta, {"chr1": chr1_seq})

    # Extract two distant 500bp regions from chr1
    seg_a = chr1_seq[500:1000]    # chr1:500-1000
    seg_b = chr1_seq[3000:3500]   # chr1:3000-3500

    # Doubled chimeric query: [A][B][A][B] = 2000bp
    doubled_seq = seg_a + seg_b + seg_a + seg_b

    query_fasta = tmp_path / "query.fasta"
    _write_fasta(query_fasta, {"read_intra|1|1000|1|circular": doubled_seq})

    # Input CSV with metadata
    input_csv = tmp_path / "input.csv"
    pd.DataFrame([{
        "reads": "read_intra",
        "query_id": "read_intra",
        "length": 1000,
        "copy_number": 1.0,
    }]).to_csv(input_csv, index=False)

    output_csv = tmp_path / "cecc_intra.csv"
    builder = CeccBuild(
        min_query_coverage=0.3,
        min_segments=2,
        min_match_degree=50.0,
        position_tolerance=200,
        min_identity=90.0,
        min_repeat_query_gap=100,
        half_query_buffer=50,
        tmp_dir=tmp_path,
    )
    result = builder.run_pipeline(
        input_csv=input_csv,
        output_csv=output_csv,
        reference_fasta=ref_fasta,
        fasta_file=query_fasta,
    )

    assert not result.empty, "Expected CeccDNA detection for doubled chimeric read"
    assert output_csv.exists()
    assert "query_id" in result.columns

    # All detected loci should be on chr1 → IntraChr
    if "CeccClass" in result.columns:
        assert result["CeccClass"].iloc[0] == "Cecc-IntraChr"


def test_detect_genomic_overlaps_allows_boundary_jitter():
    builder = CeccBuild()
    segments = pd.DataFrame(
        {
            ColumnStandard.CHR: ["chr1", "chr1"],
            ColumnStandard.START0: [100, 101],
            ColumnStandard.END0: [200, 200],
            ColumnStandard.STRAND: ["+", "+"],
        }
    )

    assert builder.detect_genomic_overlaps_sweepline(segments) is False


def test_detect_genomic_overlaps_ignores_opposite_strands():
    builder = CeccBuild()
    segments = pd.DataFrame(
        {
            ColumnStandard.CHR: ["chr1", "chr1"],
            ColumnStandard.START0: [100, 100],
            ColumnStandard.END0: [200, 200],
            ColumnStandard.STRAND: ["+", "-"],
        }
    )

    assert builder.detect_genomic_overlaps_sweepline(segments) is False


@pytest.mark.skipif(not HAS_LAST, reason="LAST aligner tools (lastal/lastdb/last-split) not installed")
def test_cecc_build_inter_chr(tmp_path):
    """Doubled chimeric read from chr1 + chr2 → should detect Cecc-InterChr."""
    # Reference genome: two chromosomes with distinct random DNA
    chr1_seq = _random_dna(5000, seed=100)
    chr2_seq = _random_dna(5000, seed=200)
    ref_fasta = tmp_path / "ref.fasta"
    _write_fasta(ref_fasta, {"chr1": chr1_seq, "chr2": chr2_seq})

    # Extract 500bp segments from different chromosomes
    seg_chr1 = chr1_seq[500:1000]   # chr1:500-1000
    seg_chr2 = chr2_seq[2000:2500]  # chr2:2000-2500

    # Doubled chimeric query: [chr1_seg][chr2_seg][chr1_seg][chr2_seg] = 2000bp
    doubled_seq = seg_chr1 + seg_chr2 + seg_chr1 + seg_chr2

    query_fasta = tmp_path / "query.fasta"
    _write_fasta(query_fasta, {"read_inter|1|1000|1|circular": doubled_seq})

    # Input CSV with metadata
    input_csv = tmp_path / "input.csv"
    pd.DataFrame([{
        "reads": "read_inter",
        "query_id": "read_inter",
        "length": 1000,
        "copy_number": 1.0,
    }]).to_csv(input_csv, index=False)

    output_csv = tmp_path / "cecc_inter.csv"
    builder = CeccBuild(
        min_query_coverage=0.3,
        min_segments=2,
        min_match_degree=50.0,
        position_tolerance=200,
        min_identity=90.0,
        min_repeat_query_gap=100,
        half_query_buffer=50,
        tmp_dir=tmp_path,
    )
    result = builder.run_pipeline(
        input_csv=input_csv,
        output_csv=output_csv,
        reference_fasta=ref_fasta,
        fasta_file=query_fasta,
    )

    assert not result.empty, "Expected CeccDNA detection for inter-chr chimeric read"
    assert output_csv.exists()

    # Segments span chr1 and chr2 → InterChr
    if "CeccClass" in result.columns:
        assert result["CeccClass"].iloc[0] == "Cecc-InterChr"
    if "chromosomes" in result.columns:
        chroms = set(result["chromosomes"].iloc[0].split(","))
        assert len(chroms) >= 2


class TestClamp01:
    def test_within_range(self):
        assert CeccBuild._clamp01(0.5) == 0.5

    def test_negative(self):
        assert CeccBuild._clamp01(-0.5) == 0.0

    def test_above_one(self):
        assert CeccBuild._clamp01(1.5) == 1.0

    def test_boundary_zero(self):
        assert CeccBuild._clamp01(0.0) == 0.0

    def test_non_numeric(self):
        assert CeccBuild._clamp01("abc") == 0.0


class TestAsFraction:
    def test_fraction(self):
        assert CeccBuild._as_fraction(0.95) == 0.95

    def test_percentage(self):
        result = CeccBuild._as_fraction(95.0)
        assert abs(result - 0.95) < 0.001

    def test_zero(self):
        assert CeccBuild._as_fraction(0.0) == 0.0

    def test_non_numeric(self):
        assert CeccBuild._as_fraction("abc") == 0.0

    def test_over_100(self):
        result = CeccBuild._as_fraction(150.0)
        assert result == 1.0


class TestNormMapq:
    def test_max_mapq(self):
        result = CeccBuild._norm_mapq(60)
        assert abs(result - 1.0) < 0.001

    def test_zero_mapq(self):
        assert CeccBuild._norm_mapq(0) == 0.0

    def test_non_numeric(self):
        assert CeccBuild._norm_mapq("abc") == 0.0


class TestNormIdentity:
    def test_100_percent(self):
        assert CeccBuild._norm_identity(100.0) == 1.0

    def test_90_percent(self):
        assert CeccBuild._norm_identity(90.0) == 0.0

    def test_below_90(self):
        assert CeccBuild._norm_identity(80.0) == 0.0


class TestGeomMean:
    def test_single_value(self):
        assert CeccBuild._geom_mean([0.5]) == 0.5

    def test_two_values(self):
        result = CeccBuild._geom_mean([4.0, 9.0])
        assert abs(result - 6.0) < 0.001

    def test_with_zero(self):
        assert CeccBuild._geom_mean([0.0, 1.0]) == 0.0

    def test_empty(self):
        assert CeccBuild._geom_mean([]) == 0.0


class TestMergeIntervals:
    def test_non_overlapping(self):
        result = CeccBuild._merge_intervals([(1, 3), (5, 7)])
        assert result == [(1, 3), (5, 7)]

    def test_overlapping(self):
        result = CeccBuild._merge_intervals([(1, 5), (3, 7)])
        assert result == [(1, 7)]

    def test_adjacent(self):
        result = CeccBuild._merge_intervals([(1, 3), (3, 5)])
        assert result == [(1, 5)]

    def test_empty(self):
        assert CeccBuild._merge_intervals([]) == []

    def test_unsorted_input(self):
        result = CeccBuild._merge_intervals([(5, 7), (1, 3)])
        assert result == [(1, 3), (5, 7)]


class TestCoverageHelpers:
    def test_coverage_length(self):
        alns = [
            LastAlignment("chr1", 0, 100, 0, 50, 200, 100, 95.0),
            LastAlignment("chr1", 0, 100, 40, 100, 200, 100, 95.0),
        ]
        result = CeccBuild._coverage_length(alns)
        assert result == 100  # 0-100 merged

    def test_full_coverage(self):
        alns = [LastAlignment("chr1", 0, 100, 0, 200, 200, 100, 95.0)]
        result = CeccBuild._coverage_fraction(alns, 200)
        assert result == 1.0

    def test_zero_query_len(self):
        alns = [LastAlignment("chr1", 0, 100, 0, 50, 0, 100, 95.0)]
        result = CeccBuild._coverage_fraction(alns, 0)
        assert result == 0.0


class TestParseSegment:
    @pytest.fixture
    def builder(self):
        return CeccBuild()

    def test_blast_columns(self, builder):
        row = pd.Series({
            "subject_id": "chr1", "q_start": 0, "q_end": 100,
            "s_start": 501, "s_end": 600, "strand": "+", "alignment_length": 100,
        })
        seg = builder._parse_segment(row)
        assert seg.chr == "chr1"
        assert seg.start0 == 500  # s_start-1
        assert seg.end0 == 600

    def test_strand_standardization(self, builder):
        row = pd.Series({
            "subject_id": "chr1", "q_start": 0, "q_end": 100,
            "s_start": 101, "s_end": 200, "strand": "plus", "alignment_length": 100,
        })
        seg = builder._parse_segment(row)
        assert seg.strand == "+"

    def test_minus_strand(self, builder):
        row = pd.Series({
            "subject_id": "chr1", "q_start": 0, "q_end": 100,
            "s_start": 101, "s_end": 200, "strand": "minus", "alignment_length": 100,
        })
        seg = builder._parse_segment(row)
        assert seg.strand == "-"

    def test_missing_coordinates_raises(self, builder):
        row = pd.Series({"subject_id": "chr1", "q_start": 0, "q_end": 100})
        with pytest.raises(ValueError):
            builder._parse_segment(row)


class TestSameLocus:
    @pytest.fixture
    def builder(self):
        return CeccBuild()

    def test_same_position(self, builder):
        a = AlignmentSegment("chr1", 100, 200, "+", 0, 100, 100)
        b = AlignmentSegment("chr1", 100, 200, "+", 200, 300, 100)
        assert builder._same_locus(a, b, 0.95) is True

    def test_different_chrom(self, builder):
        a = AlignmentSegment("chr1", 100, 200, "+", 0, 100, 100)
        b = AlignmentSegment("chr2", 100, 200, "+", 200, 300, 100)
        assert builder._same_locus(a, b, 0.95) is False

    def test_different_strand(self, builder):
        a = AlignmentSegment("chr1", 100, 200, "+", 0, 100, 100)
        b = AlignmentSegment("chr1", 100, 200, "-", 200, 300, 100)
        assert builder._same_locus(a, b, 0.95) is False

    def test_below_threshold(self, builder):
        a = AlignmentSegment("chr1", 100, 200, "+", 0, 100, 100)
        b = AlignmentSegment("chr1", 150, 300, "+", 200, 300, 100)
        assert builder._same_locus(a, b, 0.95) is False


class TestDetectDoubledRepeatPattern:
    @pytest.fixture
    def builder(self):
        return CeccBuild(min_query_coverage=0.5, min_repeat_query_gap=50, half_query_buffer=10)

    def test_valid_detection(self, builder):
        alns = [
            LastAlignment("chr1", 1000, 1500, 0, 500, 2000, 500, 95.0, "+"),
            LastAlignment("chr1", 1000, 1500, 1000, 1500, 2000, 500, 95.0, "+"),
        ]
        is_cecc, reason, pairs = builder._detect_doubled_repeat_pattern(alns)
        assert is_cecc is True
        assert reason == "doubled_repeat_pattern"

    def test_insufficient_alignments(self, builder):
        alns = [LastAlignment("chr1", 0, 100, 0, 100, 200, 100, 95.0)]
        is_cecc, reason, _ = builder._detect_doubled_repeat_pattern(alns)
        assert is_cecc is False
        assert "insufficient" in reason

    def test_low_coverage(self, builder):
        builder_strict = CeccBuild(min_query_coverage=0.99)
        alns = [
            LastAlignment("chr1", 0, 10, 0, 10, 10000, 10, 95.0, "+"),
            LastAlignment("chr1", 0, 10, 5000, 5010, 10000, 10, 95.0, "+"),
        ]
        is_cecc, reason, _ = builder_strict._detect_doubled_repeat_pattern(alns)
        assert is_cecc is False

    def test_no_matching_position(self, builder):
        alns = [
            LastAlignment("chr1", 1000, 1100, 0, 100, 2000, 100, 95.0, "+"),
            LastAlignment("chr2", 5000, 5100, 1500, 1600, 2000, 100, 95.0, "+"),
        ]
        is_cecc, _, _ = builder._detect_doubled_repeat_pattern(alns)
        assert is_cecc is False


class TestAssignLocusId:
    @pytest.fixture
    def builder(self):
        return CeccBuild()

    def test_single_chrom_multiple_loci(self, builder):
        alns = [
            LastAlignment("chr1", 100, 200, 0, 100, 1000, 100, 95.0),
            LastAlignment("chr1", 5000, 5100, 500, 600, 1000, 100, 95.0),
        ]
        aln_to_locus, locus_info = builder._assign_locus_id(alns, merge_distance=500)
        assert len(locus_info) == 2

    def test_multi_chrom(self, builder):
        alns = [
            LastAlignment("chr1", 100, 200, 0, 100, 1000, 100, 95.0),
            LastAlignment("chr2", 100, 200, 500, 600, 1000, 100, 95.0),
        ]
        aln_to_locus, locus_info = builder._assign_locus_id(alns)
        assert len(locus_info) == 2

    def test_empty(self, builder):
        aln_to_locus, locus_info = builder._assign_locus_id([])
        assert len(aln_to_locus) == 0
        assert len(locus_info) == 0


class TestCheckStrandClosure:
    @pytest.fixture
    def builder(self):
        return CeccBuild()

    def test_all_same_even(self, builder):
        assert builder._check_strand_closure(["+_+", "+_+"]) is True

    def test_two_flips_even(self, builder):
        assert builder._check_strand_closure(["+_-", "-_+"]) is True

    def test_single_flip_odd(self, builder):
        assert builder._check_strand_closure(["+_+", "+_-"]) is False


class TestFilterOverlappingQueries:
    @pytest.fixture
    def builder(self):
        return CeccBuild()

    def test_keeps_non_overlapping(self, builder):
        df = pd.DataFrame({
            "query_id": ["q1", "q1"],
            "subject_id": ["chr1", "chr1"],
            "q_start": [0, 100],
            "q_end": [100, 200],
            "s_start": [501, 1001],
            "s_end": [600, 1100],
            "strand": ["+", "+"],
            "alignment_length": [100, 100],
        })
        result = builder.filter_overlapping_queries(df)
        assert "q1" in result["query_id"].values

    def test_empty(self, builder):
        result = builder.filter_overlapping_queries(pd.DataFrame())
        assert result.empty

    def test_single_segment_kept(self, builder):
        df = pd.DataFrame({
            "query_id": ["q1"],
            "subject_id": ["chr1"],
            "q_start": [0],
            "q_end": [100],
            "s_start": [501],
            "s_end": [600],
            "strand": ["+"],
            "alignment_length": [100],
        })
        result = builder.filter_overlapping_queries(df)
        assert len(result) == 1


class TestDetectCirclesFromLast:
    @pytest.fixture
    def builder(self):
        return CeccBuild(min_query_coverage=0.5, min_segments=2, min_match_degree=50.0,
                         min_repeat_query_gap=50, half_query_buffer=10)

    def test_basic_cecc(self, builder):
        alns = {
            "q1": [
                LastAlignment("chr1", 1000, 1500, 0, 500, 2000, 500, 95.0, "+"),
                LastAlignment("chr2", 2000, 2500, 500, 1000, 2000, 500, 95.0, "+"),
                LastAlignment("chr1", 1000, 1500, 1000, 1500, 2000, 500, 95.0, "+"),
                LastAlignment("chr2", 2000, 2500, 1500, 2000, 2000, 500, 95.0, "+"),
            ]
        }
        metadata = {"q1": {"reads": "read1", "length": 1000, "copy_number": 1.0}}
        result = builder.detect_circles_from_last(alns, metadata)
        assert not result.empty

    def test_no_detection(self, builder):
        alns = {
            "q1": [
                LastAlignment("chr1", 1000, 1100, 0, 100, 200, 100, 95.0, "+"),
            ]
        }
        metadata = {"q1": {"reads": "read1", "length": 100, "copy_number": 1.0}}
        result = builder.detect_circles_from_last(alns, metadata)
        assert result.empty

    def test_intra_chr(self, builder):
        alns = {
            "q1": [
                LastAlignment("chr1", 1000, 1500, 0, 500, 2000, 500, 95.0, "+"),
                LastAlignment("chr1", 3000, 3500, 500, 1000, 2000, 500, 95.0, "+"),
                LastAlignment("chr1", 1000, 1500, 1000, 1500, 2000, 500, 95.0, "+"),
                LastAlignment("chr1", 3000, 3500, 1500, 2000, 2000, 500, 95.0, "+"),
            ]
        }
        metadata = {"q1": {"reads": "read1", "length": 1000, "copy_number": 1.0}}
        result = builder.detect_circles_from_last(alns, metadata)
        # Check for IntraChr if detected
        if not result.empty:
            assert "IntraChr" in result["CeccClass"].iloc[0]


# ========== MAF Parser Unit Tests ========== #


class TestParseLastSplitMaf:
    """Unit tests for _parse_last_split_maf coordinate parsing.

    These tests do NOT require LAST tools — they only test the MAF parser logic
    using synthetic MAF content.
    """

    @pytest.fixture
    def builder(self):
        return CeccBuild()

    def _write_maf(self, path, content: str) -> None:
        path.write_text(content)

    def test_forward_strand_basic(self, builder, tmp_path):
        """Forward strand: query_start = start, query_end = start + size."""
        maf = tmp_path / "fwd.maf"
        self._write_maf(maf, (
            "a score=1000\n"
            "s chr1 100 500 + 5000 ACGT\n"
            "s read1 0 200 + 1000 ACGT\n"
        ))
        result = builder._parse_last_split_maf(maf)
        assert "read1" in result
        aln = result["read1"][0]
        assert aln.chrom == "chr1"
        assert aln.ref_start == 100
        assert aln.ref_end == 600  # 100 + 500
        assert aln.query_start == 0
        assert aln.query_end == 200  # 0 + 200
        assert aln.query_len == 1000
        assert aln.strand == "+"
        assert aln.score == 1000

    def test_reverse_strand_coordinate_conversion(self, builder, tmp_path):
        """Reverse strand: query_start = seq_size - start - size, query_end = seq_size - start."""
        maf = tmp_path / "rev.maf"
        # Query: start=50, size=100, seq_size=500, strand=-
        # Expected: query_start = 500 - 50 - 100 = 350, query_end = 500 - 50 = 450
        self._write_maf(maf, (
            "a score=800\n"
            "s chr2 200 300 + 10000 ACGT\n"
            "s read2 50 100 - 500 ACGT\n"
        ))
        result = builder._parse_last_split_maf(maf)
        assert "read2" in result
        aln = result["read2"][0]
        assert aln.query_start == 350
        assert aln.query_end == 450
        assert aln.query_len == 500
        assert aln.strand == "-"
        assert aln.chrom == "chr2"
        assert aln.ref_start == 200
        assert aln.ref_end == 500  # 200 + 300

    def test_multiple_alignment_blocks(self, builder, tmp_path):
        """Multiple alignment blocks for the same query."""
        maf = tmp_path / "multi.maf"
        self._write_maf(maf, (
            "a score=500\n"
            "s chr1 0 100 + 5000 ACGT\n"
            "s read1 0 100 + 2000 ACGT\n"
            "\n"
            "a score=600\n"
            "s chr2 1000 200 + 8000 ACGT\n"
            "s read1 500 200 + 2000 ACGT\n"
        ))
        result = builder._parse_last_split_maf(maf)
        assert len(result["read1"]) == 2
        assert result["read1"][0].chrom == "chr1"
        assert result["read1"][0].score == 500
        assert result["read1"][1].chrom == "chr2"
        assert result["read1"][1].score == 600

    def test_multiple_queries(self, builder, tmp_path):
        """Alignments from different queries are grouped correctly."""
        maf = tmp_path / "multi_query.maf"
        self._write_maf(maf, (
            "a score=100\n"
            "s chr1 0 100 + 5000 ACGT\n"
            "s readA 0 100 + 1000 ACGT\n"
            "\n"
            "a score=200\n"
            "s chr1 500 100 + 5000 ACGT\n"
            "s readB 0 100 + 800 ACGT\n"
        ))
        result = builder._parse_last_split_maf(maf)
        assert len(result["readA"]) == 1
        assert len(result["readB"]) == 1
        assert result["readA"][0].score == 100
        assert result["readB"][0].score == 200

    def test_score_extraction(self, builder, tmp_path):
        """Score is correctly extracted from 'a score=NNN' line."""
        maf = tmp_path / "score.maf"
        self._write_maf(maf, (
            "a score=12345 mismap=1e-5\n"
            "s chr1 0 50 + 5000 ACGT\n"
            "s read1 0 50 + 200 ACGT\n"
        ))
        result = builder._parse_last_split_maf(maf)
        assert result["read1"][0].score == 12345

    def test_score_missing_defaults_to_zero(self, builder, tmp_path):
        """Missing score defaults to 0."""
        maf = tmp_path / "no_score.maf"
        self._write_maf(maf, (
            "a mismap=1e-5\n"
            "s chr1 0 50 + 5000 ACGT\n"
            "s read1 0 50 + 200 ACGT\n"
        ))
        result = builder._parse_last_split_maf(maf)
        assert result["read1"][0].score == 0

    def test_malformed_s_line_skipped(self, builder, tmp_path):
        """Lines with fewer than 7 fields are skipped."""
        maf = tmp_path / "malformed.maf"
        self._write_maf(maf, (
            "a score=100\n"
            "s chr1 0 50\n"  # Too few fields — skipped
            "s read1 0 50 + 200 ACGT\n"  # No ref_info yet, so this becomes ref
        ))
        result = builder._parse_last_split_maf(maf)
        # The malformed ref line is skipped; read1 becomes the ref line.
        # No complete alignment pair → empty result.
        assert len(result) == 0

    def test_invalid_numeric_fields_skipped(self, builder, tmp_path):
        """Non-numeric start/size fields cause the alignment to be skipped."""
        maf = tmp_path / "invalid_num.maf"
        self._write_maf(maf, (
            "a score=100\n"
            "s chr1 abc 50 + 5000 ACGT\n"  # 'abc' is not a valid start
            "s read1 0 50 + 200 ACGT\n"
        ))
        result = builder._parse_last_split_maf(maf)
        # First s-line fails parsing → ref_info never set → second line
        # becomes ref → no query → empty
        assert len(result) == 0

    def test_empty_file(self, builder, tmp_path):
        """Empty file returns empty dict."""
        maf = tmp_path / "empty.maf"
        self._write_maf(maf, "")
        result = builder._parse_last_split_maf(maf)
        assert len(result) == 0

    def test_comment_and_blank_lines_ignored(self, builder, tmp_path):
        """Comment lines and blank lines are ignored."""
        maf = tmp_path / "comments.maf"
        self._write_maf(maf, (
            "# last version 1234\n"
            "#\n"
            "\n"
            "a score=500\n"
            "s chr1 0 100 + 5000 ACGT\n"
            "\n"
            "s read1 10 100 + 2000 ACGT\n"
        ))
        result = builder._parse_last_split_maf(maf)
        assert len(result["read1"]) == 1
        assert result["read1"][0].query_start == 10
        assert result["read1"][0].query_end == 110

    def test_identity_default_95(self, builder, tmp_path):
        """Identity defaults to 95.0 since last-split doesn't provide it."""
        maf = tmp_path / "identity.maf"
        self._write_maf(maf, (
            "a score=100\n"
            "s chr1 0 100 + 5000 ACGT\n"
            "s read1 0 100 + 1000 ACGT\n"
        ))
        result = builder._parse_last_split_maf(maf)
        assert result["read1"][0].identity == 95.0

    def test_reverse_strand_boundary_zero_start(self, builder, tmp_path):
        """Reverse strand with start=0: entire sequence aligned from end."""
        maf = tmp_path / "rev_zero.maf"
        # start=0, size=100, seq_size=100: full reverse complement
        # Expected: query_start = 100 - 0 - 100 = 0, query_end = 100 - 0 = 100
        self._write_maf(maf, (
            "a score=100\n"
            "s chr1 500 100 + 5000 ACGT\n"
            "s read1 0 100 - 100 ACGT\n"
        ))
        result = builder._parse_last_split_maf(maf)
        aln = result["read1"][0]
        assert aln.query_start == 0
        assert aln.query_end == 100
        assert aln.strand == "-"

    def test_reverse_strand_partial_alignment(self, builder, tmp_path):
        """Reverse strand partial alignment at the end of the query."""
        maf = tmp_path / "rev_partial.maf"
        # start=0, size=50, seq_size=200: first 50bp of reverse complement
        # Expected: query_start = 200 - 0 - 50 = 150, query_end = 200 - 0 = 200
        self._write_maf(maf, (
            "a score=100\n"
            "s chr1 0 100 + 5000 ACGT\n"
            "s read1 0 50 - 200 ACGT\n"
        ))
        result = builder._parse_last_split_maf(maf)
        aln = result["read1"][0]
        assert aln.query_start == 150
        assert aln.query_end == 200
