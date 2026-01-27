"""Tests for ecc_unify module."""

from pathlib import Path
import sys
import pytest
import pandas as pd
import json
import tempfile

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.modules.ecc_unify import (
    merge_eccdna_tables,
    parse_region,
    parse_chimeric_regions,
    parse_chimeric_regions_with_strand,
    reciprocal_overlap_ok,
    build_chr_index,
    build_cecc_segment_index,
    find_redundant_simple,
    find_redundant_chimeric,
    _segments_match,
    _segments_match_with_strand,
    _is_subset_of,
    _is_subset_of_with_strand,
    prepare_inferred_simple,
    prepare_inferred_chimeric,
    renumber_eccdna,
    generate_overlap_report,
    generate_overlap_stats_json,
)


class TestEccUnify:
    """Test cases for ecc_unify module."""

    def test_parse_region_simple(self):
        """Test parsing simple region strings."""
        region = "chr1:1000-2000"
        chr_name, start, end = parse_region(region)

        assert chr_name == "chr1"
        assert start == 1000
        assert end == 2000

    def test_parse_region_chimeric(self):
        """Test parsing chimeric region (takes first segment)."""
        region = "chr1:1000-2000;chr2:3000-4000"
        chr_name, start, end = parse_region(region)

        assert chr_name == "chr1"
        assert start == 1000
        assert end == 2000

    def test_parse_chimeric_regions(self):
        """Test parsing multiple chimeric regions."""
        regions_str = "chr1:1000-2000;chr2:3000-4000"
        segments = parse_chimeric_regions(regions_str)

        assert len(segments) == 2
        assert segments[0] == ("chr1", 1000, 2000)
        assert segments[1] == ("chr2", 3000, 4000)

    def test_merge_eccdna_tables_confirmed_only(self, tmp_path):
        """Test merging with confirmed table only."""
        # Create test confirmed data with correct column names
        confirmed_data = pd.DataFrame({
            'eccDNA_id': ['UECC_001', 'MECC_001'],
            'Regions': ['chr1:1000-1500', 'chr2:2000-2500'],
            'Strand': ['+', '+'],
            'Length': [500, 500],
            'eccDNA_type': ['UeccDNA', 'MeccDNA'],
            'State': ['Confirmed', 'Confirmed'],
            'Seg_total': [1, 1],
            'Hit_count': [1, 1]
        })

        confirmed_file = tmp_path / "confirmed.csv"
        confirmed_data.to_csv(confirmed_file, index=False)

        overlap_report = tmp_path / "overlap_report.txt"
        overlap_stats = tmp_path / "overlap_stats.json"

        # Run merge
        merged_df, report_text, stats_dict = merge_eccdna_tables(
            confirmed_file=confirmed_file,
            overlap_report_file=overlap_report,
            overlap_stats_json=overlap_stats
        )

        # Verify results
        assert len(merged_df) == 2
        assert 'eccDNA_id' in merged_df.columns
        assert overlap_report.exists()
        assert overlap_stats.exists()

    def test_merge_eccdna_tables_with_inferred(self, tmp_path):
        """Test merging confirmed and inferred tables."""
        # Create test data with correct column names
        confirmed_data = pd.DataFrame({
            'eccDNA_id': ['UECC_001'],
            'Regions': ['chr1:1000-1500'],
            'Strand': ['+'],
            'Length': [500],
            'eccDNA_type': ['UeccDNA'],
            'State': ['Confirmed'],
            'Seg_total': [1],
            'Hit_count': [1]
        })

        # Inferred simple eccDNA format
        inferred_simple_data = pd.DataFrame({
            'eccDNA_id': ['IUECC_001'],
            'chr': ['chr2'],
            'start0': [2000],
            'end0': [2500],
            'length': [500],
            'strand': ['+']
        })

        confirmed_file = tmp_path / "confirmed.csv"
        simple_file = tmp_path / "simple.csv"
        confirmed_data.to_csv(confirmed_file, index=False)
        inferred_simple_data.to_csv(simple_file, index=False)

        # Run merge
        merged_df, report_text, stats_dict = merge_eccdna_tables(
            confirmed_file=confirmed_file,
            inferred_simple=simple_file
        )

        # Should have both confirmed and inferred
        assert len(merged_df) >= 1  # At least confirmed ones
        assert report_text is not None
        assert isinstance(stats_dict, dict)

    def test_merge_eccdna_tables_empty_input(self, tmp_path):
        """Test handling empty input files."""
        # Create empty confirmed file with expected columns
        empty_data = pd.DataFrame(columns=[
            'eccDNA_id', 'Regions', 'Strand', 'Length',
            'eccDNA_type', 'State', 'Seg_total', 'Hit_count'
        ])
        confirmed_file = tmp_path / "empty_confirmed.csv"
        empty_data.to_csv(confirmed_file, index=False)

        # Should handle gracefully
        merged_df, report_text, stats_dict = merge_eccdna_tables(
            confirmed_file=confirmed_file
        )

        assert len(merged_df) == 0
        assert isinstance(report_text, str)
        assert isinstance(stats_dict, dict)

    def test_merge_eccdna_tables_dataframe_input(self):
        """Test using DataFrame as input instead of file path."""
        confirmed_df = pd.DataFrame({
            'eccDNA_id': ['UECC_001'],
            'Regions': ['chr1:1000-1500'],
            'Strand': ['+'],
            'Length': [500],
            'eccDNA_type': ['UeccDNA'],
            'State': ['Confirmed'],
            'Seg_total': [1],
            'Hit_count': [1]
        })

        # Inferred simple eccDNA format
        inferred_df = pd.DataFrame({
            'eccDNA_id': ['IUECC_001'],
            'chr': ['chr2'],
            'start0': [2000],
            'end0': [2500],
            'length': [500],
            'strand': ['+']
        })

        # Run merge with DataFrames
        merged_df, report_text, stats_dict = merge_eccdna_tables(
            confirmed_file=confirmed_df,
            inferred_simple=inferred_df
        )

        assert len(merged_df) >= 1
        assert 'eccDNA_id' in merged_df.columns

    # ========== Chimeric eccDNA Overlap Detection Tests ========== #

    def test_chimeric_exact_match_both_methods(self, tmp_path):
        """Test that both methods detect exact matches."""
        inferred = pd.DataFrame({
            'eccDNA_id': ['IC1', 'IC1'],
            'chr': ['chr1', 'chr2'],
            'start0': [1000, 3000],
            'end0': [2000, 4000],
            'seg_index': [0, 1]
        })

        confirmed = pd.DataFrame({
            'eccDNA_id': ['CC1'],
            'Regions': ['chr1:1000-2000;chr2:3000-4000'],
            'eccDNA_type': ['CeccDNA']
        })

        redundant_exact = find_redundant_chimeric(
            inferred, confirmed, method='exact'
        )
        assert 'IC1' in redundant_exact

        redundant_overlap = find_redundant_chimeric(
            inferred, confirmed, method='overlap', thr=0.99, tol=10
        )
        assert 'IC1' in redundant_overlap

    def test_chimeric_small_difference_exact_fails(self, tmp_path):
        """Test that exact method fails with small coordinate differences."""
        inferred = pd.DataFrame({
            'eccDNA_id': ['IC1', 'IC1'],
            'chr': ['chr1', 'chr2'],
            'start0': [1000, 3000],
            'end0': [2000, 4000],
            'seg_index': [0, 1]
        })

        confirmed = pd.DataFrame({
            'eccDNA_id': ['CC1'],
            'Regions': ['chr1:1002-2001;chr2:3001-4002'],
            'eccDNA_type': ['CeccDNA']
        })

        redundant_exact = find_redundant_chimeric(
            inferred, confirmed, method='exact'
        )
        assert 'IC1' not in redundant_exact

    def test_chimeric_small_difference_overlap_succeeds(self, tmp_path):
        """Test that overlap method succeeds with small coordinate differences."""
        inferred = pd.DataFrame({
            'eccDNA_id': ['IC1', 'IC1'],
            'chr': ['chr1', 'chr2'],
            'start0': [1000, 3000],
            'end0': [2000, 4000],
            'seg_index': [0, 1]
        })

        confirmed = pd.DataFrame({
            'eccDNA_id': ['CC1'],
            'Regions': ['chr1:1002-2001;chr2:3001-4002'],
            'eccDNA_type': ['CeccDNA']
        })

        redundant_overlap = find_redundant_chimeric(
            inferred, confirmed, method='overlap', thr=0.99, tol=10
        )
        assert 'IC1' in redundant_overlap

    def test_chimeric_boundary_at_tolerance_limit(self, tmp_path):
        """Test boundary case at exactly tolerance limit (10bp)."""
        inferred = pd.DataFrame({
            'eccDNA_id': ['IC1', 'IC1'],
            'chr': ['chr1', 'chr2'],
            'start0': [1000, 3000],
            'end0': [2000, 4000],
            'seg_index': [0, 1]
        })

        confirmed = pd.DataFrame({
            'eccDNA_id': ['CC1'],
            'Regions': ['chr1:1010-2010;chr2:3010-4010'],
            'eccDNA_type': ['CeccDNA']
        })

        redundant_overlap = find_redundant_chimeric(
            inferred, confirmed, method='overlap', thr=0.99, tol=10
        )
        assert 'IC1' in redundant_overlap

    def test_chimeric_large_difference(self, tmp_path):
        """Test that large differences are rejected by both methods."""
        inferred = pd.DataFrame({
            'eccDNA_id': ['IC1', 'IC1'],
            'chr': ['chr1', 'chr2'],
            'start0': [1000, 3000],
            'end0': [2000, 4000],
            'seg_index': [0, 1]
        })

        confirmed = pd.DataFrame({
            'eccDNA_id': ['CC1'],
            'Regions': ['chr1:1050-2050;chr2:3050-4050'],
            'eccDNA_type': ['CeccDNA']
        })

        redundant_exact = find_redundant_chimeric(
            inferred, confirmed, method='exact'
        )
        redundant_overlap = find_redundant_chimeric(
            inferred, confirmed, method='overlap', thr=0.99, tol=10
        )

        assert 'IC1' not in redundant_exact
        assert 'IC1' not in redundant_overlap

    def test_chimeric_different_segment_count_subset_detected(self, tmp_path):
        """Test that inferred CeccDNA with fewer segments is detected as subset."""
        inferred = pd.DataFrame({
            'eccDNA_id': ['IC1', 'IC1'],
            'chr': ['chr1', 'chr2'],
            'start0': [1000, 3000],
            'end0': [2000, 4000],
            'seg_index': [0, 1]
        })

        confirmed = pd.DataFrame({
            'eccDNA_id': ['CC1'],
            'Regions': ['chr1:1000-2000;chr2:3000-4000;chr3:5000-6000'],
            'eccDNA_type': ['CeccDNA']
        })

        # Exact method should not match (different segment count)
        redundant_exact = find_redundant_chimeric(
            inferred, confirmed, method='exact'
        )
        assert 'IC1' not in redundant_exact

        # Overlap method with subset detection should detect it as redundant
        # because IC1's 2 segments are a subset of CC1's 3 segments
        redundant_overlap = find_redundant_chimeric(
            inferred, confirmed, method='overlap', thr=0.99, tol=10
        )
        assert 'IC1' in redundant_overlap

    def test_chimeric_different_segment_count_not_subset(self, tmp_path):
        """Test rejection when segments don't match even with different counts."""
        inferred = pd.DataFrame({
            'eccDNA_id': ['IC1', 'IC1'],
            'chr': ['chr1', 'chr4'],  # chr4 doesn't exist in confirmed
            'start0': [1000, 7000],
            'end0': [2000, 8000],
            'seg_index': [0, 1]
        })

        confirmed = pd.DataFrame({
            'eccDNA_id': ['CC1'],
            'Regions': ['chr1:1000-2000;chr2:3000-4000;chr3:5000-6000'],
            'eccDNA_type': ['CeccDNA']
        })

        redundant_overlap = find_redundant_chimeric(
            inferred, confirmed, method='overlap', thr=0.99, tol=10
        )
        # IC1 has chr4 segment which doesn't match any confirmed segment
        assert 'IC1' not in redundant_overlap

    def test_chimeric_different_chromosome(self, tmp_path):
        """Test rejection when chromosomes differ."""
        inferred = pd.DataFrame({
            'eccDNA_id': ['IC1', 'IC1'],
            'chr': ['chr1', 'chr2'],
            'start0': [1000, 3000],
            'end0': [2000, 4000],
            'seg_index': [0, 1]
        })

        confirmed = pd.DataFrame({
            'eccDNA_id': ['CC1'],
            'Regions': ['chr1:1000-2000;chr3:3000-4000'],
            'eccDNA_type': ['CeccDNA']
        })

        redundant_overlap = find_redundant_chimeric(
            inferred, confirmed, method='overlap', thr=0.99, tol=10
        )

        assert 'IC1' not in redundant_overlap

    def test_chimeric_different_order(self, tmp_path):
        """Test rejection when segment order differs."""
        inferred = pd.DataFrame({
            'eccDNA_id': ['IC1', 'IC1', 'IC1'],
            'chr': ['chr1', 'chr2', 'chr3'],
            'start0': [1000, 3000, 5000],
            'end0': [2000, 4000, 6000],
            'seg_index': [0, 1, 2]
        })

        confirmed = pd.DataFrame({
            'eccDNA_id': ['CC1'],
            'Regions': ['chr1:1000-2000;chr3:5000-6000;chr2:3000-4000'],
            'eccDNA_type': ['CeccDNA']
        })

        redundant_overlap = find_redundant_chimeric(
            inferred, confirmed, method='overlap', thr=0.99, tol=10
        )

        assert 'IC1' not in redundant_overlap

    def test_chimeric_rotated_order_matches(self, tmp_path):
        """Rotated segment order should be treated as redundant."""
        inferred = pd.DataFrame({
            'eccDNA_id': ['IC1', 'IC1', 'IC1'],
            'chr': ['chr1', 'chr2', 'chr3'],
            'start0': [1000, 3000, 5000],
            'end0': [2000, 4000, 6000],
            'seg_index': [0, 1, 2]
        })

        confirmed = pd.DataFrame({
            'eccDNA_id': ['CC1'],
            'Regions': ['chr2:3000-4000;chr3:5000-6000;chr1:1000-2000'],
            'eccDNA_type': ['CeccDNA']
        })

        redundant_overlap = find_redundant_chimeric(
            inferred, confirmed, method='overlap', thr=0.99, tol=10
        )

        assert 'IC1' in redundant_overlap

    def test_chimeric_empty_dataframes(self, tmp_path):
        """Test handling of empty DataFrames."""
        inferred = pd.DataFrame()
        confirmed = pd.DataFrame()

        redundant_overlap = find_redundant_chimeric(
            inferred, confirmed, method='overlap'
        )

        assert len(redundant_overlap) == 0

    def test_chimeric_no_cecc_in_confirmed(self, tmp_path):
        """Test when confirmed has no CeccDNA."""
        inferred = pd.DataFrame({
            'eccDNA_id': ['IC1', 'IC1'],
            'chr': ['chr1', 'chr2'],
            'start0': [1000, 3000],
            'end0': [2000, 4000],
            'seg_index': [0, 1]
        })

        confirmed = pd.DataFrame({
            'eccDNA_id': ['UC1'],
            'Regions': ['chr1:1000-2000'],
            'eccDNA_type': ['UeccDNA']
        })

        redundant_overlap = find_redundant_chimeric(
            inferred, confirmed, method='overlap'
        )

        assert len(redundant_overlap) == 0

    def test_chimeric_multiple_inferred_one_confirmed(self, tmp_path):
        """Test multiple inferred matching one confirmed."""
        inferred = pd.DataFrame({
            'eccDNA_id': ['IC1', 'IC1', 'IC2', 'IC2', 'IC3', 'IC3'],
            'chr': ['chr1', 'chr2', 'chr1', 'chr2', 'chr1', 'chr2'],
            'start0': [1000, 3000, 1002, 3001, 5000, 7000],
            'end0': [2000, 4000, 2001, 4002, 6000, 8000],
            'seg_index': [0, 1, 0, 1, 0, 1]
        })

        confirmed = pd.DataFrame({
            'eccDNA_id': ['CC1'],
            'Regions': ['chr1:1000-2000;chr2:3000-4000'],
            'eccDNA_type': ['CeccDNA']
        })

        redundant_overlap = find_redundant_chimeric(
            inferred, confirmed, method='overlap', thr=0.99, tol=10
        )

        assert 'IC1' in redundant_overlap
        assert 'IC2' in redundant_overlap
        assert 'IC3' not in redundant_overlap

    def test_segments_match_helper_function(self):
        """Test the _segments_match helper function."""
        segs_a = [('chr1', 1000, 2000), ('chr2', 3000, 4000)]
        segs_b = [('chr1', 1002, 2001), ('chr2', 3001, 4002)]

        result = _segments_match(segs_a, segs_b, thr=0.99, tol=10)
        assert result is True

        result_strict = _segments_match(segs_a, segs_b, thr=1.0, tol=0)
        assert result_strict is False

    def test_invalid_method_parameter(self):
        """Test that invalid method parameter raises ValueError."""
        inferred = pd.DataFrame({
            'eccDNA_id': ['IC1', 'IC1'],
            'chr': ['chr1', 'chr2'],
            'start0': [1000, 3000],
            'end0': [2000, 4000],
            'seg_index': [0, 1]
        })

        confirmed = pd.DataFrame({
            'eccDNA_id': ['CC1'],
            'Regions': ['chr1:1000-2000;chr2:3000-4000'],
            'eccDNA_type': ['CeccDNA']
        })

        with pytest.raises(ValueError, match="Unknown method"):
            find_redundant_chimeric(
                inferred, confirmed, method='invalid_method'
            )


@pytest.fixture
def sample_confirmed_data():
    """Sample confirmed eccDNA data."""
    return pd.DataFrame({
        'eccDNA_id': ['UECC_001', 'MECC_001', 'CECC_001'],
        'Regions': ['chr1:1000-1500', 'chr2:2000-2500', 'chr3:3000-3500'],
        'Strand': ['+', '+', '+'],
        'Length': [500, 500, 500],
        'eccDNA_type': ['UeccDNA', 'MeccDNA', 'CeccDNA'],
        'State': ['Confirmed', 'Confirmed', 'Confirmed'],
        'Seg_total': [1, 1, 1],
        'Hit_count': [1, 1, 1]
    })


@pytest.fixture
def sample_inferred_data():
    """Sample inferred simple eccDNA data."""
    return pd.DataFrame({
        'eccDNA_id': ['IUECC_001', 'ICECC_001'],
        'chr': ['chr4', 'chr5'],
        'start0': [4000, 5000],
        'end0': [4500, 5500],
        'length': [500, 500],
        'strand': ['+', '+']
    })


# ========== reciprocal_overlap_ok Tests ========== #


class TestReciprocalOverlapOk:
    """Tests for reciprocal_overlap_ok function."""

    def test_identical_regions(self):
        """Identical regions should always match."""
        assert reciprocal_overlap_ok(100, 200, 100, 200) is True

    def test_no_overlap(self):
        """Non-overlapping regions should not match."""
        assert reciprocal_overlap_ok(100, 200, 300, 400) is False

    def test_partial_overlap_below_threshold(self):
        """Partial overlap below threshold should not match."""
        # 50% overlap: [100-200] vs [150-250], overlap=50, len_a=100, len_b=100
        assert reciprocal_overlap_ok(100, 200, 150, 250, thr=0.99, tol=0) is False

    def test_high_overlap_above_threshold(self):
        """High overlap above threshold should match."""
        # [100-1100] vs [101-1101], overlap=999, frac_a=999/1000=0.999
        assert reciprocal_overlap_ok(100, 1100, 101, 1101, thr=0.99, tol=0) is True

    def test_tolerance_based_match(self):
        """Regions within tolerance should match regardless of overlap fraction."""
        # Boundaries differ by 5bp each, within tol=10
        assert reciprocal_overlap_ok(100, 200, 105, 205, thr=0.99, tol=10) is True

    def test_tolerance_exceeded(self):
        """Regions outside tolerance with low overlap should not match."""
        # Boundaries differ by 15bp, tol=10
        assert reciprocal_overlap_ok(100, 200, 115, 215, thr=0.99, tol=10) is False

    def test_zero_tolerance_strict(self):
        """Zero tolerance forces strict reciprocal overlap check."""
        # [100,10000] vs [101,10001]: ov=9899, la=9900, lb=9900
        # frac = 9899/9900 ≈ 0.99990 > 0.999
        assert reciprocal_overlap_ok(100, 10000, 101, 10001, thr=0.999, tol=0) is True
        # [100,200] vs [105,205]: ov=95, la=100, lb=100
        # frac = 0.95 < 0.99
        assert reciprocal_overlap_ok(100, 200, 105, 205, thr=0.99, tol=0) is False

    def test_adjacent_regions_no_overlap(self):
        """Adjacent (touching) regions have no overlap."""
        assert reciprocal_overlap_ok(100, 200, 200, 300) is False

    def test_contained_region(self):
        """A small region contained in a large one fails reciprocal check."""
        # [100-1000] vs [400-500], overlap=100, frac_a=100/900=0.11, frac_b=100/100=1.0
        assert reciprocal_overlap_ok(100, 1000, 400, 500, thr=0.99, tol=0) is False


# ========== parse_chimeric_regions_with_strand Tests ========== #


class TestParseChimericRegionsWithStrand:
    """Tests for parse_chimeric_regions_with_strand function."""

    def test_basic_parse(self):
        """Parse regions with matching strands."""
        result = parse_chimeric_regions_with_strand(
            "chr1:100-200;chr2:300-400", "+;-"
        )
        assert len(result) == 2
        assert result[0] == ("chr1", 100, 200, "+")
        assert result[1] == ("chr2", 300, 400, "-")

    def test_missing_strands_default_to_plus(self):
        """Missing strand info defaults to +."""
        result = parse_chimeric_regions_with_strand(
            "chr1:100-200;chr2:300-400", "+"
        )
        assert result[0][3] == "+"
        assert result[1][3] == "+"  # Defaulted

    def test_single_segment(self):
        """Single segment with single strand."""
        result = parse_chimeric_regions_with_strand("chr1:100-200", "-")
        assert len(result) == 1
        assert result[0] == ("chr1", 100, 200, "-")

    def test_empty_strand_defaults(self):
        """Empty strand string defaults all to +."""
        result = parse_chimeric_regions_with_strand(
            "chr1:100-200;chr2:300-400", ""
        )
        assert all(seg[3] == "+" for seg in result)


# ========== build_chr_index Tests ========== #


class TestBuildChrIndex:
    """Tests for build_chr_index function."""

    def test_basic_index(self):
        """Build index from simple regions."""
        df = pd.DataFrame({
            "Regions": ["chr1:100-200", "chr1:300-400", "chr2:500-600"],
            "eccDNA_type": ["UeccDNA", "UeccDNA", "UeccDNA"],
        })
        idx = build_chr_index(df)
        assert "chr1" in idx
        assert "chr2" in idx
        assert len(idx["chr1"]) == 2
        assert len(idx["chr2"]) == 1

    def test_type_filter(self):
        """Filter by eccDNA_type."""
        df = pd.DataFrame({
            "Regions": ["chr1:100-200", "chr1:300-400"],
            "eccDNA_type": ["UeccDNA", "MeccDNA"],
        })
        idx = build_chr_index(df, type_filter="UeccDNA")
        assert len(idx.get("chr1", [])) == 1

    def test_empty_dataframe(self):
        """Empty DataFrame returns empty index."""
        df = pd.DataFrame(columns=["Regions", "eccDNA_type"])
        idx = build_chr_index(df)
        assert idx == {}

    def test_chimeric_regions_takes_first(self):
        """Chimeric region takes only first segment."""
        df = pd.DataFrame({
            "Regions": ["chr1:100-200;chr2:300-400"],
            "eccDNA_type": ["CeccDNA"],
        })
        idx = build_chr_index(df)
        assert "chr1" in idx
        assert "chr2" not in idx

    def test_sorted_by_start(self):
        """Index entries sorted by start position."""
        df = pd.DataFrame({
            "Regions": ["chr1:500-600", "chr1:100-200", "chr1:300-400"],
        })
        idx = build_chr_index(df)
        starts = [entry[0] for entry in idx["chr1"]]
        assert starts == sorted(starts)


# ========== build_cecc_segment_index Tests ========== #


class TestBuildCeccSegmentIndex:
    """Tests for build_cecc_segment_index function."""

    def test_indexes_all_segments(self):
        """All segments from CeccDNA should be indexed."""
        df = pd.DataFrame({
            "Regions": ["chr1:100-200;chr2:300-400"],
            "eccDNA_type": ["CeccDNA"],
        })
        idx = build_cecc_segment_index(df)
        assert "chr1" in idx
        assert "chr2" in idx

    def test_filters_non_cecc(self):
        """Non-CeccDNA entries are excluded."""
        df = pd.DataFrame({
            "Regions": ["chr1:100-200", "chr2:300-400"],
            "eccDNA_type": ["UeccDNA", "CeccDNA"],
        })
        idx = build_cecc_segment_index(df)
        # UeccDNA chr1:100-200 should not be included
        assert "chr1" not in idx
        assert "chr2" in idx

    def test_empty_dataframe(self):
        """Empty DataFrame returns empty index."""
        df = pd.DataFrame(columns=["Regions", "eccDNA_type"])
        idx = build_cecc_segment_index(df)
        assert idx == {}


# ========== find_redundant_simple Tests ========== #


class TestFindRedundantSimple:
    """Tests for find_redundant_simple function."""

    def test_overlapping_with_confirmed_uecc(self):
        """Inferred simple overlapping with confirmed UeccDNA."""
        inferred = pd.DataFrame({
            "eccDNA_id": ["I1", "I2"],
            "chr": ["chr1", "chr3"],
            "start0": [100, 500],
            "end0": [200, 600],
        })
        confirmed = pd.DataFrame({
            "Regions": ["chr1:100-200"],
            "eccDNA_type": ["UeccDNA"],
        })
        redundant = find_redundant_simple(inferred, confirmed, thr=0.99, tol=20)
        assert "I1" in redundant
        assert "I2" not in redundant

    def test_overlapping_with_confirmed_cecc_segment(self):
        """Inferred simple overlapping with a CeccDNA segment."""
        inferred = pd.DataFrame({
            "eccDNA_id": ["I1"],
            "chr": ["chr2"],
            "start0": [300],
            "end0": [400],
        })
        confirmed = pd.DataFrame({
            "Regions": ["chr1:100-200;chr2:300-400"],
            "eccDNA_type": ["CeccDNA"],
        })
        redundant = find_redundant_simple(inferred, confirmed, thr=0.99, tol=20)
        assert "I1" in redundant

    def test_no_overlap(self):
        """Non-overlapping entries should not be flagged."""
        inferred = pd.DataFrame({
            "eccDNA_id": ["I1"],
            "chr": ["chr5"],
            "start0": [1000],
            "end0": [2000],
        })
        confirmed = pd.DataFrame({
            "Regions": ["chr1:100-200"],
            "eccDNA_type": ["UeccDNA"],
        })
        redundant = find_redundant_simple(inferred, confirmed)
        assert len(redundant) == 0

    def test_uses_regions_column_fallback(self):
        """Falls back to Regions column if chr/start0/end0 not available."""
        inferred = pd.DataFrame({
            "eccDNA_id": ["I1"],
            "regions": ["chr1:100-200"],
        })
        confirmed = pd.DataFrame({
            "Regions": ["chr1:100-200"],
            "eccDNA_type": ["UeccDNA"],
        })
        redundant = find_redundant_simple(inferred, confirmed, thr=0.99, tol=20)
        assert "I1" in redundant

    def test_empty_inferred(self):
        """Empty inferred returns empty set."""
        inferred = pd.DataFrame(columns=["eccDNA_id", "chr", "start0", "end0"])
        confirmed = pd.DataFrame({
            "Regions": ["chr1:100-200"],
            "eccDNA_type": ["UeccDNA"],
        })
        redundant = find_redundant_simple(inferred, confirmed)
        assert len(redundant) == 0


# ========== _segments_match_with_strand Tests ========== #


class TestSegmentsMatchWithStrand:
    """Tests for _segments_match_with_strand function."""

    def test_same_strand_direct_match(self):
        """Same-strand segments match directly."""
        a = [("chr1", 100, 200, "+"), ("chr2", 300, 400, "+")]
        b = [("chr1", 100, 200, "+"), ("chr2", 300, 400, "+")]
        assert _segments_match_with_strand(a, b, thr=0.99, tol=20) is True

    def test_same_strand_cyclic_rotation(self):
        """Same-strand match with cyclic rotation."""
        a = [("chr1", 100, 200, "+"), ("chr2", 300, 400, "+")]
        b = [("chr2", 300, 400, "+"), ("chr1", 100, 200, "+")]
        assert _segments_match_with_strand(a, b, thr=0.99, tol=20) is True

    def test_reverse_complement_match(self):
        """Reverse-complement match: reversed order + opposite strands."""
        a = [("chr1", 100, 200, "+"), ("chr2", 300, 400, "+")]
        b = [("chr2", 300, 400, "-"), ("chr1", 100, 200, "-")]
        assert _segments_match_with_strand(a, b, thr=0.99, tol=20) is True

    def test_mixed_strands_no_match(self):
        """Mixed strand orientations should not match."""
        a = [("chr1", 100, 200, "+"), ("chr2", 300, 400, "+")]
        b = [("chr1", 100, 200, "+"), ("chr2", 300, 400, "-")]  # One mismatch
        assert _segments_match_with_strand(a, b, thr=0.99, tol=20) is False

    def test_different_lengths_no_match(self):
        """Different number of segments should not match."""
        a = [("chr1", 100, 200, "+")]
        b = [("chr1", 100, 200, "+"), ("chr2", 300, 400, "+")]
        assert _segments_match_with_strand(a, b, thr=0.99, tol=20) is False

    def test_empty_segments(self):
        """Empty segment lists should not match."""
        assert _segments_match_with_strand([], [], thr=0.99, tol=20) is False


# ========== _is_subset_of Tests ========== #


class TestIsSubsetOf:
    """Tests for _is_subset_of function."""

    def test_true_subset(self):
        """2 segments matching within 3 confirmed segments."""
        inferred = [("chr1", 100, 200), ("chr2", 300, 400)]
        confirmed = [("chr1", 100, 200), ("chr2", 300, 400), ("chr3", 500, 600)]
        assert _is_subset_of(inferred, confirmed, thr=0.99, tol=20) is True

    def test_same_length_returns_false(self):
        """Same number of segments defers to _segments_match."""
        inferred = [("chr1", 100, 200)]
        confirmed = [("chr1", 100, 200)]
        assert _is_subset_of(inferred, confirmed, thr=0.99, tol=20) is False

    def test_more_inferred_than_confirmed(self):
        """More inferred segments than confirmed returns False."""
        inferred = [("chr1", 100, 200), ("chr2", 300, 400), ("chr3", 500, 600)]
        confirmed = [("chr1", 100, 200)]
        assert _is_subset_of(inferred, confirmed, thr=0.99, tol=20) is False

    def test_no_match_in_confirmed(self):
        """Inferred segments not found in confirmed."""
        inferred = [("chr5", 100, 200)]
        confirmed = [("chr1", 100, 200), ("chr2", 300, 400)]
        assert _is_subset_of(inferred, confirmed, thr=0.99, tol=20) is False

    def test_empty_inferred(self):
        """Empty inferred returns False."""
        assert _is_subset_of([], [("chr1", 100, 200)], thr=0.99, tol=20) is False


# ========== _is_subset_of_with_strand Tests ========== #


class TestIsSubsetOfWithStrand:
    """Tests for _is_subset_of_with_strand function."""

    def test_same_strand_subset(self):
        """Subset match with same strands."""
        inferred = [("chr1", 100, 200, "+")]
        confirmed = [("chr1", 100, 200, "+"), ("chr2", 300, 400, "+")]
        assert _is_subset_of_with_strand(inferred, confirmed, thr=0.99, tol=20) is True

    def test_opposite_strand_subset(self):
        """Subset match with opposite strands (reverse complement)."""
        inferred = [("chr1", 100, 200, "-")]
        confirmed = [("chr1", 100, 200, "+"), ("chr2", 300, 400, "+")]
        assert _is_subset_of_with_strand(inferred, confirmed, thr=0.99, tol=20) is True

    def test_same_length_returns_false(self):
        """Same number of segments defers to _segments_match_with_strand."""
        inferred = [("chr1", 100, 200, "+")]
        confirmed = [("chr1", 100, 200, "+")]
        assert _is_subset_of_with_strand(inferred, confirmed, thr=0.99, tol=20) is False

    def test_no_strand_match(self):
        """No matching segments returns False."""
        inferred = [("chr5", 100, 200, "+")]
        confirmed = [("chr1", 100, 200, "+"), ("chr2", 300, 400, "-")]
        assert _is_subset_of_with_strand(inferred, confirmed, thr=0.99, tol=20) is False


# ========== prepare_inferred_simple Tests ========== #


class TestPrepareInferredSimple:
    """Tests for prepare_inferred_simple function."""

    def test_basic_conversion(self):
        """Converts inferred format to standard format."""
        df = pd.DataFrame({
            "eccDNA_id": ["I1", "I2"],
            "chr": ["chr1", "chr2"],
            "start0": [100, 300],
            "end0": [200, 400],
            "length": [100, 100],
            "strand": ["+", "-"],
        })
        result = prepare_inferred_simple(df, redundant_ids=set())
        assert len(result) == 2
        assert "Regions" in result.columns
        assert result["eccDNA_type"].iloc[0] == "UeccDNA"
        assert result["State"].iloc[0] == "Inferred"
        assert result["Regions"].iloc[0] == "chr1:100-200"

    def test_filters_redundant(self):
        """Redundant IDs are excluded."""
        df = pd.DataFrame({
            "eccDNA_id": ["I1", "I2"],
            "chr": ["chr1", "chr2"],
            "start0": [100, 300],
            "end0": [200, 400],
            "length": [100, 100],
            "strand": ["+", "+"],
        })
        result = prepare_inferred_simple(df, redundant_ids={"I1"})
        assert len(result) == 1
        assert result["eccDNA_id"].iloc[0] == "I2"

    def test_all_redundant_returns_empty(self):
        """All entries redundant returns empty DataFrame."""
        df = pd.DataFrame({
            "eccDNA_id": ["I1"],
            "chr": ["chr1"],
            "start0": [100],
            "end0": [200],
            "length": [100],
            "strand": ["+"],
        })
        result = prepare_inferred_simple(df, redundant_ids={"I1"})
        assert result.empty

    def test_preserves_metrics(self):
        """num_split_reads, prob_present, hifi_abundance are preserved."""
        df = pd.DataFrame({
            "eccDNA_id": ["I1"],
            "chr": ["chr1"],
            "start0": [100],
            "end0": [200],
            "length": [100],
            "strand": ["+"],
            "num_split_reads": [5],
            "prob_present": [0.95],
            "hifi_abundance": [2.5],
        })
        result = prepare_inferred_simple(df, redundant_ids=set())
        assert result["reads_count"].iloc[0] == 5
        assert result["confidence_score"].iloc[0] == 0.95
        assert result["copy_number"].iloc[0] == 2.5


# ========== prepare_inferred_chimeric Tests ========== #


class TestPrepareInferredChimeric:
    """Tests for prepare_inferred_chimeric function."""

    def test_basic_conversion(self):
        """Converts chimeric segment format to unified format."""
        df = pd.DataFrame({
            "eccDNA_id": ["IC1", "IC1"],
            "chr": ["chr1", "chr2"],
            "start0": [100, 300],
            "end0": [200, 400],
            "length": [200, 200],
            "seg_total": [2, 2],
            "seg_index": [0, 1],
            "strand": ["+", "-"],
        })
        result = prepare_inferred_chimeric(df, redundant_ids=set())
        assert len(result) == 1
        assert result["Regions"].iloc[0] == "chr1:100-200;chr2:300-400"
        assert result["Strand"].iloc[0] == "+;-"
        assert result["eccDNA_type"].iloc[0] == "CeccDNA"
        assert result["State"].iloc[0] == "Inferred"

    def test_filters_redundant(self):
        """Redundant IDs excluded from chimeric conversion."""
        df = pd.DataFrame({
            "eccDNA_id": ["IC1", "IC1", "IC2", "IC2"],
            "chr": ["chr1", "chr2", "chr3", "chr4"],
            "start0": [100, 300, 500, 700],
            "end0": [200, 400, 600, 800],
            "length": [200, 200, 200, 200],
            "seg_total": [2, 2, 2, 2],
            "seg_index": [0, 1, 0, 1],
            "strand": ["+", "+", "+", "+"],
        })
        result = prepare_inferred_chimeric(df, redundant_ids={"IC1"})
        assert len(result) == 1
        assert result["eccDNA_id"].iloc[0] == "IC2"

    def test_all_redundant_returns_empty(self):
        """All entries redundant returns empty DataFrame."""
        df = pd.DataFrame({
            "eccDNA_id": ["IC1", "IC1"],
            "chr": ["chr1", "chr2"],
            "start0": [100, 300],
            "end0": [200, 400],
            "length": [200, 200],
            "seg_total": [2, 2],
            "seg_index": [0, 1],
            "strand": ["+", "+"],
        })
        result = prepare_inferred_chimeric(df, redundant_ids={"IC1"})
        assert result.empty


# ========== renumber_eccdna Tests ========== #


class TestRenumberEccdna:
    """Tests for renumber_eccdna function."""

    def test_confirmed_keep_original_ids(self):
        """Confirmed entries keep original IDs."""
        df = pd.DataFrame({
            "eccDNA_id": ["UeccDNA1", "UeccDNA2"],
            "eccDNA_type": ["UeccDNA", "UeccDNA"],
            "State": ["Confirmed", "Confirmed"],
        })
        result = renumber_eccdna(df)
        assert result["eccDNA_id"].tolist() == ["UeccDNA1", "UeccDNA2"]
        assert "original_id" in result.columns

    def test_inferred_continue_numbering(self):
        """Inferred entries get IDs continuing after confirmed."""
        df = pd.DataFrame({
            "eccDNA_id": ["UeccDNA1", "UeccDNA2", "inf_1"],
            "eccDNA_type": ["UeccDNA", "UeccDNA", "UeccDNA"],
            "State": ["Confirmed", "Confirmed", "Inferred"],
        })
        result = renumber_eccdna(df)
        # Confirmed: UeccDNA1, UeccDNA2; Inferred: UeccDNA3
        assert result.iloc[2]["eccDNA_id"] == "UeccDNA3"

    def test_empty_dataframe(self):
        """Empty DataFrame handled gracefully."""
        df = pd.DataFrame(columns=["eccDNA_id", "eccDNA_type", "State"])
        result = renumber_eccdna(df)
        assert len(result) == 0

    def test_mixed_types_sorted(self):
        """Mixed eccDNA types are sorted by type then state."""
        df = pd.DataFrame({
            "eccDNA_id": ["CeccDNA1", "UeccDNA1", "MeccDNA1"],
            "eccDNA_type": ["CeccDNA", "UeccDNA", "MeccDNA"],
            "State": ["Confirmed", "Confirmed", "Confirmed"],
        })
        result = renumber_eccdna(df)
        assert result["eccDNA_type"].tolist() == ["UeccDNA", "MeccDNA", "CeccDNA"]

    def test_multiple_inferred_per_type(self):
        """Multiple inferred per type get sequential numbering."""
        df = pd.DataFrame({
            "eccDNA_id": ["UeccDNA1", "inf_a", "inf_b"],
            "eccDNA_type": ["UeccDNA", "UeccDNA", "UeccDNA"],
            "State": ["Confirmed", "Inferred", "Inferred"],
        })
        result = renumber_eccdna(df)
        inferred = result[result["State"] == "Inferred"]
        assert inferred["eccDNA_id"].tolist() == ["UeccDNA2", "UeccDNA3"]


# ========== generate_overlap_stats_json Tests ========== #


class TestGenerateOverlapStatsJson:
    """Tests for generate_overlap_stats_json function."""

    def test_basic_stats(self):
        """Generate basic overlap statistics."""
        confirmed = pd.DataFrame({
            "eccDNA_type": ["UeccDNA", "MeccDNA", "CeccDNA"],
        })
        stats = generate_overlap_stats_json(
            confirmed, None, None, set(), set()
        )
        assert stats["confirmed"]["total"] == 3
        assert stats["confirmed"]["UeccDNA"] == 1
        assert stats["summary"]["final_total"] == 3

    def test_with_inferred(self):
        """Stats include inferred counts."""
        confirmed = pd.DataFrame({
            "eccDNA_type": ["UeccDNA"],
        })
        inferred_simple = pd.DataFrame({
            "eccDNA_id": ["I1", "I2", "I3"],
        })
        stats = generate_overlap_stats_json(
            confirmed, inferred_simple, None, {"I1"}, set()
        )
        assert stats["inferred_simple"]["total"] == 3
        assert stats["inferred_simple"]["overlapping"] == 1
        assert stats["inferred_simple"]["non_redundant"] == 2
        assert stats["summary"]["final_total"] == 3  # 1 confirmed + 2 non-redundant

    def test_save_to_file(self, tmp_path):
        """Stats saved to JSON file."""
        confirmed = pd.DataFrame({"eccDNA_type": ["UeccDNA"]})
        output = tmp_path / "stats.json"
        generate_overlap_stats_json(
            confirmed, None, None, set(), set(), output_file=output
        )
        assert output.exists()
        data = json.loads(output.read_text())
        assert "confirmed" in data
        assert "summary" in data


# ========== generate_overlap_report Tests ========== #


class TestGenerateOverlapReport:
    """Tests for generate_overlap_report function."""

    def test_basic_report(self):
        """Generate basic overlap report text."""
        confirmed = pd.DataFrame({
            "eccDNA_type": ["UeccDNA", "MeccDNA"],
        })
        report = generate_overlap_report(confirmed, None, None, set(), set())
        assert "eccDNA Overlap Analysis Report" in report
        assert "UeccDNA: 1 sequences" in report
        assert "MeccDNA: 1 sequences" in report

    def test_report_with_redundant(self):
        """Report includes redundant IDs when present."""
        confirmed = pd.DataFrame({
            "eccDNA_type": ["UeccDNA"],
        })
        inferred_simple = pd.DataFrame({
            "eccDNA_id": ["I1", "I2"],
            "chr": ["chr1", "chr2"],
            "start0": [100, 300],
            "end0": [200, 400],
        })
        report = generate_overlap_report(
            confirmed, inferred_simple, None, {"I1"}, set()
        )
        assert "Redundant IDs" in report
        assert "I1" in report

    def test_save_report_to_file(self, tmp_path):
        """Report saved to file."""
        confirmed = pd.DataFrame({"eccDNA_type": ["UeccDNA"]})
        output = tmp_path / "report.txt"
        generate_overlap_report(
            confirmed, None, None, set(), set(), output_file=output
        )
        assert output.exists()
        content = output.read_text()
        assert "eccDNA Overlap Analysis Report" in content

    def test_no_inferred(self):
        """Report handles no inferred data."""
        confirmed = pd.DataFrame({"eccDNA_type": ["UeccDNA"]})
        report = generate_overlap_report(confirmed, None, None, set(), set())
        assert "No inferred simple eccDNA provided" in report
        assert "No inferred chimeric eccDNA provided" in report
