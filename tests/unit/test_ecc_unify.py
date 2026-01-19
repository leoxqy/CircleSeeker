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
    find_redundant_chimeric,
    _segments_match
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

    def test_chimeric_different_segment_count(self, tmp_path):
        """Test rejection when segment counts differ."""
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

        redundant_exact = find_redundant_chimeric(
            inferred, confirmed, method='exact'
        )
        redundant_overlap = find_redundant_chimeric(
            inferred, confirmed, method='overlap', thr=0.99, tol=10
        )

        assert 'IC1' not in redundant_exact
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
