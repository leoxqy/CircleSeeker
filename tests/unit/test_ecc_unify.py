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

from circleseeker.modules.ecc_unify import merge_eccdna_tables, parse_region, parse_chimeric_regions


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
        # Create test confirmed data
        confirmed_data = pd.DataFrame({
            'ecc_id': ['UECC_001', 'MECC_001'],
            'chromosome': ['chr1', 'chr2'],
            'start': [1000, 2000],
            'end': [1500, 2500],
            'type': ['UECC', 'MECC']
        })

        confirmed_file = tmp_path / "confirmed.csv"
        confirmed_data.to_csv(confirmed_file, index=False)

        output_file = tmp_path / "unified.csv"
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
        assert 'ecc_id' in merged_df.columns
        assert overlap_report.exists()
        assert overlap_stats.exists()

    def test_merge_eccdna_tables_with_inferred(self, tmp_path):
        """Test merging confirmed and inferred tables."""
        # Create test data
        confirmed_data = pd.DataFrame({
            'ecc_id': ['UECC_001'],
            'chromosome': ['chr1'],
            'start': [1000],
            'end': [1500],
            'type': ['UECC']
        })

        inferred_simple_data = pd.DataFrame({
            'ecc_id': ['IUECC_001'],
            'chromosome': ['chr2'],
            'start': [2000],
            'end': [2500],
            'type': ['IUECC']
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
        # Create empty confirmed file
        empty_data = pd.DataFrame()
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
            'ecc_id': ['UECC_001'],
            'chromosome': ['chr1'],
            'start': [1000],
            'end': [1500],
            'type': ['UECC']
        })

        inferred_df = pd.DataFrame({
            'ecc_id': ['IUECC_001'],
            'chromosome': ['chr2'],
            'start': [2000],
            'end': [2500],
            'type': ['IUECC']
        })

        # Run merge with DataFrames
        merged_df, report_text, stats_dict = merge_eccdna_tables(
            confirmed_file=confirmed_df,
            inferred_simple=inferred_df
        )

        assert len(merged_df) >= 1
        assert 'ecc_id' in merged_df.columns


@pytest.fixture
def sample_confirmed_data():
    """Sample confirmed eccDNA data."""
    return pd.DataFrame({
        'ecc_id': ['UECC_001', 'MECC_001', 'CECC_001'],
        'chromosome': ['chr1', 'chr2', 'chr3'],
        'start': [1000, 2000, 3000],
        'end': [1500, 2500, 3500],
        'type': ['UECC', 'MECC', 'CECC'],
        'length': [500, 500, 500]
    })


@pytest.fixture
def sample_inferred_data():
    """Sample inferred eccDNA data."""
    return pd.DataFrame({
        'ecc_id': ['IUECC_001', 'ICECC_001'],
        'chromosome': ['chr4', 'chr5'],
        'start': [4000, 5000],
        'end': [4500, 5500],
        'type': ['IUECC', 'ICECC'],
        'length': [500, 500]
    })