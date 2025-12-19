"""Tests for ecc_summary module."""

from pathlib import Path
import sys
import pytest
import pandas as pd
import json
import tempfile
from unittest.mock import patch, MagicMock

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.modules.ecc_summary import EccSummary
from circleseeker import __version__


class TestEccSummary:
    """Test cases for ecc_summary module."""

    def test_ecc_summary_initialization(self, tmp_path):
        """Test EccSummary class initialization."""
        sample_name = "test_sample"
        output_dir = tmp_path / "summary_output"
        output_dir.mkdir()

        summary = EccSummary(
            sample_name=sample_name,
            output_dir=output_dir
        )

        assert summary.sample_name == sample_name
        assert summary.output_dir == output_dir
        assert summary.version == __version__
        assert isinstance(summary.read_stats, dict)
        assert isinstance(summary.ctcr_stats, dict)
        assert isinstance(summary.eccdna_stats, dict)

    def test_process_fasta(self, tmp_path):
        """Test FASTA file processing."""
        # Create test FASTA file
        fasta_content = """>seq1
ATCGATCGATCGATCG
>seq2
GCTAGCTAGCTAGCTA
>seq3
TTTTAAAAAGGGGCCCC"""

        fasta_file = tmp_path / "test.fasta"
        with open(fasta_file, 'w') as f:
            f.write(fasta_content)

        output_dir = tmp_path / "summary_output"
        output_dir.mkdir()

        summary = EccSummary("test_sample", output_dir)
        summary.process_fasta(fasta_file)

        # Check that read stats were updated
        assert summary.read_stats.get('total_sequences', 0) == 3
        assert summary.read_stats.get('total_length', 0) > 0

    def test_process_processed_csv(self, tmp_path):
        """Test processed CSV file processing."""
        # Create test processed CSV
        processed_data = pd.DataFrame({
            'readName': ['read1', 'read2', 'read3'],
            'repN': [1, 2, 1],
            'copyNum': [5.2, 3.1, 7.8],
            'consLen': [150, 200, 180],
            'classification': ['circular', 'linear', 'circular']
        })

        csv_file = tmp_path / "processed.csv"
        processed_data.to_csv(csv_file, index=False)

        output_dir = tmp_path / "summary_output"
        output_dir.mkdir()

        summary = EccSummary("test_sample", output_dir)
        summary.process_processed_csv(csv_file)

        # Check that CTCR stats were updated
        assert 'total_reads' in summary.ctcr_stats
        assert 'circular_reads' in summary.ctcr_stats
        assert summary.ctcr_stats['total_reads'] == 3

    def test_process_merged_csv(self, tmp_path):
        """Test merged CSV file processing."""
        # Create test merged CSV
        merged_data = pd.DataFrame({
            'ecc_id': ['UECC_001', 'MECC_001', 'CECC_001'],
            'type': ['UECC', 'MECC', 'CECC'],
            'chromosome': ['chr1', 'chr2', 'chr3'],
            'start': [1000, 2000, 3000],
            'end': [1500, 2500, 3500],
            'length': [500, 500, 500]
        })

        csv_file = tmp_path / "merged.csv"
        merged_data.to_csv(csv_file, index=False)

        output_dir = tmp_path / "summary_output"
        output_dir.mkdir()

        summary = EccSummary("test_sample", output_dir)
        summary.process_merged_csv(csv_file)

        # Check that eccDNA stats were updated
        assert 'total_eccdna' in summary.eccdna_stats
        assert 'uecc_count' in summary.eccdna_stats
        assert 'mecc_count' in summary.eccdna_stats
        assert 'cecc_count' in summary.eccdna_stats
        assert summary.eccdna_stats['total_eccdna'] == 3

    def test_process_overlap_stats(self, tmp_path):
        """Test overlap statistics processing."""
        # Create test overlap stats JSON
        overlap_stats = {
            "overlap_stats": {
                "total_confirmed": 10,
                "total_inferred": 5,
                "overlapping_pairs": 2,
                "overlap_percentage": 20.0
            }
        }

        stats_file = tmp_path / "overlap_stats.json"
        with open(stats_file, 'w') as f:
            json.dump(overlap_stats, f)

        output_dir = tmp_path / "summary_output"
        output_dir.mkdir()

        summary = EccSummary("test_sample", output_dir)
        summary.process_overlap_stats(stats_file)

        # Check that overlap stats were loaded
        assert summary.overlap_stats == overlap_stats["overlap_stats"]

    def test_generate_html_report(self, tmp_path):
        """Test HTML report generation."""
        output_dir = tmp_path / "summary_output"
        output_dir.mkdir()

        summary = EccSummary("test_sample", output_dir)

        # Add some test data
        summary.read_stats = {'total_sequences': 100, 'total_length': 50000}
        summary.ctcr_stats = {'total_reads': 100, 'circular_reads': 75}
        summary.eccdna_stats = {'total_eccdna': 10, 'uecc_count': 6, 'mecc_count': 3, 'cecc_count': 1}

        summary.generate_html_report()

        # Check that HTML file was created
        html_file = output_dir / f"{summary.sample_name}_report.html"
        assert html_file.exists()

        # Check basic content
        content = html_file.read_text()
        assert "CircleSeeker Analysis Report" in content
        assert summary.sample_name in content

    def test_generate_text_summary(self, tmp_path):
        """Test text summary generation."""
        output_dir = tmp_path / "summary_output"
        output_dir.mkdir()

        summary = EccSummary("test_sample", output_dir)

        # Add some test data
        summary.read_stats = {'total_sequences': 100, 'total_length': 50000}
        summary.ctcr_stats = {'total_reads': 100, 'circular_reads': 75}
        summary.eccdna_stats = {'total_eccdna': 10, 'uecc_count': 6, 'mecc_count': 3, 'cecc_count': 1}

        summary.generate_text_summary()

        # Check that text file was created
        text_file = output_dir / f"{summary.sample_name}_summary.txt"
        assert text_file.exists()

        # Check basic content
        content = text_file.read_text()
        assert "CircleSeeker Analysis Summary" in content
        assert summary.sample_name in content

    def test_empty_files_handling(self, tmp_path):
        """Test handling of empty input files."""
        # Create empty files
        empty_fasta = tmp_path / "empty.fasta"
        empty_csv = tmp_path / "empty.csv"
        empty_fasta.touch()
        empty_csv.touch()

        output_dir = tmp_path / "summary_output"
        output_dir.mkdir()

        summary = EccSummary("test_sample", output_dir)

        # Should handle empty files gracefully
        summary.process_fasta(empty_fasta)
        summary.process_processed_csv(empty_csv)
        summary.process_merged_csv(empty_csv)

        # Should still be able to generate reports
        summary.generate_html_report()
        summary.generate_text_summary()

        html_file = output_dir / f"{summary.sample_name}_report.html"
        text_file = output_dir / f"{summary.sample_name}_summary.txt"

        assert html_file.exists()
        assert text_file.exists()


@pytest.fixture
def sample_summary(tmp_path):
    """Create a sample EccSummary instance for testing."""
    output_dir = tmp_path / "summary_output"
    output_dir.mkdir()
    return EccSummary("test_sample", output_dir)
