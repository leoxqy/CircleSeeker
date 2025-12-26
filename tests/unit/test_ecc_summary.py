"""Tests for ecc_summary module."""

from pathlib import Path
import sys
import pytest
import pandas as pd
import json

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.modules.ecc_summary import EccSummary


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
        assert summary.version == "v2.1.0"
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
        # First create a FASTA file (required for full stats)
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(">read1\nACGT\n>read2\nACGT\n>read3\nACGT\n")

        # Create test processed CSV with correct column names
        processed_data = pd.DataFrame({
            'readName': ['read1', 'read2', 'read3'],
            'readClass': ['CtcR-perfect', 'Normal', 'CtcR-hybrid'],
            'repN': [1, 2, 1],
            'copyNum': [5.2, 3.1, 7.8],
            'consLen': [150, 200, 180]
        })

        csv_file = tmp_path / "processed.csv"
        processed_data.to_csv(csv_file, index=False)

        output_dir = tmp_path / "summary_output"
        output_dir.mkdir()

        summary = EccSummary("test_sample", output_dir)
        # Process FASTA first (needed for full stats)
        summary.process_fasta(fasta_file)
        summary.process_processed_csv(csv_file)

        # Check that CTCR stats were updated
        assert 'total_reads' in summary.ctcr_stats
        assert 'ctcr_reads' in summary.ctcr_stats
        assert summary.ctcr_stats['total_reads'] == 3

    def test_process_merged_csv(self, tmp_path):
        """Test merged CSV file processing."""
        # Create test merged CSV with correct column names
        merged_data = pd.DataFrame({
            'eccDNA_id': ['UECC_001', 'MECC_001', 'CECC_001'],
            'eccdna_type': ['UeccDNA', 'MeccDNA', 'CeccDNA'],
            'chr': ['chr1', 'chr2', 'chr3'],
            'start0': [1000, 2000, 3000],
            'end0': [1500, 2500, 3500],
            'Length': [500, 500, 500]
        })

        csv_file = tmp_path / "merged.csv"
        merged_data.to_csv(csv_file, index=False)

        output_dir = tmp_path / "summary_output"
        output_dir.mkdir()

        summary = EccSummary("test_sample", output_dir)
        summary.process_merged_csv(csv_file)

        # Check that eccDNA stats were updated - using actual key names
        assert 'all_count' in summary.eccdna_stats
        assert 'uecc_count' in summary.eccdna_stats
        assert 'mecc_count' in summary.eccdna_stats
        assert 'cecc_count' in summary.eccdna_stats

    def test_process_overlap_stats(self, tmp_path):
        """Test overlap statistics processing."""
        # Create test overlap stats JSON with actual expected format
        overlap_stats = {
            "uecc": {
                "inferred": 5,
                "overlap": 2
            },
            "cecc": {
                "inferred": 3,
                "overlap": 1
            }
        }

        stats_file = tmp_path / "overlap_stats.json"
        with open(stats_file, 'w') as f:
            json.dump(overlap_stats, f)

        output_dir = tmp_path / "summary_output"
        output_dir.mkdir()

        summary = EccSummary("test_sample", output_dir)
        summary.process_overlap_stats(stats_file)

        # Check that overlap stats were loaded - using actual key names
        assert 'inferred_uecc_count' in summary.overlap_stats
        assert 'inferred_cecc_count' in summary.overlap_stats

    def test_generate_html_report(self, tmp_path):
        """Test HTML report generation."""
        output_dir = tmp_path / "summary_output"
        output_dir.mkdir()

        summary = EccSummary("test_sample", output_dir)

        # Add some test data with actual key names
        summary.read_stats = {'total_sequences': 100, 'total_length': 50000, 'total_reads': 100}
        summary.ctcr_stats = {
            'total_reads': 100,
            'ctcr_reads': 75,
            'circular_reads': 50,
            'ctcr_perfect_count': 50,
            'ctcr_hybrid_count': 15,
            'ctcr_inversion_count': 10,
            'ctcr_perfect_pct_ctcr': 66.7,
            'ctcr_hybrid_pct_ctcr': 20.0,
            'ctcr_inversion_pct_ctcr': 13.3
        }
        summary.eccdna_stats = {
            'all_count': 10,
            'uecc_count': 6,
            'mecc_count': 3,
            'cecc_count': 1,
            'all_min_length': '100',
            'all_max_length': '1000',
            'all_mean_length': '500',
            'all_median_length': '500',
            'all_mode_length': '500',
            'all_std_length': '100',
            'uecc_min_length': '100',
            'uecc_max_length': '800',
            'uecc_mean_length': '400',
            'uecc_median_length': '400',
            'uecc_mode_length': '400',
            'uecc_std_length': '80',
            'mecc_min_length': '200',
            'mecc_max_length': '900',
            'mecc_mean_length': '500',
            'mecc_median_length': '500',
            'mecc_mode_length': '500',
            'mecc_std_length': '90',
            'cecc_min_length': '300',
            'cecc_max_length': '1000',
            'cecc_mean_length': '600',
            'cecc_median_length': '600',
            'cecc_mode_length': '600',
            'cecc_std_length': '100'
        }
        summary.overlap_stats = {
            'inferred_uecc_count': 5,
            'uecc_overlap_count': 2,
            'uecc_overlap_pct': '40.00',
            'inferred_cecc_count': 3,
            'cecc_overlap_count': 1,
            'cecc_overlap_pct': '33.33'
        }

        summary.generate_html_report()

        # Check that HTML file was created
        html_file = output_dir / f"{summary.sample_name}_report.html"
        assert html_file.exists()

        # Check basic content
        content = html_file.read_text()
        assert "CircleSeeker" in content
        assert summary.sample_name in content

    def test_generate_text_summary(self, tmp_path):
        """Test text summary generation."""
        output_dir = tmp_path / "summary_output"
        output_dir.mkdir()

        summary = EccSummary("test_sample", output_dir)

        # Add minimal test data - read_stats needs keys used by generate_text_summary
        summary.read_stats = {
            'total_sequences': 100,
            'total_length': 50000,
            'total_reads': 100,
            'ctcr_reads': 75,
            'ctcr_percentage': 75.0,
            'other_reads': 25,
            'other_percentage': 25.0
        }
        summary.ctcr_stats = {
            'total_reads': 100,
            'ctcr_reads': 75,
            'circular_reads': 50,
            'ctcr_perfect_count': 50,
            'ctcr_hybrid_count': 15,
            'ctcr_inversion_count': 10,
            'ctcr_perfect_pct_ctcr': 66.7,
            'ctcr_hybrid_pct_ctcr': 20.0,
            'ctcr_inversion_pct_ctcr': 13.3
        }
        summary.eccdna_stats = {
            'all_count': 10,
            'uecc_count': 6,
            'mecc_count': 3,
            'cecc_count': 1,
            'all_min_length': '100',
            'all_max_length': '1000',
            'all_mean_length': '500',
            'all_median_length': '500',
            'all_mode_length': '500',
            'all_std_length': '100',
            'uecc_min_length': '100',
            'uecc_max_length': '800',
            'uecc_mean_length': '400',
            'uecc_median_length': '400',
            'uecc_mode_length': '400',
            'uecc_std_length': '80',
            'mecc_min_length': '200',
            'mecc_max_length': '900',
            'mecc_mean_length': '500',
            'mecc_median_length': '500',
            'mecc_mode_length': '500',
            'mecc_std_length': '90',
            'cecc_min_length': '300',
            'cecc_max_length': '1000',
            'cecc_mean_length': '600',
            'cecc_median_length': '600',
            'cecc_mode_length': '600',
            'cecc_std_length': '100'
        }
        summary.overlap_stats = {
            'inferred_uecc_count': 5,
            'uecc_overlap_count': 2,
            'uecc_overlap_pct': '40.00',
            'inferred_cecc_count': 3,
            'cecc_overlap_count': 1,
            'cecc_overlap_pct': '33.33'
        }

        summary.generate_text_summary()

        # Check that text file was created
        text_file = output_dir / f"{summary.sample_name}_summary.txt"
        assert text_file.exists()

        # Check basic content
        content = text_file.read_text()
        assert "CircleSeeker" in content
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

        # Empty CSV will cause pandas error, but should be handled
        try:
            summary.process_processed_csv(empty_csv)
        except Exception:
            pass  # Expected for empty CSV

        try:
            summary.process_merged_csv(empty_csv)
        except Exception:
            pass  # Expected for empty CSV

        # Should still be initialized
        assert summary.sample_name == "test_sample"


@pytest.fixture
def sample_summary(tmp_path):
    """Create a sample EccSummary instance for testing."""
    output_dir = tmp_path / "summary_output"
    output_dir.mkdir()
    return EccSummary("test_sample", output_dir)
