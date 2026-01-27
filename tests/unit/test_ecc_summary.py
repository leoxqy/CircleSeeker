"""Tests for ecc_summary module."""

from pathlib import Path
from datetime import datetime
import sys
import pytest
import pandas as pd
import json

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.modules.ecc_summary import EccSummary
from circleseeker.__version__ import __version__


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
        assert summary.version == f"v{__version__}"
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


# ---------------------------------------------------------------------------
# New test classes covering previously untested areas
# ---------------------------------------------------------------------------


class TestCollectReadStatistics:
    """Tests for EccSummary.collect_read_statistics()."""

    def test_basic_fasta_and_csv_with_all_subtypes(self, tmp_path):
        """All three CtcR subtypes present in CSV; total reads from FASTA."""
        fasta = tmp_path / "reads.fasta"
        fasta.write_text(">r1\nACGT\n>r2\nTTTT\n>r3\nGGGG\n>r4\nCCCC\n>r5\nAAAA\n")

        df = pd.DataFrame({
            "readName": ["r1", "r2", "r3", "r4", "r5"],
            "readClass": [
                "CtcR-perfect", "CtcR-hybrid", "CtcR-inversion",
                "Normal", "Normal",
            ],
        })
        csv_file = tmp_path / "processed.csv"
        df.to_csv(csv_file, index=False)

        output_dir = tmp_path / "out"
        output_dir.mkdir()
        summary = EccSummary("s1", output_dir)
        result = summary.collect_read_statistics(fasta, csv_file)

        assert result["total_reads"] == 5
        assert result["ctcr_reads"] == 3
        assert result["other_reads"] == 2
        assert summary.ctcr_stats["ctcr_perfect_count"] == 1
        assert summary.ctcr_stats["ctcr_hybrid_count"] == 1
        assert summary.ctcr_stats["ctcr_inversion_count"] == 1
        # Percentages of CtcR should each be ~33.3
        assert summary.ctcr_stats["ctcr_perfect_pct_ctcr"] == pytest.approx(33.3, abs=0.1)

    def test_no_ctcr_reads(self, tmp_path):
        """FASTA has reads but CSV has no CtcR entries -> ctcr_reads=0."""
        fasta = tmp_path / "reads.fasta"
        fasta.write_text(">r1\nACGT\n>r2\nTTTT\n")

        df = pd.DataFrame({
            "readName": ["r1", "r2"],
            "readClass": ["Normal", "Normal"],
        })
        csv_file = tmp_path / "processed.csv"
        df.to_csv(csv_file, index=False)

        output_dir = tmp_path / "out"
        output_dir.mkdir()
        summary = EccSummary("s1", output_dir)
        result = summary.collect_read_statistics(fasta, csv_file)

        assert result["total_reads"] == 2
        assert result["ctcr_reads"] == 0
        assert result["other_reads"] == 2
        assert result["ctcr_percentage"] == 0

    def test_missing_fasta_file(self, tmp_path):
        """Missing FASTA file -> total_reads=0 (graceful)."""
        csv_file = tmp_path / "processed.csv"
        pd.DataFrame({"readName": ["r1"], "readClass": ["CtcR-perfect"]}).to_csv(
            csv_file, index=False
        )

        output_dir = tmp_path / "out"
        output_dir.mkdir()
        summary = EccSummary("s1", output_dir)
        result = summary.collect_read_statistics(tmp_path / "missing.fasta", csv_file)

        assert result["total_reads"] == 0
        assert result["total_length"] == 0

    def test_missing_csv_file(self, tmp_path):
        """Missing CSV file -> ctcr counts remain 0 (graceful)."""
        fasta = tmp_path / "reads.fasta"
        fasta.write_text(">r1\nACGT\n")

        output_dir = tmp_path / "out"
        output_dir.mkdir()
        summary = EccSummary("s1", output_dir)
        result = summary.collect_read_statistics(fasta, tmp_path / "missing.csv")

        assert result["total_reads"] == 1
        assert result["ctcr_reads"] == 0

    def test_empty_fasta(self, tmp_path):
        """Empty FASTA file -> total_reads=0."""
        fasta = tmp_path / "empty.fasta"
        fasta.write_text("")

        csv_file = tmp_path / "processed.csv"
        pd.DataFrame({"readName": [], "readClass": []}).to_csv(csv_file, index=False)

        output_dir = tmp_path / "out"
        output_dir.mkdir()
        summary = EccSummary("s1", output_dir)
        result = summary.collect_read_statistics(fasta, csv_file)

        assert result["total_reads"] == 0
        assert result["total_length"] == 0
        assert result["ctcr_reads"] == 0


class TestCollectEccdnaStatistics:
    """Tests for EccSummary.collect_eccdna_statistics()."""

    def test_all_three_types_present(self, tmp_path):
        """UeccDNA, MeccDNA, CeccDNA all present."""
        df = pd.DataFrame({
            "eccDNA_id": ["U1", "U2", "M1", "C1"],
            "eccDNA_type": ["UeccDNA", "UeccDNA", "MeccDNA", "CeccDNA"],
            "chr": ["chr1"] * 4,
            "start0": [100, 200, 300, 400],
            "end0": [200, 400, 600, 800],
            "Length": [100, 200, 300, 400],
        })
        csv_file = tmp_path / "merged.csv"
        df.to_csv(csv_file, index=False)

        output_dir = tmp_path / "out"
        output_dir.mkdir()
        summary = EccSummary("s1", output_dir)
        result = summary.collect_eccdna_statistics(csv_file)

        assert result["all_count"] == 4
        assert result["uecc_count"] == 2
        assert result["mecc_count"] == 1
        assert result["cecc_count"] == 1

    def test_only_uecc_present(self, tmp_path):
        """Only UeccDNA present -> mecc/cecc counts = 0."""
        df = pd.DataFrame({
            "eccDNA_id": ["U1", "U2"],
            "eccDNA_type": ["UeccDNA", "UeccDNA"],
            "chr": ["chr1", "chr2"],
            "start0": [100, 200],
            "end0": [300, 500],
            "Length": [200, 300],
        })
        csv_file = tmp_path / "merged.csv"
        df.to_csv(csv_file, index=False)

        output_dir = tmp_path / "out"
        output_dir.mkdir()
        summary = EccSummary("s1", output_dir)
        result = summary.collect_eccdna_statistics(csv_file)

        assert result["all_count"] == 2
        assert result["uecc_count"] == 2
        assert result["mecc_count"] == 0
        assert result["cecc_count"] == 0

    def test_empty_dataframe(self, tmp_path):
        """Empty CSV with correct headers -> all counts = 0."""
        df = pd.DataFrame({
            "eccDNA_id": pd.Series([], dtype=str),
            "eccDNA_type": pd.Series([], dtype=str),
            "chr": pd.Series([], dtype=str),
            "start0": pd.Series([], dtype=int),
            "end0": pd.Series([], dtype=int),
            "Length": pd.Series([], dtype=int),
        })
        csv_file = tmp_path / "merged.csv"
        df.to_csv(csv_file, index=False)

        output_dir = tmp_path / "out"
        output_dir.mkdir()
        summary = EccSummary("s1", output_dir)
        result = summary.collect_eccdna_statistics(csv_file)

        assert result["all_count"] == 0
        assert result["uecc_count"] == 0
        assert result["mecc_count"] == 0
        assert result["cecc_count"] == 0

    def test_file_not_found(self, tmp_path):
        """Missing file -> empty stats (all zeros)."""
        output_dir = tmp_path / "out"
        output_dir.mkdir()
        summary = EccSummary("s1", output_dir)
        result = summary.collect_eccdna_statistics(tmp_path / "missing.csv")

        assert result["all_count"] == 0
        assert result["uecc_count"] == 0
        assert result["mecc_count"] == 0
        assert result["cecc_count"] == 0
        # Length stats should be "0" strings
        assert result["all_min_length"] == "0"

    def test_single_entry_std_zero(self, tmp_path):
        """Single eccDNA entry -> std_length should be 0."""
        df = pd.DataFrame({
            "eccDNA_id": ["U1"],
            "eccDNA_type": ["UeccDNA"],
            "chr": ["chr1"],
            "start0": [100],
            "end0": [600],
            "Length": [500],
        })
        csv_file = tmp_path / "merged.csv"
        df.to_csv(csv_file, index=False)

        output_dir = tmp_path / "out"
        output_dir.mkdir()
        summary = EccSummary("s1", output_dir)
        result = summary.collect_eccdna_statistics(csv_file)

        assert result["all_count"] == 1
        # std with a single value should be 0; formatted as "0.0" (one decimal place)
        assert result["all_std_length"] == "0.0"
        assert result["uecc_std_length"] == "0.0"


class TestCollectOverlapStatistics:
    """Tests for EccSummary.collect_overlap_statistics()."""

    def test_valid_json_with_inferred_simple_and_chimeric(self, tmp_path):
        """Valid JSON with inferred_simple and inferred_chimeric keys."""
        data = {
            "inferred_simple": {
                "total": 10,
                "overlapping": 3,
                "overlapping_pct": 30.0,
            },
            "inferred_chimeric": {
                "total": 5,
                "overlapping": 2,
                "overlapping_pct": 40.0,
            },
        }
        json_file = tmp_path / "overlap.json"
        with open(json_file, "w") as f:
            json.dump(data, f)

        output_dir = tmp_path / "out"
        output_dir.mkdir()
        summary = EccSummary("s1", output_dir)
        result = summary.collect_overlap_statistics(json_file)

        assert result["inferred_uecc_count"] == 10
        assert result["uecc_overlap_count"] == 3
        assert result["uecc_overlap_pct"] == "30.00"
        assert result["inferred_cecc_count"] == 5
        assert result["cecc_overlap_count"] == 2
        assert result["cecc_overlap_pct"] == "40.00"

    def test_missing_json_file(self, tmp_path):
        """Missing JSON file -> all zeros."""
        output_dir = tmp_path / "out"
        output_dir.mkdir()
        summary = EccSummary("s1", output_dir)
        result = summary.collect_overlap_statistics(tmp_path / "missing.json")

        assert result["inferred_uecc_count"] == 0
        assert result["uecc_overlap_count"] == 0
        assert result["uecc_overlap_pct"] == "0.00"
        assert result["inferred_cecc_count"] == 0
        assert result["cecc_overlap_count"] == 0
        assert result["cecc_overlap_pct"] == "0.00"

    def test_malformed_json(self, tmp_path):
        """Malformed JSON -> all zeros."""
        json_file = tmp_path / "bad.json"
        json_file.write_text("{not valid json!!!")

        output_dir = tmp_path / "out"
        output_dir.mkdir()
        summary = EccSummary("s1", output_dir)
        result = summary.collect_overlap_statistics(json_file)

        assert result["inferred_uecc_count"] == 0
        assert result["inferred_cecc_count"] == 0

    def test_partial_data_missing_keys(self, tmp_path):
        """JSON present but some expected keys missing -> defaults to 0."""
        data = {
            "inferred_simple": {"total": 7},
            # inferred_chimeric missing entirely
        }
        json_file = tmp_path / "partial.json"
        with open(json_file, "w") as f:
            json.dump(data, f)

        output_dir = tmp_path / "out"
        output_dir.mkdir()
        summary = EccSummary("s1", output_dir)
        result = summary.collect_overlap_statistics(json_file)

        assert result["inferred_uecc_count"] == 7
        assert result["uecc_overlap_count"] == 0
        # chimeric section absent -> defaults
        assert result["inferred_cecc_count"] == 0
        assert result["cecc_overlap_count"] == 0


class TestGetEmptyEccdnaStats:
    """Tests for EccSummary._get_empty_eccdna_stats()."""

    def test_returns_correct_keys(self, tmp_path):
        """Should return keys for all/uecc/mecc/cecc prefixes with count and length suffixes."""
        output_dir = tmp_path / "out"
        output_dir.mkdir()
        summary = EccSummary("s1", output_dir)
        stats = summary._get_empty_eccdna_stats()

        expected_prefixes = ["all", "uecc", "mecc", "cecc"]
        expected_suffixes = [
            "min_length", "max_length", "mean_length",
            "median_length", "mode_length", "std_length",
        ]
        for prefix in expected_prefixes:
            assert f"{prefix}_count" in stats
            for suffix in expected_suffixes:
                assert f"{prefix}_{suffix}" in stats

    def test_all_values_are_zero(self, tmp_path):
        """Count values should be 0 (int) and length values should be '0' (str)."""
        output_dir = tmp_path / "out"
        output_dir.mkdir()
        summary = EccSummary("s1", output_dir)
        stats = summary._get_empty_eccdna_stats()

        for prefix in ["all", "uecc", "mecc", "cecc"]:
            assert stats[f"{prefix}_count"] == 0
            assert isinstance(stats[f"{prefix}_count"], int)
            for suffix in [
                "min_length", "max_length", "mean_length",
                "median_length", "mode_length", "std_length",
            ]:
                assert stats[f"{prefix}_{suffix}"] == "0"
                assert isinstance(stats[f"{prefix}_{suffix}"], str)


class TestGenerateHtmlReport:
    """Tests for EccSummary.generate_html_report()."""

    def _populate_summary(self, summary):
        """Populate a summary with all required stats for report generation."""
        summary.read_stats = {
            "total_reads": 100, "total_sequences": 100, "total_length": 50000,
            "ctcr_reads": 60, "ctcr_percentage": 60.0,
            "other_reads": 40, "other_percentage": 40.0,
        }
        summary.ctcr_stats = {
            "total_reads": 100, "ctcr_reads": 60, "circular_reads": 40,
            "ctcr_perfect_count": 40, "ctcr_hybrid_count": 15,
            "ctcr_inversion_count": 5,
            "ctcr_perfect_pct_ctcr": 66.7, "ctcr_hybrid_pct_ctcr": 25.0,
            "ctcr_inversion_pct_ctcr": 8.3,
            "ctcr_perfect_pct_total": 40.0, "ctcr_hybrid_pct_total": 15.0,
            "ctcr_inversion_pct_total": 5.0,
        }
        summary.eccdna_stats = summary._get_empty_eccdna_stats()
        summary.eccdna_stats.update({"all_count": 10, "uecc_count": 5, "mecc_count": 3, "cecc_count": 2})
        summary.overlap_stats = {
            "inferred_uecc_count": 3, "uecc_overlap_count": 1, "uecc_overlap_pct": "33.33",
            "inferred_cecc_count": 2, "cecc_overlap_count": 0, "cecc_overlap_pct": "0.00",
        }

    def test_custom_output_path(self, tmp_path):
        """HTML report written to custom output path."""
        output_dir = tmp_path / "out"
        output_dir.mkdir()
        summary = EccSummary("mysample", output_dir)
        self._populate_summary(summary)

        custom_path = tmp_path / "custom_report.html"
        summary.generate_html_report(output_path=custom_path)
        assert custom_path.exists()
        assert custom_path.read_text().startswith("<!DOCTYPE html>")

    def test_content_includes_metadata(self, tmp_path):
        """HTML includes sample_name, version, and analysis_date."""
        output_dir = tmp_path / "out"
        output_dir.mkdir()
        summary = EccSummary("my_sample_42", output_dir)
        self._populate_summary(summary)

        html = summary.generate_html_report()
        assert "my_sample_42" in html
        assert summary.version in html
        # analysis_date is today's date in YYYY-MM-DD format
        today_str = datetime.now().strftime("%Y-%m-%d")
        assert today_str in html

    def test_html_structure_tags(self, tmp_path):
        """HTML output has proper structure: DOCTYPE, html, head, body."""
        output_dir = tmp_path / "out"
        output_dir.mkdir()
        summary = EccSummary("s1", output_dir)
        self._populate_summary(summary)

        html = summary.generate_html_report()
        assert "<!DOCTYPE html>" in html
        assert "<html" in html
        assert "<head>" in html
        assert "<body>" in html
        assert "</html>" in html


class TestGenerateTextSummary:
    """Tests for EccSummary.generate_text_summary()."""

    def _populate_summary(self, summary):
        """Populate a summary with all required stats for text generation."""
        summary.read_stats = {
            "total_reads": 200, "total_sequences": 200, "total_length": 100000,
            "ctcr_reads": 120, "ctcr_percentage": 60.0,
            "other_reads": 80, "other_percentage": 40.0,
        }
        summary.ctcr_stats = {
            "total_reads": 200, "ctcr_reads": 120, "circular_reads": 80,
            "ctcr_perfect_count": 80, "ctcr_hybrid_count": 30,
            "ctcr_inversion_count": 10,
            "ctcr_perfect_pct_ctcr": 66.7, "ctcr_hybrid_pct_ctcr": 25.0,
            "ctcr_inversion_pct_ctcr": 8.3,
            "ctcr_perfect_pct_total": 40.0, "ctcr_hybrid_pct_total": 15.0,
            "ctcr_inversion_pct_total": 5.0,
        }
        summary.eccdna_stats = summary._get_empty_eccdna_stats()
        summary.eccdna_stats.update({"all_count": 20, "uecc_count": 10, "mecc_count": 7, "cecc_count": 3})
        summary.overlap_stats = {
            "inferred_uecc_count": 5, "uecc_overlap_count": 2, "uecc_overlap_pct": "40.00",
            "inferred_cecc_count": 3, "cecc_overlap_count": 1, "cecc_overlap_pct": "33.33",
        }

    def test_contains_section_headers(self, tmp_path):
        """Text summary should contain the three main section headers."""
        output_dir = tmp_path / "out"
        output_dir.mkdir()
        summary = EccSummary("s1", output_dir)
        self._populate_summary(summary)

        text = summary.generate_text_summary()
        assert "READ CLASSIFICATION" in text
        assert "eccDNA STATISTICS" in text
        assert "INFERENCE ANALYSIS" in text

    def test_contains_version_and_sample(self, tmp_path):
        """Text summary contains version string and sample name."""
        output_dir = tmp_path / "out"
        output_dir.mkdir()
        summary = EccSummary("demo_sample", output_dir)
        self._populate_summary(summary)

        text = summary.generate_text_summary()
        assert "demo_sample" in text
        assert summary.version in text
        assert "CircleSeeker" in text

    def test_custom_output_path(self, tmp_path):
        """Text summary written to custom output path."""
        output_dir = tmp_path / "out"
        output_dir.mkdir()
        summary = EccSummary("s1", output_dir)
        self._populate_summary(summary)

        custom_path = tmp_path / "my_summary.txt"
        summary.generate_text_summary(output_path=custom_path)
        assert custom_path.exists()
        content = custom_path.read_text()
        assert "READ CLASSIFICATION" in content


class TestProcessAll:
    """Tests for EccSummary.process_all()."""

    def _create_test_files(self, tmp_path):
        """Create a complete set of valid test input files."""
        fasta = tmp_path / "reads.fasta"
        fasta.write_text(">r1\nACGT\n>r2\nTTTT\n>r3\nGGGG\n")

        processed_csv = tmp_path / "processed.csv"
        pd.DataFrame({
            "readName": ["r1", "r2", "r3"],
            "readClass": ["CtcR-perfect", "CtcR-hybrid", "Normal"],
        }).to_csv(processed_csv, index=False)

        merged_csv = tmp_path / "merged.csv"
        pd.DataFrame({
            "eccDNA_id": ["U1", "M1"],
            "eccDNA_type": ["UeccDNA", "MeccDNA"],
            "chr": ["chr1", "chr2"],
            "start0": [100, 200],
            "end0": [300, 500],
            "Length": [200, 300],
        }).to_csv(merged_csv, index=False)

        overlap_json = tmp_path / "overlap.json"
        with open(overlap_json, "w") as f:
            json.dump({
                "inferred_simple": {"total": 4, "overlapping": 1, "overlapping_pct": 25.0},
                "inferred_chimeric": {"total": 2, "overlapping": 0, "overlapping_pct": 0.0},
            }, f)

        return fasta, processed_csv, merged_csv, overlap_json

    def test_all_valid_files_returns_tuple(self, tmp_path):
        """process_all with valid files returns a tuple of (html, text)."""
        fasta, processed_csv, merged_csv, overlap_json = self._create_test_files(tmp_path)

        output_dir = tmp_path / "out"
        output_dir.mkdir()
        summary = EccSummary("full_test", output_dir)
        result = summary.process_all(fasta, processed_csv, merged_csv, overlap_json)

        assert isinstance(result, tuple)
        assert len(result) == 2
        html_content, text_content = result
        assert isinstance(html_content, str)
        assert isinstance(text_content, str)
        # Both should be non-empty
        assert len(html_content) > 0
        assert len(text_content) > 0

    def test_all_valid_files_creates_report_files(self, tmp_path):
        """process_all creates both HTML and text files."""
        fasta, processed_csv, merged_csv, overlap_json = self._create_test_files(tmp_path)

        output_dir = tmp_path / "out"
        output_dir.mkdir()
        summary = EccSummary("full_test", output_dir)
        summary.process_all(fasta, processed_csv, merged_csv, overlap_json)

        html_file = output_dir / "full_test_report.html"
        text_file = output_dir / "full_test_summary.txt"
        assert html_file.exists()
        assert text_file.exists()

    def test_missing_files_graceful(self, tmp_path):
        """process_all with missing files should not raise, just produce reports with zeros."""
        output_dir = tmp_path / "out"
        output_dir.mkdir()
        summary = EccSummary("missing_test", output_dir)

        # All files missing
        result = summary.process_all(
            tmp_path / "missing.fasta",
            tmp_path / "missing_processed.csv",
            tmp_path / "missing_merged.csv",
            tmp_path / "missing_overlap.json",
        )

        assert isinstance(result, tuple)
        html_content, text_content = result
        assert "missing_test" in html_content
        assert "missing_test" in text_content
