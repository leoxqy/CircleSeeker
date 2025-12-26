"""Tests for external CD-HIT tools."""

from pathlib import Path
import sys
import pytest
from unittest.mock import patch, MagicMock
import csv

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.external.cd_hit import CDHitEst, CDHitRunner, CDHitConfig


class TestCDHitConfig:
    """Test cases for CDHitConfig."""

    def test_cdhit_config_defaults(self):
        """Test CDHitConfig default values."""
        config = CDHitConfig()
        assert config.similarity_threshold == 0.99
        assert config.threads == 8
        assert config.memory_limit == 8000
        assert config.word_length == 10
        assert config.length_diff_cutoff == 0.99
        assert config.alignment_coverage_short == 0.95
        assert config.alignment_coverage_long == 0.95
        assert config.include_gaps == 1
        assert config.desc_len == 0

    def test_cdhit_config_custom(self):
        """Test CDHitConfig with custom values."""
        config = CDHitConfig(
            similarity_threshold=0.95,
            threads=16,
            memory_limit=16000,
            word_length=8,
            length_diff_cutoff=0.9,
            alignment_coverage_short=0.8,
            alignment_coverage_long=0.8,
            include_gaps=0,
            desc_len=20
        )
        assert config.similarity_threshold == 0.95
        assert config.threads == 16
        assert config.memory_limit == 16000
        assert config.word_length == 8
        assert config.length_diff_cutoff == 0.9
        assert config.alignment_coverage_short == 0.8
        assert config.alignment_coverage_long == 0.8
        assert config.include_gaps == 0
        assert config.desc_len == 20


class TestCDHitEst:
    """Test cases for CDHitEst."""

    @patch.object(CDHitEst, '_check_installation')
    def test_cdhit_initialization_default(self, mock_check):
        """Test CDHitEst initialization with defaults."""
        cdhit = CDHitEst()
        assert cdhit.tool_name == "cd-hit-est"
        assert cdhit.threads == 24  # Default
        assert cdhit.c == 0.99
        assert cdhit.n == 10
        assert cdhit.s == 0.99
        assert cdhit.aS == 0.95
        assert cdhit.aL == 0.95
        assert cdhit.G == 1
        assert cdhit.d == 0
        assert cdhit.M == 0
        mock_check.assert_called_once()

    @patch.object(CDHitEst, '_check_installation')
    def test_cdhit_initialization_with_config(self, mock_check):
        """Test CDHitEst initialization with config."""
        config = CDHitConfig(threads=4, similarity_threshold=0.95)
        cdhit = CDHitEst(config=config)
        assert cdhit.threads == 4
        assert cdhit.c == 0.95

    @patch.object(CDHitEst, '_check_installation')
    def test_cdhit_initialization_custom_params(self, mock_check):
        """Test CDHitEst with custom parameters."""
        cdhit = CDHitEst(
            threads=8,
            c=0.98,
            n=8,
            s=0.95,
            aS=0.9,
            aL=0.9,
            G=0,
            d=50,
            M=4000,
            extra=["-verbose"]
        )
        assert cdhit.threads == 8
        assert cdhit.c == 0.98
        assert cdhit.n == 8
        assert cdhit.s == 0.95
        assert cdhit.aS == 0.9
        assert cdhit.aL == 0.9
        assert cdhit.G == 0
        assert cdhit.d == 50
        assert cdhit.M == 4000
        assert cdhit.extra == ["-verbose"]

    @patch.object(CDHitEst, '_check_installation')
    def test_build_command(self, mock_check, tmp_path):
        """Test _build_command method."""
        input_fasta = tmp_path / "input.fasta"
        output_prefix = tmp_path / "output"

        cdhit = CDHitEst(threads=4, c=0.95, n=8)
        command = cdhit._build_command(input_fasta, output_prefix)

        assert "cd-hit-est" in command
        assert "-i" in command and str(input_fasta) in command
        assert "-o" in command and str(output_prefix) in command
        assert "-c" in command and "0.95" in command
        assert "-n" in command and "8" in command
        assert "-T" in command and "4" in command

    @patch.object(CDHitEst, '_check_installation')
    def test_build_command_with_extra(self, mock_check, tmp_path):
        """Test _build_command with extra arguments."""
        input_fasta = tmp_path / "input.fasta"
        output_prefix = tmp_path / "output"

        cdhit = CDHitEst(extra=["-verbose", "-print_overlap"])
        command = cdhit._build_command(input_fasta, output_prefix)

        assert "-verbose" in command
        assert "-print_overlap" in command

    @patch.object(CDHitEst, '_check_installation')
    @patch.object(CDHitEst, 'run')
    def test_cluster_sequences(self, mock_run, mock_check, tmp_path):
        """Test cluster_sequences method."""
        mock_run.return_value = ("stdout output", "stderr output")

        input_fasta = tmp_path / "input.fasta"
        output_prefix = tmp_path / "clustered"
        clstr_file = tmp_path / "clustered.clstr"

        input_fasta.write_text(">seq1\nATCG\n>seq2\nGCTA\n")
        clstr_file.touch()  # Simulate cluster file creation

        cdhit = CDHitEst()
        result = cdhit.cluster_sequences(input_fasta, output_prefix)

        mock_run.assert_called_once()
        command = mock_run.call_args[0][0]
        assert "cd-hit-est" in command
        assert str(input_fasta) in command
        assert str(output_prefix) in command
        assert result == clstr_file

    @patch.object(CDHitEst, '_check_installation')
    @patch.object(CDHitEst, 'run')
    def test_cluster_sequences_failure(self, mock_run, mock_check, tmp_path):
        """Test cluster_sequences with failure."""
        mock_run.side_effect = Exception("CD-HIT-EST failed")

        input_fasta = tmp_path / "input.fasta"
        output_prefix = tmp_path / "output"

        cdhit = CDHitEst()

        with pytest.raises(Exception):
            cdhit.cluster_sequences(input_fasta, output_prefix)

    def test_parse_clstr(self, tmp_path):
        """Test _parse_clstr method."""
        # Create sample cluster file
        clstr_content = """>Cluster 0
0	100nt, >seq1... *
1	95nt, >seq2... at 95.0%
>Cluster 1
0	80nt, >seq3... *
>Cluster 2
0	120nt, >seq4... *
1	110nt, >seq5... at 90.0%
2	115nt, >seq6... at 88.0%"""

        clstr_file = tmp_path / "test.clstr"
        clstr_file.write_text(clstr_content)

        members, representatives = CDHitEst._parse_clstr(clstr_file)

        # Check members
        assert "0" in members
        assert "1" in members
        assert "2" in members
        assert len(members["0"]) == 2  # seq1, seq2
        assert len(members["1"]) == 1  # seq3
        assert len(members["2"]) == 3  # seq4, seq5, seq6

        # Check representatives
        assert representatives["0"] == "seq1"
        assert representatives["1"] == "seq3"
        assert representatives["2"] == "seq4"

    def test_parse_clstr_missing_file(self, tmp_path):
        """Test _parse_clstr with missing file."""
        missing_file = tmp_path / "missing.clstr"

        with pytest.raises(FileNotFoundError):
            CDHitEst._parse_clstr(missing_file)

    def test_write_id2cluster_csv(self, tmp_path):
        """Test _write_id2cluster_csv method."""
        members = {
            "0": ["seq1", "seq2"],
            "1": ["seq3"],
            "2": ["seq4", "seq5", "seq6"]
        }
        representatives = {
            "0": "seq1",
            "1": "seq3",
            "2": "seq4"
        }

        output_csv = tmp_path / "clusters.csv"
        result = CDHitEst._write_id2cluster_csv(members, representatives, output_csv)

        assert result == output_csv
        assert output_csv.exists()

        # Check CSV content
        with output_csv.open('r') as f:
            reader = csv.DictReader(f)
            rows = list(reader)

        assert len(rows) == 6  # Total sequences

        # Check header
        assert "id" in rows[0]
        assert "cluster" in rows[0]
        assert "is_representative" in rows[0]
        assert "cluster_size" in rows[0]

        # Check some specific rows
        seq1_row = next(r for r in rows if r["id"] == "seq1")
        assert seq1_row["cluster"] == "0"
        assert seq1_row["is_representative"] == "true"
        assert seq1_row["cluster_size"] == "2"

        seq2_row = next(r for r in rows if r["id"] == "seq2")
        assert seq2_row["cluster"] == "0"
        assert seq2_row["is_representative"] == "false"
        assert seq2_row["cluster_size"] == "2"

    @patch.object(CDHitEst, '_check_installation')
    def test_export_id2cluster(self, mock_check, tmp_path):
        """Test export_id2cluster method."""
        # Create sample cluster file
        clstr_content = """>Cluster 0
0	100nt, >seq1... *
1	95nt, >seq2... at 95.0%
>Cluster 1
0	80nt, >seq3... *"""

        clstr_file = tmp_path / "test.clstr"
        clstr_file.write_text(clstr_content)

        cdhit = CDHitEst()
        result = cdhit.export_id2cluster(clstr_file)

        expected_csv = tmp_path / "test.id2cluster.csv"
        assert result == expected_csv
        assert expected_csv.exists()

    @patch.object(CDHitEst, '_check_installation')
    def test_export_id2cluster_custom_output(self, mock_check, tmp_path):
        """Test export_id2cluster with custom output path."""
        clstr_content = """>Cluster 0
0	100nt, >seq1... *"""

        clstr_file = tmp_path / "test.clstr"
        custom_csv = tmp_path / "custom_clusters.csv"
        clstr_file.write_text(clstr_content)

        cdhit = CDHitEst()
        result = cdhit.export_id2cluster(clstr_file, custom_csv)

        assert result == custom_csv
        assert custom_csv.exists()

    @patch.object(CDHitEst, '_check_installation')
    def test_get_cluster_stats(self, mock_check, tmp_path):
        """Test get_cluster_stats method."""
        clstr_content = """>Cluster 0
0	100nt, >seq1... *
1	95nt, >seq2... at 95.0%
>Cluster 1
0	80nt, >seq3... *
>Cluster 2
0	120nt, >seq4... *
1	110nt, >seq5... at 90.0%
2	115nt, >seq6... at 88.0%"""

        clstr_file = tmp_path / "test.clstr"
        clstr_file.write_text(clstr_content)

        cdhit = CDHitEst()
        stats = cdhit.get_cluster_stats(clstr_file)

        assert stats["total_clusters"] == 3
        assert stats["total_sequences"] == 6
        assert stats["singleton_clusters"] == 1  # Cluster 1
        assert stats["multi_member_clusters"] == 2  # Clusters 0 and 2
        assert stats["avg_cluster_size"] == 2.0
        assert stats["reduction_rate"] == 50.0  # (1 - 3/6) * 100

    @patch.object(CDHitEst, '_check_installation')
    def test_get_cluster_stats_missing_file(self, mock_check, tmp_path):
        """Test get_cluster_stats with missing file."""
        missing_file = tmp_path / "missing.clstr"

        cdhit = CDHitEst()
        stats = cdhit.get_cluster_stats(missing_file)

        assert stats == {}


class TestCDHitRunner:
    """Test cases for CDHitRunner."""

    @patch.object(CDHitEst, '_check_installation')
    def test_cdhit_runner_initialization(self, mock_check):
        """Test CDHitRunner initialization."""
        runner = CDHitRunner(threads=8, c=0.95)
        assert runner.threads == 8
        assert runner.c == 0.95
        mock_check.assert_called_once()

    @patch.object(CDHitEst, '_check_installation')
    @patch.object(CDHitEst, 'cluster_sequences')
    def test_cdhit_runner_run(self, mock_cluster, mock_check, tmp_path):
        """Test CDHitRunner run method."""
        input_fasta = tmp_path / "input.fasta"
        output_prefix = tmp_path / "output"

        runner = CDHitRunner()
        runner.run(input_fasta, output_prefix)

        mock_cluster.assert_called_once_with(input_fasta, output_prefix)


@pytest.fixture
def sample_fasta_file(tmp_path):
    """Create a sample FASTA file for testing."""
    fasta_file = tmp_path / "sample.fasta"
    fasta_content = """>seq1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>seq2
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCC
>seq3
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
>seq4
TTTTAAAAGGGGCCCCTTTTAAAAGGGGCCCCTTTTAAAAGGGGCCCCTTTT"""
    fasta_file.write_text(fasta_content)
    return fasta_file


@pytest.fixture
def sample_clstr_file(tmp_path):
    """Create a sample .clstr file for testing."""
    clstr_file = tmp_path / "sample.clstr"
    clstr_content = """>Cluster 0
0	64nt, >seq1... *
1	64nt, >seq2... at 98.44%
>Cluster 1
0	64nt, >seq3... *
>Cluster 2
0	64nt, >seq4... *"""
    clstr_file.write_text(clstr_content)
    return clstr_file