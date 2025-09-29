"""Tests for external BLAST tools."""

from pathlib import Path
import sys
import pytest
from unittest.mock import patch, MagicMock

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.external.blast import MakeBlastDB, BlastN, BlastRunner


class TestMakeBlastDB:
    """Test cases for MakeBlastDB."""

    def test_makeblastdb_initialization(self):
        """Test MakeBlastDB initialization."""
        makeblastdb = MakeBlastDB()
        assert makeblastdb.tool_name == "makeblastdb"

    @patch('subprocess.run')
    def test_makeblastdb_run(self, mock_run, tmp_path):
        """Test MakeBlastDB run method."""
        mock_run.return_value = MagicMock(returncode=0, stdout="", stderr="")

        fasta_file = tmp_path / "reference.fasta"
        fasta_file.write_text(">seq1\nATCGATCG\n")

        db_prefix = tmp_path / "blast_db"

        makeblastdb = MakeBlastDB()
        makeblastdb.build_database(fasta_file, db_prefix, dbtype="nucl")

        # Verify subprocess was called
        mock_run.assert_called_once()
        call_args = mock_run.call_args[0][0]
        assert "makeblastdb" in call_args
        assert "-in" in call_args
        assert "-out" in call_args

    @patch('subprocess.run')
    def test_makeblastdb_failure(self, mock_run, tmp_path):
        """Test MakeBlastDB failure handling."""
        mock_run.return_value = MagicMock(returncode=1, stdout="", stderr="Error")

        fasta_file = tmp_path / "reference.fasta"
        fasta_file.write_text(">seq1\nATCGATCG\n")

        db_prefix = tmp_path / "blast_db"

        makeblastdb = MakeBlastDB()

        with pytest.raises(Exception):
            makeblastdb.build_database(fasta_file, db_prefix, dbtype="nucl")


class TestBlastN:
    """Test cases for BlastN."""

    def test_blastn_initialization(self):
        """Test BlastN initialization."""
        blastn = BlastN()
        assert blastn.tool_name == "blastn"

    @patch('subprocess.run')
    def test_blastn_run(self, mock_run, tmp_path):
        """Test BlastN run method."""
        # Mock successful BLAST output
        blast_output = "query1\tsubject1\t95.0\t100\t5\t0\t1\t100\t1\t100\t1e-50\t180"
        mock_run.return_value = MagicMock(returncode=0, stdout=blast_output, stderr="")

        query_file = tmp_path / "query.fasta"
        query_file.write_text(">query1\nATCGATCG\n")

        db_path = tmp_path / "blast_db"
        output_file = tmp_path / "blast_results.tsv"

        blastn = BlastN()
        blastn.run(
            database=db_path,
            query_file=query_file,
            output_file=output_file,
            evalue=1e-5,
            max_target_seqs=100
        )

        # Verify subprocess was called
        mock_run.assert_called_once()
        call_args = mock_run.call_args[0][0]
        assert "blastn" in call_args
        assert "-query" in call_args
        assert "-db" in call_args


class TestBlastRunner:
    """Test cases for BlastRunner."""

    def test_blast_runner_initialization(self):
        """Test BlastRunner initialization."""
        runner = BlastRunner(num_threads=4, evalue=1e-5, max_target_seqs=100)
        assert runner.num_threads == 4
        assert runner.evalue == 1e-5
        assert runner.max_target_seqs == 100

    @patch('circleseeker.external.blast.BlastN')
    def test_blast_runner_run(self, mock_blastn_class, tmp_path):
        """Test BlastRunner run method."""
        # Mock BlastN instance
        mock_blastn = MagicMock()
        mock_blastn_class.return_value = mock_blastn

        query_file = tmp_path / "query.fasta"
        query_file.write_text(">query1\nATCGATCG\n")

        db_path = tmp_path / "blast_db"
        output_file = tmp_path / "blast_results.tsv"

        runner = BlastRunner(num_threads=4)
        runner.run(
            database=db_path,
            query_file=query_file,
            output_file=output_file
        )

        # Verify BlastN was called
        mock_blastn_class.assert_called_once()
        mock_blastn.run.assert_called_once()

    def test_blast_runner_default_params(self):
        """Test BlastRunner default parameters."""
        runner = BlastRunner()
        assert runner.num_threads == 1
        assert runner.evalue == 1e-5
        assert runner.max_target_seqs == 100


@pytest.fixture
def sample_fasta_file(tmp_path):
    """Create a sample FASTA file for testing."""
    fasta_file = tmp_path / "sample.fasta"
    fasta_content = """>seq1
ATCGATCGATCGATCGATCG
>seq2
GCTAGCTAGCTAGCTAGCTA
>seq3
TTTTAAAAGGGGCCCCTTTT"""
    fasta_file.write_text(fasta_content)
    return fasta_file


@pytest.fixture
def sample_blast_output():
    """Sample BLAST output in TSV format."""
    return """query1\tsubject1\t95.0\t100\t5\t0\t1\t100\t1\t100\t1e-50\t180
query2\tsubject2\t90.0\t150\t15\t0\t1\t150\t1\t150\t1e-40\t160
query3\tsubject1\t85.0\t200\t30\t0\t1\t200\t50\t250\t1e-30\t140"""