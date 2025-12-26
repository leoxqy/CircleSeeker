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

    @patch.object(MakeBlastDB, '_check_installation')
    def test_makeblastdb_initialization(self, mock_check):
        """Test MakeBlastDB initialization."""
        makeblastdb = MakeBlastDB()
        assert makeblastdb.tool_name == "makeblastdb"
        mock_check.assert_called_once()

    @patch.object(MakeBlastDB, '_check_installation')
    @patch.object(MakeBlastDB, 'run')
    def test_makeblastdb_run(self, mock_run, mock_check, tmp_path):
        """Test MakeBlastDB run method."""
        mock_run.return_value = ('', '')

        fasta_file = tmp_path / "reference.fasta"
        fasta_file.write_text(">seq1\nATCGATCG\n")

        db_prefix = tmp_path / "blast_db"

        makeblastdb = MakeBlastDB()
        makeblastdb.build_database(fasta_file, db_prefix, dbtype="nucl")

        # Verify run was called
        mock_run.assert_called_once()
        call_args = mock_run.call_args[0][0]
        assert "makeblastdb" in call_args
        assert "-in" in call_args
        assert "-out" in call_args

    @patch.object(MakeBlastDB, '_check_installation')
    @patch.object(MakeBlastDB, 'run')
    def test_makeblastdb_failure(self, mock_run, mock_check, tmp_path):
        """Test MakeBlastDB failure handling."""
        mock_run.side_effect = Exception("makeblastdb failed")

        fasta_file = tmp_path / "reference.fasta"
        fasta_file.write_text(">seq1\nATCGATCG\n")

        db_prefix = tmp_path / "blast_db"

        makeblastdb = MakeBlastDB()

        with pytest.raises(Exception):
            makeblastdb.build_database(fasta_file, db_prefix, dbtype="nucl")


class TestBlastN:
    """Test cases for BlastN."""

    @patch.object(BlastN, '_check_installation')
    def test_blastn_initialization(self, mock_check):
        """Test BlastN initialization."""
        blastn = BlastN()
        assert blastn.tool_name == "blastn"
        mock_check.assert_called_once()

    @patch.object(BlastN, '_check_installation')
    @patch.object(BlastN, 'run')
    def test_blastn_run(self, mock_run, mock_check, tmp_path):
        """Test BlastN run_blast method."""
        mock_run.return_value = ('', '')

        query_file = tmp_path / "query.fasta"
        query_file.write_text(">query1\nATCGATCG\n")

        db_path = tmp_path / "blast_db"
        output_file = tmp_path / "blast_results.tsv"

        blastn = BlastN()
        blastn.run_blast(
            database=db_path,
            query_file=query_file,
            output_file=output_file,
            evalue="1e-5",
            max_target_seqs=100
        )

        # Verify run was called
        mock_run.assert_called_once()
        call_args = mock_run.call_args[0][0]
        assert "blastn" in call_args
        assert "-query" in call_args
        assert "-db" in call_args


class TestBlastRunner:
    """Test cases for BlastRunner."""

    @patch.object(MakeBlastDB, '_check_installation')
    @patch.object(BlastN, '_check_installation')
    def test_blast_runner_initialization(self, mock_blastn_check, mock_makeblastdb_check):
        """Test BlastRunner initialization."""
        runner = BlastRunner(num_threads=4, evalue="1e-5", max_target_seqs=100)
        assert runner.num_threads == 4
        assert runner.evalue == "1e-5"
        assert runner.max_target_seqs == 100

    @patch.object(MakeBlastDB, '_check_installation')
    @patch.object(BlastN, '_check_installation')
    @patch.object(BlastN, 'run')
    def test_blast_runner_run(self, mock_run, mock_blastn_check, mock_makeblastdb_check, tmp_path):
        """Test BlastRunner run method."""
        mock_run.return_value = ('', '')

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

        # Verify BlastN.run was called (via run_blast)
        mock_run.assert_called_once()

    @patch.object(MakeBlastDB, '_check_installation')
    @patch.object(BlastN, '_check_installation')
    def test_blast_runner_default_params(self, mock_blastn_check, mock_makeblastdb_check):
        """Test BlastRunner default parameters."""
        runner = BlastRunner()
        assert runner.num_threads == 8  # Default from source
        assert runner.evalue == "1e-50"  # Default from source
        assert runner.max_target_seqs == 1000  # Default from source


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