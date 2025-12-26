"""Tests for external Samtools tools."""

from pathlib import Path
import sys
import pytest
from unittest.mock import patch, MagicMock

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.external.samtools import Samtools


class TestSamtools:
    """Test cases for Samtools."""

    @patch.object(Samtools, '_check_installation')
    def test_samtools_initialization(self, mock_check):
        """Test Samtools initialization."""
        samtools = Samtools()
        assert samtools.tool_name == "samtools"
        assert samtools.threads == 1  # Default from base class
        mock_check.assert_called_once()

    @patch.object(Samtools, '_check_installation')
    def test_samtools_initialization_with_threads(self, mock_check):
        """Test Samtools initialization with threads."""
        samtools = Samtools(threads=4)
        assert samtools.threads == 4

    @patch.object(Samtools, '_check_installation')
    @patch.object(Samtools, 'run')
    def test_sort_bam(self, mock_run, mock_check, tmp_path):
        """Test sort_bam method."""
        # Configure mock to return expected tuple
        mock_run.return_value = ('', '')

        input_sam = tmp_path / "input.sam"
        output_bam = tmp_path / "output" / "sorted.bam"

        # Create input file
        input_sam.write_text("@HD\tVN:1.6\n")

        samtools = Samtools(threads=2)
        samtools.sort_bam(input_sam, output_bam)

        # Verify run was called with correct command
        mock_run.assert_called_once()
        call_args = mock_run.call_args[0][0]
        assert "samtools" in call_args
        assert "sort" in call_args
        assert "-@" in call_args
        assert "2" in call_args
        assert "-o" in call_args
        assert str(output_bam) in call_args
        assert str(input_sam) in call_args

        # Verify capture_output parameter (code uses True)
        assert mock_run.call_args[1]['capture_output'] is True

        # Verify output directory is created
        assert output_bam.parent.exists()

    @patch.object(Samtools, '_check_installation')
    @patch.object(Samtools, 'run')
    def test_index_bam(self, mock_run, mock_check, tmp_path):
        """Test index_bam method."""
        # Configure mock to return expected tuple
        mock_run.return_value = ('', '')

        bam_file = tmp_path / "test.bam"
        bam_file.touch()

        samtools = Samtools()
        samtools.index_bam(bam_file)

        # Verify run was called with correct command
        mock_run.assert_called_once()
        call_args = mock_run.call_args[0][0]
        assert "samtools" in call_args
        assert "index" in call_args
        assert str(bam_file) in call_args

        # Verify capture_output parameter (code uses True)
        assert mock_run.call_args[1]['capture_output'] is True

    @patch.object(Samtools, '_check_installation')
    @patch.object(Samtools, 'run')
    def test_faidx(self, mock_run, mock_check, tmp_path):
        """Test faidx method."""
        # Configure mock to return expected tuple
        mock_run.return_value = ('', '')

        reference_fasta = tmp_path / "reference.fasta"
        reference_fasta.write_text(">chr1\nATCGATCG\n")

        samtools = Samtools()
        samtools.faidx(reference_fasta)

        # Verify run was called with correct command
        mock_run.assert_called_once()
        call_args = mock_run.call_args[0][0]
        assert "samtools" in call_args
        assert "faidx" in call_args
        assert str(reference_fasta) in call_args

        # Verify capture_output parameter (code uses True)
        assert mock_run.call_args[1]['capture_output'] is True

    @patch.object(Samtools, '_check_installation')
    @patch.object(Samtools, 'run')
    def test_sort_bam_with_custom_threads(self, mock_run, mock_check, tmp_path):
        """Test sort_bam with custom thread count."""
        # Configure mock to return expected tuple
        mock_run.return_value = ('', '')

        input_sam = tmp_path / "input.sam"
        output_bam = tmp_path / "sorted.bam"

        input_sam.write_text("@HD\tVN:1.6\n")

        samtools = Samtools(threads=8)
        samtools.sort_bam(input_sam, output_bam)

        call_args = mock_run.call_args[0][0]
        # Find the thread count in the command
        thread_index = call_args.index("-@")
        assert call_args[thread_index + 1] == "8"

    @patch.object(Samtools, '_check_installation')
    @patch.object(Samtools, 'run')
    def test_multiple_operations(self, mock_run, mock_check, tmp_path):
        """Test multiple samtools operations in sequence."""
        # Configure mock to return expected tuple
        mock_run.return_value = ('', '')

        sam_file = tmp_path / "input.sam"
        bam_file = tmp_path / "sorted.bam"
        ref_file = tmp_path / "reference.fasta"

        sam_file.write_text("@HD\tVN:1.6\n")
        ref_file.write_text(">chr1\nATCG\n")

        samtools = Samtools(threads=4)

        # Sort, index, and create reference index
        samtools.sort_bam(sam_file, bam_file)
        samtools.index_bam(bam_file)
        samtools.faidx(ref_file)

        # Verify all three operations were called
        assert mock_run.call_count == 3

        # Check first call (sort)
        first_call = mock_run.call_args_list[0][0][0]
        assert "sort" in first_call

        # Check second call (index)
        second_call = mock_run.call_args_list[1][0][0]
        assert "index" in second_call

        # Check third call (faidx)
        third_call = mock_run.call_args_list[2][0][0]
        assert "faidx" in third_call

    @patch.object(Samtools, '_check_installation')
    @patch.object(Samtools, 'run')
    def test_error_handling(self, mock_run, mock_check, tmp_path):
        """Test error handling when samtools operations fail."""
        mock_run.side_effect = Exception("Samtools command failed")

        sam_file = tmp_path / "input.sam"
        bam_file = tmp_path / "output.bam"
        sam_file.write_text("@HD\tVN:1.6\n")

        samtools = Samtools()

        with pytest.raises(Exception):
            samtools.sort_bam(sam_file, bam_file)


@pytest.fixture
def sample_sam_file(tmp_path):
    """Create a sample SAM file for testing."""
    sam_file = tmp_path / "sample.sam"
    sam_content = """@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:248956422
@PG	ID:minimap2	PN:minimap2	VN:2.24-r1122
read1	0	chr1	100	60	50M	*	0	0	ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG	*
read2	0	chr1	200	60	50M	*	0	0	GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA	*"""
    sam_file.write_text(sam_content)
    return sam_file


@pytest.fixture
def sample_reference_fasta(tmp_path):
    """Create a sample reference FASTA file for testing."""
    ref_file = tmp_path / "reference.fasta"
    ref_content = """>chr1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>chr2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"""
    ref_file.write_text(ref_content)
    return ref_file