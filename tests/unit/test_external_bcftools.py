"""Tests for external Bcftools tools."""

from pathlib import Path
import sys
import pytest
from unittest.mock import patch, MagicMock

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.external.bcftools import Bcftools


class TestBcftools:
    """Test cases for Bcftools."""

    @patch.object(Bcftools, '_check_installation')
    def test_bcftools_initialization(self, mock_check):
        """Test Bcftools initialization."""
        bcftools = Bcftools()
        assert bcftools.tool_name == "bcftools"
        assert bcftools.threads == 1  # Default from base class
        mock_check.assert_called_once()

    @patch.object(Bcftools, '_check_installation')
    def test_bcftools_initialization_with_threads(self, mock_check):
        """Test Bcftools initialization with threads."""
        bcftools = Bcftools(threads=4)
        assert bcftools.threads == 4

    @patch.object(Bcftools, '_check_installation')
    @patch.object(Bcftools, 'run')
    def test_sort(self, mock_run, mock_check, tmp_path):
        """Test sort method."""
        # Configure mock to return expected tuple
        mock_run.return_value = ('', '')

        input_vcf = tmp_path / "input.vcf"
        output_bcf = tmp_path / "output" / "sorted.bcf"

        # Create input file
        input_vcf.write_text("##fileformat=VCFv4.2\n")

        bcftools = Bcftools()
        bcftools.sort(input_vcf, output_bcf)

        # Verify run was called with correct command
        mock_run.assert_called_once()
        call_args = mock_run.call_args[0][0]
        assert "bcftools" in call_args
        assert "sort" in call_args
        assert "-m" in call_args
        assert "4G" in call_args  # Default memory limit
        assert "-O" in call_args
        assert "b" in call_args  # BCF output format
        assert "-o" in call_args
        assert str(output_bcf) in call_args
        assert str(input_vcf) in call_args

        # Verify capture_output parameter (code uses True)
        assert mock_run.call_args[1]['capture_output'] is True

        # Verify output directory is created
        assert output_bcf.parent.exists()

    @patch.object(Bcftools, '_check_installation')
    @patch.object(Bcftools, 'run')
    def test_sort_custom_memory(self, mock_run, mock_check, tmp_path):
        """Test sort method with custom memory limit."""
        # Configure mock to return expected tuple
        mock_run.return_value = ('', '')

        input_vcf = tmp_path / "input.vcf"
        output_bcf = tmp_path / "sorted.bcf"

        input_vcf.write_text("##fileformat=VCFv4.2\n")

        bcftools = Bcftools()
        bcftools.sort(input_vcf, output_bcf, memory_limit="8G")

        call_args = mock_run.call_args[0][0]
        # Find the memory limit in the command
        memory_index = call_args.index("-m")
        assert call_args[memory_index + 1] == "8G"

    @patch.object(Bcftools, '_check_installation')
    @patch.object(Bcftools, 'run')
    def test_index(self, mock_run, mock_check, tmp_path):
        """Test index method."""
        # Configure mock to return expected tuple
        mock_run.return_value = ('', '')

        bcf_file = tmp_path / "test.bcf"
        bcf_file.touch()

        bcftools = Bcftools()
        bcftools.index(bcf_file)

        # Verify run was called with correct command
        mock_run.assert_called_once()
        call_args = mock_run.call_args[0][0]
        assert "bcftools" in call_args
        assert "index" in call_args
        assert str(bcf_file) in call_args
        assert "-f" not in call_args  # force=False by default

        # Verify capture_output parameter (code uses True)
        assert mock_run.call_args[1]['capture_output'] is True

    @patch.object(Bcftools, '_check_installation')
    @patch.object(Bcftools, 'run')
    def test_index_force(self, mock_run, mock_check, tmp_path):
        """Test index method with force option."""
        # Configure mock to return expected tuple
        mock_run.return_value = ('', '')

        bcf_file = tmp_path / "test.bcf"
        bcf_file.touch()

        bcftools = Bcftools()
        bcftools.index(bcf_file, force=True)

        call_args = mock_run.call_args[0][0]
        assert "bcftools" in call_args
        assert "index" in call_args
        assert "-f" in call_args  # force=True
        assert str(bcf_file) in call_args

    @patch.object(Bcftools, '_check_installation')
    @patch.object(Bcftools, 'run')
    def test_sort_and_index_workflow(self, mock_run, mock_check, tmp_path):
        """Test typical sort and index workflow."""
        # Configure mock to return expected tuple
        mock_run.return_value = ('', '')

        input_vcf = tmp_path / "input.vcf"
        output_bcf = tmp_path / "sorted.bcf"

        input_vcf.write_text("##fileformat=VCFv4.2\n")

        bcftools = Bcftools()

        # Sort then index
        bcftools.sort(input_vcf, output_bcf, memory_limit="2G")
        bcftools.index(output_bcf, force=True)

        # Verify both operations were called
        assert mock_run.call_count == 2

        # Check first call (sort)
        first_call = mock_run.call_args_list[0][0][0]
        assert "sort" in first_call
        assert "2G" in first_call

        # Check second call (index)
        second_call = mock_run.call_args_list[1][0][0]
        assert "index" in second_call
        assert "-f" in second_call

    @patch.object(Bcftools, '_check_installation')
    @patch.object(Bcftools, 'run')
    def test_error_handling(self, mock_run, mock_check, tmp_path):
        """Test error handling when bcftools operations fail."""
        mock_run.side_effect = Exception("Bcftools command failed")

        input_vcf = tmp_path / "input.vcf"
        output_bcf = tmp_path / "output.bcf"
        input_vcf.write_text("##fileformat=VCFv4.2\n")

        bcftools = Bcftools()

        with pytest.raises(Exception):
            bcftools.sort(input_vcf, output_bcf)

    @patch.object(Bcftools, '_check_installation')
    @patch.object(Bcftools, 'run')
    def test_sort_different_memory_units(self, mock_run, mock_check, tmp_path):
        """Test sort with different memory unit formats."""
        # Configure mock to return expected tuple
        mock_run.return_value = ('', '')

        input_vcf = tmp_path / "input.vcf"
        output_bcf = tmp_path / "output.bcf"
        input_vcf.write_text("##fileformat=VCFv4.2\n")

        bcftools = Bcftools()

        # Test different memory formats
        memory_formats = ["1G", "2048M", "500m", "16G"]
        for memory in memory_formats:
            mock_run.reset_mock()
            bcftools.sort(input_vcf, output_bcf, memory_limit=memory)
            call_args = mock_run.call_args[0][0]
            memory_index = call_args.index("-m")
            assert call_args[memory_index + 1] == memory

    @patch.object(Bcftools, '_check_installation')
    @patch.object(Bcftools, 'run')
    def test_sort_output_format(self, mock_run, mock_check, tmp_path):
        """Test that sort always outputs BCF format."""
        # Configure mock to return expected tuple
        mock_run.return_value = ('', '')

        input_vcf = tmp_path / "input.vcf"
        output_bcf = tmp_path / "output.bcf"
        input_vcf.write_text("##fileformat=VCFv4.2\n")

        bcftools = Bcftools()
        bcftools.sort(input_vcf, output_bcf)

        call_args = mock_run.call_args[0][0]
        # Verify BCF output format is specified
        output_index = call_args.index("-O")
        assert call_args[output_index + 1] == "b"


@pytest.fixture
def sample_vcf_file(tmp_path):
    """Create a sample VCF file for testing."""
    vcf_file = tmp_path / "sample.vcf"
    vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	.	A	T	60	PASS	DP=30
chr1	200	.	G	C	45	PASS	DP=25
chr1	300	.	C	G	55	PASS	DP=35"""
    vcf_file.write_text(vcf_content)
    return vcf_file


@pytest.fixture
def sample_bcf_file(tmp_path):
    """Create a sample BCF file for testing."""
    bcf_file = tmp_path / "sample.bcf"
    # Create a dummy BCF file (in real scenario this would be binary)
    bcf_file.write_bytes(b"BCF\x02\x02")  # Minimal BCF header
    return bcf_file