"""Tests for external TideHunter tools."""

from pathlib import Path
import sys
import pytest
from unittest.mock import patch, MagicMock, mock_open

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.external.tidehunter import TideHunter, TideHunterRunner


class TestTideHunter:
    """Test cases for TideHunter."""

    @patch.object(TideHunter, '_check_installation')
    def test_tidehunter_initialization(self, mock_check):
        """Test TideHunter initialization."""
        tidehunter = TideHunter()
        assert tidehunter.tool_name in {"TideHunter", "tidehunter"}
        assert tidehunter.threads == 1  # Default from base class
        mock_check.assert_called_once()

    @patch.object(TideHunter, '_check_installation')
    def test_tidehunter_initialization_with_threads(self, mock_check):
        """Test TideHunter initialization with threads."""
        tidehunter = TideHunter(threads=4)
        assert tidehunter.threads == 4

    @patch.object(TideHunter, '_check_installation')
    @patch('subprocess.run')
    @patch('builtins.open', new_callable=mock_open)
    def test_tidehunter_run_analysis(self, mock_file, mock_run, mock_check, tmp_path):
        """Test TideHunter run_analysis method."""
        # Mock successful run
        mock_run.return_value = MagicMock(returncode=0, stderr="")

        input_file = tmp_path / "input.fasta"
        input_file.write_text(">seq1\nATCGATCG\n")

        output_file = tmp_path / "output" / "tidehunter_results.txt"

        # Mock file stat for logging
        with patch.object(Path, 'stat') as mock_stat:
            mock_stat.return_value = MagicMock(st_size=1024)

            tidehunter = TideHunter(threads=2)
            tidehunter.run_analysis(
                input_file=input_file,
                output_file=output_file,
                k=16,
                w=1,
                p=100,
                P=2000000,
                e=0.1,
                f=2
            )

        # Verify subprocess was called
        mock_run.assert_called_once()
        call_args = mock_run.call_args[0][0]
        assert call_args[0] in {"TideHunter", "tidehunter"}
        assert "-f" in call_args
        assert "-t" in call_args
        assert "-k" in call_args
        assert str(input_file) in call_args

        # Verify output file was opened for writing
        assert mock_file.call_count == 2
        mock_file.assert_any_call(output_file, "w")
        mock_file.assert_any_call(output_file.parent / "tidehunter.log", "w")

    @patch.object(TideHunter, '_check_installation')
    @patch('subprocess.run')
    @patch('builtins.open', new_callable=mock_open)
    def test_tidehunter_run_analysis_custom_params(self, mock_file, mock_run, mock_check, tmp_path):
        """Test TideHunter with custom parameters."""
        mock_run.return_value = MagicMock(returncode=0, stderr="")

        input_file = tmp_path / "input.fasta"
        # Use subdirectory to avoid Path.stat mock breaking mkdir
        output_file = tmp_path / "output_dir" / "output.txt"

        # Create a mock for the output file's stat method only (after it's created)
        with patch.object(Path, 'stat') as mock_stat:
            mock_stat.return_value = MagicMock(st_size=2048)

            tidehunter = TideHunter()
            tidehunter.run_analysis(
                input_file=input_file,
                output_file=output_file,
                k=20,
                w=2,
                p=150,
                P=3000000,
                e=0.2,
                f=3
            )

        call_args = mock_run.call_args[0][0]
        assert "-k" in call_args and "20" in call_args
        assert "-w" in call_args and "2" in call_args
        assert "-p" in call_args and "150" in call_args
        assert "-P" in call_args and "3000000" in call_args
        assert "-e" in call_args and "0.2" in call_args
        assert "-f" in call_args and "3" in call_args

    @patch.object(TideHunter, '_check_installation')
    @patch('subprocess.run')
    @patch('builtins.open', new_callable=mock_open)
    def test_tidehunter_failure(self, mock_file, mock_run, mock_check, tmp_path):
        """Test TideHunter failure handling."""
        # Mock failed run
        mock_run.side_effect = Exception("Command failed")
        mock_run.return_value = MagicMock(returncode=1, stderr="Error message")

        input_file = tmp_path / "input.fasta"
        output_file = tmp_path / "output.txt"

        tidehunter = TideHunter()

        # Should handle subprocess exception
        with pytest.raises(Exception):
            tidehunter.run_analysis(
                input_file=input_file,
                output_file=output_file
            )

    @patch.object(TideHunter, '_check_installation')
    @patch('subprocess.run')
    @patch('builtins.open', new_callable=mock_open)
    def test_tidehunter_subprocess_error(self, mock_file, mock_run, mock_check, tmp_path):
        """Test TideHunter subprocess error handling."""
        from subprocess import CalledProcessError

        # Mock subprocess error
        error = CalledProcessError(1, "tidehunter", stderr="TideHunter error")
        mock_run.side_effect = error

        input_file = tmp_path / "input.fasta"
        output_file = tmp_path / "output.txt"

        tidehunter = TideHunter()

        with pytest.raises(Exception):  # Should raise ExternalToolError
            tidehunter.run_analysis(
                input_file=input_file,
                output_file=output_file
            )


class TestTideHunterRunner:
    """Test cases for TideHunterRunner."""

    @patch.object(TideHunter, '_check_installation')
    def test_tidehunter_runner_initialization(self, mock_check):
        """Test TideHunterRunner initialization."""
        runner = TideHunterRunner(num_threads=4)
        assert runner.num_threads == 4
        assert runner.threads == 4

    @patch.object(TideHunter, '_check_installation')
    def test_tidehunter_runner_default_threads(self, mock_check):
        """Test TideHunterRunner default threads."""
        runner = TideHunterRunner()
        assert runner.num_threads == 8  # Default value
        assert runner.threads == 8

    @patch.object(TideHunter, '_check_installation')
    @patch('circleseeker.external.tidehunter.TideHunter.run_analysis')
    def test_tidehunter_runner_run(self, mock_run_analysis, mock_check, tmp_path):
        """Test TideHunterRunner run method."""
        input_fasta = tmp_path / "input.fasta"
        output_path = tmp_path / "output.txt"

        runner = TideHunterRunner(num_threads=2)
        runner.run(input_fasta, output_path)

        # Verify run_analysis was called with correct parameters
        mock_run_analysis.assert_called_once_with(
            input_file=Path(input_fasta),
            output_file=Path(output_path)
        )

    @patch.object(TideHunter, '_check_installation')
    @patch('circleseeker.external.tidehunter.TideHunter.run_analysis')
    def test_tidehunter_runner_run_string_paths(self, mock_run_analysis, mock_check):
        """Test TideHunterRunner run method with string paths."""
        input_fasta = "/path/to/input.fasta"
        output_path = "/path/to/output.txt"

        runner = TideHunterRunner()
        runner.run(input_fasta, output_path)

        # Verify paths were converted to Path objects
        mock_run_analysis.assert_called_once_with(
            input_file=Path(input_fasta),
            output_file=Path(output_path)
        )


@pytest.fixture
def sample_fasta_file(tmp_path):
    """Create a sample FASTA file for testing."""
    fasta_file = tmp_path / "sample.fasta"
    fasta_content = """>seq1
ATCGATCGATCGATCGATCGATCGATCGATCGATCG
>seq2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
>seq3
TTTTAAAAGGGGCCCCTTTTAAAAGGGGCCCCTTTT"""
    fasta_file.write_text(fasta_content)
    return fasta_file
