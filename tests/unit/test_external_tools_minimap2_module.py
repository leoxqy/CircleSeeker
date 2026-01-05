"""Tests for Minimap2Module in external_tools."""

from pathlib import Path
import sys
from unittest.mock import MagicMock, patch
import subprocess

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.modules.external_tools import Minimap2Module


class TestMinimap2Module:
    """Test cases for Minimap2Module."""

    @patch.object(Minimap2Module, "check_tool_availability", return_value=True)
    @patch("circleseeker.modules.external_tools.subprocess.run")
    @patch("circleseeker.modules.external_tools.subprocess.Popen")
    def test_execute_does_not_pipe_minimap2_stderr(
        self, mock_popen, mock_run, mock_check, tmp_path
    ):
        reads_file = tmp_path / "reads.fasta"
        reference_file = tmp_path / "reference.fasta"
        output_bam = tmp_path / "aligned.bam"

        reads_file.write_text(">r1\nACGT\n")
        reference_file.write_text(">ref\nACGT\n")

        minimap2_proc = MagicMock()
        minimap2_proc.returncode = 0
        minimap2_proc.wait.return_value = 0
        minimap2_proc.stdout = MagicMock()

        samtools_proc = MagicMock()
        samtools_proc.returncode = 0
        samtools_proc.communicate.return_value = (b"", b"")

        mock_popen.side_effect = [minimap2_proc, samtools_proc]

        def run_side_effect(cmd, *args, **kwargs):
            if cmd[:2] == ["samtools", "flagstat"]:
                return subprocess.CompletedProcess(cmd, 0, stdout="10 + 0 mapped", stderr="")
            return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

        mock_run.side_effect = run_side_effect

        module = Minimap2Module()
        result = module.execute(
            reads_file=reads_file,
            reference_file=reference_file,
            output_bam=output_bam,
            threads=2,
        )

        assert result.success is True
        assert mock_popen.call_count == 2

        first_kwargs = mock_popen.call_args_list[0].kwargs
        assert first_kwargs.get("stderr") is None
        assert first_kwargs.get("stdout") == subprocess.PIPE
