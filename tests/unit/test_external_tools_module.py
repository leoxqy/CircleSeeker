"""Unit tests for external tools module wrappers."""

from __future__ import annotations

import subprocess
from pathlib import Path
from unittest.mock import MagicMock, patch, call

import pytest

from circleseeker.modules.external_tools import (
    TideHunterModule,
    CDHitModule,
    Minimap2Module,
)
from circleseeker.modules.base import ModuleResult


class TestTideHunterModule:
    """Tests for TideHunterModule."""

    @pytest.fixture
    def module(self):
        """Create TideHunterModule instance."""
        return TideHunterModule()

    def test_init_sets_tool_name(self, module):
        """TideHunterModule should have correct tool name."""
        assert module.tool_name == "TideHunter"

    def test_init_default_parameters(self, module):
        """TideHunterModule should have correct default parameters."""
        assert module.period_range == "50,10000"
        assert module.consensus_threshold == 0.8
        assert module.format_type == 2

    def test_check_tool_availability_found(self, module):
        """check_tool_availability should return True when tool exists."""
        with patch("shutil.which", return_value="/usr/bin/TideHunter"):
            assert module.check_tool_availability() is True

    def test_check_tool_availability_not_found(self, module):
        """check_tool_availability should return False when tool missing."""
        with patch("shutil.which", return_value=None):
            assert module.check_tool_availability() is False

    def test_get_tool_version(self, module):
        """get_tool_version should return version string."""
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(stdout="TideHunter 1.5.0", returncode=0)
            version = module.get_tool_version()
            assert version == "TideHunter 1.5.0"

    def test_get_tool_version_error(self, module):
        """get_tool_version should return 'unknown' on error."""
        with patch("subprocess.run", side_effect=OSError):
            version = module.get_tool_version()
            assert version == "unknown"

    def test_validate_inputs_missing_input_file(self, module):
        """validate_inputs should raise ValueError when input_file missing."""
        with pytest.raises(ValueError, match="input_file is required"):
            module.validate_inputs()

    def test_validate_inputs_nonexistent_file(self, module, tmp_path):
        """validate_inputs should raise FileNotFoundError for missing file."""
        with pytest.raises(FileNotFoundError):
            module.validate_inputs(input_file=tmp_path / "nonexistent.fa")

    def test_validate_inputs_success(self, module, tmp_path):
        """validate_inputs should return True for valid inputs."""
        input_file = tmp_path / "test.fa"
        input_file.write_text(">seq1\nACGT\n")
        assert module.validate_inputs(input_file=input_file) is True

    def test_execute_tool_not_available(self, module, tmp_path):
        """execute should return failure when tool not available."""
        input_file = tmp_path / "test.fa"
        input_file.write_text(">seq1\nACGT\n")

        with patch.object(module, "check_tool_availability", return_value=False):
            result = module.execute(input_file=input_file)

        assert result.success is False
        assert "not found" in result.error_message

    def test_execute_success(self, module, tmp_path):
        """execute should return success on successful run."""
        input_file = tmp_path / "test.fa"
        input_file.write_text(">seq1\nACGT\n")
        output_file = tmp_path / "output.txt"

        with patch.object(module, "check_tool_availability", return_value=True), \
             patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0)
            output_file.write_text("# header\ndata line\n")

            result = module.execute(
                input_file=input_file,
                output_file=output_file,
                threads=4
            )

        assert result.success is True
        assert "ecc_candidates" in result.output_files

    def test_execute_builds_correct_command(self, module, tmp_path):
        """execute should build correct TideHunter command."""
        input_file = tmp_path / "test.fa"
        input_file.write_text(">seq1\nACGT\n")
        output_file = tmp_path / "output.txt"

        with patch.object(module, "check_tool_availability", return_value=True), \
             patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0)
            output_file.write_text("")

            module.execute(
                input_file=input_file,
                output_file=output_file,
                threads=8,
                period_range="100,5000"
            )

            cmd = mock_run.call_args[0][0]
            assert cmd[0] == "TideHunter"
            assert "-f" in cmd
            assert "-p" in cmd
            assert "100,5000" in cmd
            assert "-t" in cmd
            assert "8" in cmd

    def test_execute_counts_candidates(self, module, tmp_path):
        """execute should count candidate lines in output."""
        input_file = tmp_path / "test.fa"
        input_file.write_text(">seq1\nACGT\n")
        output_file = tmp_path / "output.txt"

        with patch.object(module, "check_tool_availability", return_value=True), \
             patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0)
            output_file.write_text("# comment\ncandidate1\ncandidate2\n")

            result = module.execute(
                input_file=input_file,
                output_file=output_file
            )

        assert result.metrics.get("num_candidates") == 2

    def test_execute_handles_subprocess_error(self, module, tmp_path):
        """execute should handle subprocess errors gracefully."""
        input_file = tmp_path / "test.fa"
        input_file.write_text(">seq1\nACGT\n")

        with patch.object(module, "check_tool_availability", return_value=True), \
             patch("subprocess.run") as mock_run:
            mock_run.side_effect = subprocess.CalledProcessError(
                1, "TideHunter", stderr="Error message"
            )

            result = module.execute(input_file=input_file)

        assert result.success is False
        assert "failed" in result.error_message.lower()


class TestCDHitModule:
    """Tests for CDHitModule."""

    @pytest.fixture
    def module(self):
        """Create CDHitModule instance."""
        return CDHitModule()

    def test_init_sets_tool_name(self, module):
        """CDHitModule should have correct tool name."""
        assert module.tool_name == "cd-hit-est"

    def test_init_default_parameters(self, module):
        """CDHitModule should have correct default parameters."""
        assert module.similarity == 0.99
        assert module.memory == 0

    def test_check_tool_availability_found(self, module):
        """check_tool_availability should return True when tool exists."""
        with patch("shutil.which", return_value="/usr/bin/cd-hit-est"):
            assert module.check_tool_availability() is True

    def test_check_tool_availability_not_found(self, module):
        """check_tool_availability should return False when tool missing."""
        with patch("shutil.which", return_value=None):
            assert module.check_tool_availability() is False

    def test_get_tool_version(self, module):
        """get_tool_version should extract version from stderr."""
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(
                stderr="CD-HIT version 4.8.1\nother info",
                returncode=1  # CD-HIT returns non-zero without args
            )
            version = module.get_tool_version()
            assert "CD-HIT" in version

    def test_validate_inputs_missing_input_file(self, module):
        """validate_inputs should raise ValueError when input_file missing."""
        with pytest.raises(ValueError, match="input_file is required"):
            module.validate_inputs()

    def test_execute_tool_not_available(self, module, tmp_path):
        """execute should return failure when tool not available."""
        input_file = tmp_path / "test.fa"
        input_file.write_text(">seq1\nACGT\n")

        with patch.object(module, "check_tool_availability", return_value=False):
            result = module.execute(input_file=input_file)

        assert result.success is False
        assert "not found" in result.error_message

    def test_execute_success(self, module, tmp_path):
        """execute should return success on successful run."""
        input_file = tmp_path / "test.fa"
        input_file.write_text(">seq1\nACGT\n")
        output_prefix = tmp_path / "output"
        output_clstr = Path(f"{output_prefix}.clstr")

        with patch.object(module, "check_tool_availability", return_value=True), \
             patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0)
            output_clstr.write_text(">Cluster 0\n0\t100bp\n>Cluster 1\n0\t50bp\n")

            result = module.execute(
                input_file=input_file,
                output_prefix=str(output_prefix)
            )

        assert result.success is True

    def test_execute_builds_correct_command(self, module, tmp_path):
        """execute should build correct CD-HIT command."""
        input_file = tmp_path / "test.fa"
        input_file.write_text(">seq1\nACGT\n")
        output_prefix = tmp_path / "output"

        with patch.object(module, "check_tool_availability", return_value=True), \
             patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0)
            Path(f"{output_prefix}.clstr").write_text("")

            module.execute(
                input_file=input_file,
                output_prefix=str(output_prefix),
                threads=4,
                similarity=0.95
            )

            cmd = mock_run.call_args[0][0]
            assert cmd[0] == "cd-hit-est"
            assert "-i" in cmd
            assert "-o" in cmd
            assert "-c" in cmd
            assert "0.95" in cmd
            assert "-T" in cmd
            assert "4" in cmd

    def test_execute_counts_clusters(self, module, tmp_path):
        """execute should count clusters in output."""
        input_file = tmp_path / "test.fa"
        input_file.write_text(">seq1\nACGT\n")
        output_prefix = tmp_path / "output"
        output_clstr = Path(f"{output_prefix}.clstr")

        with patch.object(module, "check_tool_availability", return_value=True), \
             patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0)
            output_clstr.write_text(">Cluster 0\n>Cluster 1\n>Cluster 2\n")

            result = module.execute(
                input_file=input_file,
                output_prefix=str(output_prefix)
            )

        assert result.metrics.get("num_clusters") == 3


class TestMinimap2Module:
    """Tests for Minimap2Module."""

    @pytest.fixture
    def module(self):
        """Create Minimap2Module instance."""
        return Minimap2Module()

    def test_init_sets_tool_name(self, module):
        """Minimap2Module should have correct tool name."""
        assert module.tool_name == "minimap2"

    def test_check_tool_availability_both_found(self, module):
        """check_tool_availability should return True when both tools exist."""
        def which_side_effect(cmd):
            return f"/usr/bin/{cmd}" if cmd in ["minimap2", "samtools"] else None

        with patch("shutil.which", side_effect=which_side_effect):
            assert module.check_tool_availability() is True

    def test_check_tool_availability_minimap2_missing(self, module):
        """check_tool_availability should return False when minimap2 missing."""
        def which_side_effect(cmd):
            return "/usr/bin/samtools" if cmd == "samtools" else None

        with patch("shutil.which", side_effect=which_side_effect):
            assert module.check_tool_availability() is False

    def test_check_tool_availability_samtools_missing(self, module):
        """check_tool_availability should return False when samtools missing."""
        def which_side_effect(cmd):
            return "/usr/bin/minimap2" if cmd == "minimap2" else None

        with patch("shutil.which", side_effect=which_side_effect):
            assert module.check_tool_availability() is False

    def test_get_tool_version(self, module):
        """get_tool_version should return version string."""
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(stdout="2.24-r1122", returncode=0)
            version = module.get_tool_version()
            assert version == "2.24-r1122"

    def test_validate_inputs_missing_reads_file(self, module, tmp_path):
        """validate_inputs should raise ValueError when reads_file missing."""
        ref_file = tmp_path / "ref.fa"
        ref_file.write_text(">chr1\nACGT\n")

        with pytest.raises(ValueError, match="reads_file is required"):
            module.validate_inputs(reference_file=ref_file)

    def test_validate_inputs_missing_reference_file(self, module, tmp_path):
        """validate_inputs should raise ValueError when reference_file missing."""
        reads_file = tmp_path / "reads.fa"
        reads_file.write_text(">read1\nACGT\n")

        with pytest.raises(ValueError, match="reference_file is required"):
            module.validate_inputs(reads_file=reads_file)

    def test_validate_inputs_success(self, module, tmp_path):
        """validate_inputs should return True for valid inputs."""
        reads_file = tmp_path / "reads.fa"
        reads_file.write_text(">read1\nACGT\n")
        ref_file = tmp_path / "ref.fa"
        ref_file.write_text(">chr1\nACGT\n")

        assert module.validate_inputs(reads_file=reads_file, reference_file=ref_file)

    def test_execute_tool_not_available(self, module, tmp_path):
        """execute should return failure when tools not available."""
        reads_file = tmp_path / "reads.fa"
        reads_file.write_text(">read1\nACGT\n")
        ref_file = tmp_path / "ref.fa"
        ref_file.write_text(">chr1\nACGT\n")

        with patch.object(module, "check_tool_availability", return_value=False):
            result = module.execute(reads_file=reads_file, reference_file=ref_file)

        assert result.success is False
        assert "not found" in result.error_message

    def test_execute_uses_popen_pipeline(self, module, tmp_path):
        """execute should use Popen pipeline for minimap2 | samtools."""
        reads_file = tmp_path / "reads.fa"
        reads_file.write_text(">read1\nACGT\n")
        ref_file = tmp_path / "ref.fa"
        ref_file.write_text(">chr1\nACGT\n")
        output_bam = tmp_path / "output.bam"

        with patch.object(module, "check_tool_availability", return_value=True), \
             patch("subprocess.Popen") as mock_popen, \
             patch("subprocess.run") as mock_run:

            # Setup mock processes
            mock_minimap2 = MagicMock()
            mock_minimap2.stdout = MagicMock()
            mock_minimap2.returncode = 0
            mock_minimap2.wait.return_value = 0

            mock_samtools = MagicMock()
            mock_samtools.communicate.return_value = (b"", b"")
            mock_samtools.returncode = 0

            mock_popen.side_effect = [mock_minimap2, mock_samtools]

            # Mock samtools index and flagstat
            mock_run.return_value = MagicMock(
                stdout="100 + 0 mapped (100%)",
                returncode=0
            )

            result = module.execute(
                reads_file=reads_file,
                reference_file=ref_file,
                output_bam=output_bam
            )

        # Verify Popen was called twice (minimap2 and samtools)
        assert mock_popen.call_count == 2

    def test_execute_handles_pipeline_error(self, module, tmp_path):
        """execute should handle pipeline errors gracefully."""
        reads_file = tmp_path / "reads.fa"
        reads_file.write_text(">read1\nACGT\n")
        ref_file = tmp_path / "ref.fa"
        ref_file.write_text(">chr1\nACGT\n")

        with patch.object(module, "check_tool_availability", return_value=True), \
             patch("subprocess.Popen") as mock_popen:

            mock_minimap2 = MagicMock()
            mock_minimap2.stdout = MagicMock()
            mock_minimap2.returncode = 1
            mock_minimap2.wait.return_value = 1

            mock_samtools = MagicMock()
            mock_samtools.communicate.return_value = (b"", b"Error")
            mock_samtools.returncode = 1

            mock_popen.side_effect = [mock_minimap2, mock_samtools]

            result = module.execute(
                reads_file=reads_file,
                reference_file=ref_file
            )

        assert result.success is False


class TestModuleResultIntegration:
    """Tests for ModuleResult usage in external tools."""

    def test_result_has_correct_module_name(self):
        """ModuleResult should include correct module name."""
        module = TideHunterModule()
        result = ModuleResult(success=True, module_name=module.name)
        assert result.module_name == module.name

    def test_result_add_output(self):
        """ModuleResult should store outputs correctly."""
        result = ModuleResult(success=True, module_name="test")
        result.add_output("key", Path("/path/to/file"))
        assert "key" in result.output_files
        assert result.output_files["key"] == Path("/path/to/file")

    def test_result_add_metric(self):
        """ModuleResult should store metrics correctly."""
        result = ModuleResult(success=True, module_name="test")
        result.add_metric("count", 42)
        assert "count" in result.metrics
        assert result.metrics["count"] == 42

    def test_result_add_warning(self):
        """ModuleResult should store warnings correctly."""
        result = ModuleResult(success=True, module_name="test")
        result.add_warning("Warning message")
        assert "Warning message" in result.warnings
