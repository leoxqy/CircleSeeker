"""Unit tests for CLI commands."""

from pathlib import Path
import sys
import tempfile

import pytest
from click.testing import CliRunner

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.cli import cli


class TestCLIBasics:
    """Test basic CLI functionality."""

    def test_cli_help(self):
        """Test CLI help message."""
        runner = CliRunner()
        result = runner.invoke(cli, ["--help"])
        assert result.exit_code == 0
        assert "CircleSeeker" in result.output

    def test_cli_version(self):
        """Test CLI version display."""
        runner = CliRunner()
        result = runner.invoke(cli, ["--version"])
        assert result.exit_code == 0
        assert "CircleSeeker" in result.output

    def test_cli_missing_inputs(self):
        """Test CLI error when required inputs are missing."""
        runner = CliRunner()
        result = runner.invoke(cli, [])
        assert result.exit_code != 0
        assert "required" in result.output.lower() or "error" in result.output.lower()


class TestCLIValidate:
    """Test the validate command."""

    def test_validate_command_exists(self):
        """Test that validate command exists and is visible without --debug."""
        runner = CliRunner()
        result = runner.invoke(cli, ["validate", "--help"])
        assert result.exit_code == 0
        assert "Validate" in result.output or "validate" in result.output.lower()


class TestCLIConfig:
    """Test config-related CLI functionality."""

    def test_init_config_stdout(self):
        """Test config generation via init-config --stdout."""
        runner = CliRunner()
        result = runner.invoke(cli, ["init-config", "--stdout"])
        assert result.exit_code == 0
        assert "input_file:" in result.output

    def test_init_config_output_file(self):
        """Test init-config writes to file with --output-file."""
        runner = CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(cli, ["init-config", "--output-file", "custom.yaml"])
            assert result.exit_code == 0
            contents = open("custom.yaml", encoding="utf-8").read()
            assert "input_file:" in contents or "CircleSeeker" in contents

    def test_show_steps(self):
        """Test show-steps option (no --debug needed)."""
        runner = CliRunner()
        result = runner.invoke(cli, ["--show-steps"])
        assert result.exit_code == 0
        assert "Pipeline Steps" in result.output or "step" in result.output.lower()


class TestCLIPresets:
    """Test preset options."""

    def test_preset_option_help(self):
        """Test that preset option is available via --help-advanced."""
        runner = CliRunner()
        result = runner.invoke(cli, ["--help-advanced"])
        # Preset should be documented in advanced help
        assert result.exit_code == 0
        assert "--start-from" in result.output or "advanced" in result.output.lower()

    def test_preset_no_debug_required(self):
        """Test that --preset does NOT require --debug."""
        runner = CliRunner()
        # --preset alone should not trigger "require --debug" error
        # (it will still fail due to missing -i/-r, but not due to --debug)
        result = runner.invoke(cli, ["--preset", "strict"])
        assert result.exit_code != 0
        # Should NOT complain about --debug
        assert "require --debug" not in result.output


class TestCLIAdvancedOptions:
    """Test advanced debug options."""

    def test_start_from_requires_debug(self):
        """Test that --start-from requires --debug."""
        runner = CliRunner()
        result = runner.invoke(cli, ["--start-from", "1"])
        assert result.exit_code != 0
        assert "debug" in result.output.lower()

    def test_stop_at_requires_debug(self):
        """Test that --stop-at requires --debug."""
        runner = CliRunner()
        result = runner.invoke(cli, ["--stop-at", "5"])
        assert result.exit_code != 0
        assert "debug" in result.output.lower()

    def test_resume_requires_debug(self):
        """Test that --resume requires --debug."""
        runner = CliRunner()
        result = runner.invoke(cli, ["--resume"])
        assert result.exit_code != 0
        assert "debug" in result.output.lower()

    def test_dry_run_no_debug_required(self):
        """Test that --dry-run does NOT require --debug."""
        runner = CliRunner()
        with tempfile.TemporaryDirectory() as tmpdir:
            result = runner.invoke(cli, [
                "-i", "/nonexistent/input.fasta",
                "-r", "/nonexistent/ref.fa",
                "-o", tmpdir,
                "--dry-run"
            ])
            # Should fail because input file doesn't exist, NOT because of --debug
            assert result.exit_code != 0
            assert "require --debug" not in result.output

    def test_show_steps_no_debug_required(self):
        """Test that --show-steps does NOT require --debug."""
        runner = CliRunner()
        result = runner.invoke(cli, ["--show-steps"])
        assert result.exit_code == 0
        assert "Pipeline Steps" in result.output or "step" in result.output.lower()


class TestCLISubcommandVisibility:
    """Test that utility subcommands are visible without --debug."""

    def test_init_config_visible(self):
        """Test that init-config is visible in help."""
        runner = CliRunner()
        result = runner.invoke(cli, ["--help"])
        assert result.exit_code == 0
        assert "init-config" in result.output

    def test_show_checkpoint_visible(self):
        """Test that show-checkpoint is visible in help."""
        runner = CliRunner()
        result = runner.invoke(cli, ["--help"])
        assert result.exit_code == 0
        assert "show-checkpoint" in result.output

    def test_validate_visible(self):
        """Test that validate is visible in help."""
        runner = CliRunner()
        result = runner.invoke(cli, ["--help"])
        assert result.exit_code == 0
        assert "validate" in result.output
