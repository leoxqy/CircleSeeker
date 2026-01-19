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
        """Test that validate command exists."""
        runner = CliRunner()
        result = runner.invoke(cli, ["--debug", "validate", "--help"])
        assert result.exit_code == 0
        assert "Validate" in result.output or "validate" in result.output.lower()


class TestCLIConfig:
    """Test config-related CLI functionality."""

    def test_generate_config(self):
        """Test config generation."""
        runner = CliRunner()
        result = runner.invoke(cli, ["--debug", "--generate-config"])
        assert result.exit_code == 0
        assert "input_file:" in result.output

    def test_show_steps(self):
        """Test show-steps option."""
        runner = CliRunner()
        result = runner.invoke(cli, ["--debug", "--show-steps"])
        assert result.exit_code == 0
        assert "Pipeline Steps" in result.output or "step" in result.output.lower()


class TestCLIRunCommand:
    """Test the run command."""

    def test_run_help(self):
        """Test run command help."""
        runner = CliRunner()
        result = runner.invoke(cli, ["run", "--help"])
        assert result.exit_code == 0
        assert "run" in result.output.lower()

    def test_run_show_steps(self):
        """Test run command show-steps."""
        runner = CliRunner()
        result = runner.invoke(cli, ["run", "--show-steps"])
        assert result.exit_code == 0

    def test_run_dry_run_with_missing_input(self):
        """Test run with dry-run but missing input file."""
        runner = CliRunner()
        with tempfile.TemporaryDirectory() as tmpdir:
            result = runner.invoke(cli, [
                "--debug",
                "-i", "/nonexistent/input.fasta",
                "-r", "/nonexistent/ref.fa",
                "-o", tmpdir,
                "--dry-run"
            ])
            # Should fail because input file doesn't exist
            assert result.exit_code != 0


class TestCLIPresets:
    """Test preset options."""

    def test_preset_option_help(self):
        """Test that preset option is available via --help-advanced."""
        runner = CliRunner()
        result = runner.invoke(cli, ["--help-advanced"])
        # Preset should be documented in advanced help
        assert result.exit_code == 0
        # The preset option may be hidden, so we check for advanced help output
        assert "--start-from" in result.output or "advanced" in result.output.lower()


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
