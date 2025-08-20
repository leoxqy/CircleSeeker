# tests/test_cli.py
"""Tests for the command-line interface."""

from pathlib import Path
from click.testing import CliRunner
from circleseeker2.cli import cli
from circleseeker2 import __version__

def test_cli_version():
    """Test the --version flag."""
    runner = CliRunner()
    result = runner.invoke(cli, ["--version"])
    assert result.exit_code == 0
    assert __version__ in result.output

def test_cli_run_command(temp_dir: Path):
    """Test the 'run' command with basic arguments."""
    runner = CliRunner()
    input_file = temp_dir / "input.fa"
    ref_file = temp_dir / "ref.fa"
    output_dir = temp_dir / "output"
    input_file.touch()
    ref_file.touch()

    args = [
        "run",
        str(input_file),
        str(ref_file),
        "-o", str(output_dir),
        "-p", "test_run",
        "-t", "4",
        "--dry-run" # Use dry-run to avoid running the full pipeline
    ]

    result = runner.invoke(cli, args)
    
    assert result.exit_code == 0
    assert output_dir.exists()
    
    # Check that the effective config was saved
    config_path = output_dir / "config.yaml"
    assert config_path.exists()

    with open(config_path, 'r') as f:
        import yaml
        data = yaml.safe_load(f)
        assert data["prefix"] == "test_run"
        assert data["performance"]["threads"] == 4
        assert Path(data["input_file"]).name == "input.fa"
