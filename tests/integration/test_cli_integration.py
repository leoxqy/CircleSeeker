from __future__ import annotations

from click.testing import CliRunner

import pytest

from circleseeker.cli import cli
from circleseeker import __version__


@pytest.mark.integration
def test_cli_version() -> None:
    runner = CliRunner()
    result = runner.invoke(cli, ["--version"])
    assert result.exit_code == 0
    assert result.output == f"CircleSeeker {__version__}\n"


@pytest.mark.integration
def test_cli_help() -> None:
    runner = CliRunner()
    result = runner.invoke(cli, ["-h"])
    assert result.exit_code == 0
    assert "CircleSeeker: Comprehensive eccDNA detection" in result.output


@pytest.mark.integration
def test_cli_missing_required_inputs_shows_help() -> None:
    runner = CliRunner()
    result = runner.invoke(cli, [])
    assert result.exit_code == 1
    assert "Error: Both --input and --reference are required" in result.output
    assert "Usage:" in result.output


@pytest.mark.integration
def test_cli_advanced_option_requires_debug() -> None:
    runner = CliRunner()
    result = runner.invoke(cli, ["--start-from", "1"])
    assert result.exit_code == 2
    assert "require --debug" in result.output
    assert "--start-from" in result.output


@pytest.mark.integration
def test_cli_generate_config_requires_debug() -> None:
    runner = CliRunner()
    result = runner.invoke(cli, ["--generate-config"])
    assert result.exit_code == 2
    assert "--generate-config" in result.output


@pytest.mark.integration
def test_cli_generate_config_with_debug() -> None:
    runner = CliRunner()
    result = runner.invoke(cli, ["--debug", "--generate-config"])
    assert result.exit_code == 0
    assert "CircleSeeker Configuration File" in result.output
    assert "input_file:" in result.output


@pytest.mark.integration
def test_cli_show_steps_with_debug() -> None:
    runner = CliRunner()
    result = runner.invoke(cli, ["--debug", "--show-steps"])
    assert result.exit_code == 0
    assert "CircleSeeker Pipeline Steps:" in result.output


@pytest.mark.integration
def test_cli_run_show_steps() -> None:
    runner = CliRunner()
    result = runner.invoke(cli, ["run", "--show-steps"])
    assert result.exit_code == 0
    assert "CircleSeeker Pipeline Steps:" in result.output


@pytest.mark.integration
def test_cli_init_config_writes_file() -> None:
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(cli, ["--debug", "init-config", "-o", "config.yaml"])
        assert result.exit_code == 0

        contents = open("config.yaml", encoding="utf-8").read()
        assert "CircleSeeker Configuration File" in contents
