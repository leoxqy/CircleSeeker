"""Integration tests for config/CLI/default precedence."""

from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner

import pytest

from circleseeker.cli import cli

pytestmark = pytest.mark.integration


def _write_fasta(path: Path, *, name: str, sequence: str) -> None:
    path.write_text(f">{name}\n{sequence}\n", encoding="utf-8")


def test_run_uses_config_when_cli_not_provided() -> None:
    runner = CliRunner()
    with runner.isolated_filesystem():
        reads = Path("reads.fa")
        ref = Path("ref.fa")
        _write_fasta(reads, name="read1", sequence="ACGT" * 10)
        _write_fasta(ref, name="chr1", sequence="A" * 200)

        cfg_out = Path("cfg_out")
        config_path = Path("config.yaml")
        config_path.write_text(
            "\n".join(
                [
                    f"input_file: {reads}",
                    f"reference: {ref}",
                    f"output_dir: {cfg_out}",
                    "threads: 3",
                    "",
                ]
            ),
            encoding="utf-8",
        )

        result = runner.invoke(cli, ["run", "-c", str(config_path), "--dry-run"])
        assert result.exit_code == 0, result.output
        assert "CircleSeeker Pipeline Steps:" in result.output
        assert f"Would process: {reads}" in result.output
        assert f"With reference: {ref}" in result.output
        assert str(cfg_out.absolute()) in result.output
        assert "Using 3 threads" in result.output
        assert not cfg_out.exists()


def test_cli_overrides_config() -> None:
    runner = CliRunner()
    with runner.isolated_filesystem():
        reads_cfg = Path("reads_cfg.fa")
        ref_cfg = Path("ref_cfg.fa")
        _write_fasta(reads_cfg, name="read_cfg", sequence="ACGT" * 10)
        _write_fasta(ref_cfg, name="chr_cfg", sequence="A" * 200)

        reads_cli = Path("reads_cli.fa")
        ref_cli = Path("ref_cli.fa")
        _write_fasta(reads_cli, name="read_cli", sequence="TGCA" * 10)
        _write_fasta(ref_cli, name="chr_cli", sequence="C" * 200)

        cfg_out = Path("cfg_out")
        config_path = Path("config.yaml")
        config_path.write_text(
            "\n".join(
                [
                    f"input_file: {reads_cfg}",
                    f"reference: {ref_cfg}",
                    f"output_dir: {cfg_out}",
                    "threads: 3",
                    "",
                ]
            ),
            encoding="utf-8",
        )

        cli_out = Path("cli_out")
        result = runner.invoke(
            cli,
            [
                "run",
                "-c",
                str(config_path),
                "-i",
                str(reads_cli),
                "-r",
                str(ref_cli),
                "-o",
                str(cli_out),
                "-t",
                "7",
                "--dry-run",
            ],
        )
        assert result.exit_code == 0, result.output
        assert f"Would process: {reads_cli}" in result.output
        assert f"With reference: {ref_cli}" in result.output
        assert str(cli_out.absolute()) in result.output
        assert "Using 7 threads" in result.output
        assert not cfg_out.exists()
        assert not cli_out.exists()


def test_defaults_apply_without_config() -> None:
    runner = CliRunner()
    with runner.isolated_filesystem():
        reads = Path("reads.fa")
        ref = Path("ref.fa")
        _write_fasta(reads, name="read1", sequence="ACGT" * 10)
        _write_fasta(ref, name="chr1", sequence="A" * 200)

        default_out = Path("circleseeker_output")
        result = runner.invoke(cli, ["run", "-i", str(reads), "-r", str(ref), "--dry-run"])
        assert result.exit_code == 0, result.output
        assert str(default_out.absolute()) in result.output
        assert "Using 8 threads" in result.output
        assert not default_out.exists()


def test_config_missing_required_fields_errors() -> None:
    runner = CliRunner()
    with runner.isolated_filesystem():
        config_path = Path("config.yaml")
        config_path.write_text("threads: 3\n", encoding="utf-8")

        result = runner.invoke(cli, ["run", "-c", str(config_path), "--dry-run"])
        assert result.exit_code == 1
        assert "Error: Both --input and --reference are required" in result.output
