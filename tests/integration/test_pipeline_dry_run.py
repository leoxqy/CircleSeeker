"""Integration tests for the `--dry-run` execution mode."""

from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner

import pytest

from circleseeker.cli import cli
from circleseeker.core.pipeline import Pipeline

pytestmark = pytest.mark.integration


def _write_fasta(path: Path, *, name: str, sequence: str) -> None:
    path.write_text(f">{name}\n{sequence}\n", encoding="utf-8")


def test_run_dry_run_prints_plan_without_creating_output_dir() -> None:
    runner = CliRunner()
    with runner.isolated_filesystem():
        reads = Path("reads.fa")
        ref = Path("ref.fa")
        _write_fasta(reads, name="read1", sequence="ACGT" * 10)
        _write_fasta(ref, name="chr1", sequence="A" * 200)

        out_dir = Path("out")
        result = runner.invoke(
            cli,
            [
                "run",
                "-i",
                str(reads),
                "-r",
                str(ref),
                "-o",
                str(out_dir),
                "-t",
                "4",
                "--dry-run",
            ],
        )
        assert result.exit_code == 0, result.output
        assert "CircleSeeker Pipeline Steps:" in result.output
        assert f"Total: {len(Pipeline.STEPS)} steps" in result.output
        assert "check_dependencies" in result.output
        assert "tidehunter" in result.output
        assert f"Would process: {reads}" in result.output
        assert f"With reference: {ref}" in result.output
        assert str(out_dir.absolute()) in result.output
        assert "Using 4 threads" in result.output
        assert not out_dir.exists()

