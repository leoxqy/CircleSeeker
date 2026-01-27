"""Integration tests for checkpoint/resume/force behavior."""

from __future__ import annotations

import json
from pathlib import Path

from click.testing import CliRunner

import pytest

from circleseeker.cli import cli
from circleseeker.config import Config
from circleseeker.core.pipeline import Pipeline
from circleseeker.exceptions import PipelineError

pytestmark = pytest.mark.integration


def _write_fasta(path: Path, *, name: str, sequence: str) -> None:
    path.write_text(f">{name}\n{sequence}\n", encoding="utf-8")


def _make_config(tmp_path: Path, *, prefix: str = "sample") -> Config:
    reads = tmp_path / "reads.fa"
    reference = tmp_path / "ref.fa"
    _write_fasta(reads, name="read1", sequence="ACGT" * 10)
    _write_fasta(reference, name="chr1", sequence="A" * 200)

    cfg = Config()
    cfg.input_file = reads
    cfg.reference = reference
    cfg.output_dir = tmp_path / "out"
    cfg.prefix = prefix
    cfg.skip_organize = True  # keep checkpoint file for assertions
    return cfg


def _load_checkpoint(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_checkpoint_created_and_skips_completed_steps(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    cfg = _make_config(tmp_path, prefix="run1")
    pipeline = Pipeline(cfg)

    executed: list[str] = []

    def fake_execute_step(step) -> list[str]:
        executed.append(step.name)
        return []

    monkeypatch.setattr(pipeline, "_execute_step", fake_execute_step)
    pipeline.run(stop_at=2)

    assert pipeline.state_file.exists()
    checkpoint = _load_checkpoint(pipeline.state_file)
    assert checkpoint["version"] == "2.0"
    assert checkpoint["completed_steps"][:2] == ["check_dependencies", "tidehunter"]

    cfg2 = _make_config(tmp_path, prefix="run1")
    pipeline2 = Pipeline(cfg2)

    def should_not_run(_step) -> list[str]:
        raise AssertionError("Expected completed steps to be skipped when resuming")

    monkeypatch.setattr(pipeline2, "_execute_step", should_not_run)
    pipeline2.run(stop_at=2)


def test_force_reruns_completed_steps(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    cfg = _make_config(tmp_path, prefix="run2")
    pipeline = Pipeline(cfg)

    monkeypatch.setattr(pipeline, "_execute_step", lambda _step: [])
    pipeline.run(stop_at=2)

    cfg2 = _make_config(tmp_path, prefix="run2")
    pipeline2 = Pipeline(cfg2)

    executed: list[str] = []

    def fake_execute_step(step) -> list[str]:
        executed.append(step.name)
        return []

    monkeypatch.setattr(pipeline2, "_execute_step", fake_execute_step)
    pipeline2.run(stop_at=2, force=True)

    assert executed == ["check_dependencies", "tidehunter"]


def test_resume_after_failure_continues_from_failed_step(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cfg = _make_config(tmp_path, prefix="run3")
    pipeline = Pipeline(cfg)

    def fail_on_tidehunter(step) -> list[str]:
        if step.name == "tidehunter":
            raise RuntimeError("boom")
        return []

    monkeypatch.setattr(pipeline, "_execute_step", fail_on_tidehunter)
    with pytest.raises(PipelineError):
        pipeline.run(stop_at=2)

    checkpoint = _load_checkpoint(pipeline.state_file)
    assert checkpoint["completed_steps"] == ["check_dependencies"]
    assert checkpoint["failed_step"] == "tidehunter"

    cfg2 = _make_config(tmp_path, prefix="run3")
    pipeline2 = Pipeline(cfg2)

    executed: list[str] = []

    def resume_only_failed(step) -> list[str]:
        executed.append(step.name)
        if step.name == "check_dependencies":
            raise AssertionError("Completed step should be skipped on resume")
        return []

    monkeypatch.setattr(pipeline2, "_execute_step", resume_only_failed)
    pipeline2.run(stop_at=2)

    assert executed == ["tidehunter"]
    checkpoint2 = _load_checkpoint(pipeline2.state_file)
    assert checkpoint2["failed_step"] is None
    assert checkpoint2["completed_steps"][:2] == ["check_dependencies", "tidehunter"]


def test_cli_resume_ignores_force_flag(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    captured: dict[str, object] = {}

    def fake_run(
        self,
        start_from: int | None = None,
        stop_at: int | None = None,
        force: bool = False,
    ) -> dict:
        captured["force"] = force
        captured["start_from"] = start_from
        captured["stop_at"] = stop_at
        return {}

    monkeypatch.setattr(Pipeline, "run", fake_run)

    runner = CliRunner()
    with runner.isolated_filesystem():
        reads = Path("reads.fa")
        reference = Path("ref.fa")
        _write_fasta(reads, name="read1", sequence="ACGT" * 10)
        _write_fasta(reference, name="chr1", sequence="A" * 200)

        result = runner.invoke(
            cli,
            [
                "--debug",
                "-i",
                str(reads),
                "-r",
                str(reference),
                "-o",
                "out",
                "-t",
                "4",  # Use safe value for CI environments
                "--resume",
                "--force",
                "--start-from",
                "1",
                "--stop-at",
                "1",
            ],
        )

        assert result.exit_code == 0, result.output
        assert captured["force"] is False
        assert captured["start_from"] == 1
        assert captured["stop_at"] == 1
