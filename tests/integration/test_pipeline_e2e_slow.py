"""Slow end-to-end smoke tests that require external tools.

These tests are designed to be skipped by default in environments that do not
have the full CircleSeeker toolchain installed.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from circleseeker.config import Config
from circleseeker.core.pipeline import Pipeline

pytestmark = [pytest.mark.integration, pytest.mark.slow, pytest.mark.external]


def _has_inference_toolchain() -> bool:
    has_cresil = shutil.which("cresil") is not None
    has_cyrcular = shutil.which("cyrcular") is not None

    if has_cresil:
        return True

    if has_cyrcular:
        return shutil.which("bcftools") is not None and shutil.which("varlociraptor") is not None

    return False


def _has_required_tools() -> bool:
    has_tidehunter = shutil.which("TideHunter") is not None or shutil.which("tidehunter") is not None
    return (
        has_tidehunter
        and shutil.which("minimap2") is not None
        and shutil.which("samtools") is not None
        and shutil.which("cd-hit-est") is not None
        and _has_inference_toolchain()
    )


def _write_fasta(path: Path, *, name: str, sequence: str) -> None:
    path.write_text(f">{name}\n{sequence}\n", encoding="utf-8")


@pytest.mark.skipif(not _has_required_tools(), reason="External toolchain not available")
def test_pipeline_runs_through_tidehunter(tmp_path: Path) -> None:
    reads = tmp_path / "reads.fa"
    reference = tmp_path / "ref.fa"
    _write_fasta(reads, name="read1", sequence="ACGT" * 50)
    _write_fasta(reference, name="chr1", sequence="A" * 500)

    cfg = Config()
    cfg.input_file = reads
    cfg.reference = reference
    cfg.output_dir = tmp_path / "out"
    cfg.prefix = "e2e"
    cfg.skip_organize = True

    pipeline = Pipeline(cfg)
    pipeline.run(stop_at=2)

    # Checkpoint should be written in the final output dir.
    assert pipeline.state_file.exists()

    # TideHunter output should exist in the temp directory.
    tidehunter_output = pipeline.temp_dir / f"{cfg.prefix}.TH.ecc_candidates.txt"
    assert tidehunter_output.exists()
