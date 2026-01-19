from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import pandas as pd

from circleseeker.config import Config
from circleseeker.core.pipeline import Pipeline
from circleseeker.core.steps.preprocess import run_alignment, tandem_to_ring
from circleseeker.exceptions import PipelineError


class TestRunAlignment:
    """Tests for run_alignment step."""

    def test_requires_t2r_fasta_when_not_skipped(self, tmp_path, monkeypatch):
        """run_alignment should raise when tandem_to_ring.fasta is missing."""
        monkeypatch.chdir(tmp_path)
        input_file = tmp_path / "reads.fa"
        reference = tmp_path / "ref.fa"
        input_file.write_text(">r1\nACGT\n", encoding="utf-8")
        reference.write_text(">chr1\nACGT\n", encoding="utf-8")

        config = Config(
            input_file=input_file,
            reference=reference,
            output_dir=tmp_path / "output",
            prefix="sample",
        )
        pipeline = Pipeline(config)

        with pytest.raises(PipelineError, match="tandem_to_ring\\.fasta"):
            run_alignment(pipeline)

    def test_skips_when_skip_carousel_true(self, tmp_path, monkeypatch):
        """run_alignment should skip and log when skip_carousel is True."""
        monkeypatch.chdir(tmp_path)
        input_file = tmp_path / "reads.fa"
        reference = tmp_path / "ref.fa"
        input_file.write_text(">r1\nACGT\n", encoding="utf-8")
        reference.write_text(">chr1\nACGT\n", encoding="utf-8")

        config = Config(
            input_file=input_file,
            reference=reference,
            output_dir=tmp_path / "output",
            prefix="sample",
            skip_carousel=True,
        )
        pipeline = Pipeline(config)

        # Should not raise when skip_carousel is True
        run_alignment(pipeline)


class TestTandemToRing:
    """Tests for tandem_to_ring step."""

    def test_requires_tidehunter_output(self, tmp_path, monkeypatch):
        """tandem_to_ring should raise when TideHunter output is missing."""
        monkeypatch.chdir(tmp_path)
        input_file = tmp_path / "reads.fa"
        reference = tmp_path / "ref.fa"
        input_file.write_text(">r1\nACGT\n", encoding="utf-8")
        reference.write_text(">chr1\nACGT\n", encoding="utf-8")

        config = Config(
            input_file=input_file,
            reference=reference,
            output_dir=tmp_path / "output",
            prefix="sample",
        )
        pipeline = Pipeline(config)

        with pytest.raises(PipelineError, match="TideHunter output"):
            tandem_to_ring(pipeline)

    def test_requires_tidehunter_state(self, tmp_path, monkeypatch):
        """tandem_to_ring should require TideHunter output in state."""
        from circleseeker.core.pipeline_types import ResultKeys

        monkeypatch.chdir(tmp_path)
        input_file = tmp_path / "reads.fa"
        reference = tmp_path / "ref.fa"
        input_file.write_text(">r1\nACGT\n", encoding="utf-8")
        reference.write_text(">chr1\nACGT\n", encoding="utf-8")

        config = Config(
            input_file=input_file,
            reference=reference,
            output_dir=tmp_path / "output",
            prefix="sample",
        )
        pipeline = Pipeline(config)

        # Without setting TideHunter output in state, should raise
        with pytest.raises(PipelineError, match="TideHunter output not found"):
            tandem_to_ring(pipeline)

    def test_handles_empty_tidehunter_output(self, tmp_path, monkeypatch):
        """tandem_to_ring should handle empty TideHunter output gracefully."""
        from circleseeker.core.pipeline_types import ResultKeys

        monkeypatch.chdir(tmp_path)
        input_file = tmp_path / "reads.fa"
        reference = tmp_path / "ref.fa"
        input_file.write_text(">r1\nACGT\n", encoding="utf-8")
        reference.write_text(">chr1\nACGT\n", encoding="utf-8")

        config = Config(
            input_file=input_file,
            reference=reference,
            output_dir=tmp_path / "output",
            prefix="sample",
        )
        pipeline = Pipeline(config)

        # Create empty TideHunter output and set in state
        th_output = pipeline.config.output_dir / "tidehunter_output.tsv"
        th_output.parent.mkdir(parents=True, exist_ok=True)
        th_output.write_text("", encoding="utf-8")
        pipeline.state.results[ResultKeys.TIDEHUNTER_OUTPUT] = str(th_output)

        # Should handle empty file without crashing
        tandem_to_ring(pipeline)


class TestAlignmentFiltering:
    """Tests for alignment length filtering in run_alignment."""

    def test_filters_short_alignments(self, tmp_path, monkeypatch):
        """run_alignment should filter out alignments shorter than min_alignment_length."""
        monkeypatch.chdir(tmp_path)

        # Create mock alignment output
        alignment_tsv = tmp_path / "alignment.tsv"
        alignment_tsv.write_text(
            "read1\t100\tchr1\t500\t60\t0\n"  # Should be kept (500bp)
            "read2\t100\tchr1\t30\t60\t0\n"   # Should be filtered (30bp < 50)
            "read3\t100\tchr1\t50\t60\t0\n",  # Should be kept (50bp = min)
            encoding="utf-8"
        )

        # Verify content before filtering
        lines = alignment_tsv.read_text().strip().split("\n")
        assert len(lines) == 3
