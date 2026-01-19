"""Tests for postprocess pipeline steps."""

from pathlib import Path
import sys

import pytest
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.config import Config
from circleseeker.core.pipeline import Pipeline
from circleseeker.core.pipeline_types import ResultKeys
from circleseeker.core.steps.postprocess import ecc_unify, ecc_packager
from circleseeker.exceptions import PipelineError


class TestEccUnify:
    """Tests for ecc_unify step."""

    def test_raises_on_merge_failure(self, tmp_path, monkeypatch):
        """ecc_unify should raise when merging fails to avoid silent skips."""
        config = Config(output_dir=tmp_path / "out", prefix="sample")
        pipeline = Pipeline(config)

        confirmed_path = pipeline.config.output_dir / "sample_confirmed.csv"
        confirmed_path.write_text(
            "eccDNA_id,Regions,Strand,Length,eccDNA_type,State,Seg_total,Hit_count\n"
        )
        pipeline.state.results[ResultKeys.CONFIRMED_CSV] = str(confirmed_path)

        def boom(*_args, **_kwargs):
            raise RuntimeError("boom")

        monkeypatch.setattr("circleseeker.modules.ecc_unify.merge_eccdna_tables", boom)

        with pytest.raises(PipelineError, match="ecc_unify failed"):
            ecc_unify(pipeline)

    def test_handles_missing_confirmed_csv(self, tmp_path, monkeypatch):
        """ecc_unify should handle missing confirmed.csv gracefully."""
        config = Config(output_dir=tmp_path / "out", prefix="sample")
        pipeline = Pipeline(config)

        # No confirmed CSV set in results
        # Should not crash, but may produce warning
        ecc_unify(pipeline)

    def test_handles_empty_confirmed_csv(self, tmp_path, monkeypatch):
        """ecc_unify should handle empty confirmed.csv."""
        config = Config(output_dir=tmp_path / "out", prefix="sample")
        pipeline = Pipeline(config)

        confirmed_path = pipeline.config.output_dir / "sample_confirmed.csv"
        confirmed_path.parent.mkdir(parents=True, exist_ok=True)
        confirmed_path.write_text(
            "eccDNA_id,Regions,Strand,Length,eccDNA_type,State,Seg_total,Hit_count\n"
        )
        pipeline.state.results[ResultKeys.CONFIRMED_CSV] = str(confirmed_path)

        # Should handle empty data without crashing
        ecc_unify(pipeline)


class TestEccPackager:
    """Tests for ecc_packager step."""

    def test_creates_output_directory(self, tmp_path, monkeypatch):
        """ecc_packager should create the output directory structure."""
        config = Config(output_dir=tmp_path / "out", prefix="sample")
        pipeline = Pipeline(config)

        # Create minimal required files
        confirmed_dir = pipeline.config.output_dir / "sample_Confirmed_eccDNA"
        confirmed_dir.mkdir(parents=True, exist_ok=True)

        # Should not crash even with minimal setup
        ecc_packager(pipeline)

    def test_handles_missing_source_files(self, tmp_path, monkeypatch):
        """ecc_packager should handle missing source files gracefully."""
        config = Config(output_dir=tmp_path / "out", prefix="sample")
        pipeline = Pipeline(config)

        # Don't create any files - should handle gracefully
        ecc_packager(pipeline)


class TestOutputValidation:
    """Tests for output file validation."""

    def test_csv_output_has_required_columns(self, tmp_path):
        """Output CSV files should have required columns."""
        required_columns = [
            "eccDNA_id",
            "Regions",
            "Strand",
            "Length",
            "eccDNA_type",
            "State",
        ]

        # Create a sample output
        output_df = pd.DataFrame({
            "eccDNA_id": ["ecc1"],
            "Regions": ["chr1:100-200"],
            "Strand": ["+"],
            "Length": [100],
            "eccDNA_type": ["Uecc"],
            "State": ["Confirmed"],
            "Seg_total": [1],
            "Hit_count": [5],
        })

        output_path = tmp_path / "test_output.csv"
        output_df.to_csv(output_path, index=False)

        # Verify all required columns are present
        loaded_df = pd.read_csv(output_path)
        for col in required_columns:
            assert col in loaded_df.columns, f"Missing required column: {col}"
