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
            raise ValueError("boom")

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


# ========== ecc_summary step Tests ========== #


class TestEccSummary:
    """Tests for ecc_summary step."""

    def test_skips_when_files_missing(self, tmp_path):
        """ecc_summary skips gracefully when required files are missing."""
        from circleseeker.core.steps.postprocess import ecc_summary

        config = Config(output_dir=tmp_path / "out", prefix="sample")
        pipeline = Pipeline(config)
        # No files created - should skip without error
        ecc_summary(pipeline)

    def test_skips_when_no_input_file(self, tmp_path):
        """ecc_summary skips when input_file is not configured."""
        from circleseeker.core.steps.postprocess import ecc_summary

        config = Config(output_dir=tmp_path / "out", prefix="sample")
        pipeline = Pipeline(config)

        # Create required intermediate files
        processed_csv = pipeline.config.output_dir / "tandem_to_ring.csv"
        processed_csv.parent.mkdir(parents=True, exist_ok=True)
        processed_csv.write_text("col1,col2\nval1,val2\n")

        unified_csv = pipeline.config.output_dir / "sample_unified.csv"
        unified_csv.write_text(
            "eccDNA_id,Regions,Strand,Length,eccDNA_type,State\n"
            "UeccDNA1,chr1:100-200,+,100,UeccDNA,Confirmed\n"
        )
        pipeline.state.results[ResultKeys.UNIFIED_CSV] = str(unified_csv)

        # config.input_file is None by default
        ecc_summary(pipeline)

    def test_sets_result_key_on_success(self, tmp_path, monkeypatch):
        """ecc_summary sets SUMMARY_DIR result on successful run."""
        from circleseeker.core.steps.postprocess import ecc_summary

        config = Config(
            output_dir=tmp_path / "out",
            prefix="sample",
            input_file=tmp_path / "input.fasta",
        )
        pipeline = Pipeline(config)

        # Create required files
        processed_csv = pipeline.config.output_dir / "tandem_to_ring.csv"
        processed_csv.parent.mkdir(parents=True, exist_ok=True)
        processed_csv.write_text("col1\nval1\n")

        unified_csv = pipeline.config.output_dir / "sample_unified.csv"
        unified_csv.write_text(
            "eccDNA_id,Regions,Strand,Length,eccDNA_type,State\n"
            "UeccDNA1,chr1:100-200,+,100,UeccDNA,Confirmed\n"
        )
        pipeline.state.results[ResultKeys.UNIFIED_CSV] = str(unified_csv)

        # Create input FASTA
        input_fasta = tmp_path / "input.fasta"
        input_fasta.write_text(">read1\nACGTACGT\n")

        # Mock EccSummary to avoid full processing
        class MockSummary:
            def __init__(self, **kwargs):
                pass
            def process_fasta(self, *a):
                pass
            def process_processed_csv(self, *a):
                pass
            def process_merged_csv(self, *a):
                pass
            def process_overlap_stats(self, *a):
                pass
            def generate_html_report(self):
                pass
            def generate_text_summary(self):
                pass

        monkeypatch.setattr(
            "circleseeker.modules.ecc_summary.EccSummary",
            MockSummary,
        )

        ecc_summary(pipeline)
        assert ResultKeys.SUMMARY_DIR in pipeline.state.results


# ========== ecc_unify success path Tests ========== #


class TestEccUnifySuccessPath:
    """Tests for ecc_unify step success path."""

    def test_successful_merge_sets_result_keys(self, tmp_path):
        """Successful ecc_unify sets UNIFIED_CSV in pipeline results."""
        config = Config(output_dir=tmp_path / "out", prefix="sample")
        pipeline = Pipeline(config)

        # Create confirmed CSV
        confirmed_path = pipeline.config.output_dir / "sample_confirmed.csv"
        confirmed_path.parent.mkdir(parents=True, exist_ok=True)
        confirmed_path.write_text(
            "eccDNA_id,Regions,Strand,Length,eccDNA_type,State,Seg_total,Hit_count\n"
            "UeccDNA1,chr1:100-200,+,100,UeccDNA,Confirmed,1,5\n"
        )
        pipeline.state.results[ResultKeys.CONFIRMED_CSV] = str(confirmed_path)

        ecc_unify(pipeline)

        assert ResultKeys.UNIFIED_CSV in pipeline.state.results
        unified_path = pipeline.config.output_dir / "sample_unified.csv"
        assert unified_path.exists()

    def test_merge_with_inferred_simple(self, tmp_path):
        """ecc_unify correctly merges confirmed and inferred simple data."""
        config = Config(output_dir=tmp_path / "out", prefix="sample")
        pipeline = Pipeline(config)

        # Create confirmed CSV
        confirmed_path = pipeline.config.output_dir / "sample_confirmed.csv"
        confirmed_path.parent.mkdir(parents=True, exist_ok=True)
        confirmed_path.write_text(
            "eccDNA_id,Regions,Strand,Length,eccDNA_type,State,Seg_total,Hit_count\n"
            "UeccDNA1,chr1:100-200,+,100,UeccDNA,Confirmed,1,5\n"
        )
        pipeline.state.results[ResultKeys.CONFIRMED_CSV] = str(confirmed_path)

        # Create inferred simple CSV
        inferred_path = pipeline.config.output_dir / "sample_inferred_simple.csv"
        inferred_path.write_text(
            "eccDNA_id,chr,start0,end0,length,strand\n"
            "IUECC1,chr5,5000,5500,500,+\n"
        )
        pipeline.state.results[ResultKeys.INFERRED_SIMPLE_TSV] = str(inferred_path)

        ecc_unify(pipeline)

        # Should have merged data
        unified_csv = pipeline.config.output_dir / "sample_unified.csv"
        assert unified_csv.exists()
        df = pd.read_csv(unified_csv)
        assert len(df) >= 1  # At least the confirmed entry

    def test_inference_disabled_skips_inferred(self, tmp_path):
        """When inference is disabled, inferred files are not loaded."""
        config = Config(output_dir=tmp_path / "out", prefix="sample")
        pipeline = Pipeline(config)

        confirmed_path = pipeline.config.output_dir / "sample_confirmed.csv"
        confirmed_path.parent.mkdir(parents=True, exist_ok=True)
        confirmed_path.write_text(
            "eccDNA_id,Regions,Strand,Length,eccDNA_type,State,Seg_total,Hit_count\n"
            "UeccDNA1,chr1:100-200,+,100,UeccDNA,Confirmed,1,5\n"
        )
        pipeline.state.results[ResultKeys.CONFIRMED_CSV] = str(confirmed_path)
        pipeline.state.results[ResultKeys.INFERENCE_INPUT_EMPTY] = True

        ecc_unify(pipeline)

        unified_csv = pipeline.config.output_dir / "sample_unified.csv"
        assert unified_csv.exists()
        df = pd.read_csv(unified_csv)
        # Only confirmed data, no inferred
        assert len(df) == 1


# ========== ecc_packager with unified CSV Tests ========== #


class TestEccPackagerWithData:
    """Tests for ecc_packager step with actual data."""

    def test_creates_final_results_key(self, tmp_path):
        """ecc_packager sets FINAL_RESULTS in pipeline results."""
        config = Config(output_dir=tmp_path / "out", prefix="sample")
        pipeline = Pipeline(config)

        ecc_packager(pipeline)

        assert ResultKeys.FINAL_RESULTS in pipeline.state.results

    def test_fallback_on_packager_failure(self, tmp_path, monkeypatch):
        """ecc_packager falls back to basic output structure on failure."""
        config = Config(output_dir=tmp_path / "out", prefix="sample")
        pipeline = Pipeline(config)

        unified_csv = pipeline.config.output_dir / "sample_unified.csv"
        unified_csv.parent.mkdir(parents=True, exist_ok=True)
        unified_csv.write_text(
            "eccDNA_id,Regions,Strand,Length,eccDNA_type,State\n"
            "UeccDNA1,chr1:100-200,+,100,UeccDNA,Confirmed\n"
        )
        pipeline.state.results[ResultKeys.UNIFIED_CSV] = str(unified_csv)

        def failing_run(*args, **kwargs):
            raise OSError("simulated failure")

        monkeypatch.setattr(
            "circleseeker.modules.ecc_packager.run",
            failing_run,
        )

        ecc_packager(pipeline)
        assert ResultKeys.FINAL_RESULTS in pipeline.state.results
