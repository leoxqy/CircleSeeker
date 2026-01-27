"""Tests for core/steps/umc.py step executors."""

from pathlib import Path
import sys

import pandas as pd
import pytest
from unittest.mock import MagicMock, patch, PropertyMock

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.core.steps.umc import (
    um_classify,
    cecc_build,
    umc_process,
    cd_hit,
    ecc_dedup,
    _get_path_from_results,
    _collect_ecc_dedup_inputs,
    _rename_dedup_outputs,
)
from circleseeker.core.pipeline_types import ResultKeys
from circleseeker.exceptions import PipelineError


# -------------------- Fixtures -------------------- #

@pytest.fixture
def mock_pipeline(tmp_path):
    """Create a minimal mock Pipeline object."""
    pipeline = MagicMock()
    pipeline.config.output_dir = tmp_path
    pipeline.config.prefix = "sample"
    pipeline.config.enable_xecc = True
    pipeline.config.threads = 4
    pipeline.config.keep_tmp = False
    pipeline.config.performance.threads = 4
    pipeline.logger = MagicMock()
    pipeline.state.results = {}

    def _set_result(key, value):
        pipeline.state.results[key] = value

    pipeline._set_result = MagicMock(side_effect=_set_result)

    def _serialize_path(p):
        return str(p)

    pipeline._serialize_path_for_state = _serialize_path
    return pipeline


@pytest.fixture
def alignment_tsv(tmp_path):
    """Create a fake alignment TSV file."""
    output = tmp_path / "sample_alignment_results.tsv"
    # query_id has 3 pipes → valid format: read|repN|length|copy
    rows = [
        "read1|rep1|1000|2\tchr1\t99.5\t100\t1000\t+\t0\t1000\t0\t1000\t1000\t60",
        "read1|rep1|1000|2\tchr2\t98.0\t100\t500\t+\t0\t500\t0\t500\t500\t30",
    ]
    output.write_text("\n".join(rows) + "\n")
    return output


# -------------------- _get_path_from_results -------------------- #


class TestGetPathFromResults:
    def test_key_exists_file_exists(self, tmp_path):
        pipeline = MagicMock()
        f = tmp_path / "test.csv"
        f.write_text("data")
        pipeline.state.results = {"key": str(f)}
        assert _get_path_from_results(pipeline, "key") == f

    def test_key_exists_file_missing(self, tmp_path):
        pipeline = MagicMock()
        pipeline.state.results = {"key": str(tmp_path / "missing.csv")}
        assert _get_path_from_results(pipeline, "key") is None

    def test_key_missing(self):
        pipeline = MagicMock()
        pipeline.state.results = {}
        assert _get_path_from_results(pipeline, "key") is None


# -------------------- _collect_ecc_dedup_inputs -------------------- #


class TestCollectEccDedupInputs:
    def test_collects_existing_files(self, mock_pipeline, tmp_path):
        # Create some files
        uecc_proc = tmp_path / "uecc_processed.csv"
        uecc_proc.write_text("data")
        mock_pipeline.state.results[ResultKeys.UECC_PROCESSED] = str(uecc_proc)

        result = _collect_ecc_dedup_inputs(mock_pipeline)
        assert result["uecc_input"] == uecc_proc
        assert result["mecc_input"] is None  # not set

    def test_empty_results(self, mock_pipeline):
        result = _collect_ecc_dedup_inputs(mock_pipeline)
        for v in result.values():
            assert v is None


# -------------------- _rename_dedup_outputs -------------------- #


class TestRenameDedupOutputs:
    def test_renames_existing_files(self, mock_pipeline, tmp_path):
        old_file = tmp_path / "sample_Mecc.fa"
        old_file.write_text("seq")
        _rename_dedup_outputs(mock_pipeline)
        assert not old_file.exists()
        assert (tmp_path / "sample_MeccDNA_C.fasta").exists()

    def test_skips_already_renamed(self, mock_pipeline, tmp_path):
        # Both old and new exist → should not overwrite
        old_file = tmp_path / "sample_Mecc.fa"
        old_file.write_text("old")
        new_file = tmp_path / "sample_MeccDNA_C.fasta"
        new_file.write_text("new")
        _rename_dedup_outputs(mock_pipeline)
        assert old_file.exists()
        assert new_file.read_text() == "new"

    def test_no_files_to_rename(self, mock_pipeline):
        # Should not raise
        _rename_dedup_outputs(mock_pipeline)


# -------------------- um_classify -------------------- #


class TestUmClassify:
    def test_missing_alignment_output_raises(self, mock_pipeline):
        mock_pipeline._resolve_stored_path = MagicMock(return_value=None)
        with pytest.raises(PipelineError, match="Alignment output not found"):
            um_classify(mock_pipeline)

    def test_empty_alignment_file(self, mock_pipeline, tmp_path):
        empty_file = tmp_path / "sample_alignment_results.tsv"
        empty_file.touch()
        mock_pipeline._resolve_stored_path = MagicMock(return_value=empty_file)
        mock_pipeline.config.tools.um_classify = {"gap_threshold": 10.0}

        um_classify(mock_pipeline)
        # Should set counts to 0
        assert mock_pipeline.state.results.get(ResultKeys.UECC_COUNT) == 0
        assert mock_pipeline.state.results.get(ResultKeys.MECC_COUNT) == 0

    def test_invalid_query_id_format_raises(self, mock_pipeline, tmp_path):
        bad_file = tmp_path / "sample_alignment_results.tsv"
        # Query IDs without 3 pipes
        bad_file.write_text("badid\tchr1\t99\t100\n")
        mock_pipeline._resolve_stored_path = MagicMock(return_value=bad_file)
        mock_pipeline.config.tools.um_classify = {"gap_threshold": 10.0}

        with pytest.raises(PipelineError, match="query_id format is invalid"):
            um_classify(mock_pipeline)

    def test_valid_alignment_calls_classifier(self, mock_pipeline, alignment_tsv):
        mock_pipeline._resolve_stored_path = MagicMock(return_value=alignment_tsv)
        mock_pipeline.config.tools.um_classify = {
            "gap_threshold": 10.0,
            "theta_full": 0.95,
            "theta_u": 0.95,
            "theta_m": 0.95,
        }

        with patch("circleseeker.modules.um_classify.UMeccClassifier") as MockClassifier:
            mock_instance = MockClassifier.return_value
            mock_instance.classify_alignment_results = MagicMock()
            um_classify(mock_pipeline)
            mock_instance.classify_alignment_results.assert_called_once()


# -------------------- cecc_build -------------------- #


class TestCeccBuild:
    def test_empty_input_creates_empty_output(self, mock_pipeline, tmp_path):
        # No unclassified CSV → should create empty output
        cecc_build(mock_pipeline)
        assert mock_pipeline.state.results.get(ResultKeys.CECC_BUILD_COUNT) == 0

    def test_empty_unclassified_csv(self, mock_pipeline, tmp_path):
        unclass = tmp_path / "um_classify.unclassified.csv"
        unclass.write_text("header\n")  # header only, no data rows
        cecc_build(mock_pipeline)
        assert mock_pipeline.state.results.get(ResultKeys.CECC_BUILD_COUNT) == 0

    def test_calls_builder_with_data(self, mock_pipeline, tmp_path):
        unclass = tmp_path / "um_classify.unclassified.csv"
        unclass.write_text("query_id,chr,start,end\nread1,chr1,100,200\n")

        with patch("circleseeker.modules.cecc_build.CeccBuild") as MockBuilder:
            mock_instance = MockBuilder.return_value
            mock_instance.tmp_dir = None
            mock_instance.keep_tmp = False
            mock_instance.threads = 4
            mock_instance.run_pipeline = MagicMock()

            # run_pipeline creates the output
            def _create_output(**kwargs):
                output_csv = kwargs.get("output_csv", tmp_path / "cecc_build.csv")
                pd.DataFrame({"query_id": ["read1"]}).to_csv(output_csv, index=False)

            mock_instance.run_pipeline.side_effect = _create_output

            mock_pipeline.config.tools.cecc_build = {
                "overlap_threshold": 0.95,
                "min_segments": 2,
            }
            cecc_build(mock_pipeline)
            mock_instance.run_pipeline.assert_called_once()


# -------------------- umc_process -------------------- #


class TestUmcProcess:
    def test_missing_fasta_skips(self, mock_pipeline):
        # No circular FASTA file → should skip
        umc_process(mock_pipeline)
        mock_pipeline.logger.warning.assert_called()

    def test_with_fasta_calls_processor(self, mock_pipeline, tmp_path):
        fasta = tmp_path / "tandem_to_ring.fasta"
        fasta.write_text(">seq1\nATCG\n")

        with patch("circleseeker.modules.umc_process.UMCProcess") as MockProc:
            mock_instance = MockProc.return_value
            mock_instance.run = MagicMock()
            umc_process(mock_pipeline)
            mock_instance.run.assert_called_once()


# -------------------- cd_hit -------------------- #


class TestCdHit:
    def test_missing_fasta_skips(self, mock_pipeline):
        cd_hit(mock_pipeline)
        # Should warn for each missing file
        assert mock_pipeline.logger.warning.call_count >= 1

    def test_with_fasta_runs_clustering(self, mock_pipeline, tmp_path):
        fasta = tmp_path / "sample_UeccDNA_pre.fasta"
        fasta.write_text(">seq1\nATCG\n")

        # Create expected output
        output_path = tmp_path / "sample_U"

        with patch("circleseeker.external.cd_hit.CDHitEst") as MockCdHit:
            mock_instance = MockCdHit.return_value

            def _create_cluster(input_f, output_f):
                output_f.write_text(">seq1\nATCG\n")
                clstr = output_f.with_suffix(".clstr")
                clstr.write_text(">Cluster 0\n0\t4nt, >seq1... *\n")
                return str(clstr)

            mock_instance.cluster_sequences = MagicMock(side_effect=_create_cluster)
            cd_hit(mock_pipeline)
            mock_instance.cluster_sequences.assert_called_once()


# -------------------- ecc_dedup -------------------- #


class TestEccDedup:
    def test_calls_harmonizer(self, mock_pipeline, tmp_path):
        with patch("circleseeker.modules.ecc_dedup.EccDedup") as MockDedup, \
             patch("circleseeker.modules.ecc_dedup.organize_umc_files") as MockOrg:
            mock_instance = MockDedup.return_value
            mock_instance.run_deduplication = MagicMock(return_value={"uecc": pd.DataFrame()})
            MockOrg.return_value = {}

            ecc_dedup(mock_pipeline)
            mock_instance.run_deduplication.assert_called_once()

    def test_creates_confirmed_csv(self, mock_pipeline, tmp_path):
        confirmed = tmp_path / "sample_eccDNA_Confirmed.csv"
        confirmed.write_text("eccDNA_id,type\nU1,Uecc\n")

        with patch("circleseeker.modules.ecc_dedup.EccDedup") as MockDedup, \
             patch("circleseeker.modules.ecc_dedup.organize_umc_files") as MockOrg:
            mock_instance = MockDedup.return_value
            mock_instance.run_deduplication = MagicMock(return_value={})
            MockOrg.return_value = {}

            ecc_dedup(mock_pipeline)
            assert ResultKeys.CONFIRMED_CSV in mock_pipeline.state.results
