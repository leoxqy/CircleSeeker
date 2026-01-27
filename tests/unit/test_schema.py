"""Tests for schema validation module."""

from pathlib import Path
from unittest.mock import MagicMock
import sys

import pytest

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.core.steps.schema import (
    _format_template,
    _read_header,
    _validate_required_columns,
    _validate_tsv_no_header,
    resolve_artifact_path,
    validate_step_contract,
)
from circleseeker.core.steps.contracts import ArtifactSpec, StepContract
from circleseeker.exceptions import PipelineError


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _mock_pipeline(output_dir, final_dir, prefix="sample", **config_attrs):
    """Create a mock Pipeline with the given directory paths and config attributes."""
    pipeline = MagicMock()
    pipeline.config.output_dir = output_dir
    pipeline.config.prefix = prefix
    pipeline.final_output_dir = final_dir
    for attr, val in config_attrs.items():
        setattr(pipeline.config, attr, val)
    return pipeline


# ===========================================================================
# _format_template
# ===========================================================================

class TestFormatTemplate:
    """Tests for _format_template."""

    def test_basic_prefix_substitution(self):
        result = _format_template("{prefix}_output.txt", prefix="sample")
        assert result == "sample_output.txt"

    def test_no_placeholder_passes_through(self):
        result = _format_template("static_name.csv", prefix="sample")
        assert result == "static_name.csv"

    def test_missing_key_returns_template_unchanged(self):
        # Template references {other} which is not provided -- triggers KeyError
        template = "{other}_output.txt"
        result = _format_template(template, prefix="sample")
        assert result == template


# ===========================================================================
# _read_header
# ===========================================================================

class TestReadHeader:
    """Tests for _read_header."""

    def test_read_csv_header(self, tmp_path):
        csv_file = tmp_path / "data.csv"
        csv_file.write_text("col_a,col_b,col_c\n1,2,3\n", encoding="utf-8")
        header = _read_header(csv_file, delimiter=",")
        assert header == ["col_a", "col_b", "col_c"]

    def test_read_tsv_header(self, tmp_path):
        tsv_file = tmp_path / "data.tsv"
        tsv_file.write_text("col_x\tcol_y\tcol_z\nval1\tval2\tval3\n", encoding="utf-8")
        header = _read_header(tsv_file, delimiter="\t")
        assert header == ["col_x", "col_y", "col_z"]

    def test_empty_file_returns_empty_list(self, tmp_path):
        empty_file = tmp_path / "empty.csv"
        empty_file.write_text("", encoding="utf-8")
        header = _read_header(empty_file, delimiter=",")
        assert header == []


# ===========================================================================
# _validate_required_columns
# ===========================================================================

class TestValidateRequiredColumns:
    """Tests for _validate_required_columns."""

    def test_all_columns_present(self, tmp_path):
        csv_file = tmp_path / "good.csv"
        csv_file.write_text("chr,start0,end0\nchr1,100,200\n", encoding="utf-8")
        spec = ArtifactSpec(
            name="test",
            base="output",
            template="good.csv",
            kind="csv",
            required_columns=("chr", "start0", "end0"),
        )
        # Should not raise
        _validate_required_columns(csv_file, spec)

    def test_missing_columns_raises_pipeline_error(self, tmp_path):
        csv_file = tmp_path / "bad.csv"
        csv_file.write_text("chr,start0\n", encoding="utf-8")
        spec = ArtifactSpec(
            name="test",
            base="output",
            template="bad.csv",
            kind="csv",
            required_columns=("chr", "start0", "end0", "strand"),
        )
        with pytest.raises(PipelineError, match="missing columns"):
            _validate_required_columns(csv_file, spec)

    def test_missing_columns_error_contains_column_names(self, tmp_path):
        csv_file = tmp_path / "partial.csv"
        csv_file.write_text("chr,start0\n", encoding="utf-8")
        spec = ArtifactSpec(
            name="test",
            base="output",
            template="partial.csv",
            kind="csv",
            required_columns=("chr", "start0", "end0"),
        )
        with pytest.raises(PipelineError, match="end0"):
            _validate_required_columns(csv_file, spec)

    def test_no_required_columns_passes(self, tmp_path):
        csv_file = tmp_path / "any.csv"
        csv_file.write_text("whatever,columns\n", encoding="utf-8")
        spec = ArtifactSpec(
            name="test",
            base="output",
            template="any.csv",
            kind="csv",
            required_columns=(),
        )
        # Should not raise
        _validate_required_columns(csv_file, spec)

    def test_empty_file_passes(self, tmp_path):
        csv_file = tmp_path / "empty.csv"
        csv_file.write_text("", encoding="utf-8")
        spec = ArtifactSpec(
            name="test",
            base="output",
            template="empty.csv",
            kind="csv",
            required_columns=("chr", "start0"),
        )
        # Empty header -> returns early without error
        _validate_required_columns(csv_file, spec)


# ===========================================================================
# _validate_tsv_no_header
# ===========================================================================

class TestValidateTsvNoHeader:
    """Tests for _validate_tsv_no_header."""

    def test_valid_tsv_passes(self, tmp_path):
        tsv_file = tmp_path / "valid.tsv"
        fields = "\t".join([f"f{i}" for i in range(12)])
        tsv_file.write_text(fields + "\n", encoding="utf-8")
        # Should not raise
        _validate_tsv_no_header(tsv_file)

    def test_too_few_fields_raises_pipeline_error(self, tmp_path):
        tsv_file = tmp_path / "bad.tsv"
        # Only 3 fields, far fewer than the default min_fields=12
        tsv_file.write_text("a\tb\tc\n", encoding="utf-8")
        with pytest.raises(PipelineError, match="Schema validation failed"):
            _validate_tsv_no_header(tsv_file)

    def test_empty_file_passes(self, tmp_path):
        tsv_file = tmp_path / "empty.tsv"
        tsv_file.write_text("", encoding="utf-8")
        # Should not raise
        _validate_tsv_no_header(tsv_file)

    def test_custom_min_fields(self, tmp_path):
        tsv_file = tmp_path / "custom.tsv"
        tsv_file.write_text("a\tb\tc\n", encoding="utf-8")
        # min_fields=3 -> 3 fields is enough
        _validate_tsv_no_header(tsv_file, min_fields=3)
        # min_fields=4 -> 3 fields is not enough
        with pytest.raises(PipelineError):
            _validate_tsv_no_header(tsv_file, min_fields=4)

    def test_multiple_error_lines_truncated(self, tmp_path):
        tsv_file = tmp_path / "many_errors.tsv"
        # Write 10 lines, each with only 2 fields (less than default min_fields=12)
        lines = ["x\ty\n" for _ in range(10)]
        tsv_file.write_text("".join(lines), encoding="utf-8")
        with pytest.raises(PipelineError, match=r"and 5 more errors"):
            _validate_tsv_no_header(tsv_file)


# ===========================================================================
# resolve_artifact_path
# ===========================================================================

class TestResolveArtifactPath:
    """Tests for resolve_artifact_path."""

    def test_output_base_resolves_to_output_dir(self, tmp_path):
        out_dir = tmp_path / "output"
        fin_dir = tmp_path / "final"
        pipeline = _mock_pipeline(out_dir, fin_dir, prefix="demo")
        spec = ArtifactSpec(
            name="result",
            base="output",
            template="{prefix}_result.csv",
        )
        result = resolve_artifact_path(pipeline, spec)
        assert result == Path(out_dir) / "demo_result.csv"

    def test_final_base_resolves_to_final_dir(self, tmp_path):
        out_dir = tmp_path / "output"
        fin_dir = tmp_path / "final"
        pipeline = _mock_pipeline(out_dir, fin_dir, prefix="demo")
        spec = ArtifactSpec(
            name="packaged",
            base="final",
            template="{prefix}_packaged.csv",
        )
        result = resolve_artifact_path(pipeline, spec)
        assert result == Path(fin_dir) / "demo_packaged.csv"

    def test_config_base_reads_attribute(self, tmp_path):
        out_dir = tmp_path / "output"
        fin_dir = tmp_path / "final"
        ref_path = tmp_path / "ref.fasta"
        pipeline = _mock_pipeline(out_dir, fin_dir, reference=str(ref_path))
        spec = ArtifactSpec(
            name="reference",
            base="config",
            template="reference",
        )
        result = resolve_artifact_path(pipeline, spec)
        assert result == Path(str(ref_path))

    def test_config_base_returns_none_for_missing_attribute(self, tmp_path):
        out_dir = tmp_path / "output"
        fin_dir = tmp_path / "final"
        pipeline = _mock_pipeline(out_dir, fin_dir)
        # Ensure the attribute truly does not exist on the mock
        pipeline.config.configure_mock(**{})
        # getattr with default None for a missing attribute on MagicMock
        # returns a new MagicMock, so we explicitly set it to None
        pipeline.config.nonexistent_attr = None
        spec = ArtifactSpec(
            name="missing",
            base="config",
            template="nonexistent_attr",
        )
        result = resolve_artifact_path(pipeline, spec)
        assert result is None

    def test_unknown_base_returns_formatted_path(self, tmp_path):
        out_dir = tmp_path / "output"
        fin_dir = tmp_path / "final"
        pipeline = _mock_pipeline(out_dir, fin_dir, prefix="sample")
        spec = ArtifactSpec(
            name="arbitrary",
            base="custom",
            template="{prefix}_custom.txt",
        )
        result = resolve_artifact_path(pipeline, spec)
        assert result == Path("sample_custom.txt")


# ===========================================================================
# validate_step_contract
# ===========================================================================

class TestValidateStepContract:
    """Tests for validate_step_contract."""

    def test_pre_validates_inputs(self, tmp_path):
        out_dir = tmp_path / "output"
        fin_dir = tmp_path / "final"
        out_dir.mkdir()
        pipeline = _mock_pipeline(out_dir, fin_dir, prefix="s")

        # Create the input file so validation passes
        input_file = out_dir / "s_input.csv"
        input_file.write_text("chr,start0\n1,100\n", encoding="utf-8")

        contract = StepContract(
            step_name="test_step",
            inputs=(
                ArtifactSpec(
                    name="input_csv",
                    base="output",
                    template="{prefix}_input.csv",
                    kind="csv",
                    required=True,
                    required_columns=("chr", "start0"),
                ),
            ),
        )
        # Should not raise
        validate_step_contract(pipeline, contract, when="pre")

    def test_post_validates_outputs(self, tmp_path):
        out_dir = tmp_path / "output"
        fin_dir = tmp_path / "final"
        out_dir.mkdir()
        pipeline = _mock_pipeline(out_dir, fin_dir, prefix="s")

        # Required output file is missing -> should raise
        contract = StepContract(
            step_name="test_step",
            outputs=(
                ArtifactSpec(
                    name="output_csv",
                    base="output",
                    template="{prefix}_output.csv",
                    kind="csv",
                    required=True,
                ),
            ),
        )
        with pytest.raises(PipelineError, match="missing file"):
            validate_step_contract(pipeline, contract, when="post")

    def test_invalid_when_raises_value_error(self, tmp_path):
        out_dir = tmp_path / "output"
        fin_dir = tmp_path / "final"
        pipeline = _mock_pipeline(out_dir, fin_dir)
        contract = StepContract(step_name="test_step")
        with pytest.raises(ValueError, match="Invalid when"):
            validate_step_contract(pipeline, contract, when="during")
