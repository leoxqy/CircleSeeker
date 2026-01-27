"""Tests for pipeline_types module."""

from pathlib import Path
import hashlib
import sys
import time

import pytest

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.core.pipeline_types import (
    ResultKeys,
    PipelineStep,
    StepMetadata,
    PipelineState,
)


# ---------------------------------------------------------------------------
# TestResultKeys
# ---------------------------------------------------------------------------


class TestResultKeys:
    """Tests for the ResultKeys constants class."""

    def test_key_values_are_strings(self):
        """Every public class attribute on ResultKeys must be a string."""
        for attr in dir(ResultKeys):
            if attr.startswith("_"):
                continue
            value = getattr(ResultKeys, attr)
            assert isinstance(value, str), f"ResultKeys.{attr} is not a string"

    def test_no_duplicate_values(self):
        """All ResultKeys values must be unique to avoid ambiguous look-ups."""
        values = [
            getattr(ResultKeys, attr)
            for attr in dir(ResultKeys)
            if not attr.startswith("_")
        ]
        assert len(values) == len(set(values)), "Duplicate values found in ResultKeys"

    def test_known_keys_exist(self):
        """Spot-check a selection of keys that the pipeline relies on."""
        assert ResultKeys.T2R_CSV == "tandem_to_ring_csv"
        assert ResultKeys.T2R_FASTA == "tandem_to_ring_fasta"
        assert ResultKeys.ALIGNMENT_OUTPUT == "alignment_output"
        assert ResultKeys.MINIMAP2_BAM == "minimap2_bam"
        assert ResultKeys.FINAL_RESULTS == "final_results"
        assert ResultKeys.TIDEHUNTER_OUTPUT == "tidehunter_output"
        assert ResultKeys.UNIFIED_CSV == "unified_csv"


# ---------------------------------------------------------------------------
# TestPipelineStep
# ---------------------------------------------------------------------------


class TestPipelineStep:
    """Tests for the PipelineStep dataclass."""

    def test_default_values(self):
        """Only name and description are required; other fields have defaults."""
        step = PipelineStep(name="align", description="Run alignment")
        assert step.name == "align"
        assert step.description == "Run alignment"
        assert step.display_name is None
        assert step.required is True
        assert step.skip_condition is None
        assert step.depends_on == []

    def test_all_fields_set(self):
        """All fields can be provided explicitly."""
        step = PipelineStep(
            name="filter",
            description="Filter reads",
            display_name="Read Filter",
            required=False,
            skip_condition="no_reads",
            depends_on=["align", "index"],
        )
        assert step.name == "filter"
        assert step.description == "Filter reads"
        assert step.display_name == "Read Filter"
        assert step.required is False
        assert step.skip_condition == "no_reads"
        assert step.depends_on == ["align", "index"]

    def test_depends_on_default_is_independent_per_instance(self):
        """Each instance should get its own empty list (default_factory)."""
        step_a = PipelineStep(name="a", description="A")
        step_b = PipelineStep(name="b", description="B")
        step_a.depends_on.append("x")
        assert step_b.depends_on == [], "Mutable default leaked between instances"

    def test_required_default_true(self):
        """The required field defaults to True."""
        step = PipelineStep(name="s", description="d")
        assert step.required is True


# ---------------------------------------------------------------------------
# TestStepMetadata
# ---------------------------------------------------------------------------


class TestStepMetadata:
    """Tests for the StepMetadata dataclass."""

    def test_default_values(self):
        """Only step_name and start_time are required; others have defaults."""
        meta = StepMetadata(step_name="align", start_time=100.0)
        assert meta.step_name == "align"
        assert meta.start_time == 100.0
        assert meta.end_time is None
        assert meta.duration is None
        assert meta.status == "running"
        assert meta.error_message is None
        assert meta.input_files == []
        assert meta.output_files == []
        assert meta.file_checksums == {}

    def test_default_status_is_running(self):
        """Status defaults to 'running'."""
        meta = StepMetadata(step_name="x", start_time=0.0)
        assert meta.status == "running"

    def test_all_fields_set(self):
        """All fields can be explicitly provided."""
        meta = StepMetadata(
            step_name="merge",
            start_time=1000.0,
            end_time=1050.0,
            duration=50.0,
            status="completed",
            error_message=None,
            input_files=["/a.bam"],
            output_files=["/b.bam"],
            file_checksums={"/b.bam": "abc123"},
        )
        assert meta.step_name == "merge"
        assert meta.end_time == 1050.0
        assert meta.duration == 50.0
        assert meta.status == "completed"
        assert meta.input_files == ["/a.bam"]
        assert meta.output_files == ["/b.bam"]
        assert meta.file_checksums == {"/b.bam": "abc123"}


# ---------------------------------------------------------------------------
# TestPipelineStateInit
# ---------------------------------------------------------------------------


class TestPipelineStateInit:
    """Tests for PipelineState construction and defaults."""

    def test_default_construction(self):
        """Construct with only the required field (completed_steps)."""
        state = PipelineState(completed_steps=[])
        assert state.completed_steps == []
        assert state.current_step is None
        assert state.failed_step is None
        assert state.results == {}
        assert state.step_metadata == {}
        assert state.pipeline_start_time is None
        assert state.last_checkpoint_time is None
        assert state.config_hash is None

    def test_version_is_2_0(self):
        """The default version must be '2.0'."""
        state = PipelineState(completed_steps=[])
        assert state.version == "2.0"


# ---------------------------------------------------------------------------
# TestAddStepMetadata
# ---------------------------------------------------------------------------


class TestAddStepMetadata:
    """Tests for PipelineState.add_step_metadata."""

    def test_add_new_step(self):
        """Adding metadata for a step that does not yet exist creates it."""
        state = PipelineState(completed_steps=[])
        state.add_step_metadata("align")
        assert "align" in state.step_metadata
        meta = state.step_metadata["align"]
        assert meta.step_name == "align"
        assert meta.status == "running"
        assert meta.start_time > 0

    def test_add_with_input_files(self):
        """Input files are stored in the metadata."""
        state = PipelineState(completed_steps=[])
        state.add_step_metadata("align", input_files=["/reads.fq"])
        meta = state.step_metadata["align"]
        assert meta.input_files == ["/reads.fq"]

    def test_add_with_extra_kwargs(self):
        """Extra keyword arguments that match metadata attributes are set."""
        state = PipelineState(completed_steps=[])
        state.add_step_metadata("align", status="queued")
        meta = state.step_metadata["align"]
        assert meta.status == "queued"

    def test_update_existing_step(self):
        """Calling add_step_metadata again updates rather than replaces."""
        state = PipelineState(completed_steps=[])
        state.add_step_metadata("align")
        original_start = state.step_metadata["align"].start_time
        # Second call should NOT create a new StepMetadata
        state.add_step_metadata("align", input_files=["/new.fq"])
        meta = state.step_metadata["align"]
        assert meta.start_time == original_start, "start_time should not change"
        assert meta.input_files == ["/new.fq"]

    def test_ignored_kwargs_not_set(self):
        """Keywords that are not valid StepMetadata attributes are silently ignored."""
        state = PipelineState(completed_steps=[])
        state.add_step_metadata("align", nonexistent_field="ignored")
        meta = state.step_metadata["align"]
        assert not hasattr(meta, "nonexistent_field") or getattr(meta, "nonexistent_field", None) is None


# ---------------------------------------------------------------------------
# TestCompleteStep
# ---------------------------------------------------------------------------


class TestCompleteStep:
    """Tests for PipelineState.complete_step (without checksum logic)."""

    def _prepare(self) -> PipelineState:
        state = PipelineState(completed_steps=[])
        state.add_step_metadata("align")
        return state

    def test_basic_complete(self):
        """Completing a step sets status, end_time, and duration."""
        state = self._prepare()
        state.complete_step("align")
        meta = state.step_metadata["align"]
        assert meta.status == "completed"
        assert meta.end_time is not None
        assert meta.duration is not None
        assert meta.duration >= 0

    def test_complete_with_output_files(self):
        """Output files are stored on the metadata."""
        state = self._prepare()
        state.complete_step("align", output_files=["/out.bam"])
        meta = state.step_metadata["align"]
        assert meta.output_files == ["/out.bam"]

    def test_adds_to_completed_steps(self):
        """The step name is appended to completed_steps."""
        state = self._prepare()
        state.complete_step("align")
        assert "align" in state.completed_steps

    def test_no_duplicate_completed_steps(self):
        """Completing the same step twice does not produce duplicate entries."""
        state = self._prepare()
        state.complete_step("align")
        state.complete_step("align")
        assert state.completed_steps.count("align") == 1


# ---------------------------------------------------------------------------
# TestCompleteStepChecksums
# ---------------------------------------------------------------------------


class TestCompleteStepChecksums:
    """Tests for compute_checksums=True in complete_step."""

    def test_checksum_computed_for_real_file(self, tmp_path):
        """MD5 checksum is computed and stored for an existing file."""
        # Create a temporary file with known content
        content = b"hello pipeline"
        temp_file = tmp_path / "output.txt"
        temp_file.write_bytes(content)

        expected_md5 = hashlib.md5(content).hexdigest()

        state = PipelineState(completed_steps=[])
        state.add_step_metadata("align")
        state.complete_step(
            "align",
            output_files=[str(temp_file)],
            compute_checksums=True,
        )
        meta = state.step_metadata["align"]
        assert str(temp_file) in meta.file_checksums
        assert meta.file_checksums[str(temp_file)] == expected_md5

    def test_checksum_skips_missing_file(self, tmp_path):
        """A missing file does not cause an error; no checksum is stored."""
        missing = str(tmp_path / "does_not_exist.txt")

        state = PipelineState(completed_steps=[])
        state.add_step_metadata("align")
        state.complete_step(
            "align",
            output_files=[missing],
            compute_checksums=True,
        )
        meta = state.step_metadata["align"]
        assert missing not in meta.file_checksums


# ---------------------------------------------------------------------------
# TestFailStep
# ---------------------------------------------------------------------------


class TestFailStep:
    """Tests for PipelineState.fail_step."""

    def _prepare(self) -> PipelineState:
        state = PipelineState(completed_steps=[])
        state.add_step_metadata("align")
        return state

    def test_sets_status_failed(self):
        """Status is set to 'failed'."""
        state = self._prepare()
        state.fail_step("align", "disk full")
        assert state.step_metadata["align"].status == "failed"

    def test_sets_error_message(self):
        """Error message is stored in metadata."""
        state = self._prepare()
        state.fail_step("align", "disk full")
        assert state.step_metadata["align"].error_message == "disk full"

    def test_sets_failed_step_attribute(self):
        """PipelineState.failed_step is updated to the failed step name."""
        state = self._prepare()
        state.fail_step("align", "disk full")
        assert state.failed_step == "align"

    def test_fail_sets_end_time_and_duration(self):
        """Failing a step populates end_time and duration."""
        state = self._prepare()
        state.fail_step("align", "timeout")
        meta = state.step_metadata["align"]
        assert meta.end_time is not None
        assert meta.duration is not None
        assert meta.duration >= 0


# ---------------------------------------------------------------------------
# TestGetTotalRuntime
# ---------------------------------------------------------------------------


class TestGetTotalRuntime:
    """Tests for PipelineState.get_total_runtime."""

    def test_returns_none_without_start_time(self):
        """If pipeline_start_time is not set, returns None."""
        state = PipelineState(completed_steps=[])
        assert state.get_total_runtime() is None

    def test_returns_positive_float_with_start_time(self):
        """If pipeline_start_time is set, returns a non-negative float."""
        state = PipelineState(completed_steps=[], pipeline_start_time=time.time() - 1.0)
        runtime = state.get_total_runtime()
        assert runtime is not None
        assert isinstance(runtime, float)
        assert runtime >= 0.0


# ---------------------------------------------------------------------------
# TestGetStepSummary
# ---------------------------------------------------------------------------


class TestGetStepSummary:
    """Tests for PipelineState.get_step_summary."""

    def test_empty_summary(self):
        """No steps means an empty summary dict."""
        state = PipelineState(completed_steps=[])
        assert state.get_step_summary() == {}

    def test_summary_completed_step(self):
        """Summary for a completed step has correct keys and values."""
        state = PipelineState(completed_steps=[])
        state.add_step_metadata("align")
        state.complete_step("align")

        summary = state.get_step_summary()
        assert "align" in summary
        entry = summary["align"]
        assert entry["status"] == "completed"
        assert entry["duration"] is not None
        assert entry["duration"] >= 0
        assert entry["start_time"] is not None
        assert entry["error"] is None

    def test_summary_failed_step(self):
        """Summary for a failed step includes the error message."""
        state = PipelineState(completed_steps=[])
        state.add_step_metadata("index")
        state.fail_step("index", "corrupted reference")

        summary = state.get_step_summary()
        assert "index" in summary
        entry = summary["index"]
        assert entry["status"] == "failed"
        assert entry["error"] == "corrupted reference"
        assert entry["duration"] is not None

    def test_summary_start_time_is_iso_format(self):
        """The start_time in the summary is an ISO-formatted string."""
        state = PipelineState(completed_steps=[])
        state.add_step_metadata("align")

        summary = state.get_step_summary()
        start_time_str = summary["align"]["start_time"]
        assert isinstance(start_time_str, str)
        # ISO format contains 'T' separator between date and time
        assert "T" in start_time_str
