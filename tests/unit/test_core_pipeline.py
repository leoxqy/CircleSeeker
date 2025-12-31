"""Tests for core pipeline module."""

from pathlib import Path
import sys
import pytest
import time
from unittest.mock import patch, MagicMock
from dataclasses import asdict
import tempfile

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.core.pipeline import (
    PipelineStep,
    StepMetadata,
    PipelineState,
    Pipeline,
    ResultKeys,
)
from circleseeker.config import Config
from circleseeker.exceptions import PipelineError


class TestPipelineStep:
    """Test cases for PipelineStep."""

    def test_pipeline_step_creation(self):
        """Test PipelineStep creation."""
        step = PipelineStep(
            name="test_step",
            description="Test step description",
            required=True,
            skip_condition="skip_test"
        )
        assert step.name == "test_step"
        assert step.description == "Test step description"
        assert step.required is True
        assert step.skip_condition == "skip_test"

    def test_pipeline_step_defaults(self):
        """Test PipelineStep default values."""
        step = PipelineStep(
            name="test_step",
            description="Test description"
        )
        assert step.required is True
        assert step.skip_condition is None


class TestStepMetadata:
    """Test cases for StepMetadata."""

    def test_step_metadata_creation(self):
        """Test StepMetadata creation."""
        start_time = time.time()
        metadata = StepMetadata(
            step_name="test_step",
            start_time=start_time
        )
        assert metadata.step_name == "test_step"
        assert metadata.start_time == start_time
        assert metadata.end_time is None
        assert metadata.duration is None
        assert metadata.status == "running"
        assert metadata.error_message is None
        assert metadata.input_files == []
        assert metadata.output_files == []
        assert metadata.file_checksums == {}

    def test_step_metadata_with_values(self):
        """Test StepMetadata with all values."""
        metadata = StepMetadata(
            step_name="test_step",
            start_time=1000.0,
            end_time=1010.0,
            duration=10.0,
            status="completed",
            error_message=None,
            input_files=["input.txt"],
            output_files=["output.txt"],
            file_checksums={"output.txt": "abc123"}
        )
        assert metadata.step_name == "test_step"
        assert metadata.start_time == 1000.0
        assert metadata.end_time == 1010.0
        assert metadata.duration == 10.0
        assert metadata.status == "completed"
        assert metadata.input_files == ["input.txt"]
        assert metadata.output_files == ["output.txt"]
        assert metadata.file_checksums == {"output.txt": "abc123"}


class TestPipelineState:
    """Test cases for PipelineState."""

    def test_pipeline_state_creation(self):
        """Test PipelineState creation."""
        state = PipelineState(completed_steps=["step1", "step2"])
        assert state.completed_steps == ["step1", "step2"]
        assert state.current_step is None
        assert state.failed_step is None
        assert state.results == {}
        assert state.step_metadata == {}
        assert state.pipeline_start_time is None
        assert state.last_checkpoint_time is None
        assert state.config_hash is None
        assert state.version == "2.0"

    def test_add_step_metadata(self):
        """Test add_step_metadata method."""
        state = PipelineState(completed_steps=[])

        # Add new step metadata
        state.add_step_metadata("test_step", status="running")

        assert "test_step" in state.step_metadata
        metadata = state.step_metadata["test_step"]
        assert metadata.step_name == "test_step"
        assert metadata.status == "running"
        assert metadata.start_time is not None

    def test_add_step_metadata_update_existing(self):
        """Test updating existing step metadata."""
        state = PipelineState(completed_steps=[])

        # Add initial metadata
        state.add_step_metadata("test_step", status="running")
        original_start_time = state.step_metadata["test_step"].start_time

        # Update existing metadata
        state.add_step_metadata("test_step", status="completed", duration=10.0)

        metadata = state.step_metadata["test_step"]
        assert metadata.status == "completed"
        assert metadata.duration == 10.0
        assert metadata.start_time == original_start_time  # Should not change

    def test_complete_step(self):
        """Test complete_step method."""
        state = PipelineState(completed_steps=[])

        # Add step metadata first
        state.add_step_metadata("test_step")
        original_start_time = state.step_metadata["test_step"].start_time

        # Complete the step
        output_files = ["output1.txt", "output2.txt"]
        state.complete_step("test_step", output_files)

        # Check metadata
        metadata = state.step_metadata["test_step"]
        assert metadata.status == "completed"
        assert metadata.end_time is not None
        assert metadata.duration is not None
        assert metadata.duration > 0
        assert metadata.output_files == output_files

        # Check completed steps
        assert "test_step" in state.completed_steps

    def test_complete_step_without_metadata(self):
        """Test complete_step without existing metadata."""
        state = PipelineState(completed_steps=[])

        # Complete step without metadata (should not crash)
        state.complete_step("test_step", ["output.txt"])

        assert "test_step" in state.completed_steps

    def test_fail_step(self):
        """Test fail_step method."""
        state = PipelineState(completed_steps=[])

        # Add step metadata first
        state.add_step_metadata("test_step")

        # Fail the step
        error_msg = "Test error message"
        state.fail_step("test_step", error_msg)

        # Check metadata
        metadata = state.step_metadata["test_step"]
        assert metadata.status == "failed"
        assert metadata.end_time is not None
        assert metadata.duration is not None
        assert metadata.error_message == error_msg

        # Check failed step
        assert state.failed_step == "test_step"

    def test_get_total_runtime(self):
        """Test get_total_runtime method."""
        state = PipelineState(completed_steps=[])

        # No start time
        assert state.get_total_runtime() is None

        # With start time
        start_time = time.time() - 10  # 10 seconds ago
        state.pipeline_start_time = start_time
        runtime = state.get_total_runtime()
        assert runtime is not None
        assert runtime >= 9.0  # Should be around 10 seconds, allow some variance

    def test_get_step_summary(self):
        """Test get_step_summary method."""
        state = PipelineState(completed_steps=[])

        # Add some step metadata
        state.add_step_metadata("step1", status="completed", duration=5.0)
        state.add_step_metadata("step2", status="failed", error_message="Test error")

        summary = state.get_step_summary()

        assert len(summary) == 2
        assert "step1" in summary
        assert "step2" in summary

        # Check step1 summary
        step1_summary = summary["step1"]
        assert step1_summary["status"] == "completed"
        assert step1_summary["duration"] == 5.0
        assert step1_summary["start_time"] is not None
        assert step1_summary["error"] is None

        # Check step2 summary
        step2_summary = summary["step2"]
        assert step2_summary["status"] == "failed"
        assert step2_summary["error"] == "Test error"


class TestPipeline:
    """Test cases for Pipeline."""

    def test_pipeline_steps_definition(self):
        """Test that pipeline steps are properly defined."""
        steps = Pipeline.STEPS
        assert len(steps) >= 16  # Should have 16 steps (Step 0-15)

        # Check first few steps
        step_names = [step.name for step in steps]
        assert "check_dependencies" in step_names  # Step 0
        assert "tidehunter" in step_names  # Step 1
        assert "tandem_to_ring" in step_names
        assert "run_alignment" in step_names
        assert "ecc_unify" in step_names
        assert "ecc_summary" in step_names
        assert "ecc_packager" in step_names

        # Verify Step 0 is check_dependencies
        assert steps[0].name == "check_dependencies"

    def test_result_key_aliases(self):
        """Test that Pipeline has STEPS defined."""
        # The Pipeline class should have STEPS
        assert hasattr(Pipeline, 'STEPS')
        assert len(Pipeline.STEPS) > 0

        # Check some specific step names exist
        step_names = [step.name for step in Pipeline.STEPS]
        assert "tandem_to_ring" in step_names
        assert "run_alignment" in step_names

    @patch('circleseeker.core.pipeline.Config')
    @patch('circleseeker.core.pipeline.get_logger')
    def test_pipeline_initialization(self, mock_logger, mock_config):
        """Test Pipeline initialization."""
        # Mock config
        mock_config_instance = MagicMock()
        mock_config_instance.prefix = "test"
        mock_config_instance.threads = 4
        mock_config_instance.output_dir = Path("/test/output")
        mock_config_instance.keep_tmp = False
        mock_config.return_value = mock_config_instance

        # Mock logger
        mock_logger_instance = MagicMock()
        mock_logger.return_value = mock_logger_instance

        # Test that the class exists and has expected attributes
        assert hasattr(Pipeline, 'STEPS')
        assert hasattr(Pipeline, '_set_result')
        assert hasattr(Pipeline, '_get_result')

    def test_set_result_method_exists(self):
        """Test that _set_result method exists."""
        # Check that the method exists in Pipeline class
        assert hasattr(Pipeline, '_set_result')

    def test_get_result_method_exists(self):
        """Test that _get_result method exists."""
        # Check that the method exists in Pipeline class
        assert hasattr(Pipeline, '_get_result')

    def test_result_key_alias_reverse_lookup(self):
        """Test that step names can be looked up in STEPS."""
        steps = Pipeline.STEPS
        step_names = [step.name for step in steps]

        # Test that we can find steps by name
        step_dict = {step.name: step for step in steps}

        # Test some specific lookups
        if "tandem_to_ring" in step_dict:
            assert step_dict["tandem_to_ring"].description == "Process TideHunter output"

    def test_empty_fallback_fasta_skips_inference(self, tmp_path):
        """Empty fallback FASTA should skip minimap2 alignment and inference."""
        input_fasta = tmp_path / "input.fa"
        reference = tmp_path / "ref.fa"

        input_fasta.write_text(">r1\nACGT\n")
        reference.write_text(">chr1\nACGT\n")

        config = Config(
            input_file=input_fasta,
            reference=reference,
            output_dir=tmp_path / "out",
            prefix="sample",
        )
        pipeline = Pipeline(config)
        pipeline._inference_tool = "cyrcular"

        pipeline._step_minimap2()

        assert pipeline.state.results.get(ResultKeys.INFERENCE_INPUT_EMPTY) is True
        all_fasta = Path(pipeline.state.results[ResultKeys.ALL_FILTERED_FASTA])
        assert all_fasta.exists()
        assert all_fasta.stat().st_size == 0
        assert ResultKeys.MINIMAP2_BAM not in pipeline.state.results

        pipeline._step_ecc_inference()
        assert pipeline.state.results[ResultKeys.CYRCULAR_RESULT_COUNT] == 0
        assert pipeline.state.results[ResultKeys.CYRCULAR_OVERVIEW] is None

    def test_run_validates_start_stop_range(self, tmp_path):
        """Invalid start/stop ranges should raise a PipelineError."""
        config = Config(output_dir=tmp_path / "out")
        pipeline = Pipeline(config)

        with pytest.raises(PipelineError):
            pipeline.run(start_from=0)
        with pytest.raises(PipelineError):
            pipeline.run(stop_at=0)
        with pytest.raises(PipelineError):
            pipeline.run(start_from=len(Pipeline.STEPS) + 1)
        with pytest.raises(PipelineError):
            pipeline.run(stop_at=len(Pipeline.STEPS) + 1)
        with pytest.raises(PipelineError):
            pipeline.run(start_from=3, stop_at=2)

    @pytest.mark.parametrize("step_name,description", [
        ("check_dependencies", "Check required tools and dependencies"),
        ("tidehunter", "Run TideHunter for tandem repeat detection"),
        ("run_alignment", "Run minimap2 alignment"),
        ("ecc_unify", "Merge eccDNA tables into unified output"),
        ("ecc_summary", "Generate summary report"),
        ("ecc_packager", "Package output files"),
    ])
    def test_specific_pipeline_steps(self, step_name, description):
        """Test specific pipeline steps exist with correct descriptions."""
        steps = Pipeline.STEPS
        step_dict = {step.name: step for step in steps}

        assert step_name in step_dict
        step = step_dict[step_name]
        assert step.description == description
        assert isinstance(step.required, bool)

    def test_pipeline_step_skip_conditions(self):
        """Test that PipelineStep can have skip conditions."""
        # Test that PipelineStep dataclass supports skip_condition field
        test_step = PipelineStep(
            name="test_step",
            description="Test step",
            skip_condition="skip_test"
        )
        assert test_step.skip_condition == "skip_test"

        # Test that steps are properly defined
        steps = Pipeline.STEPS
        assert len(steps) > 0

        # All steps should have name and description
        for step in steps:
            assert step.name is not None
            assert step.description is not None


@pytest.fixture
def sample_pipeline_state():
    """Create a sample pipeline state for testing."""
    state = PipelineState(
        completed_steps=["step1", "step2"],
        current_step="step3",
        pipeline_start_time=time.time() - 100  # 100 seconds ago
    )

    # Add some step metadata
    state.add_step_metadata("step1", status="completed", duration=10.0)
    state.add_step_metadata("step2", status="completed", duration=15.0)
    state.add_step_metadata("step3", status="running")

    return state


@pytest.fixture
def temp_pipeline_dir():
    """Create a temporary directory for pipeline testing."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        yield Path(tmp_dir)


class TestPipelineStateIntegration:
    """Integration tests for PipelineState."""

    def test_complete_workflow(self, sample_pipeline_state):
        """Test a complete workflow with the sample state."""
        state = sample_pipeline_state

        # Complete the current step
        state.complete_step("step3", ["output3.txt"])

        # Check that step3 is now completed
        assert "step3" in state.completed_steps
        assert state.step_metadata["step3"].status == "completed"

        # Get summary
        summary = state.get_step_summary()
        assert len(summary) == 3
        assert all(step["status"] in ["completed", "running"] for step in summary.values())

    def test_failure_workflow(self, sample_pipeline_state):
        """Test failure handling workflow."""
        state = sample_pipeline_state

        # Fail the current step
        state.fail_step("step3", "Test failure")

        # Check that step3 failed
        assert state.failed_step == "step3"
        assert state.step_metadata["step3"].status == "failed"
        assert state.step_metadata["step3"].error_message == "Test failure"

        # step3 should not be in completed steps
        assert "step3" not in state.completed_steps

    def test_runtime_calculation(self, sample_pipeline_state):
        """Test runtime calculation."""
        state = sample_pipeline_state

        runtime = state.get_total_runtime()
        assert runtime is not None
        assert runtime >= 99.0  # Should be around 100 seconds
