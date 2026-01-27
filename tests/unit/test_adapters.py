"""Tests for adapters module."""

from pathlib import Path
from unittest.mock import patch, MagicMock
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.modules.adapters import CLIModuleAdapter, TandemToRingAdapter, UMClassifyAdapter
from circleseeker.modules.base import ModuleResult


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def dummy_module(tmp_path):
    """Create a minimal Python file to act as a valid module_path."""
    script = tmp_path / "dummy_module.py"
    script.write_text("print('hello')\n")
    return script


@pytest.fixture
def adapter(dummy_module):
    """Return a CLIModuleAdapter pointed at the dummy module."""
    return CLIModuleAdapter(module_path=str(dummy_module), name="test_adapter")


# ---------------------------------------------------------------------------
# TestCLIModuleAdapterInit
# ---------------------------------------------------------------------------

class TestCLIModuleAdapterInit:
    """Tests for CLIModuleAdapter.__init__."""

    def test_init_with_existing_file(self, dummy_module):
        adapter = CLIModuleAdapter(module_path=str(dummy_module), name="my_adapter")
        assert adapter.module_path == dummy_module
        assert adapter.name == "my_adapter"

    def test_init_stores_path_as_path_object(self, dummy_module):
        adapter = CLIModuleAdapter(module_path=str(dummy_module))
        assert isinstance(adapter.module_path, Path)

    def test_init_nonexistent_file_raises(self, tmp_path):
        missing = tmp_path / "nonexistent.py"
        with pytest.raises(FileNotFoundError, match="Module not found"):
            CLIModuleAdapter(module_path=str(missing))


# ---------------------------------------------------------------------------
# TestValidateInputs
# ---------------------------------------------------------------------------

class TestValidateInputs:
    """Tests for CLIModuleAdapter.validate_inputs."""

    def test_validate_existing_input_file(self, adapter, tmp_path):
        real_file = tmp_path / "data.bam"
        real_file.write_text("fake bam content")
        result = adapter.validate_inputs(input_file=str(real_file))
        assert result is True

    def test_validate_existing_alignment_file(self, adapter, tmp_path):
        real_file = tmp_path / "alignment.paf"
        real_file.write_text("fake alignment")
        result = adapter.validate_inputs(alignment_file=str(real_file))
        assert result is True

    def test_validate_missing_input_file_raises(self, adapter, tmp_path):
        missing = tmp_path / "does_not_exist.bam"
        with pytest.raises(FileNotFoundError, match="Input file not found"):
            adapter.validate_inputs(input_file=str(missing))

    def test_validate_missing_fasta_file_raises(self, adapter, tmp_path):
        missing = tmp_path / "missing.fasta"
        with pytest.raises(FileNotFoundError, match="Input file not found"):
            adapter.validate_inputs(fasta_file=str(missing))

    def test_validate_no_relevant_keys_returns_true(self, adapter):
        """When kwargs contain no recognised input keys, validation passes."""
        result = adapter.validate_inputs(output_dir="/tmp", threads=4)
        assert result is True

    def test_validate_none_value_skipped(self, adapter):
        """A key present but with a falsy/None value should not trigger a check."""
        result = adapter.validate_inputs(input_file=None)
        assert result is True


# ---------------------------------------------------------------------------
# TestBuildCommand
# ---------------------------------------------------------------------------

class TestBuildCommand:
    """Tests for CLIModuleAdapter.build_command."""

    def test_basic_param_mapping(self, adapter, tmp_path):
        cmd = adapter.build_command(
            input_file="/data/in.bam",
            output_file="/data/out.csv",
            sample_name="sample1",
        )
        assert "-i" in cmd
        assert "/data/in.bam" in cmd
        assert "-o" in cmd
        assert "/data/out.csv" in cmd
        assert "-p" in cmd
        assert "sample1" in cmd

    def test_threads_mapping(self, adapter):
        cmd = adapter.build_command(threads=8)
        assert "--threads" in cmd
        assert "8" in cmd

    def test_debug_mapping(self, adapter):
        cmd = adapter.build_command(debug=True)
        assert "--debug" in cmd
        assert "True" in cmd

    def test_none_values_skipped(self, adapter):
        cmd = adapter.build_command(input_file=None, output_file=None)
        # Only the base command: [sys.executable, module_path]
        assert len(cmd) == 2

    def test_double_dash_params_passed_directly(self, adapter):
        cmd = adapter.build_command(**{"--min-length": "500", "--max-gap": "100"})
        assert "--min-length" in cmd
        assert "500" in cmd
        assert "--max-gap" in cmd
        assert "100" in cmd

    def test_empty_kwargs_produces_minimal_command(self, adapter, dummy_module):
        cmd = adapter.build_command()
        assert cmd == [sys.executable, str(dummy_module)]

    def test_unknown_key_without_dash_prefix_ignored(self, adapter):
        """Keys not in param_mapping and not starting with '--' are ignored."""
        cmd = adapter.build_command(random_param="value")
        assert len(cmd) == 2  # only base command

    def test_output_dir_maps_to_o_flag(self, adapter):
        cmd = adapter.build_command(output_dir="/out")
        assert "-o" in cmd
        assert "/out" in cmd


# ---------------------------------------------------------------------------
# TestExecute
# ---------------------------------------------------------------------------

class TestExecute:
    """Tests for CLIModuleAdapter.execute with mocked subprocess."""

    def test_execute_success(self, adapter):
        mock_process = MagicMock()
        mock_process.returncode = 0
        mock_process.stdout = "done"
        mock_process.stderr = ""

        with patch("circleseeker.modules.adapters.subprocess.run", return_value=mock_process):
            result = adapter.execute(input_file="/data/in.bam")

        assert isinstance(result, ModuleResult)
        assert result.success is True
        assert result.module_name == "test_adapter"

    def test_execute_failure_returns_error(self, adapter):
        mock_process = MagicMock()
        mock_process.returncode = 1
        mock_process.stdout = ""
        mock_process.stderr = "Segmentation fault"

        with patch("circleseeker.modules.adapters.subprocess.run", return_value=mock_process):
            result = adapter.execute(input_file="/data/in.bam")

        assert result.success is False
        assert result.error_message == "Segmentation fault"

    def test_execute_exception_returns_error(self, adapter):
        with patch(
            "circleseeker.modules.adapters.subprocess.run",
            side_effect=OSError("No such file or directory"),
        ):
            result = adapter.execute()

        assert result.success is False
        assert "No such file or directory" in result.error_message

    def test_execute_calls_subprocess_with_correct_args(self, adapter, dummy_module):
        mock_process = MagicMock(returncode=0, stdout="", stderr="")

        with patch("circleseeker.modules.adapters.subprocess.run", return_value=mock_process) as mock_run:
            adapter.execute(input_file="/data/in.bam")

        mock_run.assert_called_once()
        call_args = mock_run.call_args
        cmd = call_args[0][0]
        assert cmd[0] == sys.executable
        assert cmd[1] == str(dummy_module)
        assert "-i" in cmd
        assert call_args[1]["capture_output"] is True
        assert call_args[1]["text"] is True
        assert call_args[1]["check"] is False


# ---------------------------------------------------------------------------
# TestTandemToRingAdapter
# ---------------------------------------------------------------------------

class TestTandemToRingAdapter:
    """Tests for TandemToRingAdapter initialisation."""

    @patch("pathlib.Path.exists", return_value=True)
    def test_init_sets_correct_name(self, mock_exists):
        adapter = TandemToRingAdapter()
        assert adapter.name == "tandem_to_ring"

    @patch("pathlib.Path.exists", return_value=True)
    def test_init_sets_correct_module_path(self, mock_exists):
        adapter = TandemToRingAdapter()
        expected = Path("src/circleseeker/modules/tandem_to_ring.py")
        assert adapter.module_path == expected

    @patch("pathlib.Path.exists", return_value=True)
    def test_is_instance_of_cli_module_adapter(self, mock_exists):
        adapter = TandemToRingAdapter()
        assert isinstance(adapter, CLIModuleAdapter)

    def test_init_raises_when_module_missing(self):
        """Without mocking, the hard-coded path is unlikely to exist in test env."""
        with pytest.raises(FileNotFoundError, match="Module not found"):
            # Use a monkeypatch to ensure the path does NOT exist
            with patch("pathlib.Path.exists", return_value=False):
                TandemToRingAdapter()


# ---------------------------------------------------------------------------
# TestUMClassifyAdapter
# ---------------------------------------------------------------------------

class TestUMClassifyAdapter:
    """Tests for UMClassifyAdapter initialisation."""

    @patch("pathlib.Path.exists", return_value=True)
    def test_init_sets_correct_name(self, mock_exists):
        adapter = UMClassifyAdapter()
        assert adapter.name == "um_classify"

    @patch("pathlib.Path.exists", return_value=True)
    def test_init_sets_correct_module_path(self, mock_exists):
        adapter = UMClassifyAdapter()
        expected = Path("src/circleseeker/modules/um_classify.py")
        assert adapter.module_path == expected

    @patch("pathlib.Path.exists", return_value=True)
    def test_is_instance_of_cli_module_adapter(self, mock_exists):
        adapter = UMClassifyAdapter()
        assert isinstance(adapter, CLIModuleAdapter)

    def test_init_raises_when_module_missing(self):
        with pytest.raises(FileNotFoundError, match="Module not found"):
            with patch("pathlib.Path.exists", return_value=False):
                UMClassifyAdapter()
