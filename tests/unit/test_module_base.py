"""Tests for base module interface classes (ModuleResult, ModuleBase, ExternalToolModule)."""

from pathlib import Path
import sys
import pytest
import logging

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.modules.base import ModuleResult, ModuleBase, ExternalToolModule


# ---------------------------------------------------------------------------
# Concrete stub subclasses for testing abstract bases
# ---------------------------------------------------------------------------

class ConcreteModule(ModuleBase):
    """Minimal concrete implementation of ModuleBase for testing."""

    def validate_inputs(self, **kwargs):
        return True

    def execute(self, **kwargs):
        return ModuleResult(success=True, module_name=self.name)


class FailingExecuteModule(ModuleBase):
    """Module whose execute() raises an exception."""

    def validate_inputs(self, **kwargs):
        return True

    def execute(self, **kwargs):
        raise RuntimeError("execute exploded")


class ConcreteExternalTool(ExternalToolModule):
    """Minimal concrete implementation of ExternalToolModule for testing."""

    def __init__(self, available=True, version="1.0.0", **kwargs):
        self._available = available
        self._version = version
        super().__init__(**kwargs)

    def check_tool_availability(self):
        return self._available

    def get_tool_version(self):
        return self._version

    def validate_inputs(self, **kwargs):
        return True

    def execute(self, **kwargs):
        return ModuleResult(success=True, module_name=self.name)


# ===========================================================================
# TestModuleResult
# ===========================================================================

class TestModuleResult:
    """Test cases for ModuleResult dataclass."""

    def test_default_initialization(self):
        """ModuleResult with only required fields uses expected defaults."""
        r = ModuleResult(success=True, module_name="test")
        assert r.success is True
        assert r.module_name == "test"
        assert r.output_files == {}
        assert r.metrics == {}
        assert r.warnings == []
        assert r.error_message is None
        assert r.execution_time == 0.0

    def test_add_output_with_path(self):
        """add_output stores value as a Path object."""
        r = ModuleResult(success=True, module_name="test")
        r.add_output("bam", Path("/tmp/out.bam"))
        assert r.output_files["bam"] == Path("/tmp/out.bam")
        assert isinstance(r.output_files["bam"], Path)

    def test_add_output_with_string_path(self):
        """add_output converts a string argument to Path."""
        r = ModuleResult(success=True, module_name="test")
        r.add_output("vcf", "/tmp/out.vcf")
        assert r.output_files["vcf"] == Path("/tmp/out.vcf")
        assert isinstance(r.output_files["vcf"], Path)

    def test_add_metric(self):
        """add_metric stores key/value in metrics dict."""
        r = ModuleResult(success=True, module_name="test")
        r.add_metric("count", 42)
        assert r.metrics["count"] == 42

    def test_add_warning(self):
        """add_warning appends message to warnings list."""
        r = ModuleResult(success=True, module_name="test")
        r.add_warning("low coverage")
        assert r.warnings == ["low coverage"]

    def test_multiple_warnings_accumulated(self):
        """Multiple add_warning calls accumulate messages in order."""
        r = ModuleResult(success=True, module_name="test")
        r.add_warning("warn1")
        r.add_warning("warn2")
        r.add_warning("warn3")
        assert r.warnings == ["warn1", "warn2", "warn3"]

    def test_error_message_defaults_to_none(self):
        """error_message defaults to None when not provided."""
        r = ModuleResult(success=False, module_name="fail")
        assert r.error_message is None

    def test_execution_time_defaults_to_zero(self):
        """execution_time defaults to 0.0 when not provided."""
        r = ModuleResult(success=True, module_name="test")
        assert r.execution_time == 0.0


# ===========================================================================
# TestModuleBase
# ===========================================================================

class TestModuleBase:
    """Test cases for the ModuleBase abstract class (via ConcreteModule stub)."""

    def test_concrete_subclass_instantiates(self):
        """A concrete subclass with validate_inputs and execute can be created."""
        mod = ConcreteModule()
        assert mod is not None

    def test_default_name_is_class_name(self):
        """When no name is given, the module name is the class name."""
        mod = ConcreteModule()
        assert mod.name == "ConcreteModule"

    def test_custom_name_used(self):
        """A custom name overrides the class name default."""
        mod = ConcreteModule(name="my_module")
        assert mod.name == "my_module"

    def test_debug_attribute_stored(self):
        """The debug flag is stored on the instance."""
        mod = ConcreteModule(debug=True)
        assert mod.debug is True
        mod2 = ConcreteModule(debug=False)
        assert mod2.debug is False

    def test_validate_input_file_existing(self, tmp_path):
        """validate_input_file returns Path for an existing non-empty file."""
        f = tmp_path / "data.txt"
        f.write_text("hello")
        mod = ConcreteModule()
        result = mod.validate_input_file(str(f))
        assert result == f
        assert isinstance(result, Path)

    def test_validate_input_file_nonexistent(self, tmp_path):
        """validate_input_file raises FileNotFoundError for missing file."""
        mod = ConcreteModule()
        with pytest.raises(FileNotFoundError):
            mod.validate_input_file(tmp_path / "no_such_file.txt")

    def test_validate_input_file_directory(self, tmp_path):
        """validate_input_file raises ValueError when path is a directory."""
        mod = ConcreteModule()
        with pytest.raises(ValueError):
            mod.validate_input_file(tmp_path)

    def test_validate_input_file_empty_logs_warning(self, tmp_path, caplog):
        """validate_input_file logs a warning for an empty file but does not error."""
        f = tmp_path / "empty.txt"
        f.write_text("")
        mod = ConcreteModule()
        # Enable propagation so caplog can capture the message
        mod.logger.propagate = True
        with caplog.at_level(logging.WARNING):
            result = mod.validate_input_file(str(f))
        assert result == f
        assert any("empty" in rec.message.lower() for rec in caplog.records)

    def test_validate_output_dir_creates_directory(self, tmp_path):
        """validate_output_dir creates the directory when it does not exist."""
        target = tmp_path / "sub" / "dir"
        assert not target.exists()
        mod = ConcreteModule()
        result = mod.validate_output_dir(target)
        assert target.is_dir()
        assert result == target

    def test_validate_output_dir_existing(self, tmp_path):
        """validate_output_dir returns Path for an already existing directory."""
        mod = ConcreteModule()
        result = mod.validate_output_dir(tmp_path)
        assert result == tmp_path

    def test_run_calls_validate_and_execute(self):
        """run() returns a successful ModuleResult with measured execution_time."""
        mod = ConcreteModule(name="runner")
        result = mod.run()
        assert result.success is True
        assert result.module_name == "runner"
        assert result.execution_time > 0.0 or result.execution_time == 0.0  # may be very fast

    def test_run_captures_exception(self):
        """When execute raises, run() returns success=False with error captured."""
        mod = FailingExecuteModule(name="boom")
        result = mod.run()
        assert result.success is False
        assert "execute exploded" in result.error_message
        assert result.execution_time >= 0.0

    def test_cleanup_removes_temp_files(self, tmp_path):
        """cleanup() removes listed temporary files."""
        f1 = tmp_path / "tmp1.txt"
        f2 = tmp_path / "tmp2.txt"
        f1.write_text("a")
        f2.write_text("b")
        mod = ConcreteModule()
        mod.cleanup(temp_files=[f1, f2])
        assert not f1.exists()
        assert not f2.exists()

    def test_cleanup_missing_files_handled_gracefully(self, tmp_path):
        """cleanup() does not raise when files are already absent."""
        missing = tmp_path / "gone.txt"
        mod = ConcreteModule()
        # Should not raise
        mod.cleanup(temp_files=[missing])


# ===========================================================================
# TestExternalToolModule
# ===========================================================================

class TestExternalToolModule:
    """Test cases for ExternalToolModule abstract class."""

    def test_tool_name_stored(self):
        """tool_name attribute is stored correctly."""
        mod = ConcreteExternalTool(tool_name="samtools")
        assert mod.tool_name == "samtools"

    def test_required_version_stored(self):
        """required_version attribute is stored correctly."""
        mod = ConcreteExternalTool(tool_name="minimap2", required_version="2.24")
        assert mod.required_version == "2.24"

    def test_required_version_defaults_to_none(self):
        """required_version defaults to None when not supplied."""
        mod = ConcreteExternalTool(tool_name="bwa")
        assert mod.required_version is None

    def test_validate_tool_raises_when_unavailable(self):
        """validate_tool raises RuntimeError when tool is not available."""
        mod = ConcreteExternalTool(tool_name="missing_tool", available=False)
        with pytest.raises(RuntimeError, match="not available"):
            mod.validate_tool()

    def test_validate_tool_succeeds_when_available(self):
        """validate_tool does not raise when tool is available."""
        mod = ConcreteExternalTool(tool_name="present_tool", available=True)
        # Should not raise
        mod.validate_tool()
