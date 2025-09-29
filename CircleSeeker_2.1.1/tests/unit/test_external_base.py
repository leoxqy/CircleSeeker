"""Tests for external tool base class."""

from pathlib import Path
import sys
import pytest
import subprocess
from unittest.mock import patch, MagicMock
import logging

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.external.base import ExternalTool
from circleseeker.exceptions import ExternalToolError


class TestExternalTool:
    """Test cases for ExternalTool base class."""

    def test_external_tool_initialization(self):
        """Test ExternalTool basic initialization."""
        # Create a concrete implementation for testing
        class TestTool(ExternalTool):
            tool_name = "test_tool"

        with patch.object(TestTool, '_check_installation'):
            tool = TestTool(threads=4)
            assert tool.threads == 4
            assert tool.tool_name == "test_tool"
            assert tool.logger is not None

    def test_external_tool_initialization_with_logger(self):
        """Test ExternalTool initialization with custom logger."""
        class TestTool(ExternalTool):
            tool_name = "test_tool"

        custom_logger = logging.getLogger("test_logger")

        with patch.object(TestTool, '_check_installation'):
            tool = TestTool(logger=custom_logger, threads=2)
            assert tool.logger == custom_logger
            assert tool.threads == 2

    def test_external_tool_class_attributes(self):
        """Test ExternalTool class attributes."""
        class TestTool(ExternalTool):
            tool_name = "test_tool"
            required_version = "1.0.0"
            version_command = "--version"
            version_regex = r"(\d+\.\d+\.\d+)"

        with patch.object(TestTool, '_check_installation'):
            tool = TestTool()
            assert tool.tool_name == "test_tool"
            assert tool.required_version == "1.0.0"
            assert tool.version_command == "--version"
            assert tool.version_regex == r"(\d+\.\d+\.\d+)"

    @patch('shutil.which')
    def test_check_tool_availability(self, mock_which):
        """Test check_tool_availability method."""
        class TestTool(ExternalTool):
            tool_name = "test_tool"

        mock_which.return_value = "/usr/bin/test_tool"

        with patch.object(TestTool, '_check_installation'):
            tool = TestTool()
            assert tool.check_tool_availability("test_tool") is True

        mock_which.return_value = None
        assert tool.check_tool_availability("nonexistent_tool") is False

    @patch('subprocess.run')
    def test_run_command(self, mock_run):
        """Test run_command method."""
        class TestTool(ExternalTool):
            tool_name = "test_tool"

        mock_process = MagicMock()
        mock_run.return_value = mock_process

        with patch.object(TestTool, '_check_installation'):
            tool = TestTool()
            result = tool.run_command(["echo", "test"])

            mock_run.assert_called_once_with(["echo", "test"])
            assert result == mock_process

    @patch('shutil.which')
    def test_check_installation_tool_not_found(self, mock_which):
        """Test _check_installation when tool is not found."""
        class TestTool(ExternalTool):
            tool_name = "nonexistent_tool"

        mock_which.return_value = None

        with pytest.raises(ExternalToolError) as exc_info:
            TestTool()

        assert "not found in PATH" in str(exc_info.value)
        assert "nonexistent_tool" in str(exc_info.value)

    @patch('shutil.which')
    @patch.object(ExternalTool, 'get_tool_version')
    @patch.object(ExternalTool, 'check_minimum_version')
    def test_check_installation_version_check(self, mock_version_check, mock_get_version, mock_which):
        """Test _check_installation with version requirements."""
        class TestTool(ExternalTool):
            tool_name = "test_tool"
            required_version = "2.0.0"

        mock_which.return_value = "/usr/bin/test_tool"
        mock_get_version.return_value = "1.5.0"
        mock_version_check.return_value = False

        with pytest.raises(ExternalToolError) as exc_info:
            TestTool()

        assert "version" in str(exc_info.value).lower()
        assert "1.5.0" in str(exc_info.value)
        assert "2.0.0" in str(exc_info.value)

    @patch('subprocess.run')
    def test_get_tool_version_success(self, mock_run):
        """Test get_tool_version with successful version detection."""
        class TestTool(ExternalTool):
            tool_name = "test_tool"

        mock_result = MagicMock()
        mock_result.stdout = "test_tool version 2.1.3\n"
        mock_result.stderr = ""
        mock_run.return_value = mock_result

        with patch.object(TestTool, '_check_installation'):
            tool = TestTool()
            version = tool.get_tool_version("test_tool")
            assert version == "2.1.3"

    @patch('subprocess.run')
    def test_get_tool_version_from_stderr(self, mock_run):
        """Test get_tool_version when version is in stderr."""
        class TestTool(ExternalTool):
            tool_name = "test_tool"

        mock_result = MagicMock()
        mock_result.stdout = ""
        mock_result.stderr = "test_tool v1.0.5 - help message\n"
        mock_run.return_value = mock_result

        with patch.object(TestTool, '_check_installation'):
            tool = TestTool()
            version = tool.get_tool_version("test_tool")
            assert version == "1.0.5"

    @patch('subprocess.run')
    def test_get_tool_version_no_version_command(self, mock_run):
        """Test get_tool_version when version_command is None."""
        class TestTool(ExternalTool):
            tool_name = "test_tool"
            version_command = None

        with patch.object(TestTool, '_check_installation'):
            tool = TestTool()
            version = tool.get_tool_version("test_tool")
            assert version is None

    @patch('subprocess.run')
    def test_get_tool_version_timeout(self, mock_run):
        """Test get_tool_version with timeout."""
        class TestTool(ExternalTool):
            tool_name = "test_tool"

        mock_run.side_effect = subprocess.TimeoutExpired(["test_tool", "--version"], 10)

        with patch.object(TestTool, '_check_installation'):
            tool = TestTool()
            version = tool.get_tool_version("test_tool")
            assert version is None

    def test_check_minimum_version_valid(self):
        """Test check_minimum_version with valid versions."""
        class TestTool(ExternalTool):
            tool_name = "test_tool"

        with patch.object(TestTool, '_check_installation'):
            tool = TestTool()

            # Current version meets requirement
            assert tool.check_minimum_version("2.1.0", "2.0.0") is True
            assert tool.check_minimum_version("2.0.0", "2.0.0") is True

            # Current version below requirement
            assert tool.check_minimum_version("1.9.0", "2.0.0") is False

    def test_check_minimum_version_invalid_format(self):
        """Test check_minimum_version with invalid version formats."""
        class TestTool(ExternalTool):
            tool_name = "test_tool"

        with patch.object(TestTool, '_check_installation'):
            tool = TestTool()

            # Should return True (assume OK) for unparseable versions
            assert tool.check_minimum_version("invalid", "2.0.0") is True
            assert tool.check_minimum_version("2.0.0", "invalid") is True

    @patch('shutil.which')
    def test_get_tool_info(self, mock_which):
        """Test get_tool_info method."""
        class TestTool(ExternalTool):
            tool_name = "test_tool"
            required_version = "1.0.0"

        mock_which.return_value = "/usr/bin/test_tool"

        with patch.object(TestTool, '_check_installation'):
            with patch.object(TestTool, 'get_tool_version', return_value="1.2.3"):
                tool = TestTool()
                info = tool.get_tool_info()

                assert info['name'] == "test_tool"
                assert info['available'] is True
                assert info['version'] == "1.2.3"
                assert info['required_version'] == "1.0.0"
                assert info['path'] == "/usr/bin/test_tool"

    @patch('subprocess.run')
    def test_run_success(self, mock_run):
        """Test run method with successful execution."""
        class TestTool(ExternalTool):
            tool_name = "test_tool"

        mock_result = MagicMock()
        mock_result.stdout = "success output"
        mock_result.stderr = ""
        mock_result.returncode = 0
        mock_run.return_value = mock_result

        with patch.object(TestTool, '_check_installation'):
            tool = TestTool()
            stdout, stderr = tool.run(["echo", "test"])

            assert stdout == "success output"
            assert stderr == ""

    @patch('subprocess.run')
    def test_run_with_stderr(self, mock_run):
        """Test run method with stderr output."""
        class TestTool(ExternalTool):
            tool_name = "test_tool"

        mock_result = MagicMock()
        mock_result.stdout = "output"
        mock_result.stderr = "warning message"
        mock_result.returncode = 0
        mock_run.return_value = mock_result

        with patch.object(TestTool, '_check_installation'):
            tool = TestTool()
            stdout, stderr = tool.run(["echo", "test"])

            assert stdout == "output"
            assert stderr == "warning message"

    @patch('subprocess.run')
    def test_run_timeout(self, mock_run):
        """Test run method with timeout."""
        class TestTool(ExternalTool):
            tool_name = "test_tool"

        mock_run.side_effect = subprocess.TimeoutExpired(["test_tool"], 30)

        with patch.object(TestTool, '_check_installation'):
            tool = TestTool()

            with pytest.raises(ExternalToolError) as exc_info:
                tool.run(["test_tool"], timeout=30)

            assert "timed out" in str(exc_info.value)

    @patch('subprocess.run')
    def test_run_called_process_error(self, mock_run):
        """Test run method with CalledProcessError."""
        class TestTool(ExternalTool):
            tool_name = "test_tool"

        error = subprocess.CalledProcessError(1, ["test_tool"], stderr="error message")
        mock_run.side_effect = error

        with patch.object(TestTool, '_check_installation'):
            tool = TestTool()

            with pytest.raises(ExternalToolError) as exc_info:
                tool.run(["test_tool"])

            assert "failed" in str(exc_info.value)

    @patch('subprocess.run')
    def test_run_os_error(self, mock_run):
        """Test run method with OSError."""
        class TestTool(ExternalTool):
            tool_name = "test_tool"

        mock_run.side_effect = OSError("No such file or directory")

        with patch.object(TestTool, '_check_installation'):
            tool = TestTool()

            with pytest.raises(ExternalToolError) as exc_info:
                tool.run(["test_tool"])

            assert "Failed to execute" in str(exc_info.value)

    @patch('subprocess.run')
    def test_run_no_capture_output(self, mock_run):
        """Test run method with capture_output=False."""
        class TestTool(ExternalTool):
            tool_name = "test_tool"

        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_run.return_value = mock_result

        with patch.object(TestTool, '_check_installation'):
            tool = TestTool()
            stdout, stderr = tool.run(["echo", "test"], capture_output=False)

            assert stdout == ""
            assert stderr == ""

    @patch('shutil.which')
    def test_check_installation_with_additional_tools(self, mock_which):
        """Test _check_installation with additional required tools."""
        class TestTool(ExternalTool):
            tool_name = "test_tool"

            def _get_required_tools(self):
                return ["helper_tool", "another_tool"]

        def which_side_effect(tool):
            if tool == "test_tool":
                return "/usr/bin/test_tool"
            elif tool == "helper_tool":
                return "/usr/bin/helper_tool"
            else:
                return None  # another_tool not found

        mock_which.side_effect = which_side_effect

        with pytest.raises(ExternalToolError) as exc_info:
            TestTool()

        assert "Required dependency 'another_tool' not found" in str(exc_info.value)


@pytest.fixture
def mock_tool():
    """Create a mock tool for testing."""
    class MockTool(ExternalTool):
        tool_name = "mock_tool"
        required_version = "1.0.0"

    with patch.object(MockTool, '_check_installation'):
        return MockTool()


class TestExternalToolIntegration:
    """Integration tests for ExternalTool."""

    @patch('shutil.which')
    @patch('subprocess.run')
    def test_complete_workflow(self, mock_run, mock_which):
        """Test complete workflow from initialization to execution."""
        class WorkflowTool(ExternalTool):
            tool_name = "workflow_tool"
            required_version = "2.0.0"

        # Mock tool availability
        mock_which.return_value = "/usr/bin/workflow_tool"

        # Mock version check
        version_result = MagicMock()
        version_result.stdout = "workflow_tool version 2.1.0"
        version_result.stderr = ""

        # Mock actual command execution
        exec_result = MagicMock()
        exec_result.stdout = "command completed"
        exec_result.stderr = ""
        exec_result.returncode = 0

        mock_run.side_effect = [version_result, exec_result]

        # Create and use tool
        tool = WorkflowTool(threads=4)
        stdout, stderr = tool.run(["workflow_tool", "--input", "test.txt"])

        assert stdout == "command completed"
        assert stderr == ""
        assert tool.threads == 4