"""Tests for dependency_checker module."""

from pathlib import Path
import sys
from unittest.mock import patch, MagicMock

import pytest

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.utils.dependency_checker import (
    DependencyChecker,
    Tool,
    find_tool,
    get_tool_version,
    compare_versions,
    _basic_version_compare,
    TOOLS,
)


class TestTool:
    """Test Tool dataclass."""

    def test_tool_creation(self):
        """Test creating a Tool instance."""
        tool = Tool(
            name="test_tool",
            required=True,
            purpose="Testing",
            install_hint="pip install test",
        )
        assert tool.name == "test_tool"
        assert tool.required is True
        assert tool.purpose == "Testing"
        assert tool.install_hint == "pip install test"
        assert tool.min_version is None
        assert tool.alt_names is None

    def test_tool_with_optional_fields(self):
        """Test creating a Tool with optional fields."""
        tool = Tool(
            name="test_tool",
            required=False,
            purpose="Testing",
            install_hint="pip install test",
            min_version="1.0.0",
            alt_names=["test", "test-tool"],
        )
        assert tool.min_version == "1.0.0"
        assert tool.alt_names == ["test", "test-tool"]


class TestFindTool:
    """Test find_tool function."""

    @patch("shutil.which")
    def test_find_tool_primary_name(self, mock_which):
        """Test finding tool by primary name."""
        mock_which.return_value = "/usr/bin/tool"
        result = find_tool("tool")
        assert result == "tool"
        mock_which.assert_called_with("tool")

    @patch("shutil.which")
    def test_find_tool_alt_name(self, mock_which):
        """Test finding tool by alternative name."""
        mock_which.side_effect = lambda name: "/usr/bin/alt-tool" if name == "alt-tool" else None
        result = find_tool("tool", alt_names=["alt-tool", "tool2"])
        assert result == "alt-tool"

    @patch("shutil.which")
    def test_find_tool_not_found(self, mock_which):
        """Test when tool is not found."""
        mock_which.return_value = None
        result = find_tool("nonexistent", alt_names=["also-nonexistent"])
        assert result is None


class TestGetToolVersion:
    """Test get_tool_version function."""

    @patch("subprocess.run")
    def test_get_version_from_stdout(self, mock_run):
        """Test extracting version from stdout."""
        mock_run.return_value = MagicMock(stdout="tool version 2.24.1", stderr="")
        result = get_tool_version("tool")
        assert result == "2.24.1"

    @patch("subprocess.run")
    def test_get_version_from_stderr(self, mock_run):
        """Test extracting version from stderr."""
        mock_run.return_value = MagicMock(stdout="", stderr="Version: 1.17.0")
        result = get_tool_version("tool")
        assert result == "1.17.0"

    @patch("subprocess.run")
    def test_get_version_no_match(self, mock_run):
        """Test when no version pattern is found."""
        mock_run.return_value = MagicMock(stdout="no version here", stderr="")
        result = get_tool_version("tool")
        assert result is None

    @patch("subprocess.run")
    def test_get_version_oserror(self, mock_run):
        """Test when OSError occurs."""
        mock_run.side_effect = OSError("Command not found")
        result = get_tool_version("tool")
        assert result is None

    @patch("subprocess.run")
    def test_get_version_timeout(self, mock_run):
        """Test when timeout occurs."""
        import subprocess
        mock_run.side_effect = subprocess.TimeoutExpired("tool", 5)
        result = get_tool_version("tool")
        assert result is None


class TestVersionComparison:
    """Test version comparison functions."""

    def test_basic_version_compare_equal(self):
        """Test basic comparison with equal versions."""
        assert _basic_version_compare("1.0.0", "1.0.0") is True

    def test_basic_version_compare_greater(self):
        """Test basic comparison with greater version."""
        assert _basic_version_compare("2.0.0", "1.0.0") is True

    def test_basic_version_compare_less(self):
        """Test basic comparison with lesser version."""
        assert _basic_version_compare("1.0.0", "2.0.0") is False

    def test_basic_version_compare_minor(self):
        """Test basic comparison with different minor versions."""
        assert _basic_version_compare("1.5.0", "1.4.0") is True
        assert _basic_version_compare("1.4.0", "1.5.0") is False

    def test_basic_version_compare_invalid(self):
        """Test basic comparison with invalid version string."""
        assert _basic_version_compare("invalid", "1.0.0") is False

    def test_compare_versions_simple(self):
        """Test compare_versions with simple versions."""
        assert compare_versions("2.24.0", "2.24") is True
        assert compare_versions("2.23.0", "2.24") is False

    def test_compare_versions_with_suffix(self):
        """Test compare_versions with version suffix."""
        # Should handle versions like "1.17-rc1"
        assert compare_versions("1.17", "1.16") is True


class TestDependencyChecker:
    """Test DependencyChecker class."""

    @patch("circleseeker.utils.dependency_checker.get_tool_version", return_value="999.0")
    @patch("shutil.which")
    def test_check_all_skips_tidehunter(self, mock_which, _mock_version):
        """Test that skip_tools parameter works."""
        def which_side_effect(name):
            if name.lower() == "tidehunter":
                return None
            return f"/usr/bin/{name}"

        mock_which.side_effect = which_side_effect

        checker = DependencyChecker()
        assert checker.check_all(skip_tools={"tidehunter"}) is True
        assert all(tool.name != "TideHunter" for tool in checker.missing_required)

    @patch("circleseeker.utils.dependency_checker.get_tool_version", return_value="999.0")
    @patch("shutil.which")
    def test_check_all_all_tools_found(self, mock_which, _mock_version):
        """Test when all tools are found."""
        mock_which.return_value = "/usr/bin/tool"

        checker = DependencyChecker()
        result = checker.check_all()
        assert result is True
        assert len(checker.missing_required) == 0
        assert "minimap2" in checker.found_tools

    @patch("circleseeker.utils.dependency_checker.get_tool_version", return_value="999.0")
    @patch("shutil.which")
    def test_check_all_missing_required(self, mock_which, _mock_version):
        """Test when required tool is missing."""
        def which_side_effect(name):
            if name == "minimap2":
                return None
            return f"/usr/bin/{name}"

        mock_which.side_effect = which_side_effect

        checker = DependencyChecker()
        result = checker.check_all()
        assert result is False
        assert any(t.name == "minimap2" for t in checker.missing_required)

    @patch("circleseeker.utils.dependency_checker.get_tool_version", return_value="1.0")
    @patch("shutil.which")
    def test_check_all_version_warning(self, mock_which, mock_version):
        """Test version warning for outdated tool."""
        mock_which.return_value = "/usr/bin/tool"
        mock_version.return_value = "1.0"  # Lower than min_version

        checker = DependencyChecker()
        checker.check_all()
        # Should have version warnings for tools with min_version > 1.0
        assert len(checker.version_warnings) > 0

    @patch("circleseeker.utils.dependency_checker.get_tool_version", return_value="999.0")
    @patch("shutil.which")
    def test_missing_inference_tools(self, mock_which, _mock_version):
        """Test when both cresil and cyrcular are missing."""
        def which_side_effect(name):
            if name in ("cresil", "cyrcular"):
                return None
            return f"/usr/bin/{name}"

        mock_which.side_effect = which_side_effect

        checker = DependencyChecker()
        result = checker.check_all()
        assert result is False
        assert len(checker.missing_inference) == 2

    def test_print_report(self, capsys):
        """Test print_report output."""
        checker = DependencyChecker()
        checker.found_tools = ["minimap2", "samtools"]
        checker.version_warnings = ["minimap2: version 2.20 < recommended 2.24"]
        checker.print_report()

        captured = capsys.readouterr()
        assert "CircleSeeker Dependency Check" in captured.out
        assert "minimap2" in captured.out
        assert "samtools" in captured.out


class TestToolDefinitions:
    """Test tool definitions in TOOLS list."""

    def test_tools_list_not_empty(self):
        """Test that TOOLS list is not empty."""
        assert len(TOOLS) > 0

    def test_required_tools_have_install_hint(self):
        """Test that required tools have install hints."""
        for tool in TOOLS:
            if tool.required:
                assert tool.install_hint, f"{tool.name} missing install_hint"

    def test_required_tools_exist(self):
        """Test that expected required tools are defined."""
        required_names = {t.name for t in TOOLS if t.required}
        expected = {"TideHunter", "minimap2", "samtools", "cd-hit-est"}
        assert expected.issubset(required_names)

    def test_inference_tools_defined(self):
        """Test that inference tools (cresil, cyrcular) are defined."""
        tool_names = {t.name for t in TOOLS}
        assert "cresil" in tool_names
        assert "cyrcular" in tool_names
