"""Tests for utils validators module."""

from pathlib import Path
import sys
import pytest
from unittest.mock import patch, MagicMock

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.utils.validators import validate_installation


class TestValidateInstallation:
    """Test cases for validate_installation function."""

    @patch("circleseeker.utils.validators.importlib")
    def test_validate_installation_basic_success(self, mock_importlib):
        """Test validate_installation with all modules available."""
        # Mock successful imports
        mock_importlib.import_module.return_value = MagicMock()

        issues = validate_installation(full_check=False)
        assert issues == []

        # Check that required modules were checked
        expected_modules = ['pandas', 'numpy', 'Bio', 'pysam', 'yaml', 'networkx', 'click']
        actual_calls = [call[0][0] for call in mock_importlib.import_module.call_args_list]

        for module in expected_modules:
            assert module in actual_calls

    @patch("circleseeker.utils.validators.importlib")
    def test_validate_installation_missing_module(self, mock_importlib):
        """Test validate_installation with missing module."""
        def import_side_effect(module_name):
            if module_name == 'pandas':
                raise ImportError(f"No module named '{module_name}'")
            return MagicMock()

        mock_importlib.import_module.side_effect = import_side_effect

        issues = validate_installation(full_check=False)
        assert len(issues) == 1
        assert "Missing Python module: pandas" in issues[0]

    @patch("circleseeker.utils.validators.importlib")
    def test_validate_installation_multiple_missing_modules(self, mock_importlib):
        """Test validate_installation with multiple missing modules."""
        def import_side_effect(module_name):
            if module_name in ['pandas', 'numpy']:
                raise ImportError(f"No module named '{module_name}'")
            return MagicMock()

        mock_importlib.import_module.side_effect = import_side_effect

        issues = validate_installation(full_check=False)
        assert len(issues) == 2
        assert any("pandas" in issue for issue in issues)
        assert any("numpy" in issue for issue in issues)

    @patch("circleseeker.utils.validators.importlib")
    def test_validate_installation_biopython_special_case(self, mock_importlib):
        """Test that biopython is imported as 'Bio'."""
        # Track what modules were imported
        imported_modules = []

        def import_side_effect(module_name):
            imported_modules.append(module_name)
            return MagicMock()

        mock_importlib.import_module.side_effect = import_side_effect

        validate_installation(full_check=False)

        # Should import 'Bio' not 'biopython'
        assert 'Bio' in imported_modules
        assert 'biopython' not in imported_modules

    @patch("circleseeker.utils.validators.importlib")
    def test_validate_installation_yaml_special_case(self, mock_importlib):
        """Test that yaml module is handled correctly."""
        imported_modules = []

        def import_side_effect(module_name):
            imported_modules.append(module_name)
            return MagicMock()

        mock_importlib.import_module.side_effect = import_side_effect

        validate_installation(full_check=False)

        # Should import 'yaml'
        assert 'yaml' in imported_modules

    @patch("circleseeker.utils.validators.shutil")
    @patch("circleseeker.utils.validators.importlib")
    def test_validate_installation_full_check_success(self, mock_importlib, mock_shutil):
        """Test validate_installation with full_check=True and all tools available."""
        # Mock successful module imports
        mock_importlib.import_module.return_value = MagicMock()

        # Mock successful tool checks
        mock_shutil.which.return_value = "/usr/bin/tool"

        issues = validate_installation(full_check=True)
        assert issues == []

        # Check that external tools were checked (actual list from validators.py)
        expected_tools = [
            'TideHunter', 'minimap2', 'cd-hit-est', 'samtools',
            'bcftools', 'cyrcular', 'varlociraptor', 'cresil'
        ]
        actual_calls = [call[0][0] for call in mock_shutil.which.call_args_list]

        for tool in expected_tools:
            assert tool in actual_calls

    @patch("circleseeker.utils.validators.shutil")
    @patch("circleseeker.utils.validators.importlib")
    def test_validate_installation_missing_external_tool(self, mock_importlib, mock_shutil):
        """Test validate_installation with missing external tool."""
        # Mock successful module imports
        mock_importlib.import_module.return_value = MagicMock()

        # Mock missing tool - use TideHunter (correct case)
        def which_side_effect(tool_name):
            if tool_name in {"TideHunter", "tidehunter"}:
                return None  # Tool not found (including alt name)
            return "/usr/bin/tool"

        mock_shutil.which.side_effect = which_side_effect

        issues = validate_installation(full_check=True)
        assert len(issues) == 1
        assert "External tool not found: TideHunter" in issues[0]

    @patch("circleseeker.utils.validators.shutil")
    @patch("circleseeker.utils.validators.importlib")
    def test_validate_installation_multiple_missing_tools(self, mock_importlib, mock_shutil):
        """Test validate_installation with multiple missing external tools."""
        # Mock successful module imports
        mock_importlib.import_module.return_value = MagicMock()

        # Mock multiple missing tools - use correct case from validators.py
        def which_side_effect(tool_name):
            if tool_name in {"TideHunter", "tidehunter", "cyrcular"}:
                return None  # Tools not found (including TideHunter alt name)
            return "/usr/bin/tool"

        mock_shutil.which.side_effect = which_side_effect

        issues = validate_installation(full_check=True)
        assert len(issues) == 2
        assert any("TideHunter" in issue for issue in issues)
        assert any("cyrcular" in issue for issue in issues)

    @patch("circleseeker.utils.validators.shutil")
    @patch("circleseeker.utils.validators.importlib")
    def test_validate_installation_full_check_false_skips_tools(self, mock_importlib, mock_shutil):
        """Test that full_check=False skips external tool checking."""
        # Mock successful module imports
        mock_importlib.import_module.return_value = MagicMock()

        issues = validate_installation(full_check=False)

        # shutil.which should not be called when full_check=False
        mock_shutil.which.assert_not_called()

    def test_validate_installation_returns_list(self):
        """Test that validate_installation always returns a list."""
        with patch("circleseeker.utils.validators.importlib") as mock_importlib:
            mock_importlib.import_module.return_value = MagicMock()
            result = validate_installation(full_check=False)
            assert isinstance(result, list)

        with patch("circleseeker.utils.validators.importlib") as mock_importlib:
            with patch("circleseeker.utils.validators.shutil") as mock_shutil:
                mock_importlib.import_module.return_value = MagicMock()
                mock_shutil.which.return_value = "/usr/bin/tool"
                result = validate_installation(full_check=True)
                assert isinstance(result, list)

    @patch("circleseeker.utils.validators.importlib")
    def test_validate_installation_empty_result_on_success(self, mock_importlib):
        """Test that successful validation returns empty list."""
        mock_importlib.import_module.return_value = MagicMock()

        issues = validate_installation(full_check=False)
        assert issues == []
        assert len(issues) == 0

    @patch("circleseeker.utils.validators.shutil")
    @patch("circleseeker.utils.validators.importlib")
    def test_external_tools_list_coverage(self, mock_importlib, mock_shutil):
        """Test that all expected external tools are checked."""
        # Mock successful imports and tools
        mock_importlib.import_module.return_value = MagicMock()
        mock_shutil.which.return_value = "/usr/bin/tool"

        validate_installation(full_check=True)

        # Get all tools that were checked
        checked_tools = [call[0][0] for call in mock_shutil.which.call_args_list]

        # Verify specific tools are included (using correct case from validators.py)
        expected_tools = [
            'TideHunter', 'minimap2', 'cd-hit-est', 'samtools',
            'bcftools', 'cyrcular', 'varlociraptor', 'cresil'
        ]

        for tool in expected_tools:
            assert tool in checked_tools

    @patch("circleseeker.utils.validators.importlib")
    def test_required_modules_list_coverage(self, mock_importlib):
        """Test that all expected Python modules are checked."""
        imported_modules = []

        def import_side_effect(module_name):
            imported_modules.append(module_name)
            return MagicMock()

        mock_importlib.import_module.side_effect = import_side_effect

        validate_installation(full_check=False)

        # Verify specific modules are included
        # Note: 'biopython' -> 'Bio', others should be direct
        expected_modules = ['pandas', 'numpy', 'Bio', 'pysam', 'yaml', 'networkx', 'click']

        for module in expected_modules:
            assert module in imported_modules


class TestValidateInstallationIntegration:
    """Integration tests for validate_installation."""

    def test_validate_installation_different_scenarios(self):
        """Test validate_installation under different scenarios."""
        # Test basic mode
        with patch("circleseeker.utils.validators.importlib") as mock_importlib:
            mock_importlib.import_module.return_value = MagicMock()
            issues_basic = validate_installation(full_check=False)
            assert isinstance(issues_basic, list)

        # Test full mode
        with patch("circleseeker.utils.validators.importlib") as mock_importlib:
            with patch("circleseeker.utils.validators.shutil") as mock_shutil:
                mock_importlib.import_module.return_value = MagicMock()
                mock_shutil.which.return_value = "/usr/bin/tool"
                issues_full = validate_installation(full_check=True)
                assert isinstance(issues_full, list)

    def test_validate_installation_real_run(self):
        """Test validate_installation without mocks (real environment check)."""
        # This test runs the actual validation
        # It may have issues if some modules are missing, which is expected behavior
        issues = validate_installation(full_check=False)
        assert isinstance(issues, list)

        # In a proper dev environment, most modules should be installed
        # We're just testing that the function runs without errors

    @patch("circleseeker.utils.validators.shutil")
    @patch("circleseeker.utils.validators.importlib")
    def test_validate_installation_partial_tool_availability(self, mock_importlib, mock_shutil):
        """Test with some tools available and some missing."""
        mock_importlib.import_module.return_value = MagicMock()

        # Simulate partial tool availability
        available_tools = {'minimap2', 'samtools', 'cd-hit-est'}

        def which_side_effect(tool_name):
            if tool_name in available_tools:
                return f"/usr/bin/{tool_name}"
            return None

        mock_shutil.which.side_effect = which_side_effect

        issues = validate_installation(full_check=True)

        # Should have issues for missing tools
        missing_tools = {'TideHunter', 'bcftools', 'cyrcular', 'varlociraptor', 'cresil'}
        for tool in missing_tools:
            assert any(tool in issue for issue in issues)

    @patch("circleseeker.utils.validators.importlib")
    def test_validate_installation_all_modules_missing(self, mock_importlib):
        """Test with all modules missing."""
        def import_side_effect(module_name):
            raise ImportError(f"No module named '{module_name}'")

        mock_importlib.import_module.side_effect = import_side_effect

        issues = validate_installation(full_check=False)

        # Should have 7 issues for 7 required modules
        # (pandas, numpy, biopython, pysam, yaml, networkx, click)
        assert len(issues) >= 7

    @patch("circleseeker.utils.validators.shutil")
    @patch("circleseeker.utils.validators.importlib")
    def test_validate_installation_all_tools_missing(self, mock_importlib, mock_shutil):
        """Test with all external tools missing."""
        mock_importlib.import_module.return_value = MagicMock()
        mock_shutil.which.return_value = None  # All tools missing

        issues = validate_installation(full_check=True)

        # Should have issues for all 8 external tools
        assert len(issues) == 8
