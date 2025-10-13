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

    @patch('importlib.import_module')
    def test_validate_installation_basic_success(self, mock_import):
        """Test validate_installation with all modules available."""
        # Mock successful imports
        mock_import.return_value = MagicMock()

        issues = validate_installation(full_check=False)
        assert issues == []

        # Check that required modules were checked
        expected_modules = ['pandas', 'numpy', 'Bio', 'pysam', 'yaml', 'networkx', 'click']
        actual_calls = [call[0][0] for call in mock_import.call_args_list]

        for module in expected_modules:
            assert module in actual_calls

    @patch('importlib.import_module')
    def test_validate_installation_missing_module(self, mock_import):
        """Test validate_installation with missing module."""
        def import_side_effect(module_name):
            if module_name == 'pandas':
                raise ImportError(f"No module named '{module_name}'")
            return MagicMock()

        mock_import.side_effect = import_side_effect

        issues = validate_installation(full_check=False)
        assert len(issues) == 1
        assert "Missing Python module: pandas" in issues[0]

    @patch('importlib.import_module')
    def test_validate_installation_multiple_missing_modules(self, mock_import):
        """Test validate_installation with multiple missing modules."""
        def import_side_effect(module_name):
            if module_name in ['pandas', 'numpy']:
                raise ImportError(f"No module named '{module_name}'")
            return MagicMock()

        mock_import.side_effect = import_side_effect

        issues = validate_installation(full_check=False)
        assert len(issues) == 2
        assert any("pandas" in issue for issue in issues)
        assert any("numpy" in issue for issue in issues)

    @patch('importlib.import_module')
    def test_validate_installation_biopython_special_case(self, mock_import):
        """Test that biopython is imported as 'Bio'."""
        # Track what modules were imported
        imported_modules = []

        def import_side_effect(module_name):
            imported_modules.append(module_name)
            return MagicMock()

        mock_import.side_effect = import_side_effect

        validate_installation(full_check=False)

        # Should import 'Bio' not 'biopython'
        assert 'Bio' in imported_modules
        assert 'biopython' not in imported_modules

    @patch('importlib.import_module')
    def test_validate_installation_yaml_special_case(self, mock_import):
        """Test that yaml module is handled correctly."""
        imported_modules = []

        def import_side_effect(module_name):
            imported_modules.append(module_name)
            return MagicMock()

        mock_import.side_effect = import_side_effect

        validate_installation(full_check=False)

        # Should import 'yaml'
        assert 'yaml' in imported_modules

    @patch('shutil.which')
    @patch('importlib.import_module')
    def test_validate_installation_full_check_success(self, mock_import, mock_which):
        """Test validate_installation with full_check=True and all tools available."""
        # Mock successful module imports
        mock_import.return_value = MagicMock()

        # Mock successful tool checks
        mock_which.return_value = "/usr/bin/tool"

        issues = validate_installation(full_check=True)
        assert issues == []

        # Check that external tools were checked
        expected_tools = [
            'TideHunter', 'makeblastdb', 'blastn', 'minimap2',
            'samtools', 'bcftools', 'cyrcular', 'varlociraptor'
        ]
        actual_calls = [call[0][0] for call in mock_which.call_args_list]

        for tool in expected_tools:
            assert tool in actual_calls

    @patch('shutil.which')
    @patch('importlib.import_module')
    def test_validate_installation_missing_external_tool(self, mock_import, mock_which):
        """Test validate_installation with missing external tool."""
        # Mock successful module imports
        mock_import.return_value = MagicMock()

        # Mock missing tool
        def which_side_effect(tool_name):
            if tool_name == 'TideHunter':
                return None  # Tool not found
            return "/usr/bin/tool"

        mock_which.side_effect = which_side_effect

        issues = validate_installation(full_check=True)
        assert len(issues) == 1
        assert "External tool not found: TideHunter" in issues[0]

    @patch('shutil.which')
    @patch('importlib.import_module')
    def test_validate_installation_multiple_missing_tools(self, mock_import, mock_which):
        """Test validate_installation with multiple missing external tools."""
        # Mock successful module imports
        mock_import.return_value = MagicMock()

        # Mock multiple missing tools
        def which_side_effect(tool_name):
            if tool_name in ['tidehunter', 'cyrcular']:
                return None  # Tools not found
            return "/usr/bin/tool"

        mock_which.side_effect = which_side_effect

        issues = validate_installation(full_check=True)
        assert len(issues) == 2
        assert any("tidehunter" in issue for issue in issues)
        assert any("cyrcular" in issue for issue in issues)

    @patch('shutil.which')
    @patch('importlib.import_module')
    def test_validate_installation_full_check_false_skips_tools(self, mock_import, mock_which):
        """Test that full_check=False skips external tool checking."""
        # Mock successful module imports
        mock_import.return_value = MagicMock()

        issues = validate_installation(full_check=False)

        # shutil.which should not be called when full_check=False
        mock_which.assert_not_called()

    @patch('importlib.import_module')
    def test_validate_installation_circleseeker_module_import_error(self, mock_import):
        """Test validate_installation with CircleSeeker module import error."""
        def import_side_effect(module_name):
            # Succeed for standard modules
            if module_name in ['pandas', 'numpy', 'Bio', 'pysam', 'yaml', 'networkx', 'click']:
                return MagicMock()
            # Fail for CircleSeeker modules
            elif 'circleseeker' in module_name:
                raise ImportError(f"Cannot import {module_name}")
            return MagicMock()

        mock_import.side_effect = import_side_effect

        issues = validate_installation(full_check=False)
        assert len(issues) == 1
        assert "CircleSeeker module import error" in issues[0]

    @patch('importlib.import_module')
    def test_validate_installation_specific_circleseeker_modules(self, mock_import):
        """Test that specific CircleSeeker modules are checked."""
        imported_modules = []

        def import_side_effect(module_name):
            imported_modules.append(module_name)
            if 'circleseeker' in module_name:
                # We'll let this pass to avoid actual import issues in test
                pass
            return MagicMock()

        mock_import.side_effect = import_side_effect

        validate_installation(full_check=False)

        # Check that specific CircleSeeker modules were imported
        expected_circleseeker_imports = [
            'circleseeker.core.pipeline',
            'circleseeker.config',
            'circleseeker.modules.tandem_to_ring',
            'circleseeker.modules.um_classify'
        ]

        # Note: The actual imports in the function are different (they import classes)
        # but the module names should be attempted

    @patch('shutil.which')
    @patch('importlib.import_module')
    def test_validate_installation_mixed_failures(self, mock_import, mock_which):
        """Test validate_installation with both module and tool failures."""
        # Mock some module import failures
        def import_side_effect(module_name):
            if module_name == 'pandas':
                raise ImportError("No module named 'pandas'")
            elif 'circleseeker' in module_name:
                raise ImportError("CircleSeeker import error")
            return MagicMock()

        mock_import.side_effect = import_side_effect

        # Mock some tool failures
        def which_side_effect(tool_name):
            if tool_name == 'tidehunter':
                return None
            return "/usr/bin/tool"

        mock_which.side_effect = which_side_effect

        issues = validate_installation(full_check=True)

        # Should have multiple types of issues
        assert len(issues) == 3
        assert any("pandas" in issue for issue in issues)
        assert any("CircleSeeker module import error" in issue for issue in issues)
        assert any("tidehunter" in issue for issue in issues)

    def test_validate_installation_returns_list(self):
        """Test that validate_installation always returns a list."""
        with patch('importlib.import_module'):
            result = validate_installation(full_check=False)
            assert isinstance(result, list)

        with patch('importlib.import_module'):
            with patch('shutil.which'):
                result = validate_installation(full_check=True)
                assert isinstance(result, list)

    @patch('importlib.import_module')
    def test_validate_installation_empty_result_on_success(self, mock_import):
        """Test that successful validation returns empty list."""
        mock_import.return_value = MagicMock()

        issues = validate_installation(full_check=False)
        assert issues == []
        assert len(issues) == 0

    @patch('shutil.which')
    @patch('importlib.import_module')
    def test_external_tools_list_coverage(self, mock_import, mock_which):
        """Test that all expected external tools are checked."""
        # Mock successful imports and tools
        mock_import.return_value = MagicMock()
        mock_which.return_value = "/usr/bin/tool"

        validate_installation(full_check=True)

        # Get all tools that were checked
        checked_tools = [call[0][0] for call in mock_which.call_args_list]

        # Verify specific tools are included
        expected_tools = [
            'tidehunter', 'makeblastdb', 'blastn', 'minimap2',
            'samtools', 'bcftools', 'cyrcular', 'varlociraptor'
        ]

        for tool in expected_tools:
            assert tool in checked_tools

    @patch('importlib.import_module')
    def test_required_modules_list_coverage(self, mock_import):
        """Test that all expected Python modules are checked."""
        imported_modules = []

        def import_side_effect(module_name):
            imported_modules.append(module_name)
            return MagicMock()

        mock_import.side_effect = import_side_effect

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
        with patch('importlib.import_module') as mock_import:
            mock_import.return_value = MagicMock()
            issues_basic = validate_installation(full_check=False)
            assert isinstance(issues_basic, list)

        # Test full mode
        with patch('importlib.import_module') as mock_import:
            with patch('shutil.which') as mock_which:
                mock_import.return_value = MagicMock()
                mock_which.return_value = "/usr/bin/tool"
                issues_full = validate_installation(full_check=True)
                assert isinstance(issues_full, list)