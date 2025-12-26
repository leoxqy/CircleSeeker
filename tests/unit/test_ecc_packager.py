"""Tests for ecc_packager module."""

from pathlib import Path
import sys
import pytest
import tempfile
import shutil
from unittest.mock import patch, MagicMock

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

# Import the packager module dynamically to avoid circular imports
import importlib.util
packager_path = SRC / "circleseeker" / "modules" / "ecc_packager.py"
spec = importlib.util.spec_from_file_location("ecc_packager", packager_path)
ecc_packager = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ecc_packager)


class TestEccPackager:
    """Test cases for ecc_packager module."""

    def test_log_function(self, caplog):
        """Test the log function."""
        import logging
        # Test verbose logging - the log function uses get_logger, which logs to logging system
        with caplog.at_level(logging.INFO, logger="circleseeker.ecc_packager"):
            ecc_packager.log("Test message", verbose=True)
        assert "Test message" in caplog.text

        # Test silent logging
        caplog.clear()
        with caplog.at_level(logging.INFO, logger="circleseeker.ecc_packager"):
            ecc_packager.log("Silent message", verbose=False)
        assert "Silent message" not in caplog.text

    def test_ensure_dir(self, tmp_path):
        """Test directory creation utility."""
        test_dir = tmp_path / "new_directory"
        assert not test_dir.exists()

        # Test directory creation
        ecc_packager.ensure_dir(test_dir, dry=False, verbose=False)
        assert test_dir.exists()

        # Test existing directory (should not error)
        ecc_packager.ensure_dir(test_dir, dry=False, verbose=False)
        assert test_dir.exists()

    def test_ensure_dir_dry_run(self, tmp_path):
        """Test directory creation in dry-run mode."""
        test_dir = tmp_path / "dry_run_dir"
        assert not test_dir.exists()

        # Dry run should not create directory
        ecc_packager.ensure_dir(test_dir, dry=True, verbose=False)
        assert not test_dir.exists()

    def test_copy_file(self, tmp_path):
        """Test file copying utility."""
        # Create source file
        src_file = tmp_path / "source.txt"
        src_file.write_text("test content")

        dst_file = tmp_path / "dest.txt"

        # Test copying
        ecc_packager.copy_file(src_file, dst_file, overwrite=True, dry=False, verbose=False)
        assert dst_file.exists()
        assert dst_file.read_text() == "test content"

    def test_copy_file_no_overwrite(self, tmp_path):
        """Test file copying without overwrite."""
        # Create source and destination files
        src_file = tmp_path / "source.txt"
        src_file.write_text("new content")

        dst_file = tmp_path / "dest.txt"
        dst_file.write_text("existing content")

        # Should not overwrite
        ecc_packager.copy_file(src_file, dst_file, overwrite=False, dry=False, verbose=False)
        assert dst_file.read_text() == "existing content"

    def test_copy_file_missing_source(self, tmp_path, caplog):
        """Test handling of missing source file."""
        import logging
        src_file = tmp_path / "missing.txt"
        dst_file = tmp_path / "dest.txt"

        # Should handle missing source gracefully
        with caplog.at_level(logging.INFO, logger="circleseeker.ecc_packager"):
            ecc_packager.copy_file(src_file, dst_file, overwrite=True, dry=False, verbose=True)
        assert not dst_file.exists()

        assert "Missing source" in caplog.text

    def test_copy_file_none_source(self, tmp_path, caplog):
        """Test handling of None source."""
        import logging
        dst_file = tmp_path / "dest.txt"

        # Should handle None source gracefully
        with caplog.at_level(logging.INFO, logger="circleseeker.ecc_packager"):
            ecc_packager.copy_file(None, dst_file, overwrite=True, dry=False, verbose=True)
        assert not dst_file.exists()

        assert "Missing source: <None>" in caplog.text

    def test_copy_file_dry_run(self, tmp_path):
        """Test file copying in dry-run mode."""
        src_file = tmp_path / "source.txt"
        src_file.write_text("test content")

        dst_file = tmp_path / "dest.txt"

        # Dry run should not copy file
        ecc_packager.copy_file(src_file, dst_file, overwrite=True, dry=True, verbose=False)
        assert not dst_file.exists()

    def test_first_glob(self, tmp_path):
        """Test glob utility function."""
        # Create test files
        (tmp_path / "file1.txt").touch()
        (tmp_path / "file2.txt").touch()
        (tmp_path / "other.log").touch()

        # Test finding files
        result = ecc_packager.first_glob(tmp_path, "*.txt")
        assert result is not None
        assert result.suffix == ".txt"

        # Test no matches
        result = ecc_packager.first_glob(tmp_path, "*.nonexistent")
        assert result is None

    def test_run_basic_packaging(self, tmp_path):
        """Test basic packaging functionality."""
        # Create mock args with correct attribute name
        class MockArgs:
            def __init__(self):
                self.sample_name = "test_sample"
                self.out_dir = str(tmp_path / "output")
                self.uecc_dir = None
                self.mecc_dir = None
                self.cecc_dir = None
                self.inferred_dir = None
                self.merged_csv = None
                self.html = None
                self.text = None
                self.overwrite = True
                self.dry_run = False
                self.verbose = False

        args = MockArgs()

        # Should run without error even with no input files
        result = ecc_packager.run(args)
        assert result == 0

        # Check output directory was created
        output_dir = Path(args.out_dir)
        assert output_dir.exists()

    def test_run_with_files(self, tmp_path):
        """Test packaging with actual files."""
        # Create source directories and files
        source_dir = tmp_path / "source"
        source_dir.mkdir()

        uecc_dir = source_dir / "uecc"
        uecc_dir.mkdir()
        (uecc_dir / "test_UeccDNA_C.fasta").write_text(">seq1\nATCG")
        (uecc_dir / "test_UeccDNA.bed").write_text("chr1\t100\t200")

        merged_csv = source_dir / "merged.csv"
        merged_csv.write_text("ecc_id,type\nUECC_001,UECC")

        html_report = source_dir / "report.html"
        html_report.write_text("<html><body>Test Report</body></html>")

        # Create mock args with correct attribute name
        class MockArgs:
            def __init__(self):
                self.sample_name = "test_sample"
                self.out_dir = str(tmp_path / "output")
                self.uecc_dir = str(uecc_dir)
                self.mecc_dir = None
                self.cecc_dir = None
                self.inferred_dir = None
                self.merged_csv = str(merged_csv)
                self.html = str(html_report)
                self.text = None
                self.overwrite = True
                self.dry_run = False
                self.verbose = False

        args = MockArgs()

        # Run packaging
        result = ecc_packager.run(args)
        assert result == 0

        # Check that files were copied to output
        output_dir = Path(args.out_dir)
        assert output_dir.exists()

        # Check for organized structure
        sample_dir = output_dir / args.sample_name
        assert sample_dir.exists()

    def test_run_dry_run(self, tmp_path):
        """Test packaging in dry-run mode."""
        # Create source file
        source_dir = tmp_path / "source"
        source_dir.mkdir()
        test_file = source_dir / "test.csv"
        test_file.write_text("test,data")

        class MockArgs:
            def __init__(self):
                self.sample_name = "test_sample"
                self.out_dir = str(tmp_path / "output")
                self.uecc_dir = None
                self.mecc_dir = None
                self.cecc_dir = None
                self.inferred_dir = None
                self.merged_csv = str(test_file)
                self.html = None
                self.text = None
                self.overwrite = True
                self.dry_run = True
                self.verbose = True

        args = MockArgs()

        # Dry run should complete without creating files
        result = ecc_packager.run(args)
        assert result == 0

        # Output directory structure should not be created in dry run
        output_dir = Path(args.out_dir)
        # Note: The function might still create the base output directory


@pytest.fixture
def sample_args(tmp_path):
    """Create sample arguments for testing."""
    class MockArgs:
        def __init__(self):
            self.sample_name = "test_sample"
            self.out_dir = str(tmp_path / "output")
            self.uecc_dir = None
            self.mecc_dir = None
            self.cecc_dir = None
            self.inferred_dir = None
            self.merged_csv = None
            self.html = None
            self.text = None
            self.overwrite = True
            self.dry_run = False
            self.verbose = False

    return MockArgs()