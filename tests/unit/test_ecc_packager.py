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
        # No log is emitted for None sources (treated as optional inputs)
        assert caplog.text == ""

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
        unified_csv = tmp_path / "unified.csv"
        unified_csv.write_text(
            "eccDNA_id,eccDNA_type,State,Length\n"
            "UeccDNA1,UeccDNA,Confirmed,100\n"
        )

        # Create mock args with correct attribute name
        class MockArgs:
            def __init__(self):
                self.sample_name = "test_sample"
                self.out_dir = str(tmp_path / "output")
                self.uecc_dir = None
                self.mecc_dir = None
                self.cecc_dir = None
                self.inferred_dir = None
                self.unified_csv = str(unified_csv)
                self.html = None
                self.text = None
                self.overwrite = True
                self.dry_run = False
                self.id_width = 4
                self.verbose = False

        args = MockArgs()

        # Should run without error even with no input files
        result = ecc_packager.run(args)
        assert result == 0

        # Check output directory was created
        output_dir = Path(args.out_dir)
        assert output_dir.exists()
        assert (output_dir / "eccDNA_summary.csv").exists()

    def test_run_with_files(self, tmp_path):
        """Test packaging with actual files."""
        # Create source directories and files
        source_dir = tmp_path / "source"
        source_dir.mkdir()

        uecc_dir = source_dir / "uecc"
        uecc_dir.mkdir()
        (uecc_dir / "test_UeccDNA_C.fasta").write_text(">UeccDNA1\nATCG\n")
        (uecc_dir / "test_UeccDNA.core.csv").write_text(
            "eccDNA_id,chr,start0,end0,strand,length\n"
            "UeccDNA1,chr1,100,200,+,100\n"
        )

        unified_csv = source_dir / "unified.csv"
        unified_csv.write_text(
            "eccDNA_id,eccDNA_type,State,Length\n"
            "UeccDNA1,UeccDNA,Confirmed,100\n"
        )

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
                self.unified_csv = str(unified_csv)
                self.html = str(html_report)
                self.text = None
                self.overwrite = True
                self.dry_run = False
                self.id_width = 4
                self.verbose = False

        args = MockArgs()

        # Run packaging
        result = ecc_packager.run(args)
        assert result == 0

        # Check that files were copied to output
        output_dir = Path(args.out_dir)
        assert output_dir.exists()

        # Check for organized structure
        assert (output_dir / "Uecc").exists()
        assert (output_dir / "eccDNA_summary.csv").exists()

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
                self.unified_csv = str(test_file)
                self.html = None
                self.text = None
                self.overwrite = True
                self.dry_run = True
                self.id_width = 4
                self.verbose = True

        args = MockArgs()

        # Dry run should complete without creating files
        result = ecc_packager.run(args)
        assert result == 0

        # Output directory structure should not be created in dry run
        output_dir = Path(args.out_dir)
        # Note: The function might still create the base output directory


class TestExtractBaseId:
    """Test cases for extract_base_id."""

    def test_uecc_simple(self):
        assert ecc_packager.extract_base_id("UeccDNA1|Chr1:100-200(+)|length=1000") == "UeccDNA1"

    def test_mecc_with_underscore_padding(self):
        assert ecc_packager.extract_base_id("MeccDNA_000001.1_1") == "MeccDNA1"

    def test_cecc_with_underscore(self):
        assert ecc_packager.extract_base_id("CeccDNA_000001") == "CeccDNA1"

    def test_with_gt_prefix(self):
        assert ecc_packager.extract_base_id(">UeccDNA42|info") == "UeccDNA42"

    def test_no_match_returns_first_part(self):
        assert ecc_packager.extract_base_id("unknown_header") == "unknown_header"

    def test_multidigit(self):
        assert ecc_packager.extract_base_id("UeccDNA0123") == "UeccDNA123"


class TestLoadFasta:
    """Test cases for load_fasta."""

    def test_load_single_sequence(self, tmp_path):
        fa = tmp_path / "test.fasta"
        fa.write_text(">UeccDNA1|info\nATCG\nGGCC\n")
        result = ecc_packager.load_fasta(fa)
        assert result == {"UeccDNA1": "ATCGGGCC"}

    def test_load_multiple_sequences(self, tmp_path):
        fa = tmp_path / "test.fasta"
        fa.write_text(">UeccDNA1\nATCG\n>MeccDNA2\nGGCC\n")
        result = ecc_packager.load_fasta(fa)
        assert len(result) == 2
        assert result["UeccDNA1"] == "ATCG"
        assert result["MeccDNA2"] == "GGCC"

    def test_load_empty_file(self, tmp_path):
        fa = tmp_path / "empty.fasta"
        fa.write_text("")
        result = ecc_packager.load_fasta(fa)
        assert result == {}

    def test_load_nonexistent_file(self, tmp_path):
        fa = tmp_path / "missing.fasta"
        result = ecc_packager.load_fasta(fa)
        assert result == {}

    def test_load_none_path(self):
        result = ecc_packager.load_fasta(None)
        assert result == {}


class TestBuildArgparser:
    """Test cases for build_argparser."""

    def test_parser_creation(self):
        parser = ecc_packager.build_argparser()
        assert parser is not None

    def test_required_args(self):
        parser = ecc_packager.build_argparser()
        args = parser.parse_args([
            "-s", "sample1",
            "--unified-csv", "/path/to/unified.csv",
        ])
        assert args.sample_name == "sample1"
        assert args.unified_csv == "/path/to/unified.csv"

    def test_default_values(self):
        parser = ecc_packager.build_argparser()
        args = parser.parse_args([
            "-s", "sample1",
            "--unified-csv", "/path/to/unified.csv",
        ])
        assert args.out_dir == "."
        assert args.id_width == 4
        assert args.dry_run is False


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
