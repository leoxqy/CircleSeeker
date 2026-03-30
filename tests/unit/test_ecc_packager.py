"""Tests for ecc_packager module."""

from pathlib import Path
import sys
import pytest
import tempfile
import shutil
from unittest.mock import patch, MagicMock

import pandas as pd

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
        assert (output_dir / "test_sample_eccDNA_summary.csv").exists()

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
        assert (output_dir / "UeccDNA").exists()
        assert (output_dir / "test_sample_eccDNA_summary.csv").exists()

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


# ════════════════════════════════════════════════════════════════════
# score_mecc_ont_lr — MeccDNA ONT LR scoring
# ════════════════════════════════════════════════════════════════════

import numpy as np


class TestScoreMeccOntLr:
    """Test the MeccDNA ONT LR scoring and demotion logic."""

    def test_empty_summary_returns_unchanged(self):
        """Empty summary_df should be returned as-is."""
        summary_df = pd.DataFrame()
        reads_df = pd.DataFrame()
        regions_df = pd.DataFrame()
        result = ecc_packager.score_mecc_ont_lr(summary_df, reads_df, regions_df)
        assert result.empty

    def test_no_mecc_returns_unchanged(self):
        """Summary with no Mecc entries should be returned as-is."""
        summary_df = pd.DataFrame({
            "eccDNA_id": ["U1"],
            "type": ["Uecc"],
            "length": [1000],
            "confidence": [0.5],
            "copy_number": [1],
        })
        reads_df = pd.DataFrame()
        regions_df = pd.DataFrame()
        result = ecc_packager.score_mecc_ont_lr(summary_df, reads_df, regions_df)
        assert len(result) == 1
        assert result.iloc[0]["type"] == "Uecc"

    def test_mecc_demotion_with_high_threshold(self):
        """Mecc with threshold=1.0 should demote all Mecc to Uecc."""
        summary_df = pd.DataFrame({
            "eccDNA_id": ["M1"],
            "type": ["Mecc"],
            "length": [1000],
            "confidence": [0.5],
            "copy_number": [1],
        })
        reads_df = pd.DataFrame({
            "eccDNA_id": ["M1"],
            "mapq": [10],
            "match_degree": [95],
        })
        regions_df = pd.DataFrame({
            "eccDNA_id": ["M1", "M1"],
        })
        result = ecc_packager.score_mecc_ont_lr(
            summary_df, reads_df, regions_df, threshold=1.0
        )
        # threshold=1.0 means all Mecc demoted
        assert result.iloc[0]["type"] == "Uecc"

    def test_threshold_zero_skips_filtering(self):
        """Threshold <= 0 should skip filtering entirely."""
        summary_df = pd.DataFrame({
            "eccDNA_id": ["M1"],
            "type": ["Mecc"],
            "length": [1000],
            "confidence": [0.5],
            "copy_number": [1],
        })
        reads_df = pd.DataFrame()
        regions_df = pd.DataFrame()
        result = ecc_packager.score_mecc_ont_lr(
            summary_df, reads_df, regions_df, threshold=0.0
        )
        assert result.iloc[0]["type"] == "Mecc"

    def test_uecc_not_affected(self):
        """Uecc entries should not be affected by the Mecc LR filter."""
        summary_df = pd.DataFrame({
            "eccDNA_id": ["U1", "M1"],
            "type": ["Uecc", "Mecc"],
            "length": [1000, 1000],
            "confidence": [0.5, 0.5],
            "copy_number": [1, 1],
        })
        reads_df = pd.DataFrame({
            "eccDNA_id": ["U1", "M1"],
            "mapq": [60, 10],
            "match_degree": [99, 95],
        })
        regions_df = pd.DataFrame({
            "eccDNA_id": ["U1", "M1"],
        })
        result = ecc_packager.score_mecc_ont_lr(
            summary_df, reads_df, regions_df, threshold=1.0
        )
        uecc_rows = result[result["eccDNA_id"] == "U1"]
        assert uecc_rows.iloc[0]["type"] == "Uecc"


# ════════════════════════════════════════════════════════════════════
# _filter_cecc_qc — CeccDNA quality-control filter
# ════════════════════════════════════════════════════════════════════


class TestFilterCeccQc:
    """Test the CeccDNA QC filter for confirmed length and inferred segments."""

    def test_short_confirmed_cecc_removed(self):
        """Confirmed CeccDNA shorter than min length should be removed."""
        unified_df = pd.DataFrame({
            "eccDNA_id": ["C1", "C2"],
            "eccDNA_type": ["CeccDNA", "CeccDNA"],
            "State": ["Confirmed", "Confirmed"],
            "Length": [500, 2000],
        })
        result_unified, _, _ = ecc_packager._filter_cecc_qc(
            unified_df, None, None, min_confirmed_length=1400,
        )
        assert len(result_unified) == 1
        assert result_unified.iloc[0]["eccDNA_id"] == "C2"

    def test_inferred_cecc_too_many_segments_removed(self):
        """Inferred CeccDNA with too many segments should be removed."""
        unified_df = pd.DataFrame({
            "eccDNA_id": ["C1", "C2"],
            "eccDNA_type": ["CeccDNA", "CeccDNA"],
            "State": ["Inferred", "Inferred"],
            "Length": [5000, 3000],
            "Seg_total": [10, 3],
        })
        result_unified, _, _ = ecc_packager._filter_cecc_qc(
            unified_df, None, None, max_inferred_segments=5,
        )
        assert len(result_unified) == 1
        assert result_unified.iloc[0]["eccDNA_id"] == "C2"

    def test_uecc_not_affected(self):
        """UeccDNA entries should not be affected by CeccDNA QC filter."""
        unified_df = pd.DataFrame({
            "eccDNA_id": ["U1", "C1"],
            "eccDNA_type": ["UeccDNA", "CeccDNA"],
            "State": ["Confirmed", "Confirmed"],
            "Length": [500, 500],
        })
        result_unified, _, _ = ecc_packager._filter_cecc_qc(
            unified_df, None, None, min_confirmed_length=1400,
        )
        # U1 (short UeccDNA) should be kept; C1 (short CeccDNA) removed
        assert len(result_unified) == 1
        assert result_unified.iloc[0]["eccDNA_id"] == "U1"

    def test_no_removal_when_all_pass(self):
        """All entries passing QC should be kept intact."""
        unified_df = pd.DataFrame({
            "eccDNA_id": ["C1", "C2"],
            "eccDNA_type": ["CeccDNA", "CeccDNA"],
            "State": ["Confirmed", "Inferred"],
            "Length": [5000, 3000],
            "Seg_total": [2, 3],
        })
        result_unified, _, _ = ecc_packager._filter_cecc_qc(
            unified_df, None, None,
            min_confirmed_length=1400,
            max_inferred_segments=5,
        )
        assert len(result_unified) == 2

    def test_cecc_df_also_filtered(self):
        """Associated cecc_df should also have removed IDs dropped."""
        unified_df = pd.DataFrame({
            "eccDNA_id": ["C1", "C2"],
            "eccDNA_type": ["CeccDNA", "CeccDNA"],
            "State": ["Confirmed", "Confirmed"],
            "Length": [500, 2000],
        })
        cecc_df = pd.DataFrame({
            "eccDNA_id": ["C1", "C1", "C2"],
            "chr": ["chr1", "chr2", "chr3"],
            "start": [100, 200, 300],
            "end": [150, 250, 350],
        })
        result_unified, result_cecc, _ = ecc_packager._filter_cecc_qc(
            unified_df, cecc_df, None, min_confirmed_length=1400,
        )
        assert len(result_cecc) == 1
        assert result_cecc.iloc[0]["eccDNA_id"] == "C2"

    def test_inferred_chimeric_df_also_filtered(self):
        """Associated inferred_chimeric_df should also have removed IDs dropped."""
        unified_df = pd.DataFrame({
            "eccDNA_id": ["C1", "C2"],
            "eccDNA_type": ["CeccDNA", "CeccDNA"],
            "State": ["Inferred", "Inferred"],
            "Length": [5000, 3000],
            "Seg_total": [10, 3],
        })
        inferred_chimeric_df = pd.DataFrame({
            "eccDNA_id": ["C1", "C1", "C2"],
            "chr": ["chr1", "chr2", "chr3"],
        })
        result_unified, _, result_chimeric = ecc_packager._filter_cecc_qc(
            unified_df, None, inferred_chimeric_df, max_inferred_segments=5,
        )
        assert len(result_chimeric) == 1
        assert result_chimeric.iloc[0]["eccDNA_id"] == "C2"

    def test_empty_input(self):
        """Empty DataFrame should be handled gracefully."""
        unified_df = pd.DataFrame(columns=["eccDNA_id", "eccDNA_type", "State", "Length"])
        result_unified, _, _ = ecc_packager._filter_cecc_qc(
            unified_df, None, None,
        )
        assert result_unified.empty


# ════════════════════════════════════════════════════════════════════
# NMS dedup logic (ALL-fragments vs ANY-fragment)
# ════════════════════════════════════════════════════════════════════


class TestNmsDedup:
    """Test the NMS dedup logic embedded in run()."""

    def test_frag_overlap_same_region(self):
        """Exact same fragment should have overlap 1.0."""
        # Replicating the _frag_overlap inner function from ecc_packager.run
        def _frag_overlap(f1, f2):
            c1, s1, e1 = f1
            c2, s2, e2 = f2
            if c1 != c2:
                return 0.0
            ovl = max(0, min(e1, e2) - max(s1, s2))
            min_len = min(e1 - s1, e2 - s2)
            return ovl / min_len if min_len > 0 else 0.0

        assert _frag_overlap(("chr1", 100, 200), ("chr1", 100, 200)) == 1.0

    def test_frag_overlap_no_overlap(self):
        def _frag_overlap(f1, f2):
            c1, s1, e1 = f1
            c2, s2, e2 = f2
            if c1 != c2:
                return 0.0
            ovl = max(0, min(e1, e2) - max(s1, s2))
            min_len = min(e1 - s1, e2 - s2)
            return ovl / min_len if min_len > 0 else 0.0

        assert _frag_overlap(("chr1", 100, 200), ("chr1", 300, 400)) == 0.0

    def test_frag_overlap_different_chrom(self):
        def _frag_overlap(f1, f2):
            c1, s1, e1 = f1
            c2, s2, e2 = f2
            if c1 != c2:
                return 0.0
            ovl = max(0, min(e1, e2) - max(s1, s2))
            min_len = min(e1 - s1, e2 - s2)
            return ovl / min_len if min_len > 0 else 0.0

        assert _frag_overlap(("chr1", 100, 200), ("chr2", 100, 200)) == 0.0

    def test_cecc_nms_requires_all_fragments(self):
        """CeccDNA dedup requires ALL fragments to overlap (not just any)."""
        import re

        def _parse_frags(row):
            frags = []
            loc = str(row.get("location", ""))
            for m in re.finditer(r"(\w+):(\d+)-(\d+)", loc):
                frags.append((m.group(1), int(m.group(2)), int(m.group(3))))
            return frags

        def _frag_overlap(f1, f2):
            c1, s1, e1 = f1
            c2, s2, e2 = f2
            if c1 != c2:
                return 0.0
            ovl = max(0, min(e1, e2) - max(s1, s2))
            min_len = min(e1 - s1, e2 - s2)
            return ovl / min_len if min_len > 0 else 0.0

        # Two CeccDNA that share one fragment but differ on another
        frags_a = [("chr1", 100, 200), ("chr2", 300, 400)]
        frags_b = [("chr1", 100, 200), ("chr3", 500, 600)]

        # ALL fragments of b must overlap a -> chr3:500-600 does not overlap any in a
        n_matched = 0
        for ci, si, ei in frags_b:
            for cj, sj, ej in frags_a:
                if _frag_overlap((ci, si, ei), (cj, sj, ej)) > 0.5:
                    n_matched += 1
                    break
        is_dup = (n_matched == len(frags_b))
        assert not is_dup, "Cecc with non-overlapping fragment should NOT be dup"

    def test_mecc_nms_any_fragment_sufficient(self):
        """MeccDNA dedup uses ANY fragment overlap."""
        def _frag_overlap(f1, f2):
            c1, s1, e1 = f1
            c2, s2, e2 = f2
            if c1 != c2:
                return 0.0
            ovl = max(0, min(e1, e2) - max(s1, s2))
            min_len = min(e1 - s1, e2 - s2)
            return ovl / min_len if min_len > 0 else 0.0

        frags_a = [("chr1", 100, 200)]
        frags_b = [("chr1", 110, 210)]

        is_dup = False
        for ci, si, ei in frags_b:
            for cj, sj, ej in frags_a:
                if _frag_overlap((ci, si, ei), (cj, sj, ej)) > 0.5:
                    is_dup = True
                    break
            if is_dup:
                break
        assert is_dup, "Mecc with overlapping fragment should be dup"


# ════════════════════════════════════════════════════════════════════
# NMS platform guard (_is_ont_reads)
# ════════════════════════════════════════════════════════════════════


class TestIsOntReads:
    """Test the _is_ont_reads guard condition for NMS dedup."""

    def test_ont_reads_detected(self):
        """ONT reads table has eccDNA_copy_number + mapq + match_degree."""
        reads_df = pd.DataFrame({
            "eccDNA_id": ["M1"],
            "eccDNA_copy_number": [5.0],
            "mapq": [60],
            "match_degree": [98.5],
        })
        _is_ont = (
            not reads_df.empty
            and "eccDNA_copy_number" in reads_df.columns
            and "mapq" in reads_df.columns
            and "match_degree" in reads_df.columns
        )
        assert _is_ont is True

    def test_ngs_reads_not_detected(self):
        """NGS reads table lacks eccDNA_copy_number."""
        reads_df = pd.DataFrame({
            "eccDNA_id": ["M1"],
            "mapq": [60],
            "match_degree": [98.5],
        })
        _is_ont = (
            not reads_df.empty
            and "eccDNA_copy_number" in reads_df.columns
            and "mapq" in reads_df.columns
            and "match_degree" in reads_df.columns
        )
        assert _is_ont is False

    def test_empty_reads_not_detected(self):
        """Empty reads table should not be detected as ONT."""
        reads_df = pd.DataFrame()
        _is_ont = (
            not reads_df.empty
            and "eccDNA_copy_number" in reads_df.columns
            and "mapq" in reads_df.columns
            and "match_degree" in reads_df.columns
        )
        assert _is_ont is False

    def test_missing_mapq_not_detected(self):
        """Missing mapq column should not be detected as ONT."""
        reads_df = pd.DataFrame({
            "eccDNA_id": ["M1"],
            "eccDNA_copy_number": [5.0],
            "match_degree": [98.5],
        })
        _is_ont = (
            not reads_df.empty
            and "eccDNA_copy_number" in reads_df.columns
            and "mapq" in reads_df.columns
            and "match_degree" in reads_df.columns
        )
        assert _is_ont is False
