from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import pytest
import pandas as pd

from circleseeker.modules.splitreads_adapter import (
    parse_splitreads_region,
    parse_splitreads_regions,
    convert_splitreads_to_overview_format,
)


# ---------------------------------------------------------------------------
# TestParseSplitreadsRegion
# ---------------------------------------------------------------------------
class TestParseSplitreadsRegion:
    """Tests for parse_splitreads_region()."""

    def test_basic_negative_strand(self):
        result = parse_splitreads_region("chr2:47242017-47243061_-")
        assert result == ("chr2", 47242017, 47243061, "-")

    def test_positive_strand(self):
        result = parse_splitreads_region("chr1:100-200_+")
        assert result == ("chr1", 100, 200, "+")

    def test_no_strand_defaults_plus(self):
        result = parse_splitreads_region("chr1:100-200")
        assert result == ("chr1", 100, 200, "+")

    def test_empty_string_raises(self):
        with pytest.raises(ValueError):
            parse_splitreads_region("")

    def test_missing_colon_raises(self):
        with pytest.raises(ValueError):
            parse_splitreads_region("chr1_100-200")


# ---------------------------------------------------------------------------
# TestParseSplitreadsRegions
# ---------------------------------------------------------------------------
class TestParseSplitreadsRegions:
    """Tests for parse_splitreads_regions()."""

    def test_single_region(self):
        result = parse_splitreads_regions("chr1:100-200_+")
        assert len(result) == 1
        assert result[0] == ("chr1", 100, 200, "+")

    def test_multiple_semicolon_separated(self):
        result = parse_splitreads_regions("chr1:100-200_+;chr2:300-400_-")
        assert len(result) == 2
        assert result[0] == ("chr1", 100, 200, "+")
        assert result[1] == ("chr2", 300, 400, "-")

    def test_nan_input_raises(self):
        with pytest.raises(ValueError):
            parse_splitreads_regions(float("nan"))


# ---------------------------------------------------------------------------
# TestConvertSplitreadsToOverviewFormat
# ---------------------------------------------------------------------------
class TestConvertSplitreadsToOverviewFormat:
    """Tests for convert_splitreads_to_overview_format()."""

    def test_valid_tsv_input(self, tmp_path):
        tsv_file = tmp_path / "eccDNA_final.txt"
        tsv_file.write_text(
            "id\tmerge_region\tmerge_len\tnumreads\tcoverage\tnum_region\ttotalbase\tctc\n"
            "ecc_1\tchr1:100-200_+\t100\t5\t50.0\t1\t500\tTrue\n"
            "ecc_2\tchr2:300-500_-\t200\t10\t80.0\t1\t1600\tTrue\n"
        )
        df = convert_splitreads_to_overview_format(tsv_file)
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 2
        assert "circle_id" in df.columns
        assert "regions" in df.columns
        assert "circle_length" in df.columns
        assert "num_split_reads" in df.columns
        assert df.iloc[0]["circle_id"] == "ecc_1"
        assert df.iloc[0]["circle_length"] == 100

    def test_nonexistent_file_raises(self, tmp_path):
        fake_path = tmp_path / "does_not_exist.txt"
        with pytest.raises(FileNotFoundError):
            convert_splitreads_to_overview_format(fake_path)
