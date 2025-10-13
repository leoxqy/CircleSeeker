from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.modules.cecc_build import CeccBuild
from circleseeker.utils.column_standards import ColumnStandard

DATA_DIR = Path(__file__).parent / "data" / "cecc_build"


def _sorted(df: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    if df.empty:
        return df.copy()
    return df.sort_values(keys).reset_index(drop=True)


def _normalise_roles(df: pd.DataFrame) -> pd.DataFrame:
    if "segment_role" in df.columns:
        df = df.copy()
        df["segment_role"] = df["segment_role"].astype(str)
    return df


def test_cecc_build_subset(tmp_path, monkeypatch):
    input_csv = DATA_DIR / "step4_unclassified_subset.csv"
    expected_filtered = pd.read_csv(DATA_DIR / "expected_cecc_filtered.csv")
    expected_no_filter = pd.read_csv(DATA_DIR / "expected_cecc_no_filter.csv")

    # Default run should drop queries with overlapping genomic segments
    output_csv = tmp_path / "step5_cecc.csv"
    builder = CeccBuild()
    result = builder.run_pipeline(
        input_csv=input_csv,
        output_csv=output_csv,
        overlap_threshold=0.8,
        min_segments=2,
        edge_tolerance=10,
        position_tolerance=100,
    )

    expected_filtered_sorted = _sorted(
        expected_filtered, ["query_id", "segment_in_circle"]
    )
    actual_filtered_sorted = _normalise_roles(
        _sorted(result[expected_filtered.columns], ["query_id", "segment_in_circle"])
    )
    expected_filtered_sorted = _normalise_roles(expected_filtered_sorted)
    pd.testing.assert_frame_equal(actual_filtered_sorted, expected_filtered_sorted)

    saved_filtered_sorted = _normalise_roles(
        _sorted(
            pd.read_csv(output_csv)[expected_filtered.columns],
            ["query_id", "segment_in_circle"],
        )
    )
    pd.testing.assert_frame_equal(saved_filtered_sorted, expected_filtered_sorted)

    # Disabling the new overlap filter should retain the original dataset
    output_csv_no_filter = tmp_path / "step5_cecc_no_filter.csv"
    builder_no_filter = CeccBuild()
    monkeypatch.setattr(builder_no_filter, "filter_overlapping_queries", lambda df: df)
    result_no_filter = builder_no_filter.run_pipeline(
        input_csv=input_csv,
        output_csv=output_csv_no_filter,
        overlap_threshold=0.8,
        min_segments=2,
        edge_tolerance=10,
        position_tolerance=100,
    )

    expected_no_filter_sorted = _sorted(
        expected_no_filter, ["query_id", "segment_in_circle"]
    )
    actual_no_filter_sorted = _normalise_roles(
        _sorted(
            result_no_filter[expected_no_filter.columns],
            ["query_id", "segment_in_circle"],
        )
    )
    expected_no_filter_sorted = _normalise_roles(expected_no_filter_sorted)
    pd.testing.assert_frame_equal(actual_no_filter_sorted, expected_no_filter_sorted)

    saved_no_filter_sorted = _normalise_roles(
        _sorted(
            pd.read_csv(output_csv_no_filter)[expected_no_filter.columns],
            ["query_id", "segment_in_circle"],
        )
    )
    pd.testing.assert_frame_equal(saved_no_filter_sorted, expected_no_filter_sorted)


def test_detect_circles_requires_two_segments(monkeypatch):
    builder = CeccBuild()
    df = pd.DataFrame(
        {
            "query_id": ["q1"],
            ColumnStandard.READS: ["read1"],
            ColumnStandard.CHR: ["chr1"],
            ColumnStandard.START0: [100],
            ColumnStandard.END0: [200],
            ColumnStandard.STRAND: ["+"],
            "q_start": [0],
            "q_end": [100],
            "alignment_length": [100],
            ColumnStandard.LENGTH: [100],
            ColumnStandard.COPY_NUMBER: [1.0],
        }
    )

    def fake_find_circular(group, edge_tol, pos_tol):
        return {
            "path": [0],
            "closing_at": 0,
            "cum_len": 100.0,
            "mat_degree": 100.0,
        }

    monkeypatch.setattr(builder, "find_circular", fake_find_circular)
    result = builder.detect_circles(df, edge_tol=20, pos_tol=50)
    assert result.empty, "Cecc candidates with fewer than two segments must be discarded"
