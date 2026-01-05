from pathlib import Path
import sys

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.modules.cecc_build import CeccBuild
from circleseeker.utils.column_standards import ColumnStandard


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
    rows = [
        # Query without genomic overlaps (should be retained)
        {
            "query_id": "q_ok",
            "subject_id": "chr1",
            "q_start": 0,
            "q_end": 100,
            "s_start": 501,
            "s_end": 600,
            "strand": "+",
            "alignment_length": 100,
            "reads": "read_ok",
            "length": 250,
            "copy_number": 1.0,
        },
        {
            "query_id": "q_ok",
            "subject_id": "chr1",
            "q_start": 100,
            "q_end": 200,
            "s_start": 1001,
            "s_end": 1100,
            "strand": "+",
            "alignment_length": 100,
            "reads": "read_ok",
            "length": 250,
            "copy_number": 1.0,
        },
        {
            "query_id": "q_ok",
            "subject_id": "chr1",
            "q_start": 200,
            "q_end": 300,
            "s_start": 1401,
            "s_end": 1500,
            "strand": "+",
            "alignment_length": 100,
            "reads": "read_ok",
            "length": 250,
            "copy_number": 1.0,
        },
        # Closing segment (matches the first rotated segment at q_start=100)
        {
            "query_id": "q_ok",
            "subject_id": "chr1",
            "q_start": 300,
            "q_end": 400,
            "s_start": 1001,
            "s_end": 1100,
            "strand": "+",
            "alignment_length": 100,
            "reads": "read_ok",
            "length": 250,
            "copy_number": 1.0,
        },
        # Query with overlapping genomic segments (should be filtered out)
        {
            "query_id": "q_overlap",
            "subject_id": "chr1",
            "q_start": 0,
            "q_end": 100,
            "s_start": 701,
            "s_end": 800,
            "strand": "+",
            "alignment_length": 100,
            "reads": "read_overlap",
            "length": 250,
            "copy_number": 1.0,
        },
        {
            "query_id": "q_overlap",
            "subject_id": "chr1",
            "q_start": 100,
            "q_end": 200,
            "s_start": 2001,
            "s_end": 2100,
            "strand": "+",
            "alignment_length": 100,
            "reads": "read_overlap",
            "length": 250,
            "copy_number": 1.0,
        },
        {
            "query_id": "q_overlap",
            "subject_id": "chr1",
            "q_start": 200,
            "q_end": 300,
            # Overlaps the previous segment in genomic space
            "s_start": 2051,
            "s_end": 2150,
            "strand": "+",
            "alignment_length": 100,
            "reads": "read_overlap",
            "length": 250,
            "copy_number": 1.0,
        },
        # Closing segment (matches the first rotated segment at q_start=100)
        {
            "query_id": "q_overlap",
            "subject_id": "chr1",
            "q_start": 300,
            "q_end": 400,
            "s_start": 2001,
            "s_end": 2100,
            "strand": "+",
            "alignment_length": 100,
            "reads": "read_overlap",
            "length": 250,
            "copy_number": 1.0,
        },
    ]

    input_csv = tmp_path / "step4_unclassified_subset.csv"
    pd.DataFrame(rows).to_csv(input_csv, index=False)

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

    assert set(result["query_id"].unique()) == {"q_ok"}
    result_sorted = _normalise_roles(_sorted(result, ["query_id", "segment_in_circle"]))
    q_ok = result_sorted[result_sorted["query_id"] == "q_ok"]
    assert q_ok["segment_role"].tolist() == ["head", "middle", "tail"]

    saved = pd.read_csv(output_csv)
    assert set(saved["query_id"].unique()) == {"q_ok"}

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

    assert set(result_no_filter["query_id"].unique()) == {"q_ok", "q_overlap"}
    saved_no_filter = pd.read_csv(output_csv_no_filter)
    assert set(saved_no_filter["query_id"].unique()) == {"q_ok", "q_overlap"}


def test_detect_genomic_overlaps_allows_boundary_jitter():
    builder = CeccBuild()
    segments = pd.DataFrame(
        {
            ColumnStandard.CHR: ["chr1", "chr1"],
            ColumnStandard.START0: [100, 101],
            ColumnStandard.END0: [200, 200],
            ColumnStandard.STRAND: ["+", "+"],
        }
    )

    assert builder.detect_genomic_overlaps_sweepline(segments) is False


def test_detect_genomic_overlaps_ignores_opposite_strands():
    builder = CeccBuild()
    segments = pd.DataFrame(
        {
            ColumnStandard.CHR: ["chr1", "chr1"],
            ColumnStandard.START0: [100, 100],
            ColumnStandard.END0: [200, 200],
            ColumnStandard.STRAND: ["+", "-"],
        }
    )

    assert builder.detect_genomic_overlaps_sweepline(segments) is False


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


def test_cecc_build_fallback_rotation_and_match_degree(tmp_path):
    rows = [
        {
            "query_id": "q_fallback",
            "subject_id": "chr1",
            "q_start": 100,
            "q_end": 300,
            "s_start": 1001,
            "s_end": 1200,
            "strand": "+",
            "alignment_length": 200,
            "reads": "read_fallback",
            "length": 1000,
            "copy_number": 1.0,
        },
        {
            "query_id": "q_fallback",
            "subject_id": "chr2",
            "q_start": 300,
            "q_end": 1000,
            "s_start": 2001,
            "s_end": 2700,
            "strand": "+",
            "alignment_length": 700,
            "reads": "read_fallback",
            "length": 1000,
            "copy_number": 1.0,
        },
    ]

    input_csv = tmp_path / "step4_unclassified_fallback.csv"
    pd.DataFrame(rows).to_csv(input_csv, index=False)

    output_csv = tmp_path / "step5_cecc_fallback.csv"
    builder = CeccBuild()
    result = builder.run_pipeline(
        input_csv=input_csv,
        output_csv=output_csv,
        overlap_threshold=0.95,
        min_segments=2,
        edge_tolerance=10,
        position_tolerance=50,
        min_match_degree=90.0,
        max_rotations=10,
    )

    assert set(result["query_id"].unique()) == {"q_fallback"}
    result_sorted = _normalise_roles(_sorted(result, ["query_id", "segment_in_circle"]))
    q_fallback = result_sorted[result_sorted["query_id"] == "q_fallback"]
    assert q_fallback["segment_role"].tolist() == ["head", "tail"]
    assert q_fallback["CeccClass"].nunique() == 1
    assert q_fallback["CeccClass"].iloc[0] == "Cecc-InterChr"
