"""Integration tests for Cresil adapter utilities."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from circleseeker.modules.cresil_adapter import (
    convert_cresil_to_cyrcular_format,
    parse_cresil_region,
    parse_cresil_regions,
)

pytestmark = pytest.mark.integration


@pytest.mark.parametrize(
    ("region_str", "expected"),
    [
        ("chr2:47242017-47243061_-", ("chr2", 47242017, 47243061, "-")),
        ("chr8:129264432-129269966_+", ("chr8", 129264432, 129269966, "+")),
        ("chr10:131725754-131727108_-", ("chr10", 131725754, 131727108, "-")),
    ],
)
def test_parse_cresil_region(region_str: str, expected: tuple[str, int, int, str]) -> None:
    assert parse_cresil_region(region_str) == expected


def test_parse_cresil_regions_multiple() -> None:
    multi_region = "chr2:47242017-47243061_-;chr2:47250000-47252000_-"
    parsed = parse_cresil_regions(multi_region)

    assert parsed == [
        ("chr2", 47242017, 47243061, "-"),
        ("chr2", 47250000, 47252000, "-"),
    ]


def test_convert_cresil_to_cyrcular_format(tmp_path: Path) -> None:
    cresil_data = {
        "id": ["ec1", "ec2", "ec3", "ec4", "ec5"],
        "merge_region": [
            "chr2:47242017-47243061_-",
            "chr2:128071984-128072363_+",
            "chr8:129264432-129269966_-",
            "chr1:1000-1500_-;chr1:2000-2500_-",
            "chr1:3000-3500_-|chr1:3600-3700_-",
        ],
        "merge_len": [1045, 380, 5535, 1000, 1200],
        "num_region": [1, 1, 1, 2, 3],
        "ctc": [True, True, True, True, True],
        "numreads": [5, 1, 14, 6, 2],
        "totalbase": [52522, 9120, 161471, 48000, 24000],
        "coverage": [50.26, 24.00, 29.17, 48.0, 20.0],
    }

    df = pd.DataFrame(cresil_data)
    test_file = tmp_path / "test_cresil_output.txt"
    df.to_csv(test_file, sep="\t", index=False)

    result_df = convert_cresil_to_cyrcular_format(test_file)

    expected_cols = {
        "circle_id",
        "regions",
        "circle_length",
        "segment_count",
        "num_split_reads",
        "prob_present",
        "prob_artifact",
        "af_nanopore",
    }
    assert expected_cols.issubset(result_df.columns)

    mismatch_row = result_df[result_df["circle_id"] == "ec5"].iloc[0]
    assert mismatch_row["segment_count"] == 2
