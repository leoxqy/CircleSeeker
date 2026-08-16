"""Regression test for the bg_end lexicographic-max bug in splitreads_core.

The coverage table is read with ``dtype=str``. Before the fix, ``bg_end`` was left
as a string while ``bg_start`` was converted, so the ``groupby().agg({"bg_end": "max"})``
compared lexicographically. Any merged region spanning a power-of-ten boundary got
its end truncated to the largest value with fewer digits::

    max("99999", "403001") == "99999"      # '9' > '4'

On the arabidopsis 10X benchmark this silently dropped Chr3:72970-403001
(truth UeccDNA_001569, 330,031 bp) -- its length collapsed to 27,029 bp, a 91.8%
loss. Long eccDNA are disproportionately affected because a longer region is more
likely to span a power-of-ten boundary.
"""

from pathlib import Path
import sys

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


COLUMNS = ["m_chrom", "m_start", "m_end", "bg_chrom", "bg_start", "bg_end", "depth", "d_ovl"]

# One merged region whose coverage intervals straddle 100,000. Written as strings
# to mirror the ``dtype=str`` read in _run_identify.
STRADDLING_REGION = [
    ("Chr3", "72970", "403001", "Chr3", "72970", "99999", "12", "27029"),
    ("Chr3", "72970", "403001", "Chr3", "99999", "403001", "11", "303002"),
]


def _numeric_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Reproduce the (fixed) conversion performed in _run_identify."""
    df = df.copy()
    df[["depth", "bg_start", "bg_end", "d_ovl"]] = df[
        ["depth", "bg_start", "bg_end", "d_ovl"]
    ].apply(pd.to_numeric, errors="coerce")
    return df


def _aggregate(df: pd.DataFrame) -> pd.DataFrame:
    """Reproduce the groupby-agg in _run_identify."""
    df = df.copy()
    df["mergeid"] = df["m_chrom"] + "_" + df["m_start"] + "_" + df["m_end"]
    out = (
        df.groupby(["mergeid"])
        .agg({"bg_chrom": "max", "bg_start": "min", "bg_end": "max", "depth": "mean"})
        .reset_index()
    )
    out[["bg_start", "bg_end"]] = out[["bg_start", "bg_end"]].apply(pd.to_numeric, errors="coerce")
    out["length"] = out["bg_end"] - out["bg_start"]
    return out


def test_string_bg_end_truncates_region_spanning_power_of_ten():
    """Guards the premise: without conversion, agg('max') picks the wrong end."""
    df = pd.DataFrame(STRADDLING_REGION, columns=COLUMNS)
    # Convert everything except bg_end -- the pre-fix behaviour.
    df[["depth", "bg_start", "d_ovl"]] = df[["depth", "bg_start", "d_ovl"]].apply(
        pd.to_numeric, errors="coerce"
    )
    broken = _aggregate(df)
    assert broken.loc[0, "bg_end"] == 99999
    assert broken.loc[0, "length"] == 27029


def test_numeric_bg_end_preserves_region_spanning_power_of_ten():
    """The fix: bg_end is numeric before the aggregation, so max() is arithmetic."""
    fixed = _aggregate(_numeric_columns(pd.DataFrame(STRADDLING_REGION, columns=COLUMNS)))
    assert fixed.loc[0, "bg_end"] == 403001
    assert fixed.loc[0, "length"] == 330031


@pytest.mark.parametrize(
    "ends, expected",
    [
        (["985", "5386"], 5386),          # 3 vs 4 digits
        (["9758", "26498"], 26498),       # 4 vs 5 digits
        (["99999", "403001"], 403001),    # 5 vs 6 digits
        (["9999999", "10198486"], 10198486),  # 7 vs 8 digits
    ],
)
def test_boundary_crossings_of_every_magnitude(ends, expected):
    rows = [("c", "0", "9" * 9, "c", "0", e, "10", "1") for e in ends]
    fixed = _aggregate(_numeric_columns(pd.DataFrame(rows, columns=COLUMNS)))
    assert fixed.loc[0, "bg_end"] == expected


def test_negative_length_would_drop_the_circle_entirely():
    """When the lexicographic max lands below bg_start the length goes negative,
    and the region is silently removed by the `length >= min_region_size` filter."""
    rows = [
        ("c", "1000", "5386", "c", "1000", "985", "10", "1"),
        ("c", "1000", "5386", "c", "1000", "5386", "10", "1"),
    ]
    df = pd.DataFrame(rows, columns=COLUMNS)
    df[["depth", "bg_start", "d_ovl"]] = df[["depth", "bg_start", "d_ovl"]].apply(
        pd.to_numeric, errors="coerce"
    )
    assert _aggregate(df).loc[0, "length"] < 0
    assert _aggregate(_numeric_columns(pd.DataFrame(rows, columns=COLUMNS))).loc[0, "length"] == 4386


def test_source_converts_bg_end_before_aggregating():
    """Pin the actual source: bg_end must appear in the to_numeric that precedes the agg."""
    src = (SRC / "circleseeker" / "modules" / "splitreads_core.py").read_text()
    marker = 'merge_genomecov_df[["depth", "bg_start", "bg_end", "d_ovl"]]'
    assert marker in src, "bg_end is not converted to numeric before the groupby-agg"
    assert src.index(marker) < src.index('"bg_end": "max"')
