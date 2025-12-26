from pathlib import Path
import sys

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.modules.um_classify import UMeccClassifier

DATA_DIR = Path(__file__).parent / "data" / "um_classify"
DATA_AVAILABLE = DATA_DIR.exists()


def _align(df: pd.DataFrame, reference_columns: pd.Index, sort_keys: list[str]) -> pd.DataFrame:
    if df.empty:
        return df
    df = df[reference_columns]
    return df.sort_values(sort_keys).reset_index(drop=True)


def _sorted(df: pd.DataFrame, sort_keys: list[str]) -> pd.DataFrame:
    if df.empty:
        return df
    return df.sort_values(sort_keys).reset_index(drop=True)


@pytest.mark.skipif(not DATA_AVAILABLE, reason="Test data not found")
def test_um_classify_subset(tmp_path):
    blast_tsv = DATA_DIR / "step3_subset.tsv"
    uecc_expected = pd.read_csv(DATA_DIR / "expected_uecc.csv")
    mecc_expected = pd.read_csv(DATA_DIR / "expected_mecc.csv")
    un_expected = pd.read_csv(DATA_DIR / "expected_unclassified.csv")

    classifier = UMeccClassifier()
    uecc_df, mecc_df, unclassified_df = classifier.run(blast_tsv)

    uecc_keys = ["query_id", "eStart", "eEnd"]
    mecc_keys = ["query_id", "eStart", "eEnd"]
    un_keys = ["query_id", "subject_id", "s_start", "s_end"]

    uecc_actual = _align(uecc_df, uecc_expected.columns, uecc_keys)
    mecc_actual = _align(mecc_df, mecc_expected.columns, mecc_keys)
    un_actual = _align(unclassified_df, un_expected.columns, un_keys)

    pd.testing.assert_frame_equal(uecc_actual, _sorted(uecc_expected, uecc_keys))
    pd.testing.assert_frame_equal(mecc_actual, _sorted(mecc_expected, mecc_keys))
    pd.testing.assert_frame_equal(un_actual, _sorted(un_expected, un_keys))
