from pathlib import Path
import sys

import pandas as pd
import pytest
from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.modules.ecc_dedup import CDHitClusters, eccDedup

DATA_DIR = Path(__file__).parent / "data" / "ecc_dedup"
DATA_AVAILABLE = DATA_DIR.exists()
EXPECTED = DATA_DIR / "expected"

LEGACY_TO_SNAKE = {
    'eccDNA_id': 'eccdna_id',
    'eChr': 'chr',
    'eStart0': 'start_0based',
    'eEnd0': 'end_0based',
    'eStrand': 'strand',
    'eLength': 'length',
    'MatDegree': 'match_degree',
    'copyNum': 'copy_number',
    'eRepeatNum': 'repeat_number',
    'eClass': 'eccdna_type',
    # readName removed - only using read_name now
    'orig_eccDNA_id': 'orig_eccdna_id',
}
SNAKE_TO_LEGACY = {v: k for k, v in LEGACY_TO_SNAKE.items()}


def _prepare(df: pd.DataFrame, sort_cols: list[str]) -> pd.DataFrame:
    if df.empty:
        return df
    return df.sort_values(sort_cols).reset_index(drop=True)


def _load_expected(name: str, sort_cols: list[str]) -> pd.DataFrame:
    df = pd.read_csv(EXPECTED / name)
    return _prepare(df, sort_cols)


def _coerce_str(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    df = df.copy()
    for col in cols:
        if col in df.columns:
            df[col] = df[col].astype(str)
    return df


def _project_to_expected(df: pd.DataFrame, expected_columns: list[str]) -> pd.DataFrame:
    df = df.copy()
    for old in expected_columns:
        if old not in df.columns:
            snake = LEGACY_TO_SNAKE.get(old)
            if snake and snake in df.columns:
                df[old] = df[snake]
    return df[expected_columns]


def _load_actual(path: Path, columns: pd.Index, sort_cols: list[str]) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = _project_to_expected(df, list(columns))
    return _prepare(df, sort_cols)


def _load_fasta(path: Path) -> list[tuple[str, str]]:
    records = list(SeqIO.parse(path, "fasta"))
    records.sort(key=lambda rec: rec.id)
    return [(rec.id, str(rec.seq)) for rec in records]


def _rename_to_legacy(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    for old, new in LEGACY_TO_SNAKE.items():
        if new in df.columns:
            if old in df.columns:
                df = df.drop(columns=[new])
            else:
                df = df.rename(columns={new: old})
    df = df.loc[:, ~df.columns.duplicated()]
    return df


@pytest.mark.skipif(not DATA_AVAILABLE, reason="Test data not found")
def test_ecc_dedup_subset(tmp_path):
    processor = eccDedup()
    results = processor.process_all_types(
        uecc_csv=DATA_DIR / "uecc_input.csv",
        uecc_clstr=DATA_DIR / "uecc.clstr",
        mecc_csv=DATA_DIR / "mecc_input.csv",
        mecc_clstr=DATA_DIR / "mecc.clstr",
        cecc_csv=DATA_DIR / "cecc_input.csv",
        cecc_clstr=DATA_DIR / "cecc.clstr",
        output_dir=tmp_path,
        prefix="step8",
        drop_seq=False,
    )

    # Uecc assertions
    expected_uecc = _load_expected("uecc_processed.csv", ["eccDNA_id"])
    actual_uecc = _prepare(results["Uecc"], ["eccdna_id"])
    assert {"eccdna_id", "chr", "start_0based", "end_0based", "strand", "match_degree", "copy_number", "eccdna_type"}.issubset(actual_uecc.columns)
    assert pd.api.types.is_integer_dtype(actual_uecc['copy_number'])
    assert actual_uecc['match_degree'].dropna().between(0, 100).all()
    actual_uecc_legacy = _rename_to_legacy(actual_uecc)
    actual_uecc_cmp = _coerce_str(actual_uecc_legacy[expected_uecc.columns], ["cluster_id"])
    expected_uecc_cmp = _coerce_str(expected_uecc, ["cluster_id"])
    pd.testing.assert_frame_equal(actual_uecc_cmp, expected_uecc_cmp, check_dtype=False)
    assert (actual_uecc["num_merged"] > 1).any()

    uecc_core = _load_actual(
        tmp_path / "step8_Uecc_C" / "step8_UeccDNA.core.csv",
        pd.read_csv(EXPECTED / "step8_UeccDNA.core.csv").columns,
        ["eccDNA_id"],
    )
    expected_uecc_core = _load_expected("step8_UeccDNA.core.csv", ["eccDNA_id"])
    pd.testing.assert_frame_equal(uecc_core, expected_uecc_core)

    # Mecc assertions
    expected_mecc = _load_expected("mecc_processed.csv", ["eccDNA_id", "q_start"])
    actual_mecc = _prepare(results["Mecc"], ["eccdna_id", "q_start"])
    assert {"eccdna_id", "chr", "start_0based", "end_0based", "strand", "copy_number", "eccdna_type"}.issubset(actual_mecc.columns)
    assert pd.api.types.is_integer_dtype(actual_mecc['copy_number'])
    assert actual_mecc['match_degree'].dropna().between(0, 100).all()
    actual_mecc_cmp = _coerce_str(_rename_to_legacy(actual_mecc)[expected_mecc.columns], ["cluster_id"])
    expected_mecc_cmp = _coerce_str(expected_mecc, ["cluster_id"])
    pd.testing.assert_frame_equal(actual_mecc_cmp, expected_mecc_cmp, check_dtype=False)
    assert (actual_mecc["num_merged"] > 1).any()
    assert (actual_mecc["num_merged"] == 1).any()

    mecc_core_cols = pd.read_csv(EXPECTED / "step8_MeccSites.core.csv").columns
    mecc_core = _load_actual(
        tmp_path / "step8_Mecc_C" / "step8_MeccSites.core.csv",
        mecc_core_cols,
        ["eccDNA_id", "hit_index"],
    )
    expected_mecc_core = _load_expected("step8_MeccSites.core.csv", ["eccDNA_id", "hit_index"])
    pd.testing.assert_frame_equal(mecc_core, expected_mecc_core, check_dtype=False)

    assert _load_fasta(tmp_path / "step8_Mecc_C" / "step8_MeccDNA_C.fasta") == _load_fasta(EXPECTED / "step8_Mecc.fa")
    assert (tmp_path / "step8_Mecc_C" / "step8_MeccBestSite.bed").read_text() == (
        EXPECTED / "step8_MeccBestSite.bed"
    ).read_text()

    # Cecc assertions
    expected_cecc = _load_expected("cecc_processed.csv", ["eccDNA_id", "eChr", "eStart0"])
    actual_cecc = _prepare(results["Cecc"], ["eccdna_id", "chr", "start_0based"])
    assert {"eccdna_id", "chr", "start_0based", "end_0based", "strand", "eccdna_type"}.issubset(actual_cecc.columns)
    assert set(actual_cecc['strand'].dropna().unique()).issubset({'+', '-'})
    actual_cecc_cmp = _coerce_str(_rename_to_legacy(actual_cecc)[expected_cecc.columns], ["cluster_id"])
    expected_cecc_cmp = _coerce_str(expected_cecc, ["cluster_id"])
    pd.testing.assert_frame_equal(actual_cecc_cmp, expected_cecc_cmp, check_dtype=False)
    assert (actual_cecc["num_merged"] > 1).any()
    assert (actual_cecc["num_merged"] == 1).any()

    cecc_core = _load_actual(
        tmp_path / "step8_Cecc_C" / "step8_CeccSegments.core.csv",
        pd.read_csv(EXPECTED / "step8_CeccSegments.core.csv").columns,
        ["eccDNA_id", "seg_index"],
    )
    expected_cecc_core = _load_expected("step8_CeccSegments.core.csv", ["eccDNA_id", "seg_index"])
    pd.testing.assert_frame_equal(cecc_core, expected_cecc_core, check_dtype=False)

    assert _load_fasta(tmp_path / "step8_Cecc_C" / "step8_CeccDNA_C.fasta") == _load_fasta(EXPECTED / "step8_Cecc.fa")
    assert (tmp_path / "step8_Cecc_C" / "step8_CeccJunctions.bedpe").read_text() == (
        EXPECTED / "step8_CeccJunctions.bedpe"
    ).read_text()

    # Uecc BED output exists
    assert (tmp_path / "step8_Uecc_C" / "step8_UeccDNA.bed").read_text() == (
        EXPECTED / "step8_UeccDNA.bed"
    ).read_text()


def test_ecc_dedup_csv_cluster_merging(tmp_path):
    """Ensure CSV cluster mappings are honoured when selecting representatives."""
    cluster_csv = tmp_path / "uecc.id2cluster.csv"
    cluster_csv.write_text(
        "id,cluster,is_representative\n"
        "U1,grp1,false\n"
        "U2,grp1,true\n"
        "U3,grp2,\n"
    )

    clusters = CDHitClusters()
    clusters.parse(cluster_csv)

    assert clusters.rep_of["grp1"] == "U2"
    assert clusters.rep_of["grp2"] == "U3"

    uecc_df = pd.DataFrame(
        [
            {
                "eccDNA_id": "U1",
                "eChr": "chr1",
                "eStart": 10,
                "eEnd": 70,
                "eLength": 60,
                "readName": "readA",
                "eRepeatNum": 1,
                "eStrand": "+",
            },
            {
                "eccDNA_id": "U2",
                "eChr": "chr1",
                "eStart": 20,
                "eEnd": 90,
                "eLength": 70,
                "readName": "readB",
                "eRepeatNum": 2,
                "eStrand": "+",
            },
            {
                "eccDNA_id": "U3",
                "eChr": "chr2",
                "eStart": 5,
                "eEnd": 45,
                "eLength": 40,
                "readName": "readC",
                "eRepeatNum": 3,
                "eStrand": "-",
            },
        ]
    )

    processor = eccDedup()
    processed = processor.process_uecc(uecc_df, clusters, "Uecc")
    processed = processed.sort_values("cluster_id").reset_index(drop=True)

    assert len(processed) == 2

    assert {"eccdna_id", "chr", "start_0based", "end_0based", "strand", "copy_number", "repeat_number", "match_degree"}.issubset(processed.columns)

    cluster1 = processed.loc[processed["cluster_id"] == "grp1"].iloc[0]
    assert cluster1["eccdna_id"] == "U2"
    assert cluster1["reads"] == "readA;readB"
    assert cluster1["repeat_number"] == 3
    assert cluster1["copy_number"] == 3
    assert cluster1["num_merged"] == 2
    assert cluster1["merged_from_ids"] == "U1;U2"
    assert cluster1["orig_eccdna_id"] == "U2"
    assert cluster1["eccdna_type"] == "Uecc"
    assert cluster1["strand"] == "+"

    cluster2 = processed.loc[processed["cluster_id"] == "grp2"].iloc[0]
    assert cluster2["eccdna_id"] == "U3"
    assert cluster2["num_merged"] == 1
    assert cluster2["merged_from_ids"] == "U3"
    assert cluster2["strand"] == "-"
    assert cluster2["copy_number"] == 3
