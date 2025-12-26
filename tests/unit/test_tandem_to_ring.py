from pathlib import Path
import sys

import pandas as pd
import pytest
from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.modules.tandem_to_ring import TandemToRing


DATA_DIR = Path(__file__).parent / "data" / "tandem_to_ring"
DATA_AVAILABLE = DATA_DIR.exists()


@pytest.mark.skipif(not DATA_AVAILABLE, reason="Test data not found")
def test_tandem_to_ring_process_subset(tmp_path):
    input_tsv = DATA_DIR / "step1_subset.tsv"
    expected_csv = DATA_DIR / "expected_classification.csv"
    expected_fasta = DATA_DIR / "expected_circular.fasta"

    output_csv = tmp_path / "step2_processed.csv"
    output_fasta = tmp_path / "step2_circular.fasta"

    module = TandemToRing(
        input_file=input_tsv,
        output_file=output_csv,
        circular_fasta=output_fasta,
    )

    df_main, df_classification = module.process()

    produced_csv = pd.read_csv(output_csv).sort_values("readName").reset_index(drop=True)
    expected_df = pd.read_csv(expected_csv).sort_values("readName").reset_index(drop=True)

    pd.testing.assert_frame_equal(produced_csv, expected_df)
    pd.testing.assert_frame_equal(
        df_classification.sort_values("readName").reset_index(drop=True), expected_df
    )

    produced_records = list(SeqIO.parse(output_fasta, "fasta"))
    expected_records = list(SeqIO.parse(expected_fasta, "fasta"))

    assert len(produced_records) == len(expected_records)
    produced_records.sort(key=lambda rec: rec.id)
    expected_records.sort(key=lambda rec: rec.id)

    for prod, exp in zip(produced_records, expected_records):
        assert prod.id == exp.id
        assert str(prod.seq) == str(exp.seq)

    assert not df_main.empty
    unique_ids = {f"{uid}|circular" for uid in df_main["unique_id"].tolist()}
    assert unique_ids == {rec.id for rec in produced_records}
