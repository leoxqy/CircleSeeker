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


def test_tandem_to_ring_process_subset(tmp_path):
    input_tsv = tmp_path / "step1_subset.tsv"

    # Input format matches TideHunter TSV (11 columns, no header)
    rows = [
        # Effective_Length = (1000 / 1000) * 100 = 100 -> CtcR-perfect
        ["readA", 1, 1.0, 1000, 1, 1000, 4, 99.0, 0, "0,1", "ATCG"],
        # Effective_Length = (800 / 1000) * 100 = 80 -> CtcR-hybrid
        ["readB", 1, 1.0, 1000, 1, 800, 4, 99.0, 0, "0,1", "GGCC"],
    ]
    pd.DataFrame(rows).to_csv(input_tsv, sep="\t", header=False, index=False)

    output_csv = tmp_path / "step2_processed.csv"
    output_fasta = tmp_path / "step2_circular.fasta"

    module = TandemToRing(
        input_file=input_tsv,
        output_file=output_csv,
        circular_fasta=output_fasta,
    )

    df_main, df_classification = module.process()

    produced_csv = pd.read_csv(output_csv).sort_values("readName").reset_index(drop=True)
    expected_df = pd.DataFrame(
        {"readName": ["readA", "readB"], "readClass": ["CtcR-perfect", "CtcR-hybrid"]}
    ).sort_values("readName").reset_index(drop=True)

    pd.testing.assert_frame_equal(produced_csv, expected_df)
    pd.testing.assert_frame_equal(
        df_classification.sort_values("readName").reset_index(drop=True), expected_df
    )

    produced_records = list(SeqIO.parse(output_fasta, "fasta"))
    produced_records.sort(key=lambda rec: rec.id)

    expected_fasta = {
        "readA|1|4|1|circular": "ATCGATCG",
        "readB|1|4|1|circular": "GGCCGGCC",
    }
    assert [rec.id for rec in produced_records] == sorted(expected_fasta.keys())
    for rec in produced_records:
        assert str(rec.seq) == expected_fasta[rec.id]

    assert not df_main.empty
    unique_ids = {f"{uid}|circular" for uid in df_main["unique_id"].tolist()}
    assert unique_ids == {rec.id for rec in produced_records}
