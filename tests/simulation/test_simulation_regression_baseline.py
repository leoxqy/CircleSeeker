from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from .synthetic_validation import run_synthetic_validation


def test_synthetic_simulation_regression_baseline(tmp_path):
    baseline_path = Path(__file__).parent / "baselines" / "synthetic_regression.json"
    baseline = json.loads(baseline_path.read_text())

    seeds = baseline["seeds"]
    num_uecc = int(baseline["dataset"]["num_uecc"])
    num_mecc = int(baseline["dataset"]["num_mecc"])
    num_cecc = int(baseline["dataset"]["num_cecc"])
    max_drop = float(baseline["tolerance"]["max_drop"])

    expected = baseline["expected"]["overall"]

    totals = {"recall": 0.0, "precision": 0.0, "f1": 0.0}
    for seed in seeds:
        _, overall = run_synthetic_validation(
            tmp_path / f"seed_{seed}",
            seed=int(seed),
            num_uecc=num_uecc,
            num_mecc=num_mecc,
            num_cecc=num_cecc,
        )
        totals["recall"] += overall.recall
        totals["precision"] += overall.precision
        totals["f1"] += overall.f1

    avg = {k: v / float(len(seeds)) for k, v in totals.items()}

    assert avg["recall"] >= float(expected["recall"]) - max_drop
    assert avg["precision"] >= float(expected["precision"]) - max_drop
    assert avg["f1"] >= float(expected["f1"]) - max_drop


@pytest.mark.parametrize("n_reads", [50])
def test_negative_control_low_mapq_has_no_high_confidence_calls(tmp_path, n_reads: int):
    """Negative control: ambiguous full-coverage alignments should not produce high-confidence U calls."""
    from circleseeker.modules.um_classify import UMeccClassifier

    # Synthetic alignment table: full coverage, but MAPQ=0 for every read.
    alignment_tsv = tmp_path / "neg_control.tsv"
    from .synthetic_validation import write_negative_control_alignment_tsv

    write_negative_control_alignment_tsv(alignment_tsv, n_reads=n_reads, mapq=0)

    classifier = UMeccClassifier(mapq_u_min=0)
    uecc_df, mecc_df, _ = classifier.run(alignment_tsv)

    assert mecc_df.empty
    assert uecc_df["query_id"].nunique() == n_reads
    assert {"confidence_score", "low_mapq"}.issubset(uecc_df.columns)

    # High-confidence false positives should be bounded (ideally zero).
    high_conf = uecc_df[pd.to_numeric(uecc_df["confidence_score"], errors="coerce") >= 0.8]
    assert high_conf.empty
    assert uecc_df["low_mapq"].all()
