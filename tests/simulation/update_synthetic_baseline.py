#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TESTS_DIR = ROOT / "tests"
if str(TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(TESTS_DIR))

from simulation.synthetic_validation import run_synthetic_validation  # noqa: E402


def main() -> int:
    baseline_path = Path(__file__).parent / "baselines" / "synthetic_regression.json"
    baseline = json.loads(baseline_path.read_text())

    seeds = baseline["seeds"]
    num_uecc = int(baseline["dataset"]["num_uecc"])
    num_mecc = int(baseline["dataset"]["num_mecc"])
    num_cecc = int(baseline["dataset"]["num_cecc"])

    totals = {"recall": 0.0, "precision": 0.0, "f1": 0.0}
    for seed in seeds:
        with tempfile.TemporaryDirectory() as d:
            _, overall = run_synthetic_validation(
                Path(d),
                seed=int(seed),
                num_uecc=num_uecc,
                num_mecc=num_mecc,
                num_cecc=num_cecc,
            )
            totals["recall"] += overall.recall
            totals["precision"] += overall.precision
            totals["f1"] += overall.f1

    avg = {k: round(v / float(len(seeds)), 6) for k, v in totals.items()}
    baseline.setdefault("expected", {}).setdefault("overall", {}).update(avg)
    baseline_path.write_text(json.dumps(baseline, indent=2) + "\n")
    print(f"Updated baseline: {baseline_path}")
    print(f"Expected overall metrics: {avg}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
