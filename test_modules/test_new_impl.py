#!/usr/bin/env python3
"""Module-local test for Trapeze (step 5).

Runs the new Trapeze implementation on test_modules/step05_trapeze/input,
writes outputs to output/, and compares against expected/ using lightweight checks.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from circleseeker2.modules.trapeze import Trapeze
from test_modules._common.compare_utils import (
    file_nonempty,
    compare_csv_row_counts,
)


def main() -> int:
    base = Path(__file__).parent
    ip = base / "input"
    op = base / "output"
    ex = base / "expected"
    op.mkdir(exist_ok=True)

    # Prefer unclassified alignments as input (as in legacy step5)
    in_csv = ip / "step4_unclassified.csv"
    if not in_csv.exists():
        # Fallback: if test data already provides a cecc-like CSV in input
        alt = ip / "step5_cecc.csv"
        if alt.exists():
            in_csv = alt
    if not in_csv.exists():
        print(f"Missing input: {ip}/step4_unclassified.csv (or step5_cecc.csv)")
        return 1

    out_csv = op / "step5_cecc.csv"
    tr = Trapeze()

    try:
        tr.run_pipeline(
            input_csv=in_csv,
            output_csv=out_csv,
            overlap_threshold=0.8,
            min_segments=2,
            edge_tolerance=10,
            position_tolerance=100,
        )
    except Exception as e:
        print(f"Trapeze failed: {e}")
        return 1

    if not file_nonempty(out_csv):
        print("ERROR: No cecc output produced")
        return 1

    exp_csv = ex / "step5_cecc.csv"
    if exp_csv.exists():
        ok, ca, cb = compare_csv_row_counts(out_csv, exp_csv, tolerance=0.2)
        print(f"Rows new/exp={ca}/{cb}, within 20%: {ok}")
        return 0 if ok else 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

