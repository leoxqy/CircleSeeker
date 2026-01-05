# Simulation Validation & Recall Benchmark (U/M/C)

This repository includes an end-to-end **synthetic** validation pipeline to regression-test the **candidate alignment → U/M classification → C detection** path, and to report read-level recall/precision/F1.

## Components

- `tests/simulation/eccdna_simulator.py`: generates a toy reference genome + U/Mecc/Cecc reads + ground-truth CSVs
- `tests/simulation/run_validation.py`: runs *simulation → minimap2 alignment → U/M classification → Cecc detection → metrics*
- `tests/simulation/validation_metrics.py`: compares ground truth with `uecc.csv/mecc.csv/cecc.csv` (primarily by read_id)
- `tests/simulation/run_recall_benchmark.py`: repeats the pipeline across multiple seeds and summarizes recall (mean/std)

## What is `unclassified.csv`?

In `tests/simulation/run_validation.py`:

1. `UMeccClassifier.run()` produces `uecc_df / mecc_df / unclassified_df`.
2. `unclassified_df` is written to `unclassified.csv` **and used as input to CeccBuild**.

So `unclassified.csv` is mainly the pool of reads that were **not called as U or M at the U/M stage**, including:

- most true Cecc reads (expected to be handled by the subsequent `cecc_build` step)
- borderline U/Mecc reads (e.g. U uniqueness veto or Mecc multi-locus coverage not strong enough)
- low-quality alignment rows (see `quality_category` / `Gap_Percentage`)

For “final” detection, check whether a read appears in `uecc.csv`, `mecc.csv`, or `cecc.csv`. Only reads absent from all three are truly “unclassified” at the end of this simplified pipeline.

## Why does recall fluctuate across runs?

With the same code + environment, results are reproducible for the same `seed`. Across different seeds, the simulator produces slightly different difficulty profiles (lengths, copy numbers, inter-chromosomal fraction, borderline cases), which changes how many reads end up near decision boundaries and get left in `Unclassified`.

In practice, report recall using **multiple runs (mean/std)** rather than chasing a single 100% run.

## How to run (recommended)

Run 3 repeats (default 1000U/1000M/1000C) and summarize:

```bash
python tests/simulation/run_recall_benchmark.py --runs 3
```

Outputs are written under `.tmp_work/simulation_benchmarks/.../benchmark_summary.json` (git-ignored by default).

Run a single end-to-end validation:

```bash
python tests/simulation/run_validation.py \
  --simulation-dir .tmp_work/sim_one/simulation \
  --results-dir .tmp_work/sim_one/results \
  --num-uecc 1000 --num-mecc 1000 --num-cecc 1000 \
  --seed 42
```

## Example baseline (illustrative)

On one 3-run benchmark (1000U/1000M/1000C; seeds 42/43/44), read-level recall was approximately:

- `Uecc`: ~0.976
- `Mecc`: ~0.987
- `Cecc`: ~0.999
- overall: ~0.987

These numbers are provided as a reference point; absolute values may change as parameters and implementations evolve.

---

## Pytest Regression Baseline (Reproducible)

To make simulation validation part of CI/local regression tests, the repo includes a **small, fixed-seed** baseline:

- Baseline file: `tests/simulation/baselines/synthetic_regression.json`
- Test: `tests/simulation/test_simulation_regression_baseline.py`
- Update script: `tests/simulation/update_synthetic_baseline.py`

Run:

```bash
pytest -q tests/simulation/test_simulation_regression_baseline.py
```

> Note: this test uses a synthetic alignment TSV (no external minimap2 dependency) to regression-test the U/M/C classification logic and assert that overall recall/precision/F1 does not regress significantly.

## One-command Baseline Update

When you intentionally change logic/parameters and want to update the expected metrics:

```bash
python tests/simulation/update_synthetic_baseline.py
```

This rewrites the `expected.overall` section in `tests/simulation/baselines/synthetic_regression.json`.

## Negative Control: False Positives (U-focused)

The negative control constrains low-MAPQ false positives typically seen in repetitive/low-complexity contexts:

- Construct full-coverage alignments with `MAPQ=0` and assert that **no high-confidence** U calls are produced (e.g. `confidence_score >= 0.8`).
- Test: `test_negative_control_low_mapq_has_no_high_confidence_calls` in `tests/simulation/test_simulation_regression_baseline.py`.
