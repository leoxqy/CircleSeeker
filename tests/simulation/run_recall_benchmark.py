#!/usr/bin/env python3
"""
Repeated U/M/C simulation benchmark for CircleSeeker.

This script runs the end-to-end simulation + pipeline + validation multiple times
with different seeds, then summarizes recall (by type + overall).
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path
from statistics import mean, stdev
from typing import Any, Dict, List, Optional


PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent


def _resolve_under_project(path: Path) -> Path:
    if path.is_absolute():
        return path
    return (PROJECT_ROOT / path).resolve()


def _safe_rmtree(path: Path) -> None:
    if path.exists():
        shutil.rmtree(path)


def _load_json(path: Path) -> Dict[str, Any]:
    with open(path) as f:
        return json.load(f)


def _agg(values: List[float]) -> Dict[str, float]:
    if not values:
        return {"mean": 0.0, "stdev": 0.0, "min": 0.0, "max": 0.0}
    if len(values) == 1:
        v = values[0]
        return {"mean": v, "stdev": 0.0, "min": v, "max": v}
    return {
        "mean": mean(values),
        "stdev": stdev(values),
        "min": min(values),
        "max": max(values),
    }


def _fmt(value: Optional[float]) -> str:
    if value is None:
        return "NA"
    return f"{value:.4f}"


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run repeated 1000U/1000M/1000C simulations and summarize recall"
    )
    parser.add_argument("--runs", type=int, default=3, help="Number of runs (default: 3)")
    parser.add_argument(
        "--seeds",
        type=int,
        nargs="*",
        default=None,
        help="Explicit seeds list (overrides --runs/--base-seed)",
    )
    parser.add_argument("--base-seed", type=int, default=42, help="Base seed (default: 42)")
    parser.add_argument("--num-uecc", type=int, default=1000, help="Number of UeccDNA")
    parser.add_argument("--num-mecc", type=int, default=1000, help="Number of MeccDNA")
    parser.add_argument("--num-cecc", type=int, default=1000, help="Number of CeccDNA")
    parser.add_argument("--threads", type=int, default=4, help="minimap2 threads (default: 4)")
    parser.add_argument(
        "--minimap2-preset",
        type=str,
        default="sr",
        help="minimap2 preset (default: sr)",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path(".tmp_work/simulation_benchmarks"),
        help="Benchmark output root (default: .tmp_work/simulation_benchmarks)",
    )
    parser.add_argument(
        "--no-clean",
        action="store_true",
        help="Do not delete the benchmark directory before running",
    )

    args = parser.parse_args()

    if args.seeds is not None and len(args.seeds) == 0:
        print("Error: --seeds provided but empty")
        return 2

    seeds = args.seeds if args.seeds is not None else [
        args.base_seed + i for i in range(max(0, args.runs))
    ]
    if not seeds:
        print("Error: no seeds to run (use --runs >= 1 or --seeds ...)")
        return 2

    output_root = _resolve_under_project(args.output_root)
    benchmark_dir = output_root / f"umc_{args.num_uecc}_{args.num_mecc}_{args.num_cecc}"
    if not args.no_clean:
        _safe_rmtree(benchmark_dir)
    benchmark_dir.mkdir(parents=True, exist_ok=True)

    validation_script = PROJECT_ROOT / "tests" / "simulation" / "run_validation.py"
    run_rows: List[Dict[str, Any]] = []

    print(f"Benchmark dir: {benchmark_dir}")
    for idx, seed in enumerate(seeds, start=1):
        run_dir = benchmark_dir / f"run_{idx:02d}_seed{seed}"
        sim_dir = run_dir / "simulation"
        results_dir = run_dir / "results"

        cmd = [
            sys.executable,
            str(validation_script),
            "--simulation-dir",
            str(sim_dir),
            "--results-dir",
            str(results_dir),
            "--threads",
            str(args.threads),
            "--minimap2-preset",
            args.minimap2_preset,
            "--num-uecc",
            str(args.num_uecc),
            "--num-mecc",
            str(args.num_mecc),
            "--num-cecc",
            str(args.num_cecc),
            "--seed",
            str(seed),
        ]

        print(f"\n=== Run {idx}/{len(seeds)} (seed={seed}) ===")
        subprocess.run(cmd, check=True)

        report_path = results_dir / "validation_report.json"
        report = _load_json(report_path)

        by_type = report.get("by_type", {})
        overall = report.get("overall", {})
        row = {
            "run": idx,
            "seed": seed,
            "paths": {
                "simulation_dir": str(sim_dir),
                "results_dir": str(results_dir),
                "validation_report": str(report_path),
            },
            "recall": {
                "Uecc": by_type.get("Uecc", {}).get("recall"),
                "Mecc": by_type.get("Mecc", {}).get("recall"),
                "Cecc": by_type.get("Cecc", {}).get("recall"),
                "overall": overall.get("recall"),
            },
        }
        run_rows.append(row)

    # Aggregate
    recalls = {
        "Uecc": [r["recall"]["Uecc"] for r in run_rows if r["recall"]["Uecc"] is not None],
        "Mecc": [r["recall"]["Mecc"] for r in run_rows if r["recall"]["Mecc"] is not None],
        "Cecc": [r["recall"]["Cecc"] for r in run_rows if r["recall"]["Cecc"] is not None],
        "overall": [r["recall"]["overall"] for r in run_rows if r["recall"]["overall"] is not None],
    }
    summary = {
        "config": {
            "num_uecc": args.num_uecc,
            "num_mecc": args.num_mecc,
            "num_cecc": args.num_cecc,
            "threads": args.threads,
            "minimap2_preset": args.minimap2_preset,
            "seeds": seeds,
        },
        "runs": run_rows,
        "recall_summary": {k: _agg(v) for k, v in recalls.items()},
    }

    summary_path = benchmark_dir / "benchmark_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    # Print quick table
    print("\n" + "=" * 72)
    print("Recall Summary (per run)")
    print("=" * 72)
    print("run  seed   Uecc    Mecc    Cecc    overall")
    for row in run_rows:
        r = row["recall"]
        print(
            f"{row['run']:>3}  {row['seed']:>4}  "
            f"{_fmt(r.get('Uecc')):>6}  "
            f"{_fmt(r.get('Mecc')):>6}  "
            f"{_fmt(r.get('Cecc')):>6}  "
            f"{_fmt(r.get('overall')):>7}"
        )

    print("\n" + "=" * 72)
    print("Recall Summary (aggregated)")
    print("=" * 72)
    for key in ["Uecc", "Mecc", "Cecc", "overall"]:
        s = summary["recall_summary"][key]
        print(
            f"{key:>7}: mean={s['mean']:.4f}  stdev={s['stdev']:.4f}  "
            f"min={s['min']:.4f}  max={s['max']:.4f}"
        )

    print(f"\nSaved: {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

