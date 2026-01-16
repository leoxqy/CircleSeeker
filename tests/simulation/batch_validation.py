#!/usr/bin/env python3
"""
Batch Validation Script for CircleSeeker

Tests multiple scales with multiple runs each, recording FP/FN rates.

Notes:
    This script reports an "accuracy" metric defined as:
        accuracy = TP / (TP + FP + FN)
    This avoids requiring true negatives (TN), which are not well-defined for
    eccDNA discovery across a genome-scale search space.

Usage:
    python batch_validation.py
"""

import os
import sys
import json
import time
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Any, Optional

# Add project root to path
PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
sys.path.insert(0, str(PROJECT_ROOT / "tests/simulation"))

from eccdna_simulator import SimulationConfig, run_simulation
from validation_metrics import run_validation


def _safe_div(numer: float, denom: float) -> float:
    return numer / denom if denom else 0.0


def _f1(precision: float, recall: float) -> float:
    return (2 * precision * recall / (precision + recall)) if (precision + recall) else 0.0


def _detection_accuracy(tp: int, fp: int, fn: int) -> float:
    """Accuracy without TN: TP / (TP + FP + FN)."""
    return _safe_div(tp, tp + fp + fn)


def _metrics_from_validation_report(
    report: Dict[str, Any],
    expected_gt: Dict[str, int],
) -> Optional[Dict[str, Any]]:
    """
    Convert `validation_report.json` (written by `tests/simulation/validation_metrics.py`)
    into this script's metrics structure.

    Returns None if the report indicates an empty/failed run (e.g. no predictions).
    """
    overall = report.get("overall", {}) or {}
    total_pred = int(overall.get("prediction_count", 0) or 0)
    if total_pred <= 0:
        return None

    by_type = report.get("by_type", {}) or {}
    type_map = {"uecc": "Uecc", "mecc": "Mecc", "cecc": "Cecc"}

    metrics: Dict[str, Any] = {}
    for out_key, in_key in type_map.items():
        t = by_type.get(in_key, {}) or {}
        tp = int(t.get("true_positives", 0) or 0)
        fp = int(t.get("false_positives", 0) or 0)
        fn = int(t.get("false_negatives", 0) or 0)
        gt = int(t.get("ground_truth_count", expected_gt.get(out_key, 0)) or 0)
        pred = int(t.get("prediction_count", 0) or 0)

        recall = float(t.get("recall", _safe_div(tp, gt)) or 0.0)
        precision = float(t.get("precision", _safe_div(tp, pred)) or 0.0)
        f1 = float(t.get("f1_score", _f1(precision, recall)) or 0.0)

        metrics[out_key] = {
            "ground_truth": gt,
            "predictions": pred,
            "tp": tp,
            "fp": fp,
            "fn": fn,
            "recall": recall,
            "precision": precision,
            "f1": f1,
            "accuracy": _detection_accuracy(tp, fp, fn),
        }

    total_tp = int(overall.get("true_positives", 0) or 0)
    total_fp = int(overall.get("false_positives", 0) or 0)
    total_fn = int(overall.get("false_negatives", 0) or 0)
    total_gt = int(overall.get("ground_truth_count", 0) or 0)

    overall_recall = float(overall.get("recall", _safe_div(total_tp, total_gt)) or 0.0)
    overall_precision = float(overall.get("precision", _safe_div(total_tp, total_pred)) or 0.0)
    overall_f1 = float(overall.get("f1_score", _f1(overall_precision, overall_recall)) or 0.0)
    overall_accuracy = _detection_accuracy(total_tp, total_fp, total_fn)

    metrics["overall"] = {
        "total_gt": total_gt,
        "total_pred": total_pred,
        "total_tp": total_tp,
        "total_fp": total_fp,
        "total_fn": total_fn,
        "recall": overall_recall,
        "precision": overall_precision,
        "f1": overall_f1,
        "accuracy": overall_accuracy,
    }

    return metrics


def run_single_test(
    scale: int,
    run_id: int,
    seed: int,
    base_output_dir: Path,
    threads: int = 4,
    minimap2_preset: str = "sr",
    minimap2_timeout: int = 1800,
    resume: bool = False,
) -> Dict[str, Any]:
    """Run a single validation test and return results."""

    # Calculate eccDNA counts (ratio: 70% U, 15% M, 15% C)
    num_uecc = int(scale * 0.7)
    num_mecc = int(scale * 0.15)
    num_cecc = scale - num_uecc - num_mecc

    # Setup directories
    test_dir = base_output_dir / f"scale_{scale}" / f"run_{run_id:02d}_seed{seed}"
    sim_dir = test_dir / "simulation"
    results_dir = test_dir / "results"

    sim_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

    report_path = results_dir / "validation_report.json"

    result = {
        "scale": scale,
        "run_id": run_id,
        "seed": seed,
        "num_uecc": num_uecc,
        "num_mecc": num_mecc,
        "num_cecc": num_cecc,
        "total": scale,
        "success": False,
        "error": None,
        "metrics": {},
    }

    try:
        expected_gt = {"uecc": num_uecc, "mecc": num_mecc, "cecc": num_cecc}
        if resume and report_path.exists():
            try:
                existing = json.loads(report_path.read_text())
            except Exception:
                existing = None
            if existing is not None:
                loaded = _metrics_from_validation_report(existing, expected_gt)
                if loaded is not None:
                    result["metrics"] = loaded
                    result["success"] = True
                    return result

        # Step 1: Run simulation
        config = SimulationConfig(
            output_dir=str(sim_dir),
            num_uecc=num_uecc,
            num_mecc=num_mecc,
            num_cecc=num_cecc,
            seed=seed,
        )
        run_simulation(config)

        # Step 2: Run pipeline
        from run_validation import run_circleseeker_pipeline

        reads_fasta = sim_dir / "simulated_reads.fa"
        reference_fasta = sim_dir / "reference.fa"

        pipeline_success = run_circleseeker_pipeline(
            reads_fasta=reads_fasta,
            reference_fasta=reference_fasta,
            output_dir=results_dir,
            threads=threads,
            minimap2_preset=minimap2_preset,
            minimap2_timeout=minimap2_timeout,
        )

        if not pipeline_success:
            result["error"] = "Pipeline failed"
            return result

        # Step 3: Run validation (suppress output)
        # Temporarily redirect stdout
        import io
        from contextlib import redirect_stdout

        f = io.StringIO()
        with redirect_stdout(f):
            by_type = run_validation(sim_dir, results_dir, report_path)

        type_map = {"uecc": "Uecc", "mecc": "Mecc", "cecc": "Cecc"}

        metrics: Dict[str, Any] = {}
        for out_key, in_key in type_map.items():
            m = by_type.get(in_key)
            if m is None:
                tp = fp = fn = pred = 0
                gt = expected_gt[out_key]
                recall = precision = f1 = 0.0
            else:
                tp = int(m.true_positives)
                fp = int(m.false_positives)
                fn = int(m.false_negatives)
                gt = int(m.total_ground_truth)
                pred = int(m.total_predictions)
                recall = float(m.recall)
                precision = float(m.precision)
                f1 = float(m.f1_score)

            metrics[out_key] = {
                "ground_truth": gt,
                "predictions": pred,
                "tp": tp,
                "fp": fp,
                "fn": fn,
                "recall": recall,
                "precision": precision,
                "f1": f1,
                "accuracy": _detection_accuracy(tp, fp, fn),
            }

        total_tp = sum(metrics[k]["tp"] for k in ("uecc", "mecc", "cecc"))
        total_fp = sum(metrics[k]["fp"] for k in ("uecc", "mecc", "cecc"))
        total_fn = sum(metrics[k]["fn"] for k in ("uecc", "mecc", "cecc"))
        total_gt = sum(metrics[k]["ground_truth"] for k in ("uecc", "mecc", "cecc"))
        total_pred = sum(metrics[k]["predictions"] for k in ("uecc", "mecc", "cecc"))

        overall_recall = _safe_div(total_tp, total_gt)
        overall_precision = _safe_div(total_tp, total_pred)
        overall_f1 = _f1(overall_precision, overall_recall)
        overall_accuracy = _detection_accuracy(total_tp, total_fp, total_fn)

        metrics["overall"] = {
            "total_gt": total_gt,
            "total_pred": total_pred,
            "total_tp": total_tp,
            "total_fp": total_fp,
            "total_fn": total_fn,
            "recall": overall_recall,
            "precision": overall_precision,
            "f1": overall_f1,
            "accuracy": overall_accuracy,
        }

        result["metrics"] = metrics
        result["success"] = True

    except Exception as e:
        result["error"] = str(e)
        import traceback
        traceback.print_exc()

    return result


def print_progress_table(results: List[Dict], scales: List[int], runs_per_scale: int):
    """Print a progress table of results."""
    print("\n" + "=" * 100)
    print("BATCH VALIDATION PROGRESS")
    print("=" * 100)

    # Header
    print(f"{'Scale':>8} | {'Run':>4} | {'Seed':>6} | {'Total':>6} | "
          f"{'TP':>6} | {'FP':>4} | {'FN':>4} | {'Recall':>8} | {'Precision':>10} | {'F1':>8} | {'Acc':>8} | Status")
    print("-" * 100)

    for r in results:
        if r["success"]:
            m = r["metrics"]["overall"]
            status = "OK"
            print(f"{r['scale']:>8} | {r['run_id']:>4} | {r['seed']:>6} | {r['total']:>6} | "
                  f"{m['total_tp']:>6} | {m['total_fp']:>4} | {m['total_fn']:>4} | "
                  f"{m['recall']*100:>7.2f}% | {m['precision']*100:>9.2f}% | {m['f1']*100:>7.2f}% | {m['accuracy']*100:>7.2f}% | {status}")
        else:
            status = f"FAIL: {r.get('error', 'Unknown')[:20]}"
            print(f"{r['scale']:>8} | {r['run_id']:>4} | {r['seed']:>6} | {r['total']:>6} | "
                  f"{'--':>6} | {'--':>4} | {'--':>4} | "
                  f"{'--':>8} | {'--':>10} | {'--':>8} | {'--':>8} | {status}")


def print_summary_table(results: List[Dict], scales: List[int]):
    """Print summary statistics by scale."""
    print("\n" + "=" * 120)
    print("SUMMARY BY SCALE")
    print("=" * 120)

    print(f"{'Scale':>8} | {'Runs':>5} | {'Avg TP':>8} | {'Avg FP':>8} | {'Avg FN':>8} | "
          f"{'Min Recall':>11} | {'Max Recall':>11} | {'Avg Recall':>11} | "
          f"{'Min Prec':>9} | {'Max Prec':>9} | {'Avg F1':>8}")
    print("-" * 120)

    for scale in scales:
        scale_results = [r for r in results if r["scale"] == scale and r["success"]]

        if not scale_results:
            print(f"{scale:>8} | {'N/A':>5} | No successful runs")
            continue

        tps = [r["metrics"]["overall"]["total_tp"] for r in scale_results]
        fps = [r["metrics"]["overall"]["total_fp"] for r in scale_results]
        fns = [r["metrics"]["overall"]["total_fn"] for r in scale_results]
        recalls = [r["metrics"]["overall"]["recall"] for r in scale_results]
        precisions = [r["metrics"]["overall"]["precision"] for r in scale_results]
        f1s = [r["metrics"]["overall"]["f1"] for r in scale_results]

        print(f"{scale:>8} | {len(scale_results):>5} | "
              f"{sum(tps)/len(tps):>8.1f} | {sum(fps)/len(fps):>8.1f} | {sum(fns)/len(fns):>8.1f} | "
              f"{min(recalls)*100:>10.2f}% | {max(recalls)*100:>10.2f}% | {sum(recalls)/len(recalls)*100:>10.2f}% | "
              f"{min(precisions)*100:>8.2f}% | {max(precisions)*100:>8.2f}% | {sum(f1s)/len(f1s)*100:>7.2f}%")


def _aggregate_runs(successful_runs: List[Dict[str, Any]]) -> Dict[str, Any]:
    total_tp = sum(r["metrics"]["overall"]["total_tp"] for r in successful_runs)
    total_fp = sum(r["metrics"]["overall"]["total_fp"] for r in successful_runs)
    total_fn = sum(r["metrics"]["overall"]["total_fn"] for r in successful_runs)
    total_gt = sum(r["metrics"]["overall"]["total_gt"] for r in successful_runs)
    total_pred = sum(r["metrics"]["overall"]["total_pred"] for r in successful_runs)

    recall = _safe_div(total_tp, total_gt)
    precision = _safe_div(total_tp, total_pred)
    f1 = _f1(precision, recall)
    accuracy = _detection_accuracy(total_tp, total_fp, total_fn)

    n = len(successful_runs)
    return {
        "runs": n,
        "tp": total_tp,
        "fp": total_fp,
        "fn": total_fn,
        "gt": total_gt,
        "pred": total_pred,
        "recall": recall,
        "precision": precision,
        "f1": f1,
        "accuracy": accuracy,
        "avg_fp": _safe_div(total_fp, n),
        "avg_fn": _safe_div(total_fn, n),
    }


def print_aggregate_by_scale_table(results: List[Dict], scales: List[int], runs_per_scale: int):
    """Print aggregate (sum) statistics by scale."""
    print("\n" + "=" * 120)
    print("AGGREGATE BY SCALE (SUM OVER RUNS)")
    print("=" * 120)

    print(
        f"{'Scale':>8} | {'OK':>4} | {'TP':>10} | {'FP':>10} | {'FN':>10} | "
        f"{'Acc':>8} | {'Avg FP':>10} | {'Avg FN':>10}"
    )
    print("-" * 120)

    for scale in scales:
        scale_results = [r for r in results if r["scale"] == scale and r["success"]]
        if not scale_results:
            print(f"{scale:>8} | {0:>4}/{runs_per_scale:<1} | {'--':>10} | {'--':>10} | {'--':>10} | "
                  f"{'--':>8} | {'--':>10} | {'--':>10}")
            continue

        agg = _aggregate_runs(scale_results)
        print(
            f"{scale:>8} | {len(scale_results):>4}/{runs_per_scale:<1} | "
            f"{agg['tp']:>10} | {agg['fp']:>10} | {agg['fn']:>10} | {agg['accuracy']*100:>7.2f}% | "
            f"{agg['avg_fp']:>10.2f} | {agg['avg_fn']:>10.2f}"
        )


def print_overall_aggregate(results: List[Dict], total_tests: int):
    """Print a single overall aggregate line across all scales/runs."""
    successful = [r for r in results if r.get("success")]

    print("\n" + "=" * 120)
    print("OVERALL AGGREGATE (ALL SCALES/RUNS)")
    print("=" * 120)

    if not successful:
        print(f"Successful runs: 0/{total_tests}")
        return

    agg = _aggregate_runs(successful)
    print(f"Successful runs: {agg['runs']}/{total_tests}")
    print(f"Total FP: {agg['fp']}, Total FN: {agg['fn']}, Accuracy: {agg['accuracy']*100:.2f}%")


def print_detailed_fp_fn_table(results: List[Dict], scales: List[int]):
    """Print detailed FP/FN breakdown by type."""
    print("\n" + "=" * 140)
    print("DETAILED FP/FN BY TYPE")
    print("=" * 140)

    print(f"{'Scale':>8} | {'Run':>4} | "
          f"{'U_GT':>5} | {'U_TP':>5} | {'U_FP':>5} | {'U_FN':>5} | "
          f"{'M_GT':>5} | {'M_TP':>5} | {'M_FP':>5} | {'M_FN':>5} | "
          f"{'C_GT':>5} | {'C_TP':>5} | {'C_FP':>5} | {'C_FN':>5}")
    print("-" * 140)

    for r in results:
        if not r["success"]:
            continue

        u = r["metrics"]["uecc"]
        m = r["metrics"]["mecc"]
        c = r["metrics"]["cecc"]

        print(f"{r['scale']:>8} | {r['run_id']:>4} | "
              f"{u['ground_truth']:>5} | {u['tp']:>5} | {u['fp']:>5} | {u['fn']:>5} | "
              f"{m['ground_truth']:>5} | {m['tp']:>5} | {m['fp']:>5} | {m['fn']:>5} | "
              f"{c['ground_truth']:>5} | {c['tp']:>5} | {c['fp']:>5} | {c['fn']:>5}")


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Batch validation for CircleSeeker")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads")
    parser.add_argument(
        "--minimap2-preset",
        type=str,
        default="sr",
        help="minimap2 preset to use for simulation validation (default: sr)",
    )
    parser.add_argument(
        "--minimap2-timeout",
        type=int,
        default=1800,
        help="Timeout (seconds) for minimap2 alignment step in each run (default: 1800)",
    )
    parser.add_argument("--output-dir", type=Path,
                        default=Path("tests/simulation/batch_validation_results"),
                        help="Output directory for results")
    parser.add_argument("--scales", type=str,
                        default="100,500,1000,3000,5000,10000,50000",
                        help="Comma-separated list of scales to test")
    parser.add_argument("--runs", type=int, default=5, help="Number of runs per scale")
    parser.add_argument("--base-seed", type=int, default=42, help="Base random seed")
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Skip runs with an existing non-empty validation_report.json",
    )

    args = parser.parse_args()

    # Parse scales
    scales = [int(s.strip()) for s in args.scales.split(",")]
    runs_per_scale = args.runs
    base_seed = args.base_seed

    # Change to project root
    os.chdir(PROJECT_ROOT)

    print("=" * 80)
    print("CircleSeeker Batch Validation")
    print(f"Started: {datetime.now().isoformat()}")
    print("=" * 80)
    print(f"Scales to test: {scales}")
    print(f"Runs per scale: {runs_per_scale}")
    print(f"Total tests: {len(scales) * runs_per_scale}")
    print(f"minimap2 preset: {args.minimap2_preset}")
    print(f"minimap2 timeout: {args.minimap2_timeout}s")
    print(f"Resume: {args.resume}")
    print(f"Output directory: {args.output_dir}")
    print("=" * 80)

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Run all tests
    all_results = []
    total_tests = len(scales) * runs_per_scale
    current_test = 0

    start_time = time.time()

    for scale in scales:
        print(f"\n{'='*60}")
        print(f"Testing scale: {scale}")
        print(f"{'='*60}")

        for run_id in range(1, runs_per_scale + 1):
            current_test += 1
            seed = base_seed + run_id - 1

            print(f"\n[{current_test}/{total_tests}] Scale={scale}, Run={run_id}, Seed={seed}")

            test_start = time.time()
            result = run_single_test(
                scale=scale,
                run_id=run_id,
                seed=seed,
                base_output_dir=args.output_dir,
                threads=args.threads,
                minimap2_preset=args.minimap2_preset,
                minimap2_timeout=args.minimap2_timeout,
                resume=args.resume,
            )
            test_time = time.time() - test_start
            result["test_time_seconds"] = test_time

            all_results.append(result)

            # Print quick summary
            if result["success"]:
                m = result["metrics"]["overall"]
                print(f"  Result: TP={m['total_tp']}, FP={m['total_fp']}, FN={m['total_fn']}, "
                      f"Acc={m['accuracy']*100:.2f}%, F1={m['f1']*100:.2f}% ({test_time:.1f}s)")
            else:
                print(f"  Result: FAILED - {result.get('error', 'Unknown error')}")

    total_time = time.time() - start_time

    # Print all tables
    print_progress_table(all_results, scales, runs_per_scale)
    print_summary_table(all_results, scales)
    print_aggregate_by_scale_table(all_results, scales, runs_per_scale)
    print_detailed_fp_fn_table(all_results, scales)
    print_overall_aggregate(all_results, total_tests)

    # Save full results to JSON
    results_file = args.output_dir / "batch_validation_results.json"
    successful = [r for r in all_results if r.get("success")]
    aggregate_overall = _aggregate_runs(successful) if successful else {}
    aggregate_by_scale = {
        str(scale): _aggregate_runs([r for r in successful if r["scale"] == scale])
        for scale in scales
    }
    with open(results_file, "w") as f:
        json.dump({
            "timestamp": datetime.now().isoformat(),
            "config": {
                "scales": scales,
                "runs_per_scale": runs_per_scale,
                "base_seed": base_seed,
                "threads": args.threads,
                "minimap2_preset": args.minimap2_preset,
                "minimap2_timeout": args.minimap2_timeout,
                "resume": args.resume,
            },
            "aggregate": {
                "overall": aggregate_overall,
                "by_scale": aggregate_by_scale,
            },
            "total_time_seconds": total_time,
            "results": all_results,
        }, f, indent=2)

    print(f"\n{'='*80}")
    print(f"Batch Validation Complete!")
    print(f"Total time: {total_time/60:.1f} minutes")
    print(f"Results saved to: {results_file}")
    print(f"{'='*80}")


if __name__ == "__main__":
    main()
