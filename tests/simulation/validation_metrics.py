#!/usr/bin/env python3
"""
Validation Metrics Calculator for CircleSeeker

Compares CircleSeeker output with ground truth to calculate:
- Recall (Sensitivity): TP / (TP + FN)
- Precision: TP / (TP + FP)
- F1 Score: 2 * (Precision * Recall) / (Precision + Recall)
- False Positive Rate
- False Negative Rate
"""

import csv
import json
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Dict, Set, Tuple, Optional
from collections import defaultdict


@dataclass
class ValidationMetrics:
    """Validation metrics for a single eccDNA type."""
    eccdna_type: str
    true_positives: int = 0
    false_positives: int = 0
    false_negatives: int = 0
    total_ground_truth: int = 0
    total_predictions: int = 0

    # Matched records for detailed analysis
    matched_gt_ids: Set[str] = field(default_factory=set)
    matched_pred_ids: Set[str] = field(default_factory=set)
    unmatched_gt_ids: Set[str] = field(default_factory=set)
    unmatched_pred_ids: Set[str] = field(default_factory=set)

    @property
    def recall(self) -> float:
        """Recall = TP / (TP + FN)"""
        if self.total_ground_truth == 0:
            return 0.0
        return self.true_positives / self.total_ground_truth

    @property
    def precision(self) -> float:
        """Precision = TP / (TP + FP)"""
        if self.total_predictions == 0:
            return 0.0
        return self.true_positives / self.total_predictions

    @property
    def f1_score(self) -> float:
        """F1 = 2 * (Precision * Recall) / (Precision + Recall)"""
        if self.precision + self.recall == 0:
            return 0.0
        return 2 * (self.precision * self.recall) / (self.precision + self.recall)

    @property
    def false_positive_rate(self) -> float:
        """FPR = FP / Total Predictions"""
        if self.total_predictions == 0:
            return 0.0
        return self.false_positives / self.total_predictions

    @property
    def false_negative_rate(self) -> float:
        """FNR = FN / Total Ground Truth"""
        if self.total_ground_truth == 0:
            return 0.0
        return self.false_negatives / self.total_ground_truth

    def to_dict(self) -> dict:
        return {
            "type": self.eccdna_type,
            "ground_truth_count": self.total_ground_truth,
            "prediction_count": self.total_predictions,
            "true_positives": self.true_positives,
            "false_positives": self.false_positives,
            "false_negatives": self.false_negatives,
            "recall": round(self.recall, 4),
            "precision": round(self.precision, 4),
            "f1_score": round(self.f1_score, 4),
            "false_positive_rate": round(self.false_positive_rate, 4),
            "false_negative_rate": round(self.false_negative_rate, 4),
        }


@dataclass
class GroundTruthRecord:
    """Parsed ground truth record."""
    eccdna_id: str
    eccdna_type: str
    read_id: str
    chrom: str
    start: int
    end: int
    strand: str
    length: int
    regions: str  # Full region string for multi-segment


@dataclass
class PredictionRecord:
    """Parsed prediction record from CircleSeeker output."""
    eccdna_id: str
    eccdna_type: str
    read_id: str
    chrom: str
    start: int
    end: int
    strand: str
    length: int
    regions: str = ""


class ValidationCalculator:
    """Calculate validation metrics by comparing predictions to ground truth."""

    def __init__(self,
                 position_tolerance: int = 50,
                 overlap_threshold: float = 0.8):
        """
        Args:
            position_tolerance: Max distance (bp) for position matching
            overlap_threshold: Min overlap fraction for region matching
        """
        self.position_tolerance = position_tolerance
        self.overlap_threshold = overlap_threshold

    def load_ground_truth(self, gt_path: Path) -> List[GroundTruthRecord]:
        """Load ground truth from CSV file."""
        records = []
        with open(gt_path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                record = GroundTruthRecord(
                    eccdna_id=row.get("eccdna_id", ""),
                    eccdna_type=row.get("eccdna_type", ""),
                    read_id=row.get("read_id", ""),
                    chrom=row.get("chr", ""),
                    start=int(row.get("start0", 0)),
                    end=int(row.get("end0", 0)),
                    strand=row.get("strand", "+"),
                    length=int(row.get("length", 0)),
                    regions=row.get("regions", ""),
                )
                records.append(record)
        return records

    def load_predictions(self, pred_path: Path,
                         eccdna_type: str) -> List[PredictionRecord]:
        """Load predictions from CircleSeeker output CSV."""
        records = []

        if not pred_path.exists():
            return records

        with open(pred_path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                # Handle different column naming conventions
                read_id = row.get("reads", row.get("read_id", row.get("query_id", "")))

                # Extract first part of read_id before '|'
                if "|" in read_id:
                    read_id = read_id.split("|")[0]

                record = PredictionRecord(
                    eccdna_id=row.get("eccDNA_id", row.get("eccdna_id", "")),
                    eccdna_type=eccdna_type,
                    read_id=read_id,
                    chrom=row.get("chr", row.get("subject_id", "")),
                    start=int(row.get("start0", row.get("s_start", 0))),
                    end=int(row.get("end0", row.get("s_end", 0))),
                    strand=row.get("strand", "+"),
                    length=int(row.get("length", row.get("circle_length", 0))),
                    regions=row.get("regions", ""),
                )
                records.append(record)

        return records

    def _regions_overlap(self, gt: GroundTruthRecord,
                         pred: PredictionRecord) -> bool:
        """Check if two genomic regions overlap sufficiently."""
        # Different chromosome = no match
        if gt.chrom != pred.chrom:
            return False

        # Calculate overlap
        overlap_start = max(gt.start, pred.start)
        overlap_end = min(gt.end, pred.end)
        overlap_len = max(0, overlap_end - overlap_start)

        # Check overlap fraction
        gt_len = gt.end - gt.start
        pred_len = pred.end - pred.start

        if gt_len == 0 or pred_len == 0:
            return False

        gt_overlap_frac = overlap_len / gt_len
        pred_overlap_frac = overlap_len / pred_len

        # Both must have sufficient overlap
        return (gt_overlap_frac >= self.overlap_threshold or
                pred_overlap_frac >= self.overlap_threshold)

    def _positions_match(self, gt: GroundTruthRecord,
                         pred: PredictionRecord) -> bool:
        """Check if positions match within tolerance."""
        if gt.chrom != pred.chrom:
            return False

        start_diff = abs(gt.start - pred.start)
        end_diff = abs(gt.end - pred.end)

        return (start_diff <= self.position_tolerance and
                end_diff <= self.position_tolerance)

    def _read_ids_match(self, gt: GroundTruthRecord,
                        pred: PredictionRecord) -> bool:
        """Check if read IDs match (for simulation validation)."""
        # Extract base read ID from ground truth (format: sim_uecc_0001|...)
        gt_base = gt.read_id.split("|")[0] if "|" in gt.read_id else gt.read_id
        return gt_base == pred.read_id

    def match_records(self, ground_truth: List[GroundTruthRecord],
                      predictions: List[PredictionRecord],
                      match_by_read_id: bool = True) -> ValidationMetrics:
        """Match predictions to ground truth and calculate metrics."""

        if match_by_read_id and predictions:
            # Some CircleSeeker outputs (e.g. Cecc segment tables) can emit multiple
            # rows per read; for simulation validation we evaluate at read-level.
            deduped: Dict[str, PredictionRecord] = {}
            for pred in predictions:
                deduped.setdefault(pred.read_id, pred)
            predictions = list(deduped.values())

        if not ground_truth:
            eccdna_type = predictions[0].eccdna_type if predictions else "Unknown"
        else:
            eccdna_type = ground_truth[0].eccdna_type

        metrics = ValidationMetrics(eccdna_type=eccdna_type)
        metrics.total_ground_truth = len(ground_truth)
        metrics.total_predictions = len(predictions)

        # Track which records have been matched
        matched_gt = set()
        matched_pred = set()

        # Create indices for faster lookup
        gt_by_chrom = defaultdict(list)
        for gt in ground_truth:
            gt_by_chrom[gt.chrom].append(gt)

        # Try to match each prediction
        for pred in predictions:
            best_match = None
            best_overlap = 0

            # If matching by read ID (for simulation)
            if match_by_read_id:
                for gt in ground_truth:
                    if gt.eccdna_id in matched_gt:
                        continue
                    if self._read_ids_match(gt, pred):
                        best_match = gt
                        break

            # Fall back to position matching
            if best_match is None:
                candidates = gt_by_chrom.get(pred.chrom, [])
                for gt in candidates:
                    if gt.eccdna_id in matched_gt:
                        continue

                    if self._positions_match(gt, pred):
                        best_match = gt
                        break
                    elif self._regions_overlap(gt, pred):
                        # Calculate overlap for ranking
                        overlap_start = max(gt.start, pred.start)
                        overlap_end = min(gt.end, pred.end)
                        overlap = overlap_end - overlap_start
                        if overlap > best_overlap:
                            best_overlap = overlap
                            best_match = gt

            if best_match:
                matched_gt.add(best_match.eccdna_id)
                matched_pred.add(pred.eccdna_id)
                metrics.true_positives += 1
                metrics.matched_gt_ids.add(best_match.eccdna_id)
                metrics.matched_pred_ids.add(pred.eccdna_id)

        # Calculate FP and FN
        metrics.false_positives = len(predictions) - metrics.true_positives
        metrics.false_negatives = len(ground_truth) - metrics.true_positives

        # Track unmatched records
        metrics.unmatched_gt_ids = {gt.eccdna_id for gt in ground_truth
                                    if gt.eccdna_id not in matched_gt}
        metrics.unmatched_pred_ids = {pred.eccdna_id for pred in predictions
                                      if pred.eccdna_id not in matched_pred}

        return metrics


def run_validation(simulation_dir: Path,
                   results_dir: Path,
                   output_path: Optional[Path] = None) -> Dict[str, ValidationMetrics]:
    """
    Run validation comparing CircleSeeker results to ground truth.

    Args:
        simulation_dir: Directory containing ground truth files
        results_dir: Directory containing CircleSeeker output
        output_path: Optional path to write validation report

    Returns:
        Dictionary of metrics by eccDNA type
    """
    calculator = ValidationCalculator()

    results = {}

    # Define file mappings
    type_mappings = [
        ("Uecc", "ground_truth_uecc.csv", ["uecc.csv", "confirmed_uecc.csv"]),
        ("Mecc", "ground_truth_mecc.csv", ["mecc.csv", "confirmed_mecc.csv"]),
        ("Cecc", "ground_truth_cecc.csv", ["cecc.csv", "confirmed_cecc.csv"]),
    ]

    print("\n" + "=" * 60)
    print("CircleSeeker Validation Report")
    print("=" * 60)

    for eccdna_type, gt_file, pred_files in type_mappings:
        gt_path = simulation_dir / gt_file

        if not gt_path.exists():
            print(f"\n[{eccdna_type}] Ground truth not found: {gt_path}")
            continue

        # Load ground truth
        ground_truth = calculator.load_ground_truth(gt_path)

        # Find and load predictions
        predictions = []
        pred_file_used = None
        for pred_file in pred_files:
            pred_path = results_dir / pred_file
            if pred_path.exists():
                predictions = calculator.load_predictions(pred_path, eccdna_type)
                pred_file_used = pred_file
                break

        # Calculate metrics
        metrics = calculator.match_records(ground_truth, predictions)
        results[eccdna_type] = metrics

        # Print results
        print(f"\n[{eccdna_type}]")
        print(f"  Ground Truth: {metrics.total_ground_truth}")
        print(f"  Predictions:  {metrics.total_predictions} ({pred_file_used or 'not found'})")
        print(f"  True Positives:  {metrics.true_positives}")
        print(f"  False Positives: {metrics.false_positives}")
        print(f"  False Negatives: {metrics.false_negatives}")
        print(f"  ─────────────────────────")
        print(f"  Recall:    {metrics.recall:.2%}")
        print(f"  Precision: {metrics.precision:.2%}")
        print(f"  F1 Score:  {metrics.f1_score:.2%}")

        if metrics.unmatched_gt_ids:
            print(f"  Missed: {', '.join(sorted(metrics.unmatched_gt_ids)[:5])}")
            if len(metrics.unmatched_gt_ids) > 5:
                print(f"          ... and {len(metrics.unmatched_gt_ids) - 5} more")

    # Calculate overall metrics
    total_tp = sum(m.true_positives for m in results.values())
    total_fp = sum(m.false_positives for m in results.values())
    total_fn = sum(m.false_negatives for m in results.values())
    total_gt = sum(m.total_ground_truth for m in results.values())
    total_pred = sum(m.total_predictions for m in results.values())

    overall_recall = total_tp / total_gt if total_gt > 0 else 0
    overall_precision = total_tp / total_pred if total_pred > 0 else 0
    overall_f1 = (2 * overall_precision * overall_recall /
                  (overall_precision + overall_recall)
                  if (overall_precision + overall_recall) > 0 else 0)

    print("\n" + "=" * 60)
    print("Overall Metrics")
    print("=" * 60)
    print(f"  Total Ground Truth: {total_gt}")
    print(f"  Total Predictions:  {total_pred}")
    print(f"  Total TP: {total_tp}, FP: {total_fp}, FN: {total_fn}")
    print(f"  ─────────────────────────")
    print(f"  Overall Recall:    {overall_recall:.2%}")
    print(f"  Overall Precision: {overall_precision:.2%}")
    print(f"  Overall F1 Score:  {overall_f1:.2%}")

    # Write report if output path specified
    if output_path:
        report = {
            "by_type": {k: v.to_dict() for k, v in results.items()},
            "overall": {
                "ground_truth_count": total_gt,
                "prediction_count": total_pred,
                "true_positives": total_tp,
                "false_positives": total_fp,
                "false_negatives": total_fn,
                "recall": round(overall_recall, 4),
                "precision": round(overall_precision, 4),
                "f1_score": round(overall_f1, 4),
            }
        }

        with open(output_path, 'w') as f:
            json.dump(report, f, indent=2)
        print(f"\nReport saved to: {output_path}")

    return results


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Validate CircleSeeker results against ground truth"
    )
    parser.add_argument(
        "-s", "--simulation-dir",
        type=Path,
        required=True,
        help="Directory containing ground truth files"
    )
    parser.add_argument(
        "-r", "--results-dir",
        type=Path,
        required=True,
        help="Directory containing CircleSeeker output files"
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        default=None,
        help="Output path for validation report (JSON)"
    )

    args = parser.parse_args()

    run_validation(args.simulation_dir, args.results_dir, args.output)
