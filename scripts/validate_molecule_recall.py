#!/usr/bin/env python3
"""
Molecule-level recall/precision evaluation with multi-segment matching.

This script matches predicted eccDNA molecules to ground truth molecules
using segment-wise reciprocal overlap. It supports multi-segment CeccDNA
matching with a configurable coverage threshold (default: 0.90).
"""

from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import pandas as pd

REGION_RE = re.compile(r"([^:]+):(\d+)-(\d+)")


@dataclass(frozen=True)
class Segment:
    chrom: str
    start: int
    end: int

    @property
    def length(self) -> int:
        return max(1, int(self.end) - int(self.start))


def parse_regions(value: str) -> list[Segment]:
    if value is None:
        return []
    raw = str(value)
    if not raw or raw.lower() == "nan":
        return []
    raw = raw.replace("|", ";")
    segments: list[Segment] = []
    for part in raw.split(";"):
        part = part.strip()
        if not part:
            continue
        match = REGION_RE.match(part)
        if not match:
            continue
        chrom = match.group(1)
        start = int(match.group(2))
        end = int(match.group(3))
        if end < start:
            start, end = end, start
        segments.append(Segment(chrom, start, end))
    return segments


def dedupe_segments(segments: Iterable[Segment]) -> list[Segment]:
    seen = set()
    uniq: list[Segment] = []
    for seg in segments:
        key = (seg.chrom, seg.start, seg.end)
        if key in seen:
            continue
        seen.add(key)
        uniq.append(seg)
    return uniq


def index_by_chrom(segments: Iterable[Segment]) -> dict[str, list[Segment]]:
    by_chrom: dict[str, list[Segment]] = defaultdict(list)
    for seg in segments:
        by_chrom[seg.chrom].append(seg)
    return dict(by_chrom)


def merge_segments(segments: list[Segment]) -> dict[str, list[Segment]]:
    merged: dict[str, list[Segment]] = defaultdict(list)
    by_chrom = index_by_chrom(segments)
    for chrom, segs in by_chrom.items():
        segs = sorted(segs, key=lambda s: (s.start, s.end))
        if not segs:
            continue
        cur_start = segs[0].start
        cur_end = segs[0].end
        for seg in segs[1:]:
            if seg.start <= cur_end:
                cur_end = max(cur_end, seg.end)
            else:
                merged[chrom].append(Segment(chrom, cur_start, cur_end))
                cur_start = seg.start
                cur_end = seg.end
        merged[chrom].append(Segment(chrom, cur_start, cur_end))
    return dict(merged)


def union_length(segments_by_chrom: dict[str, list[Segment]]) -> int:
    total = 0
    for segs in segments_by_chrom.values():
        for seg in segs:
            total += max(0, seg.end - seg.start)
    return total


def overlap_length(a: list[Segment], b: list[Segment]) -> int:
    total = 0
    i = 0
    j = 0
    while i < len(a) and j < len(b):
        s = max(a[i].start, b[j].start)
        e = min(a[i].end, b[j].end)
        if e > s:
            total += e - s
        if a[i].end <= b[j].end:
            i += 1
        else:
            j += 1
    return total


def segment_coverage(
    gt_segments: list[Segment],
    pred_segments: list[Segment],
) -> tuple[float, float]:
    if not gt_segments or not pred_segments:
        return 0.0, 0.0

    gt_merged = merge_segments(gt_segments)
    pred_merged = merge_segments(pred_segments)

    gt_total = union_length(gt_merged)
    pred_total = union_length(pred_merged)
    if gt_total <= 0 or pred_total <= 0:
        return 0.0, 0.0

    overlap = 0
    for chrom, gt_list in gt_merged.items():
        pred_list = pred_merged.get(chrom, [])
        if not pred_list:
            continue
        overlap += overlap_length(gt_list, pred_list)

    return overlap / gt_total, overlap / pred_total


def any_site_coverage(
    gt_segments: list[Segment],
    pred_segments: list[Segment],
) -> float:
    if not gt_segments or not pred_segments:
        return 0.0
    gt_merged = merge_segments(gt_segments)
    pred_merged = merge_segments(pred_segments)
    best = 0.0
    for chrom, gt_list in gt_merged.items():
        pred_list = pred_merged.get(chrom, [])
        if not pred_list:
            continue
        for gt_seg in gt_list:
            overlap = overlap_length([gt_seg], pred_list)
            cov = overlap / gt_seg.length
            if cov > best:
                best = cov
    return best


def match_molecules(
    gt_map: dict[str, list[Segment]],
    pred_map: dict[str, list[Segment]],
    threshold: float,
    require_pred_coverage: bool,
    match_mode: str,
) -> dict[str, int]:
    matched_gt: set[str] = set()
    matched_pred: set[str] = set()

    for pred_id, pred_segments in pred_map.items():
        best_gt: Optional[str] = None
        best_score = 0.0
        for gt_id, gt_segments in gt_map.items():
            if gt_id in matched_gt:
                continue
            if match_mode == "any_site":
                cov = any_site_coverage(gt_segments, pred_segments)
                if cov < threshold:
                    continue
                score = cov
                if score > best_score:
                    best_score = score
                    best_gt = gt_id
            else:
                cov_gt, cov_pred = segment_coverage(gt_segments, pred_segments)
                if cov_gt < threshold:
                    continue
                if require_pred_coverage and cov_pred < threshold:
                    continue
                score = cov_gt + (cov_pred if require_pred_coverage else 0.0)
                if score > best_score:
                    best_score = score
                    best_gt = gt_id
        if best_gt is not None:
            matched_gt.add(best_gt)
            matched_pred.add(pred_id)

    tp = len(matched_gt)
    fp = len(pred_map) - tp
    fn = len(gt_map) - tp
    return {"tp": tp, "fp": fp, "fn": fn}


def load_ground_truth(truth_path: Path, lib_path: Path) -> dict[str, dict[str, list[Segment]]]:
    truth = pd.read_csv(truth_path, sep="\t")
    truth = truth[truth["is_background"] == False].copy()  # noqa: E712
    truth_ids = set(truth["ecc_ids"].astype(str))

    lib = pd.read_csv(lib_path, sep="\t", usecols=["id", "region"])
    lib = lib[lib["id"].isin(truth_ids)].copy()
    lib["eccdna_type"] = lib["id"].astype(str).str.extract(r"^(\w+DNA)")[0]

    gt_by_type: dict[str, dict[str, list[Segment]]] = defaultdict(dict)
    for _, row in lib.iterrows():
        ecc_id = str(row["id"])
        ecc_type = str(row["eccdna_type"])
        segments = dedupe_segments(parse_regions(row.get("region", "")))
        if segments:
            gt_by_type[ecc_type][ecc_id] = segments
    return dict(gt_by_type)


def load_predictions(pred_path: Path) -> dict[str, dict[str, list[Segment]]]:
    df = pd.read_csv(pred_path)
    if "Regions" not in df.columns or "eccDNA_type" not in df.columns:
        raise ValueError(f"Prediction file missing Regions/eccDNA_type: {pred_path}")

    pred_by_type: dict[str, dict[str, list[Segment]]] = defaultdict(dict)
    for _, row in df.iterrows():
        ecc_id = str(row.get("eccDNA_id", row.get("eccdna_id", "")))
        ecc_type = str(row["eccDNA_type"])
        segments = dedupe_segments(parse_regions(row.get("Regions", "")))
        if segments:
            pred_by_type[ecc_type][ecc_id] = segments
    return dict(pred_by_type)


def summarize_metrics(tp: int, fp: int, fn: int) -> dict[str, float]:
    total_gt = tp + fn
    total_pred = tp + fp
    recall = tp / total_gt if total_gt > 0 else 0.0
    precision = tp / total_pred if total_pred > 0 else 0.0
    f1 = (2 * precision * recall / (precision + recall)) if (precision + recall) > 0 else 0.0
    return {
        "ground_truth": total_gt,
        "predictions": total_pred,
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "recall": round(recall, 4),
        "precision": round(precision, 4),
        "f1": round(f1, 4),
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Molecule-level recall/precision with multi-segment matching"
    )
    parser.add_argument("--truth", required=True, type=Path, help="Truth TSV (HiFi.truth.tsv)")
    parser.add_argument("--lib", required=True, type=Path, help="Truth library CSV (lib.csv)")
    parser.add_argument("--pred", required=True, type=Path, help="Prediction CSV with Regions")
    parser.add_argument("--overlap", type=float, default=0.9, help="Reciprocal overlap threshold")
    parser.add_argument(
        "--require-pred-coverage",
        action="store_true",
        help="Require predicted coverage >= threshold (default: off)",
    )
    parser.add_argument(
        "--mecc-any-site",
        action="store_true",
        help="For MeccDNA, accept match if any single site meets threshold",
    )
    parser.add_argument("--output", type=Path, default=None, help="Output JSON path")
    args = parser.parse_args()

    gt_by_type = load_ground_truth(args.truth, args.lib)
    pred_by_type = load_predictions(args.pred)

    results = {}
    totals = {"tp": 0, "fp": 0, "fn": 0}

    for ecc_type in sorted(set(gt_by_type) | set(pred_by_type)):
        gt_map = gt_by_type.get(ecc_type, {})
        pred_map = pred_by_type.get(ecc_type, {})
        match_mode = "any_site" if (ecc_type == "MeccDNA" and args.mecc_any_site) else "coverage"
        counts = match_molecules(
            gt_map,
            pred_map,
            args.overlap,
            args.require_pred_coverage,
            match_mode,
        )
        metrics = summarize_metrics(counts["tp"], counts["fp"], counts["fn"])
        metrics["match_mode"] = match_mode
        results[ecc_type] = metrics
        totals["tp"] += counts["tp"]
        totals["fp"] += counts["fp"]
        totals["fn"] += counts["fn"]

    overall = summarize_metrics(totals["tp"], totals["fp"], totals["fn"])

    report = {
        "overlap_threshold": args.overlap,
        "require_pred_coverage": bool(args.require_pred_coverage),
        "by_type": results,
        "overall": overall,
    }

    print(json.dumps(report, indent=2))
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2))
        print(f"\nSaved: {args.output}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
