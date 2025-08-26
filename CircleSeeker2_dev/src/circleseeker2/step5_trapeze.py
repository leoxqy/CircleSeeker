#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Circular pipeline — minimal, parameterized CLI.

Reads an alignment CSV grouped by `query_id`, then:
  1) Rotates segments within each query (by q_start), adding `segment_order`.
  2) Filters queries: drop if segments <= min_segments or any pair overlaps ≥ overlap_threshold (relative to the shorter segment).
  3) Detects circular paths with continuity (|gap| ≤ edge_tolerance) and position closure (same chr/strand and s_start/s_end within tolerance).
  4) Labels segments in each circle as head/body/tail.

Writes ONLY the final CSV (labeled circular results). No intermediate files.

Dependencies: pandas, numpy, optional numba.
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

# Optional acceleration
try:
    from numba import jit  # type: ignore
    NUMBA_AVAILABLE = True
except Exception:  # pragma: no cover
    NUMBA_AVAILABLE = False

LOGGER = logging.getLogger("circular_pipeline")

# Columns required for the algorithm
REQUIRED_COLS = [
    "query_id", "subject_id", "q_start", "q_end", "s_start", "s_end",
    "strand", "alignment_length", "consLen", "readName", "copyNum",
]


def ensure_required_columns(df: pd.DataFrame) -> pd.DataFrame:
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"Input CSV missing required columns: {missing}")
    # Cast numeric columns
    for col in ["q_start", "q_end", "s_start", "s_end", "alignment_length", "consLen", "copyNum"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=["q_start", "q_end", "s_start", "s_end", "alignment_length", "consLen"]).copy()
    return df


def rotate_group(g: pd.DataFrame) -> pd.DataFrame:
    """Sort by q_start, move first row to the end, add 1-based `segment_order`."""
    s = g.sort_values("q_start").reset_index(drop=True)
    if len(s) > 1:
        s = pd.concat([s.iloc[1:], s.iloc[:1]], ignore_index=True)
    s = s.copy()
    s["segment_order"] = np.arange(1, len(s) + 1)
    return s


def rotate_all(df: pd.DataFrame) -> pd.DataFrame:
    parts = [rotate_group(g) for _, g in df.groupby("query_id", sort=False)]
    return pd.concat(parts, ignore_index=True)


# Overlap check (pairwise; relative to shorter segment)
if NUMBA_AVAILABLE:
    @jit(nopython=True)  # type: ignore
    def _overlap_numba(qs, qe, thr=0.8):
        n = len(qs)
        max_r = 0.0
        for i in range(n):
            for j in range(i + 1, n):
                os = qs[i] if qs[i] > qs[j] else qs[j]
                oe = qe[i] if qe[i] < qe[j] else qe[j]
                if os < oe:
                    ov = oe - os
                    l1 = qe[i] - qs[i]
                    l2 = qe[j] - qs[j]
                    m = l1 if l1 < l2 else l2
                    if m > 0:
                        r = ov / m
                        if r > max_r:
                            max_r = r
                        if r >= thr:
                            return True, max_r
        return False, max_r
else:
    def _overlap_numba(qs, qe, thr=0.8):
        n = len(qs)
        max_r = 0.0
        for i in range(n):
            for j in range(i + 1, n):
                os = max(qs[i], qs[j])
                oe = min(qe[i], qe[j])
                if os < oe:
                    ov = oe - os
                    l1 = qe[i] - qs[i]
                    l2 = qe[j] - qs[j]
                    m = min(l1, l2)
                    if m > 0:
                        r = ov / m
                        max_r = max(max_r, r)
                        if r >= thr:
                            return True, max_r
        return False, max_r


def vectorized_overlap_filter(df_rot: pd.DataFrame, overlap_threshold: float, min_segments: int) -> pd.DataFrame:
    keep_idx = []
    for qid, g in df_rot.groupby("query_id", sort=False):
        if len(g) <= min_segments:
            continue
        qs = g["q_start"].to_numpy(np.float64)
        qe = g["q_end"].to_numpy(np.float64)
        has, _ = _overlap_numba(qs, qe, overlap_threshold)
        if not has:
            keep_idx.append(g.index)
    if not keep_idx:
        return df_rot.iloc[0:0].copy()
    return df_rot.loc[np.concatenate(keep_idx)].reset_index(drop=True)


def position_match(a: pd.Series, b: pd.Series, tol: int) -> bool:
    if a["subject_id"] != b["subject_id"]:
        return False
    if str(a["strand"]) != str(b["strand"]):
        return False
    return (abs(float(a["s_start"]) - float(b["s_start"])) <= tol and
            abs(float(a["s_end"]) - float(b["s_end"])) <= tol)


def analyze_gaps(path_df: pd.DataFrame) -> Tuple[float, float, int, bool]:
    gaps = []
    for i in range(1, len(path_df)):
        prev = path_df.iloc[i - 1]
        cur = path_df.iloc[i]
        gaps.append(float(cur["q_start"]) - float(prev["q_end"]))
    if not gaps:
        return 0.0, 0.0, 0, True
    max_gap = max(abs(g) for g in gaps)
    avg_gap = sum(abs(g) for g in gaps) / len(gaps)
    return max_gap, avg_gap, len(gaps), max_gap <= 20


def find_circular(group_df: pd.DataFrame, edge_tol: int, pos_tol: int):
    g = group_df.sort_values("segment_order").reset_index(drop=True)
    if len(g) < 2:
        return None
    cons_len = float(g.loc[0, "consLen"])

    path = [0]
    cum_len = float(g.loc[0, "alignment_length"])

    for idx in range(1, len(g)):
        prev = g.iloc[path[-1]]
        cur = g.iloc[idx]
        gap = float(cur["q_start"]) - float(prev["q_end"]) 
        if abs(gap) > edge_tol:
            return None  # break in continuity
        path.append(idx)
        cum_len += float(cur["alignment_length"])

        if cum_len >= cons_len:
            first = g.iloc[0]
            # close at current
            if position_match(first, cur, pos_tol):
                return {
                    "path": path[:-1],
                    "closing_at": idx,
                    "cum_len": cum_len - float(cur["alignment_length"]),
                    "mat_degree": round(((cum_len - float(cur["alignment_length"])) / cons_len) * 100, 2),
                }
            # close at next (if still continuous)
            if idx < len(g) - 1:
                nxt = g.iloc[idx + 1]
                next_gap = float(nxt["q_start"]) - float(cur["q_end"]) 
                if abs(next_gap) <= edge_tol and position_match(first, nxt, pos_tol):
                    return {
                        "path": path,
                        "closing_at": idx + 1,
                        "cum_len": cum_len,
                        "mat_degree": round((cum_len / cons_len) * 100, 2),
                    }
    return None


def detect_circles(df_filt: pd.DataFrame, edge_tol: int, pos_tol: int) -> pd.DataFrame:
    rows: List[Dict] = []
    for qid, g in df_filt.groupby("query_id", sort=False):
        res = find_circular(g, edge_tol, pos_tol)
        if not res:
            continue
        path_df = g.iloc[res["path"]].reset_index(drop=True)
        chroms = path_df["subject_id"].unique()
        eclass = "Cecc-InterChr" if len(chroms) > 1 else "Cecc-IntraChr"
        max_gap, avg_gap, _, _ = analyze_gaps(path_df)
        base = {
            "query_id": qid,
            "readName": g.iloc[0]["readName"],
            "eClass": eclass,
            "consLen": g.iloc[0]["consLen"],
            "copyNum": g.iloc[0]["copyNum"],
            "num_segments": len(path_df),
            "cumulative_length": res["cum_len"],
            "MatDegree": res["mat_degree"],
            "max_gap": max_gap,
            "avg_gap": round(avg_gap, 2),
            "chromosomes": ",".join(chroms.astype(str)),
        }
        for i in range(len(path_df)):
            seg = path_df.iloc[i]
            row = dict(base)
            row.update({
                "segment_in_circle": i + 1,
                "chr": seg["subject_id"],
                "start": seg["s_start"],
                "end": seg["s_end"],
                "strand": seg["strand"],
                "q_start": seg["q_start"],
                "q_end": seg["q_end"],
                "alignment_length": seg["alignment_length"],
            })
            rows.append(row)
    return pd.DataFrame(rows)


def label_roles(df_circ: pd.DataFrame) -> pd.DataFrame:
    if df_circ.empty:
        return df_circ.copy()
    # Stable sort; avoid pandas groupby.apply deprecation warnings
    df = df_circ.sort_values(["query_id", "segment_in_circle"], kind="mergesort").copy()

    parts: List[pd.DataFrame] = []
    for _, g in df.groupby("query_id", sort=False):
        n = len(g)
        roles = ["body"] * n
        if n >= 1:
            roles[0] = "head"
        if n >= 2:
            roles[-1] = "tail"
        g2 = g.copy()
        g2["segment_role"] = pd.Categorical(roles, categories=["head", "body", "tail"], ordered=True)
        parts.append(g2)

    df = pd.concat(parts, ignore_index=True)
    cols = df.columns.tolist()
    if "segment_role" in cols:
        cols.remove("segment_role")
        i = cols.index("segment_in_circle") + 1
        cols = cols[:i] + ["segment_role"] + cols[i:]
        df = df[cols]
    return df


def run_pipeline(
    input_csv: Path,
    output_csv: Path | None,
    overlap_threshold: float,
    min_segments: int,
    edge_tolerance: int,
    position_tolerance: int,
    verbose: bool,
) -> Path:
    # Logging
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.WARNING,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

    # Read
    df = pd.read_csv(input_csv, sep=None, engine="python")
    df = ensure_required_columns(df)

    # Rotate and filter in-memory
    df_rot = rotate_all(df)
    df_filt = vectorized_overlap_filter(df_rot, overlap_threshold, min_segments)

    # Detect circles and label
    df_circ = detect_circles(df_filt, edge_tolerance, position_tolerance)
    df_labeled = label_roles(df_circ)

    # Output
    if output_csv is None:
        output_csv = input_csv.with_name(f"{input_csv.stem}_circular_labeled.csv")
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    df_labeled.to_csv(output_csv, index=False)

    # Single, clean stdout line
    print(f"Saved: {output_csv}")
    return output_csv


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Detect circular structures from per-segment alignments. Writes only the final labeled CSV.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-i", "--input", type=Path, required=True, help="Input CSV with per-segment alignments")
    p.add_argument("-o", "--output", type=Path, default=None, help="Output CSV path (final labeled results)")
    p.add_argument("--overlap-threshold", type=float, default=0.8, help="Overlap threshold to discard a query (relative to shorter segment)")
    p.add_argument("--min-segments", type=int, default=3, help="Discard queries with segments <= this value")
    p.add_argument("--edge-tolerance", type=int, default=20, help="Max allowed |q-gap| between adjacent segments")
    p.add_argument("--position-tolerance", type=int, default=50, help="Genome position tolerance for closure (bp)")
    p.add_argument("-v", "--verbose", action="store_true", help="Enable debug logging")
    return p


def main(argv: List[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    if not args.input.exists():
        LOGGER.error("Input file not found: %s", args.input)
        return 2
    try:
        run_pipeline(
            input_csv=args.input,
            output_csv=args.output,
            overlap_threshold=args.overlap_threshold,
            min_segments=args.min_segments,
            edge_tolerance=args.edge_tolerance,
            position_tolerance=args.position_tolerance,
            verbose=args.verbose,
        )
    except Exception as e:  # pragma: no cover
        LOGGER.exception("Failed: %s", e)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

