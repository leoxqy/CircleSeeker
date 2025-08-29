#!/usr/bin/env python3
"""Lightweight helpers to compare outputs between new modules and legacy/expected files.

These functions are intentionally permissive to account for minor numeric or ordering
differences while still catching regressions.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable, Tuple


def file_nonempty(path: Path) -> bool:
    return path.exists() and path.stat().st_size > 0


def compare_csv_row_counts(a: Path, b: Path, tolerance: float = 0.1) -> Tuple[bool, int, int]:
    """Compare CSV row counts within tolerance fraction (default 10%)."""
    def count_rows(p: Path) -> int:
        with open(p, newline="") as f:
            return sum(1 for _ in csv.reader(f)) - 1  # minus header

    if not (a.exists() and b.exists()):
        return False, 0, 0

    ca, cb = count_rows(a), count_rows(b)
    if cb == 0:
        return (ca == 0, ca, cb)
    ok = abs(ca - cb) / max(cb, 1) <= tolerance
    return ok, ca, cb


def compare_csv_columns(a: Path, required: Iterable[str]) -> Tuple[bool, list[str]]:
    """Verify CSV has required columns."""
    with open(a, newline="") as f:
        reader = csv.DictReader(f)
        cols = list(reader.fieldnames or [])
    missing = [c for c in required if c not in cols]
    return len(missing) == 0, missing


def count_fasta_records(path: Path) -> int:
    count = 0
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count


def compare_fasta_counts(a: Path, b: Path, tolerance: float = 0.1) -> Tuple[bool, int, int]:
    if not (a.exists() and b.exists()):
        return False, 0, 0
    ca, cb = count_fasta_records(a), count_fasta_records(b)
    if cb == 0:
        return (ca == 0, ca, cb)
    ok = abs(ca - cb) / max(cb, 1) <= tolerance
    return ok, ca, cb

