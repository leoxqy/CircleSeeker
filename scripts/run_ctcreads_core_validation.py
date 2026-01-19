#!/usr/bin/env python3
"""
Server-side CtcReads-Caller core validation runner for CircleSeeker simulation tests.

Runs the existing batch validator with server-friendly defaults:
  - runs per scale: 10
  - scales: 100..100000 (configurable)
  - writes all outputs to a timestamped directory
  - captures full stdout/stderr to a log file
  - optionally archives results to a .tar.gz for transfer

Scope note:
  This validates the CtcReads-Caller *core* logic (minimap2 alignment -> U/M classification
  -> Cecc build) and intentionally skips external CtcReads detection tools
  (TideHunter + tandem_to_ring).

Usage (from repo root or unpacked sdist):
  python scripts/run_ctcreads_core_validation.py --threads 16
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import subprocess
import sys
import tarfile
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional


def _run_capture(cmd: list[str], *, cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, cwd=str(cwd), capture_output=True, text=True)


def _tool_version(cmd: list[str], *, cwd: Path) -> Optional[str]:
    try:
        proc = _run_capture(cmd, cwd=cwd)
    except FileNotFoundError:
        return None
    text = (proc.stdout or proc.stderr or "").strip()
    return text.splitlines()[0].strip() if text else None


def _git_head(repo_root: Path) -> Optional[str]:
    proc = _run_capture(["git", "rev-parse", "HEAD"], cwd=repo_root)
    return proc.stdout.strip() if proc.returncode == 0 else None


def _circleseeker_version(repo_root: Path) -> Optional[str]:
    try:
        sys.path.insert(0, str(repo_root / "src"))
        import circleseeker  # type: ignore

        return getattr(circleseeker, "__version__", None) or getattr(circleseeker, "version", None)
    except Exception:
        return None


def _archive_directory(src_dir: Path, archive_path: Path) -> Path:
    archive_path.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive_path, "w:gz") as tf:
        tf.add(src_dir, arcname=src_dir.name)
    return archive_path


def main() -> int:
    repo_root = Path(__file__).resolve().parents[1]

    parser = argparse.ArgumentParser(
        description="Run CircleSeeker CtcReads-Caller core validation on a server"
    )
    parser.add_argument("--threads", type=int, default=os.cpu_count() or 8, help="Threads for minimap2")
    parser.add_argument("--runs", type=int, default=10, help="Runs per scale (default: 10)")
    parser.add_argument("--base-seed", type=int, default=42, help="Base random seed (default: 42)")
    parser.add_argument(
        "--scales",
        type=str,
        default="100,500,1000,3000,5000,10000,50000,100000",
        help="Comma-separated scales (default includes 100000)",
    )
    parser.add_argument(
        "--minimap2-timeout",
        type=int,
        default=7200,
        help="minimap2 timeout per run in seconds (default: 7200)",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("server_validation_results"),
        help="Root folder to store results (default: server_validation_results)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Optional explicit output directory (enables resume across re-runs)",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Skip runs with an existing non-empty validation report",
    )
    parser.add_argument(
        "--archive",
        action="store_true",
        help="Create a .tar.gz archive after completion",
    )
    parser.add_argument(
        "--archive-path",
        type=Path,
        default=None,
        help="Optional archive output path (default: <output_dir>.tar.gz)",
    )

    args = parser.parse_args()

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = args.output_dir or (args.output_root / f"ctcreads_core_validation_{ts}")
    output_dir.mkdir(parents=True, exist_ok=True)

    meta: Dict[str, Any] = {
        "started_at": datetime.now().isoformat(),
        "repo_root": str(repo_root),
        "hostname": platform.node(),
        "platform": platform.platform(),
        "python": sys.version,
        "git_head": _git_head(repo_root),
        "circleseeker_version": _circleseeker_version(repo_root),
        "minimap2_version": _tool_version(["minimap2", "--version"], cwd=repo_root),
        "config": {
            "threads": args.threads,
            "runs": args.runs,
            "base_seed": args.base_seed,
            "scales": args.scales,
            "minimap2_timeout": args.minimap2_timeout,
            "resume": args.resume,
        },
    }

    meta_path = output_dir / f"run_metadata_{ts}.json"
    meta_path.write_text(json.dumps(meta, indent=2) + "\n")

    log_path = output_dir / f"ctcreads_core_validation_{ts}.log"
    batch_cmd = [
        sys.executable,
        "-u",
        str(repo_root / "tests" / "simulation" / "batch_validation.py"),
        "--threads",
        str(args.threads),
        "--output-dir",
        str(output_dir),
        "--scales",
        args.scales,
        "--runs",
        str(args.runs),
        "--base-seed",
        str(args.base_seed),
        "--minimap2-timeout",
        str(args.minimap2_timeout),
    ]
    if args.resume:
        batch_cmd.append("--resume")

    with log_path.open("w") as log_file:
        log_file.write(f"$ {' '.join(batch_cmd)}\n\n")
        log_file.flush()
        proc = subprocess.run(batch_cmd, cwd=str(repo_root), stdout=log_file, stderr=subprocess.STDOUT, text=True)

    meta["finished_at"] = datetime.now().isoformat()
    meta["returncode"] = proc.returncode
    meta["output_dir"] = str(output_dir)
    meta["log_path"] = str(log_path)
    meta_path.write_text(json.dumps(meta, indent=2) + "\n")

    if args.archive:
        archive_path = args.archive_path or output_dir.with_suffix(".tar.gz")
        _archive_directory(output_dir, archive_path)
        print(f"Archived: {archive_path}")

    print(f"Done. Results: {output_dir}")
    print(f"Log: {log_path}")
    return int(proc.returncode)


if __name__ == "__main__":
    raise SystemExit(main())
