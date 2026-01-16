from __future__ import annotations

import re
from pathlib import Path

import circleseeker


def _read_pyproject_version(pyproject_path: Path) -> str:
    in_project = False
    for raw_line in pyproject_path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if line == "[project]":
            in_project = True
            continue
        if in_project and line.startswith("[") and line.endswith("]"):
            break
        if not in_project:
            continue
        match = re.match(r'^version\s*=\s*"([^"]+)"\s*$', line)
        if match:
            return match.group(1)
    raise AssertionError(f"Could not find [project].version in {pyproject_path}")


def test_version_matches_pyproject() -> None:
    repo_root = Path(__file__).resolve().parents[2]
    pyproject_version = _read_pyproject_version(repo_root / "pyproject.toml")
    assert circleseeker.__version__ == pyproject_version
