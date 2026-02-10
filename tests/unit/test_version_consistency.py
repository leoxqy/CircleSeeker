from __future__ import annotations

import re
from pathlib import Path

import circleseeker


def _read_pyproject_version(pyproject_path: Path) -> str | None:
    """Read version from pyproject.toml.

    Returns the static version string if present, or None if version is dynamic
    (i.e. read from __version__.py at build time via setuptools).
    """
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
        # Check for dynamic = ["version"]
        if re.match(r'^dynamic\s*=\s*\[.*"version".*\]', line):
            return None
        match = re.match(r'^version\s*=\s*"([^"]+)"\s*$', line)
        if match:
            return match.group(1)
    raise AssertionError(f"Could not find [project].version or dynamic version in {pyproject_path}")


def test_version_matches_pyproject() -> None:
    repo_root = Path(__file__).resolve().parents[2]
    pyproject_version = _read_pyproject_version(repo_root / "pyproject.toml")
    if pyproject_version is None:
        # Version is dynamic — just verify __version__ is a valid PEP 440 string
        from packaging.version import Version

        Version(circleseeker.__version__)  # raises InvalidVersion if malformed
    else:
        assert circleseeker.__version__ == pyproject_version


def test_version_is_valid_pep440() -> None:
    """__version__ must be a valid PEP 440 version string."""
    from packaging.version import Version

    v = Version(circleseeker.__version__)
    assert v.release  # must have at least a major.minor.micro tuple
