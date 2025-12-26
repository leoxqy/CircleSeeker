"""Validation utilities for CircleSeeker."""

from __future__ import annotations

import shutil
import importlib
from typing import List


def validate_installation(full_check: bool = False) -> List[str]:
    """
    Validate CircleSeeker installation and dependencies.

    Args:
        full_check: If True, perform comprehensive validation

    Returns:
        List of validation issues (empty if all good)
    """
    issues = []

    # Check Python modules
    required_modules = ["pandas", "numpy", "biopython", "pysam", "yaml", "networkx", "click"]

    for module in required_modules:
        try:
            if module == "biopython":
                importlib.import_module("Bio")
            elif module == "yaml":
                importlib.import_module("yaml")
            else:
                importlib.import_module(module)
        except ImportError:
            issues.append(f"Missing Python module: {module}")

    # Check external tools (basic check)
    if full_check:
        # CircleSeeker external tools with alternative names
        external_tools = [
            ("TideHunter", ["tidehunter"]),  # bioconda may install as lowercase
            ("makeblastdb", None),
            ("blastn", None),
            ("cd-hit-est", None),  # Required for sequence clustering
            ("minimap2", None),
            ("samtools", None),
            ("bcftools", None),
            ("cyrcular", None),
            ("varlociraptor", None),
        ]

        for tool_name, alt_names in external_tools:
            found = shutil.which(tool_name)
            if not found and alt_names:
                for alt in alt_names:
                    found = shutil.which(alt)
                    if found:
                        break
            if not found:
                issues.append(f"External tool not found: {tool_name}")

    # Check CircleSeeker modules
    try:
        from circleseeker.core.pipeline import Pipeline  # noqa: F401
        from circleseeker.config import Config  # noqa: F401
        from circleseeker.modules.tandem_to_ring import TandemToRing  # noqa: F401
        from circleseeker.modules.um_classify import UMeccClassifier  # noqa: F401
    except ImportError as e:
        issues.append(f"CircleSeeker module import error: {e}")

    return issues
