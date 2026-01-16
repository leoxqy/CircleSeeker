"""Validation utilities for CircleSeeker."""

from __future__ import annotations

import shutil
import importlib
def validate_installation(full_check: bool = False) -> list[str]:
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

    # If core Python dependencies are missing, importing CircleSeeker modules will be noisy
    # and redundant. Report missing modules first and skip deeper checks.
    if any(issue.startswith("Missing Python module:") for issue in issues):
        return issues

    # Check external tools (basic check)
    if full_check:
        # CircleSeeker external tools with alternative names
        external_tools = [
            ("TideHunter", ["tidehunter"]),  # bioconda may install as lowercase
            ("minimap2", None),  # Required for sequence alignment
            ("cd-hit-est", None),  # Required for sequence clustering
            ("samtools", None),
            ("bcftools", None),  # Optional: required for Cyrcular
            ("cyrcular", None),  # Optional: inference engine (fallback)
            ("varlociraptor", None),  # Optional: required for Cyrcular
            ("cresil", None),  # Optional: inference engine (preferred)
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
        importlib.import_module("circleseeker.core.pipeline")
        importlib.import_module("circleseeker.config")
        importlib.import_module("circleseeker.modules.tandem_to_ring")
        importlib.import_module("circleseeker.modules.um_classify")
    except ImportError as e:
        issues.append(f"CircleSeeker module import error: {e}")

    return issues
