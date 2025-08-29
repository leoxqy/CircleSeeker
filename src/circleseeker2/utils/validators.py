"""Validation utilities for CircleSeeker2."""

from __future__ import annotations

import shutil
import importlib
from typing import List


def validate_installation(full_check: bool = False) -> List[str]:
    """
    Validate CircleSeeker2 installation and dependencies.
    
    Args:
        full_check: If True, perform comprehensive validation
        
    Returns:
        List of validation issues (empty if all good)
    """
    issues = []
    
    # Check Python modules
    required_modules = [
        'pandas', 'numpy', 'biopython', 'pysam', 'yaml', 'networkx', 'click'
    ]
    
    for module in required_modules:
        try:
            if module == 'biopython':
                importlib.import_module('Bio')
            elif module == 'yaml':
                importlib.import_module('yaml')
            else:
                importlib.import_module(module)
        except ImportError:
            issues.append(f"Missing Python module: {module}")
    
    # Check external tools (basic check)
    if full_check:
        external_tools = ['TideHunter', 'makeblastdb', 'blastn', 'minimap2', 'samtools', 'mosdepth']
        
        for tool in external_tools:
            if not shutil.which(tool):
                issues.append(f"External tool not found: {tool}")
    
    # Check CircleSeeker2 modules
    try:
        from circleseeker2.core.pipeline import Pipeline
        from circleseeker2.config import Config
        from circleseeker2.modules.carousel import Carousel
        from circleseeker2.modules.gatekeeper import GatekeeperClassifier
    except ImportError as e:
        issues.append(f"CircleSeeker2 module import error: {e}")
    
    return issues