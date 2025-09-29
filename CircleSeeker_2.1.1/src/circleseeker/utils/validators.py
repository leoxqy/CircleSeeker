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
        # CircleSeeker 2.1.0 external tools (post-alignment pipeline)
        external_tools = [
            'TideHunter', 'makeblastdb', 'blastn', 'minimap2', 'samtools',
            'bcftools', 'cyrcular', 'varlociraptor'
        ]
        
        for tool in external_tools:
            if not shutil.which(tool):
                issues.append(f"External tool not found: {tool}")
    
    # Check CircleSeeker modules
    try:
        from circleseeker.core.pipeline import Pipeline
        from circleseeker.config import Config
        from circleseeker.modules.tandem_to_ring import TandemToRing
        from circleseeker.modules.um_classify import UMeccClassifier
    except ImportError as e:
        issues.append(f"CircleSeeker module import error: {e}")
    
    return issues
