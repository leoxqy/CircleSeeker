"""CircleSeeker2: Comprehensive eccDNA detection from HiFi sequencing data.

This package provides tools for detecting extrachromosomal circular DNA (eccDNA) 
from PacBio HiFi long-read sequencing data using a multi-round iterative analysis approach.
"""

from circleseeker2.__version__ import (
    __version__,
    __author__,
    __email__,
    __license__,
    __description__,
)

# Import main classes when needed
try:
    from circleseeker2.core.pipeline import Pipeline
    from circleseeker2.config import Config
    from circleseeker2.exceptions import CircleSeekerError
    
    __all__ = [
        "__version__",
        "__author__", 
        "__email__",
        "__license__",
        "__description__",
        "Pipeline",
        "Config", 
        "CircleSeekerError",
    ]
except ImportError:
    # Allow partial imports during development
    __all__ = [
        "__version__",
        "__author__",
        "__email__", 
        "__license__",
        "__description__",
    ]