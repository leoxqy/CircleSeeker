"""CircleSeeker2: Comprehensive eccDNA detection from HiFi sequencing data."""

from circleseeker2.__version__ import __version__, __author__, __email__, __license__
from circleseeker2.config import Config
from circleseeker2.exceptions import CircleSeekerError

__all__ = [
    "__version__",
    "__author__",
    "__email__",
    "__license__",
    "Config",
    "CircleSeekerError",
]
