"""External tool wrappers (CircleSeeker).

This package provides Python wrappers for external bioinformatics tools:
- TideHunter: Tandem repeat detection
- Minimap2: Sequence alignment
- Samtools: BAM/SAM manipulation
- CDHitEst: Sequence clustering

Note: SplitReads-Core is built-in and doesn't require external tools.
"""

from circleseeker.external.base import ExternalTool
from circleseeker.external.tidehunter import TideHunter
from circleseeker.external.minimap2 import Minimap2
from circleseeker.external.minimap2_align import Minimap2Aligner
from circleseeker.external.samtools import Samtools
from circleseeker.external.cd_hit import CDHitEst

__all__ = [
    "ExternalTool",
    "TideHunter",
    "Minimap2",
    "Minimap2Aligner",
    "Samtools",
    "CDHitEst",
]
