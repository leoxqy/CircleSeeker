"""External tool wrappers (CircleSeeker).

This package provides Python wrappers for external bioinformatics tools:
- TideHunter: Tandem repeat detection
- Minimap2: Sequence alignment
- Samtools: BAM/SAM manipulation
- CDHitEst: Sequence clustering
- Bcftools: Variant calling utilities
- Varlociraptor: Variant calling
- Cyrcular: Circular DNA detection
- Cresil: Circular DNA detection (alternative)
"""

from circleseeker.external.base import ExternalTool
from circleseeker.external.tidehunter import TideHunter
from circleseeker.external.minimap2 import Minimap2
from circleseeker.external.minimap2_align import Minimap2Aligner
from circleseeker.external.samtools import Samtools
from circleseeker.external.cd_hit import CDHitEst
from circleseeker.external.bcftools import Bcftools
from circleseeker.external.varlociraptor import Varlociraptor
from circleseeker.external.cyrcular import Cyrcular
from circleseeker.external.cresil import Cresil

__all__ = [
    "ExternalTool",
    "TideHunter",
    "Minimap2",
    "Minimap2Aligner",
    "Samtools",
    "CDHitEst",
    "Bcftools",
    "Varlociraptor",
    "Cyrcular",
    "Cresil",
]
