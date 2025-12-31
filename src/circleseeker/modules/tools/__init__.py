"""
CircleSeeker tool modules

This package contains external tool integration, circular DNA detection
and inferred eccDNA curation modules.
"""

# Tool modules - import from parent package (circleseeker.modules)
from circleseeker.modules.cyrcular_calling import CyrcularCallingPipeline
from circleseeker.modules.iecc_curator import (
    curate_ecc_tables,
    generate_fasta_sequences,
    write_curated_tables,
    write_curated_tables_with_fasta,
    process_eccDNA,
)
from circleseeker.modules.external_tools import (
    TideHunterModule,
    BlastModule,
    CDHitModule,
    Minimap2Module,
)
from circleseeker.modules.adapters import (
    CLIModuleAdapter,
    TandemToRingAdapter,
    UMClassifyAdapter,
)

__all__ = [
    # Cyrcular pipeline
    "CyrcularCallingPipeline",
    # Inferred eccDNA curation functions
    "curate_ecc_tables",
    "generate_fasta_sequences",
    "write_curated_tables",
    "write_curated_tables_with_fasta",
    "process_eccDNA",
    # External tool modules
    "TideHunterModule",
    "BlastModule",
    "CDHitModule",
    "Minimap2Module",
    # CLI adapters
    "CLIModuleAdapter",
    "TandemToRingAdapter",
    "UMClassifyAdapter",
]
