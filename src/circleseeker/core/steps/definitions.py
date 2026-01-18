"""Canonical step ordering and user-facing metadata."""

from __future__ import annotations

from circleseeker.core.pipeline_types import PipelineStep


# Final step order; deprecated/unused steps removed.
# `skip_condition` maps to Config attributes (canonical names via properties).
PIPELINE_STEPS: list[PipelineStep] = [
    PipelineStep("check_dependencies", "Check required tools and dependencies"),
    PipelineStep(
        "tidehunter",
        "Run TideHunter for tandem repeat detection",
        skip_condition="skip_tidehunter",
    ),
    PipelineStep(
        "tandem_to_ring",
        "Process TideHunter output",
        skip_condition="skip_tandem_to_ring",
    ),
    PipelineStep("run_alignment", "Run minimap2 alignment"),
    PipelineStep("um_classify", "Classify eccDNA types"),
    PipelineStep("cecc_build", "Process Cecc candidates"),
    PipelineStep("umc_process", "Generate FASTA files"),
    PipelineStep("cd_hit", "Remove redundant sequences"),
    PipelineStep("ecc_dedup", "Coordinate results"),
    PipelineStep("read_filter", "Filter sequences"),
    PipelineStep("minimap2", "Prepare alignments / reference index"),
    PipelineStep(
        "ecc_inference",
        "Detect circular DNA (Cresil preferred, Cyrcular fallback)",
        display_name="ecc_inference",
    ),
    PipelineStep("curate_inferred_ecc", "Curate inferred eccDNA tables"),
    PipelineStep("ecc_unify", "Merge eccDNA tables into unified output"),
    PipelineStep("ecc_summary", "Generate summary report"),
    PipelineStep("ecc_packager", "Package output files", skip_condition="skip_organize"),
]

