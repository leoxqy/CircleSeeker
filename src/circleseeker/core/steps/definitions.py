"""Canonical step ordering and user-facing metadata."""

from __future__ import annotations

from circleseeker.core.pipeline_types import PipelineStep


# Final step order; deprecated/unused steps removed.
# `skip_condition` maps to Config attributes (canonical names via properties).
# `depends_on` declares explicit step dependencies for validation and documentation.
PIPELINE_STEPS: list[PipelineStep] = [
    PipelineStep(
        "check_dependencies",
        "Check required tools and dependencies",
        depends_on=[],
    ),
    PipelineStep(
        "tidehunter",
        "Run TideHunter for tandem repeat detection",
        skip_condition="skip_tidehunter",
        depends_on=["check_dependencies"],
    ),
    PipelineStep(
        "tandem_to_ring",
        "Process TideHunter output",
        skip_condition="skip_tandem_to_ring",
        depends_on=["tidehunter"],
    ),
    PipelineStep(
        "run_alignment",
        "Run minimap2 alignment",
        depends_on=["tandem_to_ring"],
    ),
    PipelineStep(
        "um_classify",
        "Classify eccDNA types",
        depends_on=["run_alignment"],
    ),
    PipelineStep(
        "cecc_build",
        "Process Cecc candidates",
        depends_on=["um_classify"],
    ),
    PipelineStep(
        "umc_process",
        "Generate FASTA files",
        depends_on=["um_classify", "cecc_build"],
    ),
    PipelineStep(
        "cd_hit",
        "Remove redundant sequences",
        depends_on=["umc_process"],
    ),
    PipelineStep(
        "ecc_dedup",
        "Coordinate results",
        depends_on=["cd_hit"],
    ),
    PipelineStep(
        "read_filter",
        "Filter sequences",
        depends_on=["tandem_to_ring"],
    ),
    PipelineStep(
        "minimap2",
        "Prepare alignments / reference index",
        depends_on=["read_filter"],
    ),
    PipelineStep(
        "ecc_inference",
        "Detect circular DNA using SplitReads-Core (built-in)",
        display_name="ecc_inference",
        depends_on=["minimap2"],
    ),
    PipelineStep(
        "curate_inferred_ecc",
        "Curate inferred eccDNA tables",
        depends_on=["ecc_inference"],
    ),
    PipelineStep(
        "ecc_unify",
        "Merge eccDNA tables into unified output",
        depends_on=["ecc_dedup", "curate_inferred_ecc"],
    ),
    PipelineStep(
        "ecc_summary",
        "Generate summary report",
        depends_on=["ecc_unify"],
    ),
    PipelineStep(
        "ecc_packager",
        "Package output files",
        skip_condition="skip_organize",
        depends_on=["ecc_summary"],
    ),
]

