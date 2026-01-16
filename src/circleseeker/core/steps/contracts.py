"""Step IO contracts.

Contracts are intentionally lightweight and declarative:
- File naming (relative to output/final dirs)
- Table column requirements
- Coordinate system notes (documented; minimally validated)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional


@dataclass(frozen=True)
class ArtifactSpec:
    """Declarative description of an input/output artifact."""

    name: str
    base: str
    template: str
    kind: str = "file"  # file/dir/csv/tsv/tsv_no_header/fasta/bam/json/txt
    required: bool = False
    sep: str = ","
    required_columns: tuple[str, ...] = field(default_factory=tuple)
    coordinate_system: Optional[str] = None


@dataclass(frozen=True)
class StepContract:
    """Input/output contract for a pipeline step."""

    step_name: str
    inputs: tuple[ArtifactSpec, ...] = field(default_factory=tuple)
    outputs: tuple[ArtifactSpec, ...] = field(default_factory=tuple)


COORD_0_BASED_HALF_OPEN = "0-based, half-open [start0, end0)"


def _out(
    name: str,
    template: str,
    *,
    kind: str = "file",
    required: bool = False,
    **kwargs: Any,
) -> ArtifactSpec:
    return ArtifactSpec(name=name, base="output", template=template, kind=kind, required=required, **kwargs)


def _final(
    name: str,
    template: str,
    *,
    kind: str = "file",
    required: bool = False,
    **kwargs: Any,
) -> ArtifactSpec:
    return ArtifactSpec(name=name, base="final", template=template, kind=kind, required=required, **kwargs)


def _cfg(name: str, attr: str, *, required: bool = False, **kwargs: Any) -> ArtifactSpec:
    return ArtifactSpec(name=name, base="config", template=attr, required=required, **kwargs)


# Minimal, stable column requirements for key tables.
_UECC_MECC_REQUIRED = (
    "query_id",
    "eccdna_type",
    "chr",
    "start0",
    "end0",
    "strand",
    "confidence_score",
    "mapq_best",
    "identity_best",
    "query_cov_best",
    "query_cov_2nd",
    "low_mapq",
    "low_identity",
)

_CECC_REQUIRED = (
    "query_id",
    "eccdna_type",
    "chr",
    "start0",
    "end0",
    "strand",
    "segment_in_circle",
    "confidence_score",
    "mapq_best",
    "identity_best",
    "query_cov_best",
    "query_cov_2nd",
    "low_mapq",
    "low_identity",
)

_CORE_REQUIRED = (
    "eccDNA_id",
    "chr",
    "start0",
    "end0",
    "strand",
    "eccdna_type",
    "confidence_score",
    "mapq_best",
    "identity_best",
    "query_cov_best",
    "query_cov_2nd",
    "low_mapq",
    "low_identity",
)


STEP_CONTRACTS: dict[str, StepContract] = {
    "check_dependencies": StepContract(step_name="check_dependencies"),
    "tidehunter": StepContract(
        step_name="tidehunter",
        inputs=(
            _cfg("input_fasta", "input_file", required=True),
        ),
        outputs=(
            _out("tidehunter_output", "{prefix}.TH.ecc_candidates.txt", kind="txt", required=True),
        ),
    ),
    "tandem_to_ring": StepContract(
        step_name="tandem_to_ring",
        inputs=(
            _out("tidehunter_output", "{prefix}.TH.ecc_candidates.txt", kind="txt"),
        ),
        outputs=(
            _out("tandem_to_ring_csv", "tandem_to_ring.csv", kind="csv", required=True),
            _out("tandem_to_ring_fasta", "tandem_to_ring.fasta", kind="fasta", required=True),
        ),
    ),
    "run_alignment": StepContract(
        step_name="run_alignment",
        inputs=(
            _cfg("reference_fasta", "reference", required=True),
            _out("tandem_to_ring_fasta", "tandem_to_ring.fasta", kind="fasta"),
        ),
        outputs=(
            _out(
                "alignment_results",
                "{prefix}_alignment_results.tsv",
                kind="tsv_no_header",
                required=True,
                sep="\t",
            ),
        ),
    ),
    "um_classify": StepContract(
        step_name="um_classify",
        inputs=(
            _out(
                "alignment_results",
                "{prefix}_alignment_results.tsv",
                kind="tsv_no_header",
                required=True,
                sep="\t",
            ),
        ),
        outputs=(
            _out("um_all", "um_classify.all.csv", kind="csv", required=False),
            _out(
                "um_uecc",
                "um_classify.uecc.csv",
                kind="csv",
                required=False,
                required_columns=_UECC_MECC_REQUIRED,
                coordinate_system=COORD_0_BASED_HALF_OPEN,
            ),
            _out(
                "um_mecc",
                "um_classify.mecc.csv",
                kind="csv",
                required=False,
                required_columns=_UECC_MECC_REQUIRED,
                coordinate_system=COORD_0_BASED_HALF_OPEN,
            ),
            _out("um_unclassified", "um_classify.unclassified.csv", kind="csv", required=False),
        ),
    ),
    "cecc_build": StepContract(
        step_name="cecc_build",
        inputs=(
            _out("um_all", "um_classify.all.csv", kind="csv"),
            _out("um_unclassified", "um_classify.unclassified.csv", kind="csv"),
        ),
        outputs=(
            _out(
                "cecc_build",
                "cecc_build.csv",
                kind="csv",
                required=True,
                required_columns=_CECC_REQUIRED,
                coordinate_system=COORD_0_BASED_HALF_OPEN,
            ),
            _final("ambiguous_uc", "{prefix}_ambiguous_uc.csv", kind="csv", required=False),
            _final("ambiguous_mc", "{prefix}_ambiguous_mc.csv", kind="csv", required=False),
            _final("read_classification", "{prefix}_read_classification.csv", kind="csv", required=False),
        ),
    ),
    "umc_process": StepContract(
        step_name="umc_process",
        outputs=(
            _out("uecc_pre_fasta", "{prefix}_UeccDNA_pre.fasta", kind="fasta", required=False),
            _out("mecc_pre_fasta", "{prefix}_MeccDNA_pre.fasta", kind="fasta", required=False),
            _out("cecc_pre_fasta", "{prefix}_CeccDNA_pre.fasta", kind="fasta", required=False),
            _out("uecc_processed_csv", "{prefix}_UeccDNA_processed.csv", kind="csv", required=False),
            _out("mecc_processed_csv", "{prefix}_MeccDNA_processed.csv", kind="csv", required=False),
            _out("cecc_processed_csv", "{prefix}_CeccDNA_processed.csv", kind="csv", required=False),
        ),
    ),
    "cd_hit": StepContract(
        step_name="cd_hit",
        outputs=(
            _out("uecc_nr99", "{prefix}_U", kind="file", required=False),
            _out("mecc_nr99", "{prefix}_M", kind="file", required=False),
            _out("cecc_nr99", "{prefix}_C", kind="file", required=False),
            _out("uecc_clusters", "{prefix}_U.clstr", kind="txt", required=False),
            _out("mecc_clusters", "{prefix}_M.clstr", kind="txt", required=False),
            _out("cecc_clusters", "{prefix}_C.clstr", kind="txt", required=False),
        ),
    ),
    "ecc_dedup": StepContract(
        step_name="ecc_dedup",
        outputs=(
            _out(
                "uecc_core",
                "{prefix}_UeccDNA.core.csv",
                kind="csv",
                required=False,
                required_columns=_CORE_REQUIRED,
                coordinate_system=COORD_0_BASED_HALF_OPEN,
            ),
            _out(
                "mecc_core",
                "{prefix}_MeccSites.core.csv",
                kind="csv",
                required=False,
                required_columns=_CORE_REQUIRED,
                coordinate_system=COORD_0_BASED_HALF_OPEN,
            ),
            _out(
                "cecc_core",
                "{prefix}_CeccSegments.core.csv",
                kind="csv",
                required=False,
                required_columns=_CORE_REQUIRED,
                coordinate_system=COORD_0_BASED_HALF_OPEN,
            ),
            _out("confirmed", "{prefix}_eccDNA_Confirmed.csv", kind="csv", required=False),
        ),
    ),
    "read_filter": StepContract(
        step_name="read_filter",
        inputs=(
            _cfg("input_fasta", "input_file", required=True),
        ),
        outputs=(
            _out("all_filtered_fasta", "{prefix}_all_filtered.fasta", kind="fasta", required=True),
        ),
    ),
    "minimap2": StepContract(
        step_name="minimap2",
        outputs=(
            _out("sorted_bam", "{prefix}_sorted.bam", kind="bam", required=False),
            _out("all_filtered_fasta", "{prefix}_all_filtered.fasta", kind="fasta", required=False),
        ),
    ),
    "ecc_inference": StepContract(
        step_name="ecc_inference",
        outputs=(
            _out("overview_tsv", "{prefix}_overview.tsv", kind="tsv", required=False, sep="\t"),
            _out("inference_failed_marker", "{prefix}_inference_failed.txt", kind="txt", required=False),
        ),
    ),
    "curate_inferred_ecc": StepContract(
        step_name="curate_inferred_ecc",
        outputs=(
            _out("inferred_dir", "{prefix}_Inferred_eccDNA", kind="dir", required=False),
        ),
    ),
    "ecc_unify": StepContract(
        step_name="ecc_unify",
        outputs=(
            _out("unified_csv", "{prefix}_unified.csv", kind="csv", required=False),
            _out("overlap_report", "{prefix}_overlap_report.txt", kind="txt", required=False),
            _out("overlap_stats", "{prefix}_overlap_stats.json", kind="json", required=False),
        ),
    ),
    "ecc_summary": StepContract(
        step_name="ecc_summary",
        outputs=(
            _out("summary_dir", "summary_output", kind="dir", required=False),
        ),
    ),
    "ecc_packager": StepContract(step_name="ecc_packager"),
}


def get_step_contract(step_name: str) -> Optional[StepContract]:
    return STEP_CONTRACTS.get(step_name)
