"""Step executors for early preprocessing stages."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from circleseeker.exceptions import PipelineError
from circleseeker.core.pipeline_types import ResultKeys

if TYPE_CHECKING:
    from circleseeker.core.pipeline import Pipeline


def check_dependencies(pipeline: Pipeline) -> None:
    """Step 0: Check all required tools and dependencies."""
    from circleseeker.utils.dependency_checker import DependencyChecker

    checker = DependencyChecker(logger=pipeline.logger.getChild("dependency_checker"))
    all_ok = checker.check_all()

    if not all_ok:
        checker.print_report()
        raise PipelineError(
            "Missing required dependencies. Please install missing tools and try again."
        )

    pipeline.logger.info("All required dependencies are available")


def tidehunter(pipeline: Pipeline) -> None:
    """Step 1: Run TideHunter."""
    from circleseeker.external.tidehunter import TideHunter

    input_file = pipeline.config.input_file
    if input_file is None:
        raise PipelineError("Input file is required")

    output_file = pipeline.config.output_dir / f"{pipeline.config.prefix}.TH.ecc_candidates.txt"

    tidehunter_tool = TideHunter(threads=pipeline.config.threads)
    tidehunter_tool.run_analysis(
        input_file=input_file,
        output_file=output_file,
        **pipeline.config.tools.tidehunter,
    )

    pipeline.state.results[ResultKeys.TIDEHUNTER_OUTPUT] = str(output_file)


def tandem_to_ring(pipeline: Pipeline) -> None:
    """Step 2: Process TideHunter output using tandem_to_ring module."""
    from circleseeker.modules.tandem_to_ring import TandemToRing

    tidehunter_output = pipeline.state.results.get(ResultKeys.TIDEHUNTER_OUTPUT)
    if not tidehunter_output:
        raise PipelineError("TideHunter output not found in state")

    csv_output = pipeline.config.output_dir / "tandem_to_ring.csv"
    fasta_output = pipeline.config.output_dir / "tandem_to_ring.fasta"

    processor = TandemToRing(
        input_file=tidehunter_output,
        output_file=csv_output,
        circular_fasta=fasta_output,
        logger=pipeline.logger.getChild("tandem_to_ring"),
    )

    pipeline.logger.info("Running tandem_to_ring module")
    processor.run()

    pipeline._set_result(ResultKeys.T2R_CSV, str(csv_output))
    pipeline._set_result(ResultKeys.T2R_FASTA, str(fasta_output))


def run_alignment(pipeline: Pipeline) -> None:
    """Step 3: Run alignment (minimap2 or LAST).

    The aligner can be configured via tools.alignment.aligner:
    - "minimap2" (default): Use minimap2 for alignment
    - "last": Use LAST for alignment (better for complex eccDNA)
    """
    alignment_output = (
        pipeline.config.output_dir / f"{pipeline.config.prefix}_alignment_results.tsv"
    )

    default_query = pipeline.config.input_file
    if default_query is None:
        raise PipelineError("Input file is required")

    query_file = Path(
        pipeline._get_result(ResultKeys.T2R_FASTA, default=default_query)
    )

    reference = pipeline.config.reference
    if reference is None:
        raise PipelineError("Reference genome is required")

    # Get alignment config (supports both old minimap2_align and new alignment format)
    align_cfg = pipeline.config.tools.alignment or pipeline.config.tools.minimap2_align or {}
    if not isinstance(align_cfg, dict):
        align_cfg = {}

    # Determine which aligner to use
    aligner = align_cfg.get("aligner", "minimap2")

    if aligner == "last":
        _run_last_alignment(pipeline, query_file, reference, alignment_output, align_cfg)
    else:
        _run_minimap2_alignment(pipeline, query_file, reference, alignment_output, align_cfg)

    pipeline.state.results[ResultKeys.ALIGNMENT_OUTPUT] = str(alignment_output)


def _run_minimap2_alignment(
    pipeline: Pipeline,
    query_file: Path,
    reference: Path,
    alignment_output: Path,
    cfg: dict,
) -> None:
    """Run minimap2 alignment."""
    from circleseeker.external.minimap2_align import Minimap2Aligner

    preset = cfg.get("preset", "sr")
    additional_args = cfg.get("additional_args", "")
    max_target_seqs = cfg.get("max_target_seqs", 200)
    # Identity filter with length-based compensation for HiFi data
    min_identity = float(cfg.get("min_identity", 99.0))
    identity_decay_per_10kb = float(cfg.get("identity_decay_per_10kb", 0.5))
    min_identity_floor = float(cfg.get("min_identity_floor", 97.0))
    # Length-based preset splitting
    split_by_length = bool(cfg.get("split_by_length", False))
    split_length = int(cfg.get("split_length", 5000))
    preset_short = cfg.get("preset_short", "sr")
    preset_long = cfg.get("preset_long", "sr")

    runner = Minimap2Aligner(
        threads=pipeline.config.threads, logger=pipeline.logger.getChild("minimap2_align")
    )
    runner.run_alignment(
        reference=reference,
        query=query_file,
        output_tsv=alignment_output,
        preset=preset,
        max_target_seqs=int(max_target_seqs),
        additional_args=additional_args,
        min_identity=min_identity,
        identity_decay_per_10kb=identity_decay_per_10kb,
        min_identity_floor=min_identity_floor,
        split_by_length=split_by_length,
        split_length=split_length,
        preset_short=preset_short,
        preset_long=preset_long,
    )


def _run_last_alignment(
    pipeline: Pipeline,
    query_file: Path,
    reference: Path,
    alignment_output: Path,
    cfg: dict,
) -> None:
    """Run LAST alignment."""
    from circleseeker.external.last import LastAligner

    min_identity = float(cfg.get("min_identity", 90.0))
    min_alignment_length = int(cfg.get("min_alignment_length", 50))

    # LAST database prefix (can be pre-built or auto-generated)
    db_prefix = cfg.get("db_prefix")
    if db_prefix:
        db_prefix = Path(db_prefix)

    runner = LastAligner(
        threads=pipeline.config.threads,
        logger=pipeline.logger.getChild("last"),
    )
    runner.run_alignment(
        query_fasta=query_file,
        reference_fasta=reference,
        output_tsv=alignment_output,
        db_prefix=db_prefix,
        min_identity=min_identity,
        min_alignment_length=min_alignment_length,
    )
