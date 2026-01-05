"""Step executors for early preprocessing stages."""

from __future__ import annotations

from pathlib import Path

from circleseeker.exceptions import PipelineError
from circleseeker.core.pipeline_types import ResultKeys


def check_dependencies(pipeline) -> None:
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


def tidehunter(pipeline) -> None:
    """Step 1: Run TideHunter."""
    from circleseeker.external.tidehunter import TideHunter

    output_file = pipeline.config.output_dir / f"{pipeline.config.prefix}.TH.ecc_candidates.txt"

    tidehunter_tool = TideHunter(threads=pipeline.config.threads)
    tidehunter_tool.run_analysis(
        input_file=pipeline.config.input_file,
        output_file=output_file,
        **pipeline.config.tools.tidehunter,
    )

    pipeline.state.results[ResultKeys.TIDEHUNTER_OUTPUT] = str(output_file)


def tandem_to_ring(pipeline) -> None:
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


def run_alignment(pipeline) -> None:
    """Step 3: Run minimap2 alignment."""
    from circleseeker.external.minimap2_align import Minimap2Aligner

    alignment_output = (
        pipeline.config.output_dir / f"{pipeline.config.prefix}_alignment_results.tsv"
    )

    query_file = Path(
        pipeline._get_result(ResultKeys.T2R_FASTA, default=pipeline.config.input_file)
    )

    minimap_cfg = pipeline.config.tools.minimap2_align or {}
    if not isinstance(minimap_cfg, dict):
        raise PipelineError(
            "Invalid minimap2_align config; expected mapping, "
            f"got {type(minimap_cfg).__name__}"
        )
    preset = minimap_cfg.get("preset", "sr")
    additional_args = minimap_cfg.get("additional_args", "")
    max_target_seqs = minimap_cfg.get("max_target_seqs", 200)

    runner = Minimap2Aligner(
        threads=pipeline.config.threads, logger=pipeline.logger.getChild("minimap2_align")
    )
    runner.run(
        reference=pipeline.config.reference,
        query=query_file,
        output_tsv=alignment_output,
        preset=preset,
        max_target_seqs=int(max_target_seqs),
        additional_args=additional_args,
    )
    pipeline.state.results[ResultKeys.ALIGNMENT_OUTPUT] = str(alignment_output)

