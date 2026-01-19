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
    skip_tools: set[str] = set()
    if getattr(pipeline.config, "skip_tidehunter", False):
        skip_tools.add("tidehunter")

    all_ok = checker.check_all(skip_tools=skip_tools or None)

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

    pipeline._set_result(ResultKeys.TIDEHUNTER_OUTPUT, str(output_file))


def tandem_to_ring(pipeline: Pipeline) -> None:
    """Step 2: Process TideHunter output using tandem_to_ring module."""
    from circleseeker.modules.tandem_to_ring import TandemToRing

    tidehunter_output = pipeline.state.results.get(ResultKeys.TIDEHUNTER_OUTPUT)
    if not tidehunter_output:
        raise PipelineError("TideHunter output not found in state")

    csv_output = pipeline.config.output_dir / "tandem_to_ring.csv"
    fasta_output = pipeline.config.output_dir / "tandem_to_ring.fasta"

    t2r_cfg = getattr(pipeline.config.tools, "tandem_to_ring", {}) or {}
    if not isinstance(t2r_cfg, dict) and not hasattr(t2r_cfg, 'get'):
        raise PipelineError(
            "Invalid tandem_to_ring config; expected mapping, "
            f"got {type(t2r_cfg).__name__}"
        )
    min_ave_match = float(t2r_cfg.get("min_ave_match", 99.0))

    processor = TandemToRing(
        input_file=tidehunter_output,
        output_file=csv_output,
        circular_fasta=fasta_output,
        min_ave_match=min_ave_match,
        logger=pipeline.logger.getChild("tandem_to_ring"),
    )

    pipeline.logger.info("Running tandem_to_ring module")
    processor.run()

    pipeline._set_result(ResultKeys.T2R_CSV, str(csv_output))
    pipeline._set_result(ResultKeys.T2R_FASTA, str(fasta_output))


def run_alignment(pipeline: Pipeline) -> None:
    """Step 3: Run minimap2 alignment."""
    from circleseeker.external.minimap2_align import Minimap2Aligner

    alignment_output = (
        pipeline.config.output_dir / f"{pipeline.config.prefix}_alignment_results.tsv"
    )

    default_query = pipeline.config.input_file
    if default_query is None:
        raise PipelineError("Input file is required")

    skip_t2r = bool(getattr(pipeline.config, "skip_tandem_to_ring", False))
    t2r_fasta = pipeline._resolve_stored_path(
        pipeline.state.results.get(ResultKeys.T2R_FASTA),
        ["tandem_to_ring.fasta", f"{pipeline.config.prefix}_circular.fasta"],
    )
    if t2r_fasta is not None:
        query_file = t2r_fasta
    elif skip_t2r:
        query_file = Path(default_query)
        pipeline.logger.warning(
            "tandem_to_ring output not found; using input reads for alignment. "
            "Ensure query IDs match 'read|repN|length|copy' for um_classify."
        )
    else:
        raise PipelineError(
            "tandem_to_ring.fasta not found. Run tandem_to_ring or provide a "
            "precomputed output and set skip_carousel: true."
        )

    reference = pipeline.config.reference
    if reference is None:
        raise PipelineError("Reference genome is required")

    alignment_cfg = getattr(pipeline.config.tools, "alignment", {}) or {}
    if not isinstance(alignment_cfg, dict) and not hasattr(alignment_cfg, 'get'):
        raise PipelineError(
            "Invalid alignment config; expected mapping, "
            f"got {type(alignment_cfg).__name__}"
        )
    aligner_name = str(alignment_cfg.get("aligner", "minimap2")).lower()
    min_alignment_length = int(alignment_cfg.get("min_alignment_length", 0) or 0)

    minimap_cfg = pipeline.config.tools.minimap2_align or {}
    if not isinstance(minimap_cfg, dict) and not hasattr(minimap_cfg, 'get'):
        raise PipelineError(
            "Invalid minimap2_align config; expected mapping, "
            f"got {type(minimap_cfg).__name__}"
        )
    preset = minimap_cfg.get("preset", "sr")
    additional_args = minimap_cfg.get("additional_args", "")
    max_target_seqs = minimap_cfg.get("max_target_seqs", 200)
    # Identity filter with length-based compensation for HiFi data:
    # - Base threshold: 99.0% (HiFi error rate is ~1%)
    # - Decay: 0.5% per 10kb (longer sequences accumulate more errors)
    # - Floor: 97.0% (never go below this)
    min_identity = float(minimap_cfg.get("min_identity", 99.0))
    if alignment_cfg.get("min_identity") is not None:
        min_identity = float(alignment_cfg.get("min_identity"))
    identity_decay_per_10kb = float(minimap_cfg.get("identity_decay_per_10kb", 0.5))
    min_identity_floor = float(minimap_cfg.get("min_identity_floor", 97.0))
    # Length-based preset splitting: use different presets for short/long sequences
    split_by_length = bool(minimap_cfg.get("split_by_length", False))
    split_length = int(minimap_cfg.get("split_length", 5000))
    preset_short = minimap_cfg.get("preset_short", "sr")
    preset_long = minimap_cfg.get("preset_long", "sr")

    if aligner_name == "last":
        pipeline.logger.warning(
            "LAST aligner is not yet implemented; falling back to minimap2. "
            "Set alignment.aligner='minimap2' to silence this warning."
        )
        aligner_name = "minimap2"

    if aligner_name not in ("minimap2", ""):
        pipeline.logger.warning("Unknown aligner '%s', falling back to minimap2", aligner_name)

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

    if min_alignment_length > 0:
        tmp_output = alignment_output.with_suffix(".tmp.tsv")
        with open(alignment_output, "r") as fin, open(tmp_output, "w") as fout:
            for line in fin:
                if not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4:
                    continue
                try:
                    aln_len = int(float(parts[3]))
                except ValueError:
                    continue
                if aln_len >= min_alignment_length:
                    fout.write(line)
        tmp_output.replace(alignment_output)
    pipeline._set_result(ResultKeys.ALIGNMENT_OUTPUT, str(alignment_output))
