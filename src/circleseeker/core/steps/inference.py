"""Step executors for inference and read filtering stages."""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import TYPE_CHECKING, Optional

import pandas as pd

from circleseeker.core.pipeline_types import ResultKeys
from circleseeker.exceptions import PipelineError

if TYPE_CHECKING:
    from circleseeker.core.pipeline import Pipeline
    from circleseeker.modules.splitreads_core import SplitReadsConfig


def _get_xecc_source_read_names(xecc_fasta: Path, logger: logging.Logger) -> set[str]:
    """
    Extract original read names from XeccDNA FASTA.

    XeccDNA FASTA IDs have format: {readName}|{repN}|{consLen}|{copyNum}|circular
    or with suffix: {readName}|...|circular__X{n}

    Returns set of original read names.
    """
    if not xecc_fasta.exists() or xecc_fasta.stat().st_size == 0:
        return set()

    read_names: set[str] = set()
    try:
        with open(xecc_fasta, "rb") as f:
            for line in f:
                if line.startswith(b">"):
                    # Parse ID: >readName|repN|consLen|copyNum|circular or similar
                    header_parts = line[1:].split()
                    if not header_parts:
                        continue  # Skip malformed headers
                    seq_id = header_parts[0].decode("utf-8", errors="replace")
                    # Original read name is the first part before |
                    original_read_name = seq_id.split("|")[0]
                    if original_read_name:
                        read_names.add(original_read_name)
    except OSError as e:
        logger.warning(f"Failed to parse XeccDNA FASTA (I/O error): {e}")
        return set()

    logger.info(f"Found {len(read_names)} unique source reads from XeccDNA")
    return read_names


def _append_reads_from_fasta(
    source_fasta: Path,
    target_fasta: Path,
    read_names: set[str],
    logger: logging.Logger,
) -> int:
    """
    Append specific reads from source FASTA to target FASTA.

    Args:
        source_fasta: Original input FASTA file
        target_fasta: Target file to append reads to
        read_names: Set of read names to extract
        logger: Logger instance

    Returns:
        Number of reads appended
    """
    if not source_fasta.exists() or not read_names:
        return 0

    appended_count = 0
    try:
        if target_fasta.exists() and target_fasta.stat().st_size > 0:
            with open(target_fasta, "rb+") as out_f:
                out_f.seek(-1, os.SEEK_END)
                if out_f.read(1) != b"\n":
                    out_f.write(b"\n")
        with open(target_fasta, "ab") as out_f:
            with open(source_fasta, "rb") as in_f:
                current_header: Optional[bytes] = None
                current_lines: list[bytes] = []
                write_current = False

                for line in in_f:
                    if line.startswith(b">"):
                        # Write previous sequence if needed
                        if current_header is not None and write_current:
                            out_f.write(current_header)
                            for seq_line in current_lines:
                                out_f.write(seq_line)
                            appended_count += 1

                        # Start new sequence
                        current_header = line
                        current_lines = []

                        # Check if this read should be included
                        header_parts = line[1:].decode("utf-8", errors="replace").split()
                        if not header_parts:
                            write_current = False
                            continue
                        header_str = header_parts[0]
                        write_current = header_str in read_names
                    else:
                        current_lines.append(line)

                # Write last sequence
                if current_header is not None and write_current:
                    out_f.write(current_header)
                    for seq_line in current_lines:
                        out_f.write(seq_line)
                    appended_count += 1

    except OSError as e:
        logger.warning(f"Failed to append XeccDNA source reads (I/O error): {e}")
        return 0

    if appended_count > 0:
        logger.info(f"Appended {appended_count} XeccDNA source reads to inference input")
    return appended_count


def read_filter(pipeline: Pipeline) -> None:
    """Step 9: Filter sequences using read_filter module (Sieve)."""
    from circleseeker.modules.read_filter import Sieve

    tandem_to_ring_csv = pipeline.config.output_dir / "tandem_to_ring.csv"

    original_fasta = pipeline.config.input_file
    if original_fasta is None:
        raise PipelineError("Input file is required")
    if not original_fasta.exists():
        pipeline.logger.error("Original input FASTA file not found; cannot run sieve step")
        raise PipelineError(f"Missing input FASTA: {original_fasta}")

    output_fasta = pipeline.config.output_dir / f"{pipeline.config.prefix}_all_filtered.fasta"

    sieve = Sieve(pipeline.logger.getChild("read_filter"))

    if not tandem_to_ring_csv.exists():
        pipeline.logger.warning(
            "TandemToRing CSV not found; skipping CtcR filtering and retaining all reads"
        )
        stats = sieve.filter_fasta_files(
            input_fastas=[original_fasta],
            output_fasta=output_fasta,
        )
    else:
        stats = sieve.run_sieve(
            tandem_to_ring_csv=tandem_to_ring_csv,
            input_fastas=[original_fasta],
            output_fasta=output_fasta,
            ctcr_classes=None,
        )

    pipeline._set_result(ResultKeys.READ_FILTER_OUTPUT_FASTA, str(output_fasta))
    pipeline._set_result(ResultKeys.READ_FILTER_TOTAL, stats.total_reads)
    pipeline._set_result(ResultKeys.READ_FILTER_FILTERED, stats.filtered_reads)
    pipeline._set_result(ResultKeys.READ_FILTER_RETAINED, stats.retained_reads)
    pipeline._set_result(ResultKeys.READ_FILTER_CTCR, stats.csv_ctcr_reads)

    pipeline.logger.info(
        f"Sieve complete: {stats.filtered_reads} CtcR reads filtered, "
        f"{stats.retained_reads} retained ({stats.filtered_percentage:.2f}% filtered)"
    )


def minimap2(pipeline: Pipeline, *, force_alignment: bool = False) -> None:
    """Step 10: Prepare mapping artifacts (alignment optional for SplitReads-Core)."""
    from circleseeker.external.minimap2 import Minimap2, MiniMapConfig

    sieve_output = pipeline._get_result(ResultKeys.READ_FILTER_OUTPUT_FASTA)
    if sieve_output and Path(sieve_output).exists() and Path(sieve_output).stat().st_size > 0:
        input_fasta = Path(sieve_output)
    else:
        input_fasta = pipeline.config.output_dir / f"{pipeline.config.prefix}_all_filtered.fasta"
        try:
            seen_fastas: set[Path] = set()
            with open(input_fasta, "wb") as out_f:
                for ecc_type in ["uecc", "mecc", "cecc"]:
                    candidate_keys = [
                        f"ecc_dedup_{ecc_type}_fasta",
                        f"{ecc_type}_fasta",
                        f"harmonizer_{ecc_type}_fasta",
                    ]
                    for key in candidate_keys:
                        fasta_path_str = pipeline.state.results.get(key)
                        if not fasta_path_str:
                            continue
                        fasta_path = Path(fasta_path_str)
                        if fasta_path in seen_fastas:
                            continue
                        if fasta_path.exists() and fasta_path.stat().st_size > 0:
                            last_byte: Optional[int] = None
                            with open(fasta_path, "rb") as in_f:
                                while True:
                                    chunk = in_f.read(1024 * 1024)
                                    if not chunk:
                                        break
                                    out_f.write(chunk)
                                    last_byte = chunk[-1]
                            if last_byte is not None and last_byte != 10:
                                out_f.write(b"\n")
                            seen_fastas.add(fasta_path)
        except OSError as e:
            pipeline.logger.error(f"Failed to create combined FASTA: {e}")
            raise

    # Append XeccDNA source reads to inference input
    # XeccDNA are rings that couldn't be classified as U/M/C, their source reads
    # should also be included in inference to see if they can be detected
    xecc_fasta = pipeline.config.output_dir / f"{pipeline.config.prefix}_XeccDNA.fasta"
    xecc_appended = 0
    if xecc_fasta.exists() and xecc_fasta.stat().st_size > 0:
        xecc_read_names = _get_xecc_source_read_names(xecc_fasta, pipeline.logger)
        if xecc_read_names:
            original_input = pipeline.config.input_file
            if original_input and Path(original_input).exists():
                xecc_appended = _append_reads_from_fasta(
                    source_fasta=Path(original_input),
                    target_fasta=input_fasta,
                    read_names=xecc_read_names,
                    logger=pipeline.logger,
                )
    else:
        pipeline.logger.debug("No XeccDNA FASTA found or empty; skipping XeccDNA source reads")

    pipeline._set_result(ResultKeys.ALL_FILTERED_FASTA, str(input_fasta))
    pipeline._set_result(ResultKeys.XECC_SOURCE_READS_ADDED, xecc_appended)

    if not input_fasta.exists() or input_fasta.stat().st_size == 0:
        pipeline.logger.warning(
            "Inference input FASTA is empty; skipping minimap2 alignment and inference steps"
        )
        pipeline._set_result(ResultKeys.INFERENCE_INPUT_EMPTY, True)
        return

    reference_mmi = pipeline._ensure_reference_mmi()
    pipeline._set_result(ResultKeys.REFERENCE_MMI, str(reference_mmi))

    # SplitReads-Core is built-in and doesn't need BAM alignment
    # It uses mappy directly for alignment
    if not force_alignment:
        pipeline.logger.info("SplitReads-Core selected; skipping minimap2 alignment (BAM not required)")
        return

    # Read minimap2 config from pipeline configuration
    minimap2_cfg = pipeline.config.tools.minimap2
    config = MiniMapConfig(
        preset=minimap2_cfg.get("preset", "map-hifi"),
        threads=pipeline.config.threads,
        output_format="bam",
        allow_secondary=True,
        build_index=True,
        sort_bam=True,
        index_bam=True,
        additional_args=minimap2_cfg.get("additional_args", ""),
    )

    minimap2_tool = Minimap2(config=config, logger=pipeline.logger.getChild("minimap2"))
    output_bam = pipeline.config.output_dir / f"{pipeline.config.prefix}_sorted.bam"

    minimap2_tool.align(
        reference=reference_mmi,
        query=input_fasta,
        output_file=output_bam,
        build_index_first=False,
    )

    pipeline._set_result(ResultKeys.MINIMAP2_BAM, str(output_bam))


def ecc_inference(pipeline: Pipeline) -> None:
    """Step 11: Circular DNA detection using SplitReads-Core (built-in, HiFi optimized)."""
    pipeline._del_result(ResultKeys.INFERENCE_FAILED)
    pipeline._del_result(ResultKeys.INFERENCE_ERROR)
    pipeline._del_result(ResultKeys.INFERENCE_ERROR_FILE)

    if pipeline._get_result(ResultKeys.INFERENCE_INPUT_EMPTY):
        pipeline.logger.warning("Inference input FASTA empty; skipping eccDNA inference")
        pipeline._set_result(ResultKeys.INFERENCE_RESULT_COUNT, 0)
        pipeline._set_result(ResultKeys.INFERENCE_OVERVIEW, None)
        return

    all_filtered = pipeline._get_result(ResultKeys.ALL_FILTERED_FASTA)
    if all_filtered:
        all_filtered_path = Path(all_filtered)
        if all_filtered_path.exists() and all_filtered_path.stat().st_size == 0:
            pipeline.logger.warning("Inference input FASTA empty; skipping eccDNA inference")
            pipeline._set_result(ResultKeys.INFERENCE_RESULT_COUNT, 0)
            pipeline._set_result(ResultKeys.INFERENCE_OVERVIEW, None)
            return

    def _mark_inference_failed(exc: Exception) -> None:
        import traceback

        # Inference is a best-effort step; keep the pipeline running but make debugging easy.
        pipeline.logger.exception(
            "eccDNA inference failed; continuing without inferred eccDNA",
        )
        error_text = f"{type(exc).__name__}: {exc}"
        pipeline._set_result(ResultKeys.INFERENCE_FAILED, True)
        pipeline._set_result(ResultKeys.INFERENCE_ERROR, error_text)
        pipeline._set_result(ResultKeys.INFERENCE_RESULT_COUNT, 0)
        pipeline._set_result(ResultKeys.INFERENCE_OVERVIEW, None)
        pipeline._set_result(ResultKeys.INFERRED_SIMPLE_TSV, None)
        pipeline._set_result(ResultKeys.INFERRED_CHIMERIC_TSV, None)
        pipeline._set_result(ResultKeys.INFERRED_SIMPLE_FASTA, None)
        pipeline._set_result(ResultKeys.INFERRED_CHIMERIC_FASTA, None)
        pipeline._set_result(ResultKeys.INFERRED_DIR, None)

        try:
            pipeline.config.output_dir.mkdir(parents=True, exist_ok=True)
            marker = pipeline.config.output_dir / f"{pipeline.config.prefix}_inference_failed.txt"
            marker.write_text(
                "\n".join(
                    [
                        "eccDNA inference failed; CircleSeeker continued without inferred eccDNA.",
                        "step: ecc_inference",
                        "tool: SplitReads-Core",
                        f"error: {error_text}",
                        "",
                        "traceback:",
                        traceback.format_exc().rstrip(),
                        "",
                        "Confirmed eccDNA outputs generated before inference should still be present.",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            pipeline._set_result(
                ResultKeys.INFERENCE_ERROR_FILE,
                pipeline._serialize_path_for_state(marker),
            )
        except OSError as write_exc:  # pragma: no cover
            pipeline.logger.debug("Failed to write inference failure marker: %s", write_exc)

    try:
        pipeline.logger.info("Using SplitReads-Core for circular DNA detection (built-in, HiFi optimized)")
        pipeline._run_splitreads_core_inference()
    except (PipelineError, OSError, ValueError, RuntimeError) as exc:
        _mark_inference_failed(exc=exc)


# Pipeline config → module config field name mapping.
# The pipeline-level SplitReadsConfig (config.py) uses descriptive names
# while the module-level SplitReadsConfig (splitreads_core.py) uses short names.
_SPLITREADS_FIELD_MAP: dict[str, str] = {
    "mapq_threshold": "mapq",
    "gap_tolerance": "allow_gap",
    "overlap_tolerance": "allow_overlap",
    "min_breakpoint_depth": "breakpoint_depth",
    "min_avg_depth": "average_depth",
}


def _build_splitreads_config(raw_cfg: object, config_cls: type[SplitReadsConfig]) -> SplitReadsConfig:
    """Create module SplitReadsConfig from pipeline config, bridging field name differences."""
    if raw_cfg is None:
        return config_cls()

    if isinstance(raw_cfg, config_cls):
        return raw_cfg

    # Extract raw dict from whatever object we have
    if isinstance(raw_cfg, dict):
        raw_dict = dict(raw_cfg)
    else:
        raw_dict = {
            k: getattr(raw_cfg, k)
            for k in dir(raw_cfg)
            if not k.startswith("_") and not callable(getattr(raw_cfg, k, None))
        }

    # Apply field name mapping (pipeline names → module names)
    mapped: dict[str, object] = {}
    for key, value in raw_dict.items():
        mapped_key = _SPLITREADS_FIELD_MAP.get(key, key)
        mapped[mapped_key] = value

    try:
        return config_cls.from_dict(mapped)
    except (ValueError, TypeError):
        return config_cls()


def run_splitreads_core_inference(pipeline: Pipeline) -> None:
    """Run SplitReads-Core inference pipeline (built-in, HiFi optimized)."""
    from circleseeker.modules.splitreads_core import SplitReadsCore, SplitReadsConfig
    from circleseeker.modules.splitreads_adapter import write_overview_tsv

    sieve_output = pipeline._get_result(ResultKeys.READ_FILTER_OUTPUT_FASTA)

    if sieve_output and Path(sieve_output).exists():
        input_fasta = Path(sieve_output)
    else:
        input_fasta = pipeline.config.output_dir / f"{pipeline.config.prefix}_all_filtered.fasta"
        if not input_fasta.exists() or input_fasta.stat().st_size == 0:
            pipeline.logger.warning(
                "No filtered FASTA found; skipping SplitReads-Core inference. "
                "Inference requires filtered reads from earlier pipeline steps."
            )
            pipeline._set_result(ResultKeys.INFERENCE_RESULT_COUNT, 0)
            pipeline._set_result(ResultKeys.INFERENCE_OVERVIEW, None)
            return

    reference_mmi_str = pipeline._get_result(ResultKeys.REFERENCE_MMI)
    reference_mmi = Path(reference_mmi_str) if reference_mmi_str else None
    if not (reference_mmi and reference_mmi.exists()):
        reference_mmi = pipeline._ensure_reference_mmi()
        pipeline._set_result(ResultKeys.REFERENCE_MMI, str(reference_mmi))
    if reference_mmi is None:
        raise PipelineError("Failed to create or locate reference MMI index")

    reference_fasta = pipeline.config.reference
    if reference_fasta is None:
        raise PipelineError("Reference genome is required")

    # Ensure .fai index exists for reference
    reference_fai = Path(str(reference_fasta) + ".fai")
    if not reference_fai.exists():
        from circleseeker.external.samtools import Samtools
        samtools = Samtools(logger=pipeline.logger.getChild("samtools"))
        samtools.faidx(reference_fasta)

    output_dir = pipeline.config.output_dir / "splitreads_core_output"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Get configuration from pipeline config
    splitreads_cfg = getattr(pipeline.config.tools, "splitreads", None)
    splitreads_config = _build_splitreads_config(splitreads_cfg, SplitReadsConfig)

    # Initialize SplitReadsCore (built-in, HiFi optimized)
    splitreads_core = SplitReadsCore(
        reference_mmi=reference_mmi,
        config=splitreads_config,
        logger=pipeline.logger.getChild("splitreads_core"),
        threads=pipeline.config.threads,
    )

    pipeline.logger.info("Running SplitReads-Core pipeline (built-in, HiFi optimized)")
    eccDNA_final = splitreads_core.run(
        input_fasta=input_fasta,
        chrom_sizes_file=reference_fai,
        output_dir=output_dir,
    )

    if not eccDNA_final.exists():
        raise FileNotFoundError(
            f"SplitReads-Core did not produce expected output: {eccDNA_final}"
        )

    pipeline.logger.info("Converting SplitReads-Core output to internal format")
    overview_tsv = pipeline.config.output_dir / f"{pipeline.config.prefix}_overview.tsv"
    write_overview_tsv(eccDNA_final, overview_tsv)

    if overview_tsv.exists():
        pipeline._set_result(
            ResultKeys.INFERENCE_OVERVIEW,
            pipeline._serialize_path_for_state(overview_tsv),
        )
        df = pd.read_csv(overview_tsv, sep="\t")
        pipeline._set_result(ResultKeys.INFERENCE_RESULT_COUNT, len(df))
        pipeline.logger.info(f"SplitReads-Core detected {len(df)} circular DNA candidates")
    else:
        pipeline.logger.warning("SplitReads-Core output conversion failed")
        pipeline._set_result(ResultKeys.INFERENCE_RESULT_COUNT, 0)


def curate_inferred_ecc(pipeline: Pipeline) -> list[str] | None:
    """Step 12: Curate inferred eccDNA tables from SplitReads-Core overview TSV."""
    from circleseeker.modules.iecc_curator import (
        curate_ecc_tables,
        write_curated_tables,
        write_curated_tables_with_fasta,
    )

    if pipeline._get_result(ResultKeys.INFERENCE_INPUT_EMPTY) or pipeline._get_result(
        ResultKeys.INFERENCE_FAILED
    ):
        pipeline.logger.warning("Inference was skipped/failed; skipping inferred eccDNA curation step")
        pipeline._set_result(ResultKeys.INFERRED_SIMPLE_TSV, None)
        pipeline._set_result(ResultKeys.INFERRED_CHIMERIC_TSV, None)
        pipeline._set_result(ResultKeys.INFERRED_SIMPLE_FASTA, None)
        pipeline._set_result(ResultKeys.INFERRED_CHIMERIC_FASTA, None)
        pipeline._set_result(ResultKeys.INFERRED_DIR, None)
        return None

    overview_path_str = pipeline._get_result(ResultKeys.INFERENCE_OVERVIEW)

    overview_path_obj = pipeline._resolve_stored_path(
        overview_path_str,
        [f"{pipeline.config.prefix}_overview.tsv"],
    )
    if not overview_path_obj:
        pipeline.logger.warning("SplitReads-Core overview table missing; skipping inferred eccDNA curation")
        pipeline._set_result(ResultKeys.INFERRED_SIMPLE_TSV, None)
        pipeline._set_result(ResultKeys.INFERRED_CHIMERIC_TSV, None)
        pipeline._set_result(ResultKeys.INFERRED_SIMPLE_FASTA, None)
        pipeline._set_result(ResultKeys.INFERRED_CHIMERIC_FASTA, None)
        pipeline._set_result(ResultKeys.INFERRED_DIR, None)
        return None

    pipeline._set_result(
        ResultKeys.INFERENCE_OVERVIEW,
        pipeline._serialize_path_for_state(overview_path_obj),
    )

    simple_df, chimeric_df = curate_ecc_tables(overview_path_obj)

    # Optional (opt-in) filter to suppress false positives from inferred
    # CeccDNA.  Empirically, 2-segment chimeric circles are the dominant FP
    # mode for split-read graph inference.  The filter activates only when
    # the user explicitly raises min_inferred_chimeric_segments above 2 or
    # sets min_inferred_two_segment_split_reads > 0 in the config.
    try:
        splitreads_cfg = pipeline.config.tools.splitreads
        min_seg = int(getattr(splitreads_cfg, "min_inferred_chimeric_segments", 2))
        min_two_seg_splits = int(getattr(splitreads_cfg, "min_inferred_two_segment_split_reads", 0))
    except (AttributeError, ValueError, TypeError):
        min_seg = 2
        min_two_seg_splits = 0

    if (
        not chimeric_df.empty
        and {"eccDNA_id", "seg_total", "num_split_reads"}.issubset(chimeric_df.columns)
        and (min_seg > 2 or min_two_seg_splits > 0)
    ):
        before = int(chimeric_df["eccDNA_id"].nunique())
        keep_mask = chimeric_df["seg_total"] >= min_seg
        if min_two_seg_splits > 0:
            keep_mask |= (chimeric_df["seg_total"] == 2) & (
                chimeric_df["num_split_reads"] >= min_two_seg_splits
            )
        chimeric_df = chimeric_df.loc[keep_mask].copy()
        after = int(chimeric_df["eccDNA_id"].nunique())
        if after != before:
            pipeline.logger.info(
                "Filtered inferred chimeric eccDNA: %d -> %d circles "
                "(min_segments=%d, min_two_segment_split_reads=%d)",
                before,
                after,
                min_seg,
                min_two_seg_splits,
            )

    inferred_prefix = pipeline.config.output_dir / pipeline.config.prefix
    reference_for_inferred: Optional[Path] = None
    if pipeline.config.reference:
        ref_path = Path(pipeline.config.reference)
        if ref_path.exists():
            reference_for_inferred = ref_path
        else:
            pipeline.logger.warning(
                f"Reference FASTA not found for inferred eccDNA FASTA generation: {ref_path}"
            )

    simple_csv_path: Optional[Path] = None
    chimeric_csv_path: Optional[Path] = None
    simple_fasta_path: Optional[Path] = None
    chimeric_fasta_path: Optional[Path] = None

    try:
        simple_csv_path, chimeric_csv_path, simple_fasta_path, chimeric_fasta_path = (
            write_curated_tables_with_fasta(
                simple_df,
                chimeric_df,
                inferred_prefix,
                reference_fasta=reference_for_inferred,
                organize_files=True,
            )
        )
    except (OSError, ValueError, KeyError) as exc:
        pipeline.logger.warning(f"Failed to generate inferred eccDNA tables with FASTA: {exc}")
        simple_csv_path, chimeric_csv_path = write_curated_tables(
            simple_df,
            chimeric_df,
            inferred_prefix,
        )

    written: list[str] = []
    if simple_csv_path and simple_csv_path.exists() and simple_csv_path.stat().st_size > 0:
        pipeline._set_result(
            ResultKeys.INFERRED_SIMPLE_TSV,
            pipeline._serialize_path_for_state(simple_csv_path),
        )
        written.append(str(simple_csv_path))
    else:
        pipeline.logger.info("No inferred simple eccDNA records; simple TSV not created")

    if chimeric_csv_path and chimeric_csv_path.exists() and chimeric_csv_path.stat().st_size > 0:
        pipeline._set_result(
            ResultKeys.INFERRED_CHIMERIC_TSV,
            pipeline._serialize_path_for_state(chimeric_csv_path),
        )
        written.append(str(chimeric_csv_path))
    else:
        pipeline.logger.info("No inferred chimeric eccDNA records; chimeric TSV not created")

    if simple_fasta_path and simple_fasta_path.exists() and simple_fasta_path.stat().st_size > 0:
        pipeline._set_result(
            ResultKeys.INFERRED_SIMPLE_FASTA,
            pipeline._serialize_path_for_state(simple_fasta_path),
        )

    if chimeric_fasta_path and chimeric_fasta_path.exists() and chimeric_fasta_path.stat().st_size > 0:
        pipeline._set_result(
            ResultKeys.INFERRED_CHIMERIC_FASTA,
            pipeline._serialize_path_for_state(chimeric_fasta_path),
        )

    inferred_dir_path = None
    if simple_csv_path:
        inferred_dir_path = simple_csv_path.parent
    elif chimeric_csv_path:
        inferred_dir_path = chimeric_csv_path.parent
    if inferred_dir_path and inferred_dir_path.exists():
        pipeline._set_result(
            ResultKeys.INFERRED_DIR,
            pipeline._serialize_path_for_state(inferred_dir_path),
        )

    if not simple_df.empty and not (
        simple_fasta_path and simple_fasta_path.exists() and simple_fasta_path.stat().st_size > 0
    ):
        pipeline.logger.warning(
            "Inferred simple eccDNA FASTA was not generated; continuing with tables only. "
            "Install pysam and ensure reference access to enable FASTA output."
        )

    if not chimeric_df.empty and not (
        chimeric_fasta_path
        and chimeric_fasta_path.exists()
        and chimeric_fasta_path.stat().st_size > 0
    ):
        pipeline.logger.warning(
            "Inferred chimeric eccDNA FASTA was not generated; continuing with tables only. "
            "Install pysam and ensure reference access to enable FASTA output."
        )

    return written or None
