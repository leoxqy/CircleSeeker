"""Step executors for inference and read filtering stages."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Optional

from circleseeker.core.pipeline_types import ResultKeys
from circleseeker.exceptions import PipelineError

if TYPE_CHECKING:
    from circleseeker.core.pipeline import Pipeline


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
        with open(xecc_fasta, "r") as f:
            for line in f:
                if line.startswith(">"):
                    # Parse ID: >readName|repN|consLen|copyNum|circular or similar
                    seq_id = line[1:].split()[0]  # Remove > and take first part
                    # Original read name is the first part before |
                    original_read_name = seq_id.split("|")[0]
                    if original_read_name:
                        read_names.add(original_read_name)
    except Exception as e:
        logger.warning(f"Failed to parse XeccDNA FASTA: {e}")
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
                        header_str = line[1:].decode("utf-8", errors="replace").split()[0]
                        write_current = header_str in read_names
                    else:
                        current_lines.append(line)

                # Write last sequence
                if current_header is not None and write_current:
                    out_f.write(current_header)
                    for seq_line in current_lines:
                        out_f.write(seq_line)
                    appended_count += 1

    except Exception as e:
        logger.warning(f"Failed to append XeccDNA source reads: {e}")
        return 0

    if appended_count > 0:
        logger.info(f"Appended {appended_count} XeccDNA source reads to inference input")
    return appended_count


def read_filter(pipeline: Pipeline) -> None:
    """Step 9: Filter sequences using read_filter module (Sieve)."""
    from circleseeker.modules.read_filter import Sieve

    carousel_csv = pipeline.config.output_dir / "tandem_to_ring.csv"

    original_fasta = pipeline.config.input_file
    if original_fasta is None:
        raise PipelineError("Input file is required")
    if not original_fasta.exists():
        pipeline.logger.error("Original input FASTA file not found; cannot run sieve step")
        raise PipelineError(f"Missing input FASTA: {original_fasta}")

    output_fasta = pipeline.config.output_dir / f"{pipeline.config.prefix}_all_filtered.fasta"

    sieve = Sieve(pipeline.logger.getChild("read_filter"))

    if not carousel_csv.exists():
        pipeline.logger.warning(
            "TandemToRing CSV not found; skipping CtcR filtering and retaining all reads"
        )
        stats = sieve.filter_fasta_files(
            input_fastas=[original_fasta],
            output_fasta=output_fasta,
        )
    else:
        stats = sieve.run_sieve(
            carousel_csv=carousel_csv,
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


def minimap2(pipeline: Pipeline) -> None:
    """Step 10: Prepare mapping artifacts (alignment optional for Cresil)."""
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
        except Exception as e:
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
        pipeline.state.results[ResultKeys.INFERENCE_INPUT_EMPTY] = True
        return

    reference_mmi = pipeline._ensure_reference_mmi()
    pipeline._set_result(ResultKeys.REFERENCE_MMI, str(reference_mmi))

    inference_tool = pipeline._select_inference_tool()
    if inference_tool == "cresil":
        pipeline.logger.info("Cresil selected; skipping minimap2 alignment (BAM not required)")
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
    """Step 11: Circular DNA detection - Cresil preferred, Cyrcular fallback."""
    pipeline.state.results.pop(ResultKeys.INFERENCE_FAILED, None)
    pipeline.state.results.pop(ResultKeys.INFERENCE_ERROR, None)
    pipeline.state.results.pop(ResultKeys.INFERENCE_ERROR_FILE, None)

    if pipeline.state.results.get(ResultKeys.INFERENCE_INPUT_EMPTY):
        pipeline.logger.warning("Inference input FASTA empty; skipping eccDNA inference")
        pipeline.state.results[ResultKeys.CYRCULAR_RESULT_COUNT] = 0
        pipeline.state.results[ResultKeys.CYRCULAR_OVERVIEW] = None
        return

    all_filtered = pipeline.state.results.get(ResultKeys.ALL_FILTERED_FASTA)
    if all_filtered:
        all_filtered_path = Path(all_filtered)
        if all_filtered_path.exists() and all_filtered_path.stat().st_size == 0:
            pipeline.logger.warning("Inference input FASTA empty; skipping eccDNA inference")
            pipeline.state.results[ResultKeys.CYRCULAR_RESULT_COUNT] = 0
            pipeline.state.results[ResultKeys.CYRCULAR_OVERVIEW] = None
            return

    def _mark_inference_failed(tool_name: str, exc: Exception) -> None:
        pipeline.logger.warning(
            "eccDNA inference failed (%s); continuing without inferred eccDNA: %s",
            tool_name,
            exc,
        )
        pipeline.state.results[ResultKeys.INFERENCE_FAILED] = True
        pipeline.state.results[ResultKeys.INFERENCE_ERROR] = str(exc)
        pipeline.state.results[ResultKeys.CYRCULAR_RESULT_COUNT] = 0
        pipeline.state.results[ResultKeys.CYRCULAR_OVERVIEW] = None
        pipeline.state.results[ResultKeys.INFERRED_SIMPLE_TSV] = None
        pipeline.state.results[ResultKeys.INFERRED_CHIMERIC_TSV] = None
        pipeline.state.results[ResultKeys.INFERRED_SIMPLE_FASTA] = None
        pipeline.state.results[ResultKeys.INFERRED_CHIMERIC_FASTA] = None
        pipeline.state.results[ResultKeys.INFERRED_DIR] = None

        try:
            pipeline.config.output_dir.mkdir(parents=True, exist_ok=True)
            marker = pipeline.config.output_dir / f"{pipeline.config.prefix}_inference_failed.txt"
            marker.write_text(
                "\n".join(
                    [
                        "eccDNA inference failed; CircleSeeker continued without inferred eccDNA.",
                        "step: ecc_inference",
                        f"tool: {tool_name}",
                        f"error: {exc}",
                        "",
                        "Confirmed eccDNA outputs generated before inference should still be present.",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            pipeline.state.results[ResultKeys.INFERENCE_ERROR_FILE] = pipeline._serialize_path_for_state(
                marker
            )
        except Exception as write_exc:  # pragma: no cover
            pipeline.logger.debug("Failed to write inference failure marker: %s", write_exc)

    try:
        tool = pipeline._select_inference_tool()
        if tool == "cresil":
            pipeline.logger.info("Using Cresil for circular DNA detection (preferred)")
            pipeline._run_cresil_inference()
        else:
            pipeline.logger.info("Using Cyrcular for circular DNA detection (fallback)")
            pipeline._run_cyrcular_inference()
    except Exception as exc:
        _mark_inference_failed(tool_name=locals().get("tool", "unknown"), exc=exc)


def run_cresil_inference(pipeline: Pipeline) -> None:
    """Run Cresil inference pipeline (replaces Steps 10-11 when available)."""
    from circleseeker.external.cresil import Cresil
    from circleseeker.modules.cresil_adapter import write_cyrcular_compatible_tsv

    sieve_output = pipeline._get_result(ResultKeys.READ_FILTER_OUTPUT_FASTA)

    if sieve_output and Path(sieve_output).exists():
        input_fasta = Path(sieve_output)
    else:
        input_fasta = pipeline.config.output_dir / f"{pipeline.config.prefix}_all_filtered.fasta"
        if not input_fasta.exists() or input_fasta.stat().st_size == 0:
            pipeline.logger.warning(
                "No filtered FASTA found; skipping Cresil inference. "
                "Inference requires filtered reads from earlier pipeline steps."
            )
            pipeline.state.results[ResultKeys.CYRCULAR_RESULT_COUNT] = 0
            pipeline.state.results[ResultKeys.CYRCULAR_OVERVIEW] = None
            return

    reference_mmi_str = pipeline._get_result(ResultKeys.REFERENCE_MMI)
    reference_mmi = Path(reference_mmi_str) if reference_mmi_str else None
    if not (reference_mmi and reference_mmi.exists()):
        reference_mmi = pipeline._ensure_reference_mmi()
        pipeline._set_result(ResultKeys.REFERENCE_MMI, str(reference_mmi))
    if reference_mmi is None:
        raise PipelineError("Failed to create or locate reference MMI index")

    cresil = Cresil(logger=pipeline.logger.getChild("cresil"))

    output_dir = pipeline.config.output_dir / "cresil_output"
    output_dir.mkdir(parents=True, exist_ok=True)

    reference_fasta = pipeline.config.reference
    if reference_fasta is None:
        raise PipelineError("Reference genome is required")

    pipeline.logger.info("Running Cresil pipeline (trim + identify)")
    eccDNA_final = cresil.run_full_pipeline(
        fasta_query=input_fasta,
        reference_fasta=reference_fasta,
        reference_mmi=reference_mmi,
        output_dir=output_dir,
        threads=pipeline.config.threads,
        split_reads=True,
    )

    pipeline.logger.info("Converting Cresil output to Cyrcular-compatible format")
    overview_tsv = pipeline.config.output_dir / f"{pipeline.config.prefix}_overview.tsv"
    write_cyrcular_compatible_tsv(eccDNA_final, overview_tsv)

    if overview_tsv.exists():
        pipeline.state.results[ResultKeys.CYRCULAR_OVERVIEW] = pipeline._serialize_path_for_state(
            overview_tsv
        )
        import pandas as pd

        df = pd.read_csv(overview_tsv, sep="\t")
        pipeline.state.results[ResultKeys.CYRCULAR_RESULT_COUNT] = len(df)
        pipeline.logger.info(f"Cresil detected {len(df)} circular DNA candidates")
    else:
        pipeline.logger.warning("Cresil output conversion failed")
        pipeline.state.results[ResultKeys.CYRCULAR_RESULT_COUNT] = 0


def run_cyrcular_inference(pipeline: Pipeline) -> None:
    """Run Cyrcular inference pipeline (original implementation)."""
    from circleseeker.modules.cyrcular_calling import PipelineConfig as CCConfig
    from circleseeker.modules.cyrcular_calling import CyrcularCallingPipeline

    bam_file = pipeline.state.results.get(ResultKeys.MINIMAP2_BAM)
    if not bam_file or not Path(bam_file).exists():
        pipeline.logger.warning("No BAM file from minimap2, skipping Cyrcular-Calling step")
        return

    reference = pipeline.config.reference
    if reference is None:
        raise PipelineError("Reference genome is required")

    cfg = CCConfig(
        bam_file=Path(bam_file),
        reference=reference,
        output_dir=pipeline.config.output_dir,
        sample_name=pipeline.config.prefix,
        threads=pipeline.config.threads,
    )

    cyrcular_pipeline = CyrcularCallingPipeline(cfg, logger=pipeline.logger.getChild("cyrcular_calling"))
    results = cyrcular_pipeline.run()

    overview = cyrcular_pipeline.file_paths.get("overview_table")
    if overview and overview.exists():
        pipeline.state.results[ResultKeys.CYRCULAR_OVERVIEW] = pipeline._serialize_path_for_state(
            overview
        )
    else:
        pipeline.logger.warning(
            "Cyrcular overview table not found; inferred eccDNA curation may be skipped"
        )
    pipeline.state.results[ResultKeys.CYRCULAR_RESULT_COUNT] = len(results) if results else 0


def curate_inferred_ecc(pipeline: Pipeline) -> list[str] | None:
    """Step 12: Curate inferred eccDNA tables from Cyrcular overview TSV."""
    from circleseeker.modules.iecc_curator import (
        curate_ecc_tables,
        write_curated_tables,
        write_curated_tables_with_fasta,
    )

    if pipeline.state.results.get(ResultKeys.INFERENCE_INPUT_EMPTY) or pipeline.state.results.get(
        ResultKeys.INFERENCE_FAILED
    ):
        pipeline.logger.warning("Inference was skipped/failed; skipping inferred eccDNA curation step")
        pipeline.state.results[ResultKeys.INFERRED_SIMPLE_TSV] = None
        pipeline.state.results[ResultKeys.INFERRED_CHIMERIC_TSV] = None
        pipeline.state.results[ResultKeys.INFERRED_SIMPLE_FASTA] = None
        pipeline.state.results[ResultKeys.INFERRED_CHIMERIC_FASTA] = None
        pipeline.state.results[ResultKeys.INFERRED_DIR] = None
        return None

    overview_path_str = pipeline.state.results.get(ResultKeys.CYRCULAR_OVERVIEW)

    overview_path_obj = pipeline._resolve_stored_path(
        overview_path_str,
        [f"{pipeline.config.prefix}_overview.tsv"],
    )
    if not overview_path_obj:
        pipeline.logger.warning("Cyrcular overview table missing; skipping inferred eccDNA curation")
        pipeline.state.results[ResultKeys.INFERRED_SIMPLE_TSV] = None
        pipeline.state.results[ResultKeys.INFERRED_CHIMERIC_TSV] = None
        pipeline.state.results[ResultKeys.INFERRED_SIMPLE_FASTA] = None
        pipeline.state.results[ResultKeys.INFERRED_CHIMERIC_FASTA] = None
        pipeline.state.results[ResultKeys.INFERRED_DIR] = None
        return None

    pipeline.state.results[ResultKeys.CYRCULAR_OVERVIEW] = pipeline._serialize_path_for_state(
        overview_path_obj
    )

    simple_df, chimeric_df = curate_ecc_tables(overview_path_obj)

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
    except Exception as exc:
        pipeline.logger.warning(f"Failed to generate inferred eccDNA tables with FASTA: {exc}")
        simple_csv_path, chimeric_csv_path = write_curated_tables(
            simple_df,
            chimeric_df,
            inferred_prefix,
        )

    written: list[str] = []
    if simple_csv_path and simple_csv_path.exists() and simple_csv_path.stat().st_size > 0:
        pipeline.state.results[ResultKeys.INFERRED_SIMPLE_TSV] = pipeline._serialize_path_for_state(
            simple_csv_path
        )
        written.append(str(simple_csv_path))
    else:
        pipeline.logger.info("No inferred simple eccDNA records; simple TSV not created")

    if chimeric_csv_path and chimeric_csv_path.exists() and chimeric_csv_path.stat().st_size > 0:
        pipeline.state.results[ResultKeys.INFERRED_CHIMERIC_TSV] = pipeline._serialize_path_for_state(
            chimeric_csv_path
        )
        written.append(str(chimeric_csv_path))
    else:
        pipeline.logger.info("No inferred chimeric eccDNA records; chimeric TSV not created")

    if simple_fasta_path and simple_fasta_path.exists() and simple_fasta_path.stat().st_size > 0:
        pipeline.state.results[ResultKeys.INFERRED_SIMPLE_FASTA] = pipeline._serialize_path_for_state(
            simple_fasta_path
        )

    if chimeric_fasta_path and chimeric_fasta_path.exists() and chimeric_fasta_path.stat().st_size > 0:
        pipeline.state.results[ResultKeys.INFERRED_CHIMERIC_FASTA] = pipeline._serialize_path_for_state(
            chimeric_fasta_path
        )

    inferred_dir_path = None
    if simple_csv_path:
        inferred_dir_path = simple_csv_path.parent
    elif chimeric_csv_path:
        inferred_dir_path = chimeric_csv_path.parent
    if inferred_dir_path and inferred_dir_path.exists():
        pipeline.state.results[ResultKeys.INFERRED_DIR] = pipeline._serialize_path_for_state(
            inferred_dir_path
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
