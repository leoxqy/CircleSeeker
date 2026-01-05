"""Step executors for inference and read filtering stages."""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional

from circleseeker.core.pipeline_types import ResultKeys
from circleseeker.exceptions import PipelineError


def read_filter(pipeline) -> None:
    """Step 9: Filter sequences using read_filter module (Sieve)."""
    from circleseeker.modules.read_filter import Sieve

    carousel_csv = pipeline.config.output_dir / "tandem_to_ring.csv"

    original_fasta = pipeline.config.input_file
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


def minimap2(pipeline) -> None:
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

    pipeline._set_result(ResultKeys.ALL_FILTERED_FASTA, str(input_fasta))

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

    config = MiniMapConfig(
        preset="map-hifi",
        threads=pipeline.config.threads,
        output_format="bam",
        allow_secondary=True,
        build_index=True,
        sort_bam=True,
        index_bam=True,
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


def ecc_inference(pipeline) -> None:
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


def run_cresil_inference(pipeline) -> None:
    """Run Cresil inference pipeline (replaces Steps 10-11 when available)."""
    from circleseeker.external.cresil import Cresil
    from circleseeker.modules.cresil_adapter import write_cyrcular_compatible_tsv

    sieve_output = pipeline._get_result(ResultKeys.READ_FILTER_OUTPUT_FASTA)

    if sieve_output and Path(sieve_output).exists():
        input_fasta = Path(sieve_output)
    else:
        input_fasta = pipeline.config.output_dir / f"{pipeline.config.prefix}_all_filtered.fasta"
        if not input_fasta.exists():
            pipeline.logger.warning("No filtered FASTA found, using original input")
            input_fasta = pipeline.config.input_file

    reference_mmi_str = pipeline._get_result(ResultKeys.REFERENCE_MMI)
    reference_mmi = Path(reference_mmi_str) if reference_mmi_str else None
    if not (reference_mmi and reference_mmi.exists()):
        reference_mmi = pipeline._ensure_reference_mmi()
        pipeline._set_result(ResultKeys.REFERENCE_MMI, str(reference_mmi))

    cresil = Cresil(logger=pipeline.logger.getChild("cresil"))

    output_dir = pipeline.config.output_dir / "cresil_output"
    output_dir.mkdir(parents=True, exist_ok=True)

    pipeline.logger.info("Running Cresil pipeline (trim + identify)")
    eccDNA_final = cresil.run_full_pipeline(
        fasta_query=input_fasta,
        reference_fasta=pipeline.config.reference,
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


def run_cyrcular_inference(pipeline) -> None:
    """Run Cyrcular inference pipeline (original implementation)."""
    from circleseeker.modules.cyrcular_calling import PipelineConfig as CCConfig
    from circleseeker.modules.cyrcular_calling import CyrcularCallingPipeline

    bam_file = pipeline.state.results.get(ResultKeys.MINIMAP2_BAM)
    if not bam_file or not Path(bam_file).exists():
        pipeline.logger.warning("No BAM file from minimap2, skipping Cyrcular-Calling step")
        return

    cfg = CCConfig(
        bam_file=Path(bam_file),
        reference=pipeline.config.reference,
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


def curate_inferred_ecc(pipeline) -> list[str] | None:
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

    written: List[str] = []
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

