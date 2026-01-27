"""Step executors for final merging, reporting, and packaging."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import TYPE_CHECKING, Optional

import pandas as pd

from circleseeker.core.pipeline_types import ResultKeys
from circleseeker.exceptions import PipelineError

if TYPE_CHECKING:
    from circleseeker.core.pipeline import Pipeline


def ecc_unify(pipeline: Pipeline) -> None:
    """Step 13: Merge eccDNA tables into unified output."""
    from circleseeker.modules.ecc_unify import merge_eccdna_tables

    confirmed_csv_path = pipeline._resolve_stored_path(
        pipeline.state.results.get(ResultKeys.CONFIRMED_CSV),
        [
            f"{pipeline.config.prefix}_confirmed.csv",
            f"{pipeline.config.prefix}_eccDNA_Confirmed.csv",
        ],
    )

    if not confirmed_csv_path:
        pipeline.logger.info("No confirmed eccDNA file found, skipping ecc_unify step")
        return

    pipeline._set_result(
        ResultKeys.CONFIRMED_CSV,
        pipeline._serialize_path_for_state(confirmed_csv_path),
    )

    inference_disabled = bool(
        pipeline._get_result(ResultKeys.INFERENCE_INPUT_EMPTY)
        or pipeline._get_result(ResultKeys.INFERENCE_FAILED)
    )

    simple_csv_path = None
    chimeric_csv_path = None
    if not inference_disabled:
        simple_csv_path = pipeline._resolve_stored_path(
            pipeline.state.results.get(ResultKeys.INFERRED_SIMPLE_TSV),
            [
                f"{pipeline.config.prefix}_inferred_simple.csv",
                Path(f"{pipeline.config.prefix}_Inferred_eccDNA") / f"{pipeline.config.prefix}_simple.csv",
            ],
        )

        chimeric_csv_path = pipeline._resolve_stored_path(
            pipeline.state.results.get(ResultKeys.INFERRED_CHIMERIC_TSV),
            [
                f"{pipeline.config.prefix}_inferred_chimeric.csv",
                Path(f"{pipeline.config.prefix}_Inferred_eccDNA") / f"{pipeline.config.prefix}_chimeric.csv",
            ],
        )

    unified_csv = pipeline.config.output_dir / f"{pipeline.config.prefix}_unified.csv"
    overlap_report = pipeline.config.output_dir / f"{pipeline.config.prefix}_overlap_report.txt"
    overlap_stats = pipeline.config.output_dir / f"{pipeline.config.prefix}_overlap_stats.json"

    pipeline.logger.info("Running ecc_unify module")
    try:
        merged_df, report_text, stats_dict = merge_eccdna_tables(
            confirmed_file=str(confirmed_csv_path),
            inferred_simple=str(simple_csv_path) if simple_csv_path and simple_csv_path.exists() else None,
            inferred_chimeric=str(chimeric_csv_path) if chimeric_csv_path and chimeric_csv_path.exists() else None,
            overlap_report_file=overlap_report,
            overlap_stats_json=overlap_stats,
        )

        merged_df.to_csv(unified_csv, index=False)
        pipeline.logger.info(f"Unified {len(merged_df)} eccDNA entries")
    except (OSError, pd.errors.ParserError, ValueError, KeyError) as e:
        pipeline.logger.error(f"ecc_unify failed: {e}")
        raise PipelineError(f"ecc_unify failed: {e}") from e

    if unified_csv.exists():
        pipeline._set_result(ResultKeys.UNIFIED_CSV, pipeline._serialize_path_for_state(unified_csv))
    if overlap_report.exists():
        pipeline._set_result(ResultKeys.OVERLAP_REPORT, pipeline._serialize_path_for_state(overlap_report))
    if overlap_stats.exists():
        pipeline._set_result(ResultKeys.OVERLAP_STATS, pipeline._serialize_path_for_state(overlap_stats))


def ecc_summary(pipeline: Pipeline) -> None:
    """Step 14: Generate summary report."""
    from circleseeker.modules.ecc_summary import EccSummary

    processed_csv = pipeline.config.output_dir / "tandem_to_ring.csv"
    unified_csv_path = pipeline._resolve_stored_path(
        pipeline.state.results.get(ResultKeys.UNIFIED_CSV),
        [f"{pipeline.config.prefix}_unified.csv"],
    )
    confirmed_csv_path = pipeline._resolve_stored_path(
        pipeline.state.results.get(ResultKeys.CONFIRMED_CSV),
        [
            f"{pipeline.config.prefix}_confirmed.csv",
            f"{pipeline.config.prefix}_eccDNA_Confirmed.csv",
        ],
    )
    stats_json_path = pipeline._resolve_stored_path(
        pipeline.state.results.get(ResultKeys.OVERLAP_STATS),
        [f"{pipeline.config.prefix}_overlap_stats.json"],
    )

    main_csv = unified_csv_path if unified_csv_path else confirmed_csv_path

    if not (processed_csv.exists() and main_csv and main_csv.exists()):
        pipeline.logger.info("Missing required files for summary, skipping")
        return

    all_fasta = pipeline.config.output_dir / f"{pipeline.config.prefix}_all.fasta"
    pipeline._create_combined_fasta(all_fasta)

    stats_json = stats_json_path if stats_json_path else (pipeline.config.output_dir / f"{pipeline.config.prefix}_overlap_stats.json")
    if not stats_json.exists():
        try:
            stats_json.parent.mkdir(parents=True, exist_ok=True)
            import json as _json

            with open(stats_json, "w") as f:
                _json.dump({"overlap_stats": {}}, f)
        except OSError as e:
            pipeline.logger.warning(f"Failed to create stats file {stats_json}: {e}")

    output_dir = pipeline.config.output_dir / "summary_output"
    output_dir.mkdir(exist_ok=True)

    summary = EccSummary(
        sample_name=pipeline.config.prefix,
        output_dir=output_dir,
        logger=pipeline.logger.getChild("ecc_summary"),
    )

    input_file = pipeline.config.input_file
    if input_file is None:
        pipeline.logger.info("No input FASTA configured for summary, skipping")
        return
    summary.process_fasta(input_file)
    summary.process_processed_csv(processed_csv)
    summary.process_merged_csv(main_csv)
    if stats_json.exists():
        summary.process_overlap_stats(stats_json)

    summary.generate_html_report()
    summary.generate_text_summary()

    pipeline._set_result(ResultKeys.SUMMARY_DIR, str(output_dir))


def ecc_packager(pipeline: Pipeline) -> None:
    """Step 15: Package output files using ecc_packager module.

    Generates standardized output format:
    - eccDNA_summary.csv: One row per eccDNA with key metrics
    - eccDNA_regions.csv: All genomic regions
    - eccDNA_reads.csv: Read-level support
    - Type-specific BED and FASTA files in Uecc/, Mecc/, Cecc/ subdirs
    """
    from circleseeker.modules import ecc_packager as ecc_packager_module

    final_output_root = pipeline.final_output_dir
    final_output_root.mkdir(parents=True, exist_ok=True)

    def resolve_dir(key: str, default_name: str) -> Optional[Path]:
        stored = pipeline.state.results.get(key)
        resolved = pipeline._resolve_stored_path(stored, [default_name])
        if resolved and resolved.exists():
            return resolved
        fallback_path = pipeline.config.output_dir / default_name
        return fallback_path if fallback_path.exists() else None

    uecc_dir = resolve_dir("ecc_dedup_uecc_dir", f"{pipeline.config.prefix}_Uecc_C")
    mecc_dir = resolve_dir("ecc_dedup_mecc_dir", f"{pipeline.config.prefix}_Mecc_C")
    cecc_dir = resolve_dir("ecc_dedup_cecc_dir", f"{pipeline.config.prefix}_Cecc_C")

    inference_disabled = bool(
        pipeline._get_result(ResultKeys.INFERENCE_INPUT_EMPTY)
        or pipeline._get_result(ResultKeys.INFERENCE_FAILED)
    )
    inferred_dir = None
    if not inference_disabled:
        inferred_dir = pipeline._resolve_stored_path(
            pipeline.state.results.get(ResultKeys.INFERRED_DIR),
            [f"{pipeline.config.prefix}_Inferred_eccDNA"],
        )
        if inferred_dir and not inferred_dir.exists():
            inferred_dir = None

    # Get unified CSV (required for new output format)
    unified_csv_path = pipeline._resolve_stored_path(
        pipeline.state.results.get(ResultKeys.UNIFIED_CSV),
        [f"{pipeline.config.prefix}_unified.csv"],
    )
    unified_csv = unified_csv_path if unified_csv_path and unified_csv_path.exists() else None

    confirmed_csv_path = pipeline._resolve_stored_path(
        pipeline.state.results.get(ResultKeys.CONFIRMED_CSV),
        [
            f"{pipeline.config.prefix}_confirmed.csv",
            f"{pipeline.config.prefix}_eccDNA_Confirmed.csv",
        ],
    )
    confirmed_csv = confirmed_csv_path if confirmed_csv_path and confirmed_csv_path.exists() else None

    summary_dir = pipeline.config.output_dir / "summary_output"
    html_report_path = summary_dir / f"{pipeline.config.prefix}_report.html"
    html_report: Optional[Path] = html_report_path if html_report_path.exists() else None
    text_summary_path = summary_dir / f"{pipeline.config.prefix}_summary.txt"
    text_summary: Optional[Path] = text_summary_path if text_summary_path.exists() else None

    # Check if we have the minimum required input (unified CSV)
    if unified_csv is None:
        pipeline.logger.warning("No unified CSV found, cannot run ecc_packager")
        pipeline._create_basic_output_structure(
            final_output_root,
            uecc_dir=uecc_dir,
            mecc_dir=mecc_dir,
            cecc_dir=cecc_dir,
            inferred_dir=inferred_dir,
            merged_csv=confirmed_csv,
            html_report=html_report,
            text_summary=text_summary,
        )
        pipeline._set_result(ResultKeys.FINAL_RESULTS, str(final_output_root))
        return

    self_config = pipeline.config
    runtime = getattr(self_config, "runtime", None)
    debug_verbose = bool(getattr(runtime, "debug_verbose", False)) if runtime else False

    mock_args = argparse.Namespace(
        sample_name=self_config.prefix,
        out_dir=str(final_output_root),
        unified_csv=str(unified_csv),
        uecc_dir=str(uecc_dir) if uecc_dir else None,
        mecc_dir=str(mecc_dir) if mecc_dir else None,
        cecc_dir=str(cecc_dir) if cecc_dir else None,
        inferred_dir=str(inferred_dir) if inferred_dir else None,
        html=str(html_report) if html_report else None,
        text=str(text_summary) if text_summary else None,
        id_width=4,  # Zero-padding width for IDs
        overwrite=True,
        dry_run=False,
        verbose=debug_verbose,
    )

    pipeline.logger.info("Running ecc_packager module to organize final results")
    packager_success = False
    try:
        result_code = ecc_packager_module.run(mock_args)
        if result_code == 0:
            packager_success = True
            pipeline.logger.info(f"Results successfully packaged in: {final_output_root}")
        else:
            pipeline.logger.warning("ecc_packager completed with warnings")
    except (OSError, ValueError, RuntimeError) as exc:
        pipeline.logger.warning(f"ecc_packager encountered issues: {exc}")

    if not packager_success:
        pipeline._create_basic_output_structure(
            final_output_root,
            uecc_dir=uecc_dir,
            mecc_dir=mecc_dir,
            cecc_dir=cecc_dir,
            inferred_dir=inferred_dir,
            merged_csv=unified_csv,
            html_report=html_report,
            text_summary=text_summary,
        )

    pipeline._set_result(ResultKeys.FINAL_RESULTS, str(final_output_root))
