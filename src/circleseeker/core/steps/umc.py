"""Step executors for U/M/C confirmed eccDNA path."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any, Optional, cast

from circleseeker.core.pipeline_types import ResultKeys
from circleseeker.exceptions import PipelineError

if TYPE_CHECKING:
    from circleseeker.core.pipeline import Pipeline


def um_classify(pipeline: Pipeline) -> None:
    """Step 4: Classify eccDNA types using um_classify module."""
    import pandas as pd

    from circleseeker.modules.um_classify import UMeccClassifier

    alignment_output = pipeline._resolve_stored_path(
        pipeline.state.results.get(ResultKeys.ALIGNMENT_OUTPUT),
        [f"{pipeline.config.prefix}_alignment_results.tsv"],
    )
    if not alignment_output:
        raise PipelineError(
            "Alignment output not found; run the alignment step before um_classify."
        )
    pipeline.state.results[ResultKeys.ALIGNMENT_OUTPUT] = pipeline._serialize_path_for_state(
        alignment_output
    )
    output_prefix = pipeline.config.output_dir / "um_classify"

    um_cfg = pipeline.config.tools.um_classify if hasattr(pipeline.config.tools, "um_classify") else {}
    if not isinstance(um_cfg, dict):
        raise PipelineError(
            "Invalid um_classify config; expected mapping, "
            f"got {type(um_cfg).__name__}"
        )
    # Get min_full_length_coverage with fallback to theta_full
    min_full_len_cov = um_cfg.get("min_full_length_coverage")
    if min_full_len_cov is None:
        min_full_len_cov = um_cfg.get("theta_full", 0.95)
    classifier = UMeccClassifier(
        gap_threshold=float(um_cfg.get("gap_threshold", 10.0)),
        min_full_length_coverage=float(min_full_len_cov),
        max_identity_gap_for_mecc=um_cfg.get("max_identity_gap_for_mecc", 5.0),
        theta_full=um_cfg.get("theta_full"),
        theta_u=um_cfg.get("theta_u"),
        theta_m=um_cfg.get("theta_m"),
        theta_u2_max=um_cfg.get("theta_u2_max"),
        mapq_u_min=um_cfg.get("mapq_u_min", 0),
        u_secondary_min_frac=um_cfg.get("u_secondary_min_frac", 0.01),
        u_secondary_min_bp=um_cfg.get("u_secondary_min_bp", 50),
        u_contig_gap_bp=um_cfg.get("u_contig_gap_bp", 1000),
        u_secondary_max_ratio=um_cfg.get("u_secondary_max_ratio", 0.05),
        u_high_coverage_threshold=um_cfg.get("u_high_coverage_threshold", 0.98),
        u_high_mapq_threshold=um_cfg.get("u_high_mapq_threshold", 50),
        theta_locus=um_cfg.get("theta_locus", 0.95),
        pos_tol_bp=um_cfg.get("pos_tol_bp", 50),
        logger=pipeline.logger.getChild("um_classify"),
    )

    pipeline.logger.info("Running um_classify module")

    if not alignment_output.exists() or alignment_output.stat().st_size == 0:
        pipeline.logger.warning(
            "Alignment output is empty - no alignments found. "
            "This may occur when reads have no similarity to the reference."
        )
        pipeline.state.results[ResultKeys.UECC_COUNT] = 0
        pipeline.state.results[ResultKeys.MECC_COUNT] = 0
        return

    alignment_df = pd.read_csv(alignment_output, sep="\t", header=None)

    try:
        query_ids = alignment_df.iloc[:, 0].astype(str)
        valid_mask = query_ids.str.count(r"\|") >= 3
    except Exception as exc:
        pipeline.logger.error(f"Failed to validate alignment query IDs: {exc}")
        raise PipelineError("Alignment query_id validation failed") from exc

    if not valid_mask.any():
        pipeline.logger.error(
            "Alignment query_id format is invalid for um_classify. "
            "Expected 'read|repN|length|copy' (at least 4 '|' fields). "
            "Skipping classification."
        )
        uecc_output = pipeline.config.output_dir / "um_classify.uecc.csv"
        mecc_output = pipeline.config.output_dir / "um_classify.mecc.csv"
        unclass_output = pipeline.config.output_dir / "um_classify.unclassified.csv"
        pd.DataFrame().to_csv(uecc_output, index=False)
        pd.DataFrame().to_csv(mecc_output, index=False)
        pd.DataFrame().to_csv(unclass_output, index=False)
        pipeline.state.results[ResultKeys.UECC_CSV] = str(uecc_output)
        pipeline.state.results[ResultKeys.MECC_CSV] = str(mecc_output)
        pipeline.state.results[ResultKeys.UNCLASSIFIED_CSV] = str(unclass_output)
        pipeline.state.results[ResultKeys.UECC_COUNT] = 0
        pipeline.state.results[ResultKeys.MECC_COUNT] = 0
        return
    if not valid_mask.all():
        pipeline.logger.warning(
            "Some alignment query IDs lack the expected 'read|repN|length|copy' format; "
            "those entries will be ignored during classification."
        )
        alignment_df = alignment_df[valid_mask].copy()

    classifier.classify_alignment_results(alignment_df, output_prefix)

    uecc_output = pipeline.config.output_dir / "um_classify.uecc.csv"
    mecc_output = pipeline.config.output_dir / "um_classify.mecc.csv"
    unclass_output = pipeline.config.output_dir / "um_classify.unclassified.csv"

    if uecc_output.exists():
        uecc_df = pd.read_csv(uecc_output)
        pipeline.state.results[ResultKeys.UECC_CSV] = str(uecc_output)
        pipeline.state.results[ResultKeys.UECC_COUNT] = len(uecc_df)

    if mecc_output.exists():
        mecc_df = pd.read_csv(mecc_output)
        pipeline.state.results[ResultKeys.MECC_CSV] = str(mecc_output)
        pipeline.state.results[ResultKeys.MECC_COUNT] = len(mecc_df)

    if unclass_output.exists():
        pipeline.state.results[ResultKeys.UNCLASSIFIED_CSV] = str(unclass_output)


def cecc_build(pipeline: Pipeline) -> None:
    """Step 5: Process Cecc candidates using cecc_build module."""
    import pandas as pd
    import circleseeker.core.pipeline as pipeline_module

    unclass_input = pipeline.config.output_dir / "um_classify.unclassified.csv"
    input_csv = unclass_input

    def _has_data_rows(path: Path) -> bool:
        try:
            with open(path, "r", encoding="utf-8", errors="replace") as handle:
                header = handle.readline()
                if not header.strip():
                    return False
                for line in handle:
                    if line.strip():
                        return True
        except OSError:
            return False
        return False

    if not input_csv.exists() or not _has_data_rows(input_csv):
        pipeline.logger.warning(
            "No candidates found for CECC analysis (missing/empty %s).",
            input_csv.name,
        )
        output_file = pipeline.config.output_dir / "cecc_build.csv"
        pd.DataFrame().to_csv(output_file, index=False)
        pipeline._set_result(ResultKeys.CECC_BUILD_OUTPUT, str(output_file))
        pipeline._set_result(ResultKeys.CECC_BUILD_COUNT, 0)
        return

    output_file = pipeline.config.output_dir / "cecc_build.csv"

    from circleseeker.modules.cecc_build import CeccBuild as DefaultCeccBuild

    builder_cls = getattr(pipeline_module, "CeccBuild", None)
    if builder_cls is None:
        builder_cls = DefaultCeccBuild

    builder = cast(type[Any], builder_cls)(logger=pipeline.logger.getChild("cecc_build"))

    pipeline.logger.info("Running cecc_build module")
    try:
        cecc_cfg = getattr(pipeline.config.tools, "cecc_build", {}) or {}
        if not isinstance(cecc_cfg, dict):
            raise PipelineError(
                "Invalid cecc_build config; expected mapping, "
                f"got {type(cecc_cfg).__name__}"
            )

        tau_gap = cecc_cfg.get("tau_gap") or cecc_cfg.get("edge_tolerance") or 20
        min_match_degree = cecc_cfg.get("min_match_degree", 95.0)
        theta_chain = cecc_cfg.get("theta_chain")
        if theta_chain is not None:
            try:
                theta_chain_val = float(theta_chain)
                min_match_degree = (
                    theta_chain_val * 100.0 if theta_chain_val <= 1.0 else theta_chain_val
                )
            except (TypeError, ValueError):
                pass

        # Prepare LAST-based detection parameters
        reference_fasta = getattr(pipeline.config, "reference", None)
        if reference_fasta:
            reference_fasta = Path(reference_fasta)
            if not reference_fasta.exists():
                reference_fasta = None

        fasta_file = pipeline.config.output_dir / "tandem_to_ring.fasta"
        if not fasta_file.exists():
            fasta_file = pipeline.config.output_dir / f"{pipeline.config.prefix}_circular.fasta"
            if not fasta_file.exists():
                fasta_file = None

        builder.run_pipeline(
            input_csv=input_csv,
            output_csv=output_file,
            overlap_threshold=cecc_cfg.get("overlap_threshold", 0.95),
            min_segments=int(cecc_cfg.get("min_segments", 2)),
            edge_tolerance=int(tau_gap),
            position_tolerance=int(cecc_cfg.get("position_tolerance", 50)),
            min_match_degree=float(min_match_degree),
            max_rotations=int(cecc_cfg.get("max_rotations", 20)),
            locus_overlap_threshold=cecc_cfg.get("locus_overlap_threshold", 0.95),
            reference_fasta=reference_fasta,
            fasta_file=fasta_file,
        )

        df_cecc = (
            pd.read_csv(output_file)
            if output_file.exists() and output_file.stat().st_size > 1
            else pd.DataFrame()
        )

        # Safety guard: even though we only analyze unclassified reads, ensure we never emit
        # Cecc calls that overlap with existing U/Mecc classifications.
        if not df_cecc.empty and "query_id" in df_cecc.columns:
            uecc_path = pipeline.config.output_dir / "um_classify.uecc.csv"
            mecc_path = pipeline.config.output_dir / "um_classify.mecc.csv"
            uecc_df = (
                pd.read_csv(uecc_path)
                if uecc_path.exists() and uecc_path.stat().st_size > 1
                else pd.DataFrame()
            )
            mecc_df = (
                pd.read_csv(mecc_path)
                if mecc_path.exists() and mecc_path.stat().st_size > 1
                else pd.DataFrame()
            )
            uecc_ids = set(uecc_df.get("query_id", pd.Series([], dtype=str)).astype(str))
            mecc_ids = set(mecc_df.get("query_id", pd.Series([], dtype=str)).astype(str))
            exclude_from_c = uecc_ids | mecc_ids
            if exclude_from_c:
                df_cecc = df_cecc[~df_cecc["query_id"].astype(str).isin(exclude_from_c)].copy()
                df_cecc.to_csv(output_file, index=False)

        pipeline._set_result(ResultKeys.CECC_BUILD_OUTPUT, str(output_file))
        pipeline._set_result(
            ResultKeys.CECC_BUILD_COUNT,
            int(df_cecc["query_id"].nunique())
            if (not df_cecc.empty and "query_id" in df_cecc.columns)
            else 0,
        )
    except Exception as e:
        pipeline.logger.warning(f"cecc_build failed: {e}")
        pd.DataFrame().to_csv(output_file, index=False)
        pipeline._set_result(ResultKeys.CECC_BUILD_OUTPUT, str(output_file))
        pipeline._set_result(ResultKeys.CECC_BUILD_COUNT, 0)


def umc_process(pipeline: Pipeline) -> None:
    """Step 6: Generate FASTA files using umc_process module."""
    from circleseeker.modules.umc_process import UMCProcess

    circular_fasta = pipeline.config.output_dir / "tandem_to_ring.fasta"
    if not circular_fasta.exists():
        circular_fasta = pipeline.config.output_dir / f"{pipeline.config.prefix}_circular.fasta"
        if not circular_fasta.exists():
            pipeline.logger.warning("Circular FASTA file not found, skipping umc_process step")
            return

    uecc_csv = pipeline.config.output_dir / "um_classify.uecc.csv"
    mecc_csv = pipeline.config.output_dir / "um_classify.mecc.csv"
    cecc_csv = pipeline.config.output_dir / "cecc_build.csv"

    processor = UMCProcess(logger=pipeline.logger.getChild("umc_process"))

    pipeline.logger.info("Running umc_process module")
    processor.run(
        fasta_file=circular_fasta,
        uecc_csv=uecc_csv if uecc_csv.exists() else None,
        mecc_csv=mecc_csv if mecc_csv.exists() else None,
        cecc_csv=cecc_csv if cecc_csv.exists() else None,
        output_dir=pipeline.config.output_dir,
        prefix=pipeline.config.prefix,
    )

    type_mapping = {"Uecc": "uecc", "Mecc": "mecc", "Cecc": "cecc"}

    for ecc_name, ecc_type in type_mapping.items():
        fasta_file = pipeline.config.output_dir / f"{pipeline.config.prefix}_{ecc_name}DNA_pre.fasta"
        if not fasta_file.exists():
            continue

        pipeline.state.results[f"{ecc_type}_fasta"] = str(fasta_file)

        processed_csv = pipeline.config.output_dir / f"{pipeline.config.prefix}_{ecc_name}DNA_processed.csv"
        if not processed_csv.exists():
            fallback_csv = fasta_file.with_suffix(".csv")
            processed_csv = fallback_csv if fallback_csv.exists() else processed_csv

        if processed_csv.exists():
            pipeline.state.results[f"{ecc_type}_processed"] = str(processed_csv)

        pipeline.logger.debug(
            "Found %s outputs: fasta=%s csv=%s",
            ecc_type,
            fasta_file,
            processed_csv if processed_csv.exists() else "missing",
        )


def cd_hit(pipeline: Pipeline) -> None:
    """Step 7: Remove redundant sequences using CD-HIT-EST wrapper."""
    from circleseeker.external.cd_hit import CDHitEst

    fasta_specs = [
        ("uecc", f"{pipeline.config.prefix}_UeccDNA_pre.fasta", f"{pipeline.config.prefix}_U"),
        ("mecc", f"{pipeline.config.prefix}_MeccDNA_pre.fasta", f"{pipeline.config.prefix}_M"),
        ("cecc", f"{pipeline.config.prefix}_CeccDNA_pre.fasta", f"{pipeline.config.prefix}_C"),
    ]

    cd_hit_tool = CDHitEst(logger=pipeline.logger.getChild("cd_hit"), threads=pipeline.config.threads)

    for ecc_type, input_name, output_prefix in fasta_specs:
        input_fasta = pipeline.config.output_dir / input_name
        if not input_fasta.exists():
            pipeline.logger.warning(f"Skipping cd-hit for {input_name}: file not found")
            continue

        output_path = pipeline.config.output_dir / output_prefix

        try:
            clstr_path = cd_hit_tool.cluster_sequences(input_fasta, output_path)
        except Exception as e:
            pipeline.logger.warning(f"cd-hit failed for {input_name}: {e}")
            continue

        rep_candidates = [output_path, output_path.with_suffix(".fasta")]
        output_fasta = next((p for p in rep_candidates if p.exists()), None)
        cluster_file = (
            clstr_path
            if clstr_path and Path(clstr_path).exists()
            else output_path.with_suffix(".clstr")
        )

        if output_fasta is not None:
            pipeline.state.results[f"{ecc_type}_nr99"] = str(output_fasta)
        if cluster_file and Path(cluster_file).exists():
            pipeline.state.results[f"{ecc_type}_clusters"] = str(cluster_file)


def _get_path_from_results(pipeline: Pipeline, key: str) -> Optional[Path]:
    """Retrieve a Path from pipeline results if the key exists.

    Args:
        pipeline: The pipeline instance
        key: The ResultKeys key to look up

    Returns:
        Path if key exists and file exists, None otherwise
    """
    if key in pipeline.state.results:
        path = Path(pipeline.state.results[key])
        return path if path.exists() else None
    return None


def _collect_ecc_dedup_inputs(pipeline: Pipeline) -> dict[str, Optional[Path]]:
    """Collect input files for the ecc_dedup step.

    Args:
        pipeline: The pipeline instance

    Returns:
        Dictionary with input and cluster paths for each eccDNA type
    """
    return {
        "uecc_input": _get_path_from_results(pipeline, ResultKeys.UECC_PROCESSED),
        "uecc_cluster": _get_path_from_results(pipeline, ResultKeys.UECC_CLUSTERS),
        "mecc_input": _get_path_from_results(pipeline, ResultKeys.MECC_PROCESSED),
        "mecc_cluster": _get_path_from_results(pipeline, ResultKeys.MECC_CLUSTERS),
        "cecc_input": _get_path_from_results(pipeline, ResultKeys.CECC_PROCESSED),
        "cecc_cluster": _get_path_from_results(pipeline, ResultKeys.CECC_CLUSTERS),
    }


def _rename_dedup_outputs(pipeline: Pipeline) -> None:
    """Rename eccDedup output files to standard names.

    Args:
        pipeline: The pipeline instance
    """
    rename_mappings = [
        (f"{pipeline.config.prefix}_Mecc.fa", f"{pipeline.config.prefix}_MeccDNA_C.fasta"),
        (f"{pipeline.config.prefix}_Cecc.fa", f"{pipeline.config.prefix}_CeccDNA_C.fasta"),
        (f"{pipeline.config.prefix}_UeccDNA.fa", f"{pipeline.config.prefix}_UeccDNA_C.fasta"),
    ]

    for old_name, new_name in rename_mappings:
        old_path = pipeline.config.output_dir / old_name
        new_path = pipeline.config.output_dir / new_name
        if old_path.exists() and not new_path.exists():
            old_path.rename(new_path)
            pipeline.logger.debug(f"Renamed {old_name} to {new_name}")


def ecc_dedup(pipeline: Pipeline) -> None:
    """Step 8: Coordinate and deduplicate results using ecc_dedup module."""
    from circleseeker.modules.ecc_dedup import eccDedup, organize_umc_files

    harmonizer = eccDedup(pipeline.logger.getChild("ecc_dedup"))
    output_dir = Path(pipeline.config.output_dir)

    # Collect input files
    inputs = _collect_ecc_dedup_inputs(pipeline)

    # Run deduplication
    results = harmonizer.run_deduplication(
        output_dir=pipeline.config.output_dir,
        prefix=pipeline.config.prefix,
        uecc_input=inputs["uecc_input"],
        uecc_cluster=inputs["uecc_cluster"],
        mecc_input=inputs["mecc_input"],
        mecc_cluster=inputs["mecc_cluster"],
        cecc_input=inputs["cecc_input"],
        cecc_cluster=inputs["cecc_cluster"],
    )

    # Rename output files to standard naming convention
    _rename_dedup_outputs(pipeline)

    umc_dirs: dict[str, Path] = {}
    if results:
        prefix_for_dirs = pipeline.config.prefix or "sample"
        try:
            base_lists = {
                "U": [
                    "UeccDNA_C.fasta",
                    "UeccDNA.bed",
                    "UeccDNA.core.csv",
                ],
                "M": [
                    "MeccDNA_C.fasta",
                    "MeccSites.bed",
                    "MeccBestSite.bed",
                    "MeccSites.core.csv",
                ],
                "C": [
                    "CeccDNA_C.fasta",
                    "CeccSegments.bed",
                    "CeccSegments.core.csv",
                    "CeccJunctions.bedpe",
                ],
            }

            def with_prefix(name: str) -> str:
                return f"{prefix_for_dirs}_{name}" if prefix_for_dirs else name

            files_to_move: dict[str, list[Path]] = {}
            for ecc_code, file_list in base_lists.items():
                collected: list[Path] = []
                for base_name in file_list:
                    candidate = pipeline.config.output_dir / with_prefix(base_name)
                    if candidate.exists():
                        collected.append(candidate)
                    elif not prefix_for_dirs:
                        fallback = pipeline.config.output_dir / base_name
                        if fallback.exists():
                            collected.append(fallback)
                if collected:
                    files_to_move[ecc_code] = collected

            umc_dirs = organize_umc_files(
                pipeline.config.output_dir,
                prefix_for_dirs,
                umc_files=files_to_move if files_to_move else None,
                auto_detect=False,
            )
            dir_key_map = {"U": "uecc", "M": "mecc", "C": "cecc"}
            for key, directory in umc_dirs.items():
                mapped = dir_key_map.get(key.upper())
                if mapped and directory.exists():
                    pipeline.state.results[f"ecc_dedup_{mapped}_dir"] = str(directory)
        except Exception as exc:
            pipeline.logger.warning(f"Failed to organize eccDedup outputs: {exc}")
            umc_dirs = {}

    def _locate_output(filename: str, ecc_code: str | None = None) -> Optional[Path]:
        candidate = output_dir / filename
        if candidate.exists():
            return candidate
        if ecc_code:
            code = ecc_code.upper()
            target_dir = umc_dirs.get(code)
            if target_dir:
                fallback = target_dir / filename
                if fallback.exists():
                    return fallback
        for fallback in output_dir.rglob(filename):
            if fallback.exists():
                return fallback
        return None

    folder_key_map = {"uecc": "U", "mecc": "M", "cecc": "C"}

    fasta_mappings = {
        "uecc": f"{pipeline.config.prefix}_UeccDNA_C.fasta",
        "mecc": f"{pipeline.config.prefix}_MeccDNA_C.fasta",
        "cecc": f"{pipeline.config.prefix}_CeccDNA_C.fasta",
    }

    for ecc_type, filename in fasta_mappings.items():
        fasta_file = _locate_output(filename, folder_key_map.get(ecc_type))
        if fasta_file:
            pipeline._set_result(f"ecc_dedup_{ecc_type}_fasta", str(fasta_file))
            pipeline.logger.debug(f"Found harmonizer output: {fasta_file}")
        else:
            pipeline.logger.warning(f"eccDedup output not found: {filename}")

    for ecc_type in ["Uecc", "Mecc", "Cecc"]:
        harmonized_csv = pipeline.config.output_dir / f"{pipeline.config.prefix}_{ecc_type}_harmonized.csv"
        if harmonized_csv.exists():
            key = {
                "Uecc": ResultKeys.UECC_HARMONIZED,
                "Mecc": ResultKeys.MECC_HARMONIZED,
                "Cecc": ResultKeys.CECC_HARMONIZED,
            }[ecc_type]
            pipeline.state.results[key] = str(harmonized_csv)

    additional_files = {
        "uecc": [
            (f"{pipeline.config.prefix}_UeccDNA.bed", "uecc_bed"),
            (f"{pipeline.config.prefix}_UeccDNA.core.csv", "uecc_core_csv"),
        ],
        "mecc": [
            (f"{pipeline.config.prefix}_MeccSites.bed", "mecc_sites_bed"),
            (f"{pipeline.config.prefix}_MeccBestSite.bed", "mecc_bestsite_bed"),
            (f"{pipeline.config.prefix}_MeccSites.core.csv", "mecc_core_csv"),
        ],
        "cecc": [
            (f"{pipeline.config.prefix}_CeccJunctions.bedpe", "cecc_junctions_bedpe"),
            (f"{pipeline.config.prefix}_CeccSegments.bed", "cecc_segments_bed"),
            (f"{pipeline.config.prefix}_CeccSegments.core.csv", "cecc_core_csv"),
        ],
    }

    for ecc_type, file_specs in additional_files.items():
        for filename, key in file_specs:
            file_path = _locate_output(filename, folder_key_map.get(ecc_type))
            if file_path:
                pipeline.state.results[key] = str(file_path)
                pipeline.logger.debug(f"Found {key}: {file_path}")

    confirmed_file = pipeline.config.output_dir / f"{pipeline.config.prefix}_eccDNA_Confirmed.csv"
    if confirmed_file.exists():
        pipeline.state.results[ResultKeys.CONFIRMED_CSV] = pipeline._serialize_path_for_state(
            confirmed_file
        )
