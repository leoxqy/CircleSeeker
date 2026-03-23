"""Step executors for U/M/C confirmed eccDNA path."""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import TYPE_CHECKING, Optional

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
    pipeline._set_result(
        ResultKeys.ALIGNMENT_OUTPUT,
        pipeline._serialize_path_for_state(alignment_output)
    )
    output_prefix = pipeline.config.output_dir / f"{pipeline.config.prefix}_um_classify"

    um_cfg = pipeline.config.tools.um_classify if hasattr(pipeline.config.tools, "um_classify") else {}
    if not isinstance(um_cfg, dict) and not hasattr(um_cfg, 'get'):
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
        # Additional filters for ambiguous classifications
        mapq_m_ambiguous_threshold=um_cfg.get("mapq_m_ambiguous_threshold", 0),
        mecc_identity_gap_threshold=um_cfg.get("mecc_identity_gap_threshold", 1.0),
        u_secondary_min_frac=um_cfg.get("u_secondary_min_frac", 0.01),
        u_secondary_min_bp=um_cfg.get("u_secondary_min_bp", 50),
        u_contig_gap_bp=um_cfg.get("u_contig_gap_bp", 1000),
        u_secondary_max_ratio=um_cfg.get("u_secondary_max_ratio", 0.05),
        u_high_coverage_threshold=um_cfg.get("u_high_coverage_threshold", 0.98),
        u_high_mapq_threshold=um_cfg.get("u_high_mapq_threshold", 50),
        theta_locus=um_cfg.get("theta_locus", 0.95),
        pos_tol_bp=um_cfg.get("pos_tol_bp", 50),
        span_ratio_min=um_cfg.get("span_ratio_min", 0.95),
        logger=pipeline.logger.getChild("um_classify"),
    )

    pipeline.logger.debug("Running um_classify module")

    if not alignment_output.exists() or alignment_output.stat().st_size == 0:
        pipeline.logger.warning(
            "Alignment output is empty - no alignments found. "
            "This may occur when reads have no similarity to the reference."
        )
        pipeline._set_result(ResultKeys.UECC_COUNT, 0)
        pipeline._set_result(ResultKeys.MECC_COUNT, 0)
        return

    alignment_df = pd.read_csv(alignment_output, sep="\t", header=None)

    try:
        query_ids = alignment_df.iloc[:, 0].astype(str)
        valid_mask = query_ids.str.count(r"\|") >= 3
    except (pd.errors.ParserError, KeyError, TypeError) as exc:
        pipeline.logger.error(f"Failed to validate alignment query IDs: {exc}")
        raise PipelineError("Alignment query_id validation failed") from exc

    if not valid_mask.any():
        message = (
            "Alignment query_id format is invalid for um_classify. "
            "Expected 'read|repN|length|copy' from tandem_to_ring output. "
            "Check that tandem_to_ring.fasta was generated and used for alignment."
        )
        pipeline.logger.error(message)
        raise PipelineError(message)
    if not valid_mask.all():
        pipeline.logger.warning(
            "Some alignment query IDs lack the expected 'read|repN|length|copy' format; "
            "those entries will be ignored during classification."
        )
        alignment_df = alignment_df[valid_mask].copy()

    classifier.classify_alignment_results(alignment_df, output_prefix)

    uecc_output = pipeline.config.output_dir / f"{pipeline.config.prefix}_um_classify.uecc.csv"
    mecc_output = pipeline.config.output_dir / f"{pipeline.config.prefix}_um_classify.mecc.csv"
    unclass_output = pipeline.config.output_dir / f"{pipeline.config.prefix}_um_classify.unclassified.csv"

    if uecc_output.exists():
        try:
            uecc_df = pd.read_csv(uecc_output)
        except pd.errors.EmptyDataError:
            uecc_df = pd.DataFrame()
        pipeline._set_result(ResultKeys.UECC_CSV, str(uecc_output))
        pipeline._set_result(ResultKeys.UECC_COUNT, len(uecc_df))

    if mecc_output.exists():
        try:
            mecc_df = pd.read_csv(mecc_output)
        except pd.errors.EmptyDataError:
            mecc_df = pd.DataFrame()
        pipeline._set_result(ResultKeys.MECC_CSV, str(mecc_output))
        pipeline._set_result(ResultKeys.MECC_COUNT, len(mecc_df))

    if unclass_output.exists():
        pipeline._set_result(ResultKeys.UNCLASSIFIED_CSV, str(unclass_output))


def cecc_build(pipeline: Pipeline) -> None:
    """Step 5: Process Cecc candidates using cecc_build module."""
    import pandas as pd

    unclass_input = pipeline.config.output_dir / f"{pipeline.config.prefix}_um_classify.unclassified.csv"
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
        output_file = pipeline.config.output_dir / f"{pipeline.config.prefix}_cecc_build.csv"
        pd.DataFrame().to_csv(output_file, index=False)
        pipeline._set_result(ResultKeys.CECC_BUILD_OUTPUT, str(output_file))
        pipeline._set_result(ResultKeys.CECC_BUILD_COUNT, 0)
        return

    output_file = pipeline.config.output_dir / f"{pipeline.config.prefix}_cecc_build.csv"

    from circleseeker.modules.cecc_build import CeccBuild

    # Get threads from config (use performance.threads as the primary source)
    threads = int(getattr(pipeline.config, "threads", pipeline.config.performance.threads))

    cecc_cfg = getattr(pipeline.config.tools, "cecc_build", {}) or {}

    builder = CeccBuild(
        logger=pipeline.logger.getChild("cecc_build"),
        threads=threads,
        tmp_dir=pipeline.config.output_dir,
        keep_tmp=pipeline.config.keep_tmp,
        prefix=pipeline.config.prefix,
        fast_last=bool(cecc_cfg.get("fast_last", False)),
    )

    pipeline.logger.debug("Running cecc_build module")
    try:
        if not isinstance(cecc_cfg, dict) and not hasattr(cecc_cfg, 'get'):
            raise PipelineError(
                "Invalid cecc_build config; expected mapping, "
                f"got {type(cecc_cfg).__name__}"
            )

        raw_tau_gap = cecc_cfg.get("tau_gap")
        if raw_tau_gap is None:
            raw_tau_gap = cecc_cfg.get("edge_tolerance")
        if raw_tau_gap is None:
            raw_tau_gap = 20
        try:
            tau_gap = int(raw_tau_gap)
        except (TypeError, ValueError):
            tau_gap = 20
        if tau_gap < 1:
            pipeline.logger.warning(
                "cecc_build tau_gap must be >= 1; using 1 instead of %s",
                tau_gap,
            )
            tau_gap = 1
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

        fasta_candidate = pipeline.config.output_dir / f"{pipeline.config.prefix}_tandem_to_ring.fasta"
        if not fasta_candidate.exists():
            fasta_candidate = pipeline.config.output_dir / f"{pipeline.config.prefix}_circular.fasta"
        fasta_file: Optional[Path] = fasta_candidate if fasta_candidate.exists() else None

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
            half_query_buffer=int(cecc_cfg.get("half_query_buffer", 50)),
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
            uecc_path = pipeline.config.output_dir / f"{pipeline.config.prefix}_um_classify.uecc.csv"
            mecc_path = pipeline.config.output_dir / f"{pipeline.config.prefix}_um_classify.mecc.csv"
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

        # ONT platform: filter out intra-chromosomal CeccDNA.
        # ONT's higher error rate causes LAST to generate spurious intra-chr
        # chimeric calls.  HiFi data is unaffected (high-accuracy LAST).
        if (
            not df_cecc.empty
            and "CeccClass" in df_cecc.columns
            and getattr(pipeline.config, "platform", "") == "ont"
        ):
            before_cnt = int(df_cecc["query_id"].nunique())
            df_cecc = df_cecc[df_cecc["CeccClass"] != "Cecc-IntraChr"].copy()
            after_cnt = int(df_cecc["query_id"].nunique())
            if after_cnt != before_cnt:
                pipeline.logger.info(
                    "ONT: removed intra-chromosomal CeccDNA: %d -> %d reads",
                    before_cnt, after_cnt,
                )
                df_cecc.to_csv(output_file, index=False)

        pipeline._set_result(ResultKeys.CECC_BUILD_OUTPUT, str(output_file))
        pipeline._set_result(
            ResultKeys.CECC_BUILD_COUNT,
            int(df_cecc["query_id"].nunique())
            if (not df_cecc.empty and "query_id" in df_cecc.columns)
            else 0,
        )
    except (OSError, ValueError, KeyError) as e:
        pipeline.logger.warning(f"cecc_build failed: {e}")
        pd.DataFrame().to_csv(output_file, index=False)
        pipeline._set_result(ResultKeys.CECC_BUILD_OUTPUT, str(output_file))
        pipeline._set_result(ResultKeys.CECC_BUILD_COUNT, 0)
        pipeline._set_result(ResultKeys.CECC_BUILD_FAILED, True)
        pipeline._set_result(ResultKeys.CECC_BUILD_ERROR, str(e))

    # ── Supplementary: ONT cross-chr split read CeccDNA detection ──
    # Only run for ONT platform. This detects additional CeccDNA by looking at
    # reads that split-align across chromosomes in the raw minimap2 PAF.
    if getattr(pipeline.config, "platform", "") == "ont":
        _run_crosschr_cecc_supplement(pipeline, output_file)


def _run_crosschr_cecc_supplement(pipeline: "Pipeline", cecc_build_csv: Path) -> None:
    """Supplementary CeccDNA detection via cross-chr split reads in raw ONT PAF.

    Aligns raw ONT reads to reference (minimap2), extracts cross-chromosome
    split alignments, clusters junctions, and appends novel CeccDNA to the
    cecc_build output CSV.

    Only runs for ONT platform.
    """
    import pandas as pd
    import numpy as np

    log = pipeline.logger.getChild("crosschr_cecc")

    input_file = getattr(pipeline.config, "input_file", None)
    reference = getattr(pipeline.config, "reference", None)
    if not input_file or not reference:
        log.info("Cross-chr CeccDNA: skipped (no input_file or reference)")
        return

    input_file = Path(input_file)
    reference = Path(reference)
    if not input_file.exists() or not reference.exists():
        log.info("Cross-chr CeccDNA: skipped (input/reference not found)")
        return

    try:
        from circleseeker.modules.ont_cecc_crosschr import (
            parse_paf,
            extract_crosschr_junctions,
            detect_cecc_from_junctions,
        )
    except ImportError:
        log.warning("Cross-chr CeccDNA: module not available, skipped")
        return

    threads = int(getattr(pipeline.config, "threads", pipeline.config.performance.threads))
    output_dir = pipeline.config.output_dir
    prefix = pipeline.config.prefix

    # Step 1: Run minimap2 on raw reads (reuse if cached)
    paf_file = output_dir / f"{prefix}_crosschr_raw.paf"
    if not paf_file.exists() or paf_file.stat().st_size == 0:
        log.info("Cross-chr CeccDNA: aligning raw reads with minimap2...")
        ref_mmi = reference.parent / (reference.name + ".ont.mmi")
        ref_to_use = str(ref_mmi) if ref_mmi.exists() else str(reference)

        cmd = [
            "minimap2", "-cx", "map-ont",
            "--secondary=yes", "-N", "10", "-p", "0.5",
            "-t", str(threads),
            ref_to_use, str(input_file),
        ]
        log.info(f"Cross-chr CeccDNA: {' '.join(cmd)}")

        import subprocess as _sp
        with open(paf_file, "w") as fout:
            proc = _sp.run(cmd, stdout=fout, stderr=_sp.PIPE, text=True)
        if proc.returncode != 0:
            log.warning(f"Cross-chr CeccDNA: minimap2 failed: {proc.stderr[:300]}")
            return
        log.info(f"Cross-chr CeccDNA: minimap2 done ({paf_file})")
    else:
        log.info(f"Cross-chr CeccDNA: reusing cached PAF ({paf_file})")

    # Step 2: Parse PAF and extract cross-chr junctions
    reads = parse_paf(paf_file)
    log.info(f"Cross-chr CeccDNA: {len(reads)} reads parsed from PAF")

    junctions, n_crosschr = extract_crosschr_junctions(reads, min_mapq=1, min_align_len=200)
    log.info(f"Cross-chr CeccDNA: {n_crosschr} cross-chr reads, {len(junctions)} junctions")

    if not junctions:
        log.info("Cross-chr CeccDNA: no cross-chr junctions found")
        return

    # Step 3: Detect CeccDNA candidates
    # Use min_junction_reads=1 for union-find connectivity,
    # min_total_reads=2 at component level
    results = detect_cecc_from_junctions(
        junctions,
        min_junction_reads=1,
        min_total_reads=2,
        min_mean_mapq=20.0,
        cluster_dist=1000,
        padding=500,
        logger=log,
    )

    if not results:
        log.info("Cross-chr CeccDNA: no candidates after filtering")
        return

    log.info(f"Cross-chr CeccDNA: {len(results)} candidates detected")

    # Step 4: Remove duplicates with existing CeccBuild results
    existing_df = pd.DataFrame()
    if cecc_build_csv.exists() and cecc_build_csv.stat().st_size > 1:
        try:
            existing_df = pd.read_csv(cecc_build_csv)
        except Exception:
            existing_df = pd.DataFrame()

    # Build existing CeccBuild location signatures for dedup
    existing_sigs = set()
    if not existing_df.empty and "query_id" in existing_df.columns:
        for qid, grp in existing_df.groupby("query_id"):
            if "chr" in grp.columns and "start0" in grp.columns and "end0" in grp.columns:
                frags = []
                for _, row in grp.iterrows():
                    frags.append((str(row["chr"]), int(row["start0"]), int(row["end0"])))
                frags.sort()
                existing_sigs.add(tuple(frags))

    # Filter out cross-chr candidates that overlap with existing CeccBuild
    novel_results = []
    for r in results:
        frags = []
        for frag in r["fragments"].split(";"):
            chrom, coords = frag.split(":")
            start, end = coords.split("-")
            frags.append((chrom, int(start), int(end)))
        frags.sort()

        # Check overlap with existing
        is_dup = False
        for existing_frags in existing_sigs:
            if _fragments_overlap(frags, list(existing_frags), threshold=0.5):
                is_dup = True
                break

        if not is_dup:
            novel_results.append(r)

    if not novel_results:
        log.info("Cross-chr CeccDNA: all candidates overlap with existing CeccBuild")
        return

    log.info(f"Cross-chr CeccDNA: {len(novel_results)} novel candidates (not in CeccBuild)")

    # Step 5: Convert to cecc_build.csv format and append
    # Each CeccDNA has multiple segments (one row per segment)
    new_rows = []
    existing_qids = set(existing_df["query_id"].unique()) if not existing_df.empty and "query_id" in existing_df.columns else set()

    # Generate unique query_ids that won't conflict
    max_cecc_num = 0
    for qid in existing_qids:
        # Extract number from CeccDNA_NNNNNN format
        parts = str(qid).split("_")
        for p in parts:
            try:
                num = int(p.lstrip("0") or "0")
                max_cecc_num = max(max_cecc_num, num)
            except ValueError:
                continue

    for idx, r in enumerate(novel_results):
        cecc_num = max_cecc_num + idx + 1
        query_id = f"CrossChr_{cecc_num:06d}"

        frags = []
        for frag in r["fragments"].split(";"):
            chrom, coords = frag.split(":")
            start, end = coords.split("-")
            frags.append((chrom, int(start), int(end)))

        chroms_str = ",".join(sorted(set(c for c, _, _ in frags)))

        for seg_idx, (chrom, start, end) in enumerate(frags, 1):
            new_rows.append({
                "query_id": query_id,
                "reads": query_id,
                "eccdna_type": "Cecc",
                "CeccClass": "Cecc-InterChr",
                "length": r["length"],
                "copy_number": float(r["n_reads"]),
                "num_segments": r["n_fragments"],
                "num_distinct_loci": r["n_fragments"],
                "cumulative_length": r["length"],
                "match_degree": 95.0,
                "mapq_best": int(r["mean_junction_mapq"]),
                "mapq_min": int(r["mean_junction_mapq"]),
                "identity_best": 90.0,
                "identity_min": 85.0,
                "query_cov_best": 0.9,
                "query_cov_2nd": 0.0,
                "confidence_score": r["score"],
                "low_mapq": r["mean_junction_mapq"] < 5,
                "low_identity": False,
                "C_cov_best": 0.9,
                "C_cov_2nd": 0.0,
                "max_gap": 0,
                "avg_gap": 0.0,
                "chromosomes": chroms_str,
                "strand_closure_valid": bool(r["strand_closure"]),
                "segment_in_circle": seg_idx,
                "chr": chrom,
                "start0": start,
                "end0": end,
                "strand": "+",
                "q_start": 0,
                "q_end": end - start,
                "alignment_length": end - start,
            })

    new_df = pd.DataFrame(new_rows)

    # Append to existing cecc_build output
    if not existing_df.empty:
        combined_df = pd.concat([existing_df, new_df], ignore_index=True)
    else:
        combined_df = new_df

    combined_df.to_csv(cecc_build_csv, index=False)

    # Update result counts
    total_count = int(combined_df["query_id"].nunique()) if "query_id" in combined_df.columns else 0
    pipeline._set_result(ResultKeys.CECC_BUILD_COUNT, total_count)

    log.info(
        f"Cross-chr CeccDNA: appended {len(novel_results)} CeccDNA "
        f"({len(new_rows)} rows) to cecc_build output. "
        f"Total CeccDNA now: {total_count}"
    )


def _fragments_overlap(frags_a: list, frags_b: list, threshold: float = 0.5) -> bool:
    """Check if two fragment lists significantly overlap."""
    matched = 0
    for ca, sa, ea in frags_a:
        for cb, sb, eb in frags_b:
            if ca != cb:
                continue
            ovl = max(0, min(ea, eb) - max(sa, sb))
            min_len = min(ea - sa, eb - sb)
            if min_len > 0 and ovl / min_len >= threshold:
                matched += 1
                break
    # At least half of the fragments must overlap
    return matched >= max(1, min(len(frags_a), len(frags_b)) // 2)


def umc_process(pipeline: Pipeline) -> None:
    """Step 6: Generate FASTA files using umc_process module."""
    from circleseeker.modules.umc_process import UMCProcess, UMCProcessConfig

    circular_fasta = pipeline.config.output_dir / f"{pipeline.config.prefix}_tandem_to_ring.fasta"
    if not circular_fasta.exists():
        circular_fasta = pipeline.config.output_dir / f"{pipeline.config.prefix}_circular.fasta"
        if not circular_fasta.exists():
            pipeline.logger.warning("Circular FASTA file not found, skipping umc_process step")
            return

    uecc_csv = pipeline.config.output_dir / f"{pipeline.config.prefix}_um_classify.uecc.csv"
    mecc_csv = pipeline.config.output_dir / f"{pipeline.config.prefix}_um_classify.mecc.csv"
    cecc_csv = pipeline.config.output_dir / f"{pipeline.config.prefix}_cecc_build.csv"

    umc_cfg = UMCProcessConfig(process_xecc=bool(pipeline.config.enable_xecc))
    processor = UMCProcess(config=umc_cfg, logger=pipeline.logger.getChild("umc_process"))

    pipeline.logger.debug("Running umc_process module")
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

        pipeline._set_result(f"{ecc_type}_fasta", str(fasta_file))

        processed_csv = pipeline.config.output_dir / f"{pipeline.config.prefix}_{ecc_name}DNA_processed.csv"
        if not processed_csv.exists():
            fallback_csv = fasta_file.with_suffix(".csv")
            processed_csv = fallback_csv if fallback_csv.exists() else processed_csv

        if processed_csv.exists():
            pipeline._set_result(f"{ecc_type}_processed", str(processed_csv))

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
            clstr_path = cd_hit_tool.cluster_sequences(input_fasta, output_path, log_prefix=pipeline.config.prefix)
        except (subprocess.CalledProcessError, OSError) as e:
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
            pipeline._set_result(f"{ecc_type}_nr99", str(output_fasta))
        if cluster_file and Path(cluster_file).exists():
            pipeline._set_result(f"{ecc_type}_clusters", str(cluster_file))


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
    """Rename EccDedup output files to standard names.

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
    from circleseeker.modules.ecc_dedup import EccDedup, organize_umc_files

    harmonizer = EccDedup(pipeline.logger.getChild("ecc_dedup"))
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
                    pipeline._set_result(f"ecc_dedup_{mapped}_dir", str(directory))
        except (OSError, KeyError, ValueError) as exc:
            pipeline.logger.warning(f"Failed to organize EccDedup outputs: {exc}")
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
            pipeline.logger.warning(f"EccDedup output not found: {filename}")

    for ecc_type in ["Uecc", "Mecc", "Cecc"]:
        harmonized_csv = pipeline.config.output_dir / f"{pipeline.config.prefix}_{ecc_type}_harmonized.csv"
        if harmonized_csv.exists():
            key = {
                "Uecc": ResultKeys.UECC_HARMONIZED,
                "Mecc": ResultKeys.MECC_HARMONIZED,
                "Cecc": ResultKeys.CECC_HARMONIZED,
            }[ecc_type]
            pipeline._set_result(key, str(harmonized_csv))

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
                pipeline._set_result(key, str(file_path))
                pipeline.logger.debug(f"Found {key}: {file_path}")

    confirmed_file = pipeline.config.output_dir / f"{pipeline.config.prefix}_eccDNA_Confirmed.csv"
    if confirmed_file.exists():
        pipeline._set_result(
            ResultKeys.CONFIRMED_CSV,
            pipeline._serialize_path_for_state(confirmed_file),
        )
