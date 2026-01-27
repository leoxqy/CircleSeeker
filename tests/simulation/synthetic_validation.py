from __future__ import annotations

import random
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd

import sys

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.modules.cecc_build import CeccBuild, LastAlignment
from circleseeker.modules.um_classify import UMeccClassifier

from .eccdna_simulator import (
    CeccDNASimulator,
    EccDNARecord,
    MeccDNASimulator,
    ReferenceGenomeGenerator,
    SimulationConfig,
    SimulationOutputWriter,
    UeccDNASimulator,
)
from .validation_metrics import ValidationCalculator, ValidationMetrics


@dataclass(frozen=True)
class OverallMetrics:
    recall: float
    precision: float
    f1: float
    tp: int
    fp: int
    fn: int
    gt: int
    pred: int


def _overall_from_by_type(by_type: Dict[str, ValidationMetrics]) -> OverallMetrics:
    total_tp = sum(m.true_positives for m in by_type.values())
    total_fp = sum(m.false_positives for m in by_type.values())
    total_fn = sum(m.false_negatives for m in by_type.values())
    total_gt = sum(m.total_ground_truth for m in by_type.values())
    total_pred = sum(m.total_predictions for m in by_type.values())

    recall = total_tp / total_gt if total_gt else 0.0
    precision = total_tp / total_pred if total_pred else 0.0
    f1 = (2 * precision * recall / (precision + recall)) if (precision + recall) else 0.0
    return OverallMetrics(
        recall=recall,
        precision=precision,
        f1=f1,
        tp=total_tp,
        fp=total_fp,
        fn=total_fn,
        gt=total_gt,
        pred=total_pred,
    )


def generate_simulation_dir(
    base_dir: Path,
    *,
    seed: int,
    num_uecc: int,
    num_mecc: int,
    num_cecc: int,
) -> Tuple[SimulationConfig, List[EccDNARecord], Path]:
    """Generate a small, deterministic simulation dataset under base_dir."""
    sim_dir = base_dir / "simulation"
    sim_dir.mkdir(parents=True, exist_ok=True)

    config = SimulationConfig(
        output_dir=str(sim_dir),
        seed=seed,
        error_rate=0.0,
        num_chromosomes=3,
        chr_length=20_000,
        num_uecc=num_uecc,
        num_mecc=num_mecc,
        num_cecc=num_cecc,
        uecc_length_range=(200, 800),
        mecc_unit_range=(120, 120),
        mecc_copy_range=(2.0, 2.0),  # integer copy number for determinism
        mecc_repeat_mutation_rate=0.0,
        cecc_segments_range=(2, 3),
        cecc_segment_length_range=(120, 260),
        cecc_inter_chr_ratio=0.5,
    )

    random.seed(config.seed)
    genome_gen = ReferenceGenomeGenerator(config)
    genome = genome_gen.generate()

    uecc_sim = UeccDNASimulator(
        genome,
        config,
        excluded_intervals_by_chrom=genome_gen.repeat_intervals_by_chrom,
    )
    mecc_sim = MeccDNASimulator(genome, config, repeat_families=genome_gen.repeat_families)
    cecc_sim = CeccDNASimulator(
        genome,
        config,
        excluded_intervals_by_chrom=genome_gen.repeat_intervals_by_chrom,
    )

    uecc_records = uecc_sim.generate(config.num_uecc)
    mecc_records = mecc_sim.generate(config.num_mecc)
    cecc_records = cecc_sim.generate(config.num_cecc)

    all_records = uecc_records + mecc_records + cecc_records
    writer = SimulationOutputWriter(sim_dir)
    writer.write_ground_truth(all_records, "ground_truth.csv")
    writer.write_ground_truth(uecc_records, "ground_truth_uecc.csv")
    writer.write_ground_truth(mecc_records, "ground_truth_mecc.csv")
    writer.write_ground_truth(cecc_records, "ground_truth_cecc.csv")

    return config, all_records, sim_dir


def write_alignment_tsv_for_records(
    records: List[EccDNARecord],
    output_path: Path,
    *,
    identity: float = 99.0,
    mapq: int = 60,
) -> Path:
    """Create a synthetic BLAST-like TSV (plus trailing MAPQ) from ground truth records."""
    rows: List[List[object]] = []

    for record in records:
        query_id = record.read_id

        if record.eccdna_type == "Uecc":
            region = record.regions[0]
            cons_len = region.length
            q_start, q_end = 1, cons_len
            if region.strand == "+":
                s_start, s_end, sstrand = region.start + 1, region.end, "plus"
            else:
                s_start, s_end, sstrand = region.end, region.start + 1, "minus"

            rows.append(
                [
                    query_id,
                    region.chrom,
                    identity,
                    cons_len,
                    0,
                    0,
                    q_start,
                    q_end,
                    s_start,
                    s_end,
                    0.0,
                    200.0,
                    sstrand,
                    mapq,
                ]
            )
            continue

        if record.eccdna_type == "Mecc":
            # Mecc uses a consensus unit length encoded in query_id.
            unit_len = record.regions[0].length
            for i, region in enumerate(record.regions):
                offset = i * unit_len
                q_start = offset + 1
                q_end = offset + unit_len
                if region.strand == "+":
                    s_start, s_end, sstrand = region.start + 1, region.end, "plus"
                else:
                    s_start, s_end, sstrand = region.end, region.start + 1, "minus"
                rows.append(
                    [
                        query_id,
                        region.chrom,
                        identity,
                        unit_len,
                        0,
                        0,
                        q_start,
                        q_end,
                        s_start,
                        s_end,
                        0.0,
                        200.0,
                        sstrand,
                        mapq,
                    ]
                )
            continue

        if record.eccdna_type == "Cecc":
            offset = 0
            for region in record.regions:
                seg_len = region.length
                q_start = offset + 1
                q_end = offset + seg_len
                offset += seg_len
                if region.strand == "+":
                    s_start, s_end, sstrand = region.start + 1, region.end, "plus"
                else:
                    s_start, s_end, sstrand = region.end, region.start + 1, "minus"

                rows.append(
                    [
                        query_id,
                        region.chrom,
                        identity,
                        seg_len,
                        0,
                        0,
                        q_start,
                        q_end,
                        s_start,
                        s_end,
                        0.0,
                        200.0,
                        sstrand,
                        mapq,
                    ]
                )
            continue

    pd.DataFrame(rows).to_csv(output_path, sep="\t", header=False, index=False)
    return output_path


def write_negative_control_alignment_tsv(
    output_path: Path,
    *,
    n_reads: int,
    cons_len: int = 200,
    mapq: int = 0,
    identity: float = 99.0,
) -> Path:
    """Synthetic negative control: full-coverage alignments with low MAPQ (ambiguous/repetitive)."""
    rows: List[List[object]] = []
    for i in range(1, n_reads + 1):
        query_id = f"neg_{i:04d}|unit_1|{cons_len}|1.0"
        rows.append(
            [
                query_id,
                "chrR",
                identity,
                cons_len,
                0,
                0,
                1,
                cons_len,
                100,
                100 + cons_len - 1,
                0.0,
                200.0,
                "plus",
                mapq,
            ]
        )
    pd.DataFrame(rows).to_csv(output_path, sep="\t", header=False, index=False)
    return output_path


def run_synthetic_validation(
    work_dir: Path,
    *,
    seed: int,
    num_uecc: int,
    num_mecc: int,
    num_cecc: int,
    um_kwargs: Dict | None = None,
    cecc_kwargs: Dict | None = None,
) -> Tuple[Dict[str, ValidationMetrics], OverallMetrics]:
    """Run a synthetic end-to-end validation (no external aligner) and return metrics."""
    _, records, sim_dir = generate_simulation_dir(
        work_dir,
        seed=seed,
        num_uecc=num_uecc,
        num_mecc=num_mecc,
        num_cecc=num_cecc,
    )

    results_dir = work_dir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    alignment_tsv = work_dir / "alignment.tsv"
    write_alignment_tsv_for_records(records, alignment_tsv)

    classifier = UMeccClassifier(**(um_kwargs or {}))
    uecc_df, mecc_df, unclassified_df = classifier.run(alignment_tsv)

    def _write_pred(df: pd.DataFrame, path: Path) -> None:
        if df is None or df.empty:
            pd.DataFrame(
                columns=[
                    "eccDNA_id",
                    "eccdna_type",
                    "reads",
                    "chr",
                    "start0",
                    "end0",
                    "strand",
                    "length",
                    "copy_number",
                ]
            ).to_csv(path, index=False)
            return
        out = df.copy()
        if "eccDNA_id" not in out.columns:
            out["eccDNA_id"] = out.get("query_id", "")
        out.to_csv(path, index=False)

    _write_pred(uecc_df, results_dir / "uecc.csv")
    _write_pred(mecc_df, results_dir / "mecc.csv")

    # CeccBuild requires headers even for empty inputs.
    unclassified_csv = results_dir / "unclassified.csv"
    if unclassified_df is None or unclassified_df.empty:
        pd.DataFrame(
            columns=CeccBuild.BASE_REQUIRED_COLS + CeccBuild.STANDARD_EXTRA_COLS
        ).to_csv(unclassified_csv, index=False)
    else:
        unclassified_df.to_csv(unclassified_csv, index=False)

    cecc_csv = results_dir / "cecc.csv"
    builder = CeccBuild()

    def _double_cecc_alignments(df: pd.DataFrame) -> dict[str, list[LastAlignment]]:
        """Build LAST-like alignments with doubled-sequence query coordinates.

        CeccBuild's graph-based detector expects the query to be a doubled sequence
        (seq+seq) so the locus pattern repeats. For simulation we don't run LAST;
        instead, we expand the synthetic BLAST-like alignments into a doubled query.
        """
        if df.empty:
            return {}

        alignments: dict[str, list[LastAlignment]] = {}
        for query_id, group in df.groupby("query_id"):
            q_end = pd.to_numeric(group.get("q_end"), errors="coerce").dropna()
            if q_end.empty:
                continue
            cons_len = int(q_end.max())
            if cons_len <= 0:
                continue

            query_len = cons_len * 2
            for _, row in group.iterrows():
                try:
                    chrom = str(row.get("chr", row.get("subject_id", "")))
                    ref_start = int(row.get("start0", row.get("s_start", 0)))
                    ref_end = int(row.get("end0", row.get("s_end", 0)))
                    # Synthetic alignments use BLAST-like query coords: 1-based start, end inclusive.
                    # Convert to 0-based half-open for CeccBuild's internal logic.
                    q_start0 = int(row.get("q_start", 1)) - 1
                    q_end0 = int(row.get("q_end", 0))
                    strand = str(row.get("strand", "+"))
                    identity = float(row.get("identity", 99.0))
                except (TypeError, ValueError):
                    continue

                base = LastAlignment(
                    chrom=chrom,
                    ref_start=ref_start,
                    ref_end=ref_end,
                    query_start=q_start0,
                    query_end=q_end0,
                    query_len=query_len,
                    score=0,
                    identity=identity,
                    strand=strand,
                )
                dup = LastAlignment(
                    chrom=chrom,
                    ref_start=ref_start,
                    ref_end=ref_end,
                    query_start=q_start0 + cons_len,
                    query_end=q_end0 + cons_len,
                    query_len=query_len,
                    score=0,
                    identity=identity,
                    strand=strand,
                )
                alignments.setdefault(str(query_id), []).extend([base, dup])

        return alignments

    # Synthetic Cecc validation: avoid requiring LAST by using the graph-based detector
    # directly on expanded (doubled) synthetic alignments.
    if unclassified_df is None or unclassified_df.empty:
        cecc_df = pd.DataFrame()
        _write_pred(cecc_df, cecc_csv)
    else:
        cecc_alignments = _double_cecc_alignments(unclassified_df)
        cecc_metadata = builder._get_metadata_from_csv(unclassified_df)
        try:
            cecc_df = builder.detect_cecc_with_closure_and_loci(
                cecc_alignments,
                cecc_metadata,
                min_distinct_loci=int((cecc_kwargs or {}).get("min_segments", 2)),
            )
        except Exception:
            cecc_df = pd.DataFrame()
        _write_pred(cecc_df, cecc_csv)

    if not cecc_df.empty and "eccDNA_id" not in cecc_df.columns:
        cecc_df = cecc_df.copy()
        cecc_df["eccDNA_id"] = cecc_df.get("query_id", "")
        cecc_df.to_csv(cecc_csv, index=False)

    calc = ValidationCalculator()
    by_type: Dict[str, ValidationMetrics] = {}

    gt_uecc = calc.load_ground_truth(sim_dir / "ground_truth_uecc.csv")
    pred_uecc = calc.load_predictions(results_dir / "uecc.csv", "Uecc")
    by_type["Uecc"] = calc.match_records(gt_uecc, pred_uecc)

    gt_mecc = calc.load_ground_truth(sim_dir / "ground_truth_mecc.csv")
    pred_mecc = calc.load_predictions(results_dir / "mecc.csv", "Mecc")
    by_type["Mecc"] = calc.match_records(gt_mecc, pred_mecc)

    gt_cecc = calc.load_ground_truth(sim_dir / "ground_truth_cecc.csv")
    pred_cecc = calc.load_predictions(results_dir / "cecc.csv", "Cecc")
    by_type["Cecc"] = calc.match_records(gt_cecc, pred_cecc)

    return by_type, _overall_from_by_type(by_type)
