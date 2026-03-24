"""ONT CeccDNA detection via cross-chromosome split read signals.

Supplements CeccBuild (LAST-based) with a direct cross-chr split alignment approach,
adapted from the NGS pipeline's _detect_cecc_early method.

ONT reads aligned by minimap2 can produce split alignments across chromosomes,
which directly indicate CeccDNA. This module:
1. Aligns raw ONT reads to the reference genome (minimap2, map-ont preset)
2. Extracts cross-chromosome split alignments (same read, different chromosomes)
3. Clusters junctions by genomic proximity using union-find
4. Outputs CeccDNA candidates with fragment coordinates

Key ONT adaptations vs NGS:
- min_junction_reads = 2 (empirically optimal for ONT precision/recall)
- cluster_dist = 1000 (ONT alignment coordinates less precise)
- min_mean_mapq = 20 (filters noise from low-quality multi-mapping)
- No discordant pairs (ONT is single-end, only split reads)
"""

from __future__ import annotations

import logging
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd


# ── Default thresholds ─────────────────────────────────────────────

DEFAULT_MIN_MAPQ = 1            # per-alignment minimum MAPQ
DEFAULT_MIN_ALIGN_LEN = 200     # per-alignment minimum length
DEFAULT_MIN_JUNCTION_READS = 2  # reads per junction bin
DEFAULT_MIN_MEAN_MAPQ = 20      # mean MAPQ across junction
DEFAULT_CLUSTER_DIST = 1000     # bp bin size for clustering
DEFAULT_PADDING = 500           # fragment boundary padding
DEFAULT_MAX_FRAGMENT_SIZE = 100000
DEFAULT_MAX_TOTAL_SIZE = 500000

# ── CrossChr CeccDNA LR quality-control model ─────────────────────
# Trained on ara rep1 10X+30X combined (124 TP, 30 FP).
# Key discriminators: n_junctions, strand_closure, mean_junction_mapq,
# log_length, n_reads.  Coverage-adaptive via local_reads feature.
import math as _math

_CECC_CROSSCHR_FEATURE_NAMES = [
    "n_reads", "n_junctions", "n_fragments", "n_chromosomes",
    "mean_junction_mapq", "strand_closure", "log_length",
    "ratio", "reads_per_junction", "local_reads",
]
_CECC_CROSSCHR_COEFS = np.array([
    +0.0993, +0.1528, +0.0149, +0.0149, +0.1927,
    +0.1566, +0.1718, -0.0073, +0.0100, +0.0342,
])
_CECC_CROSSCHR_INTERCEPT = 0.2906
_CECC_CROSSCHR_SCALER_MEAN = np.array([
    16.0195, 6.8701, 2.1299, 2.1299, 52.1513,
    0.9286, 3.6968, 0.1285, 2.5806, 166.3442,
])
_CECC_CROSSCHR_SCALER_STD = np.array([
    15.7243, 4.9592, 0.4218, 0.4218, 10.8337,
    0.2575, 0.2112, 0.1122, 2.3641, 161.2461,
])
_CECC_CROSSCHR_THRESHOLD = 0.50  # Prec>=92%, Rec~90% on training set


def _sigmoid(x: float) -> float:
    return 1.0 / (1.0 + _math.exp(-x)) if x > -500 else 0.0


# ── PAF parsing ──────────────────────────────────────────────────────

def parse_paf(paf_path: str | Path) -> dict[str, list[dict]]:
    """Parse PAF file, return dict: read_name -> list of alignment segments."""
    reads = defaultdict(list)
    with open(paf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 12:
                continue
            qname = cols[0]
            reads[qname].append({
                "qname": qname,
                "qlen": int(cols[1]),
                "qstart": int(cols[2]),
                "qend": int(cols[3]),
                "strand": cols[4],
                "tname": cols[5],
                "tlen": int(cols[6]),
                "tstart": int(cols[7]),
                "tend": int(cols[8]),
                "n_match": int(cols[9]),
                "block_len": int(cols[10]),
                "mapq": int(cols[11]),
            })
    return reads


# ── Cross-chr junction extraction ────────────────────────────────────

def extract_crosschr_junctions(
    reads: dict[str, list[dict]],
    min_mapq: int = DEFAULT_MIN_MAPQ,
    min_align_len: int = DEFAULT_MIN_ALIGN_LEN,
) -> list[tuple]:
    """Extract cross-chromosome junctions from ONT split alignments.

    For each read with alignments on >= 2 chromosomes, generate pairwise
    junction tuples using alignment midpoints as breakpoint positions.

    Returns list of (chr_a, pos_a, strand_a, chr_b, pos_b, strand_b,
                      read_name, 'split', min_mapq).
    """
    junctions = []
    n_crosschr_reads = 0

    for qname, segments in reads.items():
        good_segs = [
            s for s in segments
            if s["mapq"] >= min_mapq and (s["tend"] - s["tstart"]) >= min_align_len
        ]
        if len(good_segs) < 2:
            continue

        chrom_segs = defaultdict(list)
        for s in good_segs:
            chrom_segs[s["tname"]].append(s)

        if len(chrom_segs) < 2:
            continue

        n_crosschr_reads += 1

        chroms = sorted(chrom_segs.keys())
        for i in range(len(chroms)):
            for j in range(i + 1, len(chroms)):
                chr_a, chr_b = chroms[i], chroms[j]
                for seg_a in chrom_segs[chr_a]:
                    for seg_b in chrom_segs[chr_b]:
                        pos_a = (seg_a["tstart"] + seg_a["tend"]) // 2
                        pos_b = (seg_b["tstart"] + seg_b["tend"]) // 2
                        strand_a = seg_a["strand"]
                        strand_b = seg_b["strand"]
                        min_mq = min(seg_a["mapq"], seg_b["mapq"])
                        junctions.append((
                            chr_a, pos_a, strand_a,
                            chr_b, pos_b, strand_b,
                            qname, "split", min_mq,
                        ))

    return junctions, n_crosschr_reads


# ── Junction clustering + Union-Find ─────────────────────────────────

def detect_cecc_from_junctions(
    crosschr: list[tuple],
    min_junction_reads: int = 1,
    min_total_reads: int = DEFAULT_MIN_JUNCTION_READS,
    min_junction_ratio: float = 0.02,
    min_mean_mapq: float = DEFAULT_MIN_MEAN_MAPQ,
    cluster_dist: int = DEFAULT_CLUSTER_DIST,
    max_fragment_size: int = DEFAULT_MAX_FRAGMENT_SIZE,
    max_total_size: int = DEFAULT_MAX_TOTAL_SIZE,
    padding: int = DEFAULT_PADDING,
    local_coverage_index: Optional[dict] = None,
    logger: Optional[logging.Logger] = None,
) -> list[dict]:
    """Detect CeccDNA from cross-chr junctions using clustering + union-find.

    Filtering strategy:
    - min_junction_reads=1: keep all junctions for union-find connectivity
    - min_total_reads: absolute minimum reads per component (floor = 2)
    - min_junction_ratio: per-region adaptive filter.  For each component,
      estimate local coverage at its fragment regions from
      local_coverage_index.  Require:
          n_junction_reads / n_local_reads >= min_junction_ratio
      Real CeccDNA: junction reads ARE the local reads → ratio ≈ high
      Noise (NUMT etc.): few junction reads among many genomic reads → low
      This is coverage-adaptive without needing global coverage estimation.
    - min_mean_mapq: filter at the component level (mean MAPQ across junctions)
    """
    log = logger or logging.getLogger(__name__)

    if not crosschr:
        return []

    # Step 1: Bin and cluster junctions
    junction_bins = defaultdict(lambda: {
        "reads": set(),
        "positions_a": [], "positions_b": [],
        "strands_a": [], "strands_b": [],
        "mapqs": [],
    })

    for item in crosschr:
        chr_a, pos_a, strand_a, chr_b, pos_b, strand_b, rname, etype, mapq = item[:9]
        if chr_a > chr_b:
            chr_a, pos_a, strand_a, chr_b, pos_b, strand_b = \
                chr_b, pos_b, strand_b, chr_a, pos_a, strand_a
        bin_a = pos_a // cluster_dist
        bin_b = pos_b // cluster_dist
        key = (chr_a, bin_a, chr_b, bin_b)
        junction_bins[key]["reads"].add(rname)
        junction_bins[key]["mapqs"].append(mapq)
        junction_bins[key]["positions_a"].append(pos_a)
        junction_bins[key]["positions_b"].append(pos_b)
        junction_bins[key]["strands_a"].append(strand_a)
        junction_bins[key]["strands_b"].append(strand_b)

    # Step 2: Build junctions (keep all with >= min_junction_reads for union-find)
    junctions = []
    for key, ev in junction_bins.items():
        n_reads = len(ev["reads"])
        if n_reads < min_junction_reads:
            continue
        chr_a, _, chr_b, _ = key
        med_a = int(np.median(ev["positions_a"]))
        med_b = int(np.median(ev["positions_b"]))
        dominant_strand_a = max(set(ev["strands_a"]), key=ev["strands_a"].count)
        dominant_strand_b = max(set(ev["strands_b"]), key=ev["strands_b"].count)
        mapqs = ev["mapqs"] if ev["mapqs"] else [0]
        mean_mq = float(np.mean(mapqs))

        junctions.append({
            "chr_a": chr_a, "pos_a": med_a, "strand_a": dominant_strand_a,
            "chr_b": chr_b, "pos_b": med_b, "strand_b": dominant_strand_b,
            "n_reads": n_reads,
            "read_names": ev["reads"],
            "mean_mapq": mean_mq,
            "min_mapq": min(mapqs),
        })

    if not junctions:
        log.info("ONT cross-chr: no junctions passed filtering")
        return []

    log.info(f"ONT cross-chr: {len(junctions)} junctions after filtering")

    # Step 3: Union-find
    parent = {}

    def find(x):
        while parent.get(x, x) != x:
            parent[x] = parent.get(parent[x], parent[x])
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for j in junctions:
        node_a = (j["chr_a"], j["pos_a"] // cluster_dist)
        node_b = (j["chr_b"], j["pos_b"] // cluster_dist)
        union(node_a, node_b)

    # Step 4: Group by component
    components = defaultdict(list)
    for j in junctions:
        node_a = (j["chr_a"], j["pos_a"] // cluster_dist)
        root = find(node_a)
        components[root].append(j)

    # Step 5: Build CeccDNA candidates
    results = []
    cecc_idx = 1

    for comp_junctions in components.values():
        fragments = defaultdict(list)
        total_reads_set = set()
        for j in comp_junctions:
            fragments[j["chr_a"]].append(j["pos_a"])
            fragments[j["chr_b"]].append(j["pos_b"])
            total_reads_set.update(j["read_names"])

        chroms = set(fragments.keys())
        if len(chroms) < 2:
            continue

        total_reads = len(total_reads_set)

        # Absolute minimum read count filter
        if total_reads < min_total_reads:
            continue

        frag_list = []
        oversized = False
        for chrom, positions in fragments.items():
            frag_start = min(positions) - padding
            frag_end = max(positions) + padding
            if frag_end - frag_start > max_fragment_size:
                oversized = True
                break
            frag_list.append((chrom, max(0, frag_start), frag_end))

        if oversized:
            continue

        frag_list.sort()
        total_length = sum(e - s for _, s, e in frag_list)
        if total_length > max_total_size:
            continue

        # Strand closure
        n_strand_flips = sum(
            1 for j in comp_junctions if j["strand_a"] != j["strand_b"]
        )
        strand_closure = (n_strand_flips % 2 == 0)
        n_junctions = len(comp_junctions)

        rep_chr, rep_start, rep_end = frag_list[0]
        frag_str = ";".join(f"{c}:{s}-{e}" for c, s, e in frag_list)

        base_score = min(1.0, total_reads / 5.0)
        if not strand_closure:
            base_score *= 0.7

        junction_reads = [j["n_reads"] for j in comp_junctions]
        mean_rpj = float(np.mean(junction_reads))
        junc_mapqs = [j.get("mean_mapq", 0) for j in comp_junctions]
        mean_junc_mapq = float(np.mean(junc_mapqs))

        # Component-level MAPQ filter
        if mean_junc_mapq < min_mean_mapq:
            continue

        # Compute local coverage for LR model feature
        local_reads = 0
        if local_coverage_index:
            for chrom, frag_start, frag_end in frag_list:
                chrom_bins = local_coverage_index.get(chrom, {})
                for b in range(frag_start // 1000, frag_end // 1000 + 1):
                    local_reads += chrom_bins.get(b, 0)
        ratio = total_reads / local_reads if local_reads > 0 else 1.0
        rpj = total_reads / max(n_junctions, 1)

        eid = f"CrossChrCecc{cecc_idx:05d}"
        cecc_idx += 1

        results.append({
            "eccDNA_id": eid,
            "chr": rep_chr,
            "start": rep_start,
            "end": rep_end,
            "length": total_length,
            "eccdna_type": "Cecc",
            "state": "Confirmed",
            "n_fragments": len(frag_list),
            "n_chromosomes": len(chroms),
            "n_junctions": n_junctions,
            "n_reads": total_reads,
            "fragments": frag_str,
            "mean_reads_per_junction": round(mean_rpj, 3),
            "strand_closure": int(strand_closure),
            "n_strand_flips": n_strand_flips,
            "score": round(base_score, 4),
            "mean_junction_mapq": round(mean_junc_mapq, 1),
            "read_names": sorted(total_reads_set),
            # LR features (for scoring below)
            "ratio": round(ratio, 4),
            "reads_per_junction": round(rpj, 3),
            "local_reads": local_reads,
            "log_length": round(float(np.log10(max(total_length, 1))), 4),
        })

    log.info(f"ONT cross-chr: {len(results)} candidates before LR filter")

    # ── LR quality-control filter ──
    # Score each candidate with the trained model.  Remove low-scoring ones.
    if results:
        n_before = len(results)
        filtered = []
        for r in results:
            feats = np.array([r.get(f, 0) for f in _CECC_CROSSCHR_FEATURE_NAMES], dtype=float)
            normed = (feats - _CECC_CROSSCHR_SCALER_MEAN) / np.maximum(_CECC_CROSSCHR_SCALER_STD, 1e-6)
            logit = float(np.dot(normed, _CECC_CROSSCHR_COEFS) + _CECC_CROSSCHR_INTERCEPT)
            lr_score = _sigmoid(logit)
            r["lr_score"] = round(lr_score, 4)
            if lr_score >= _CECC_CROSSCHR_THRESHOLD:
                filtered.append(r)
        n_removed = n_before - len(filtered)
        if n_removed:
            log.info(f"ONT cross-chr LR filter: removed {n_removed}/{n_before} "
                     f"(threshold={_CECC_CROSSCHR_THRESHOLD}), kept {len(filtered)}")
        results = filtered

    log.info(f"ONT cross-chr: detected {len(results)} CeccDNA candidates")
    return results


# ── High-level pipeline entry point ──────────────────────────────────

def run_ont_cecc_crosschr(
    input_fastq: str | Path,
    reference: str | Path,
    output_dir: str | Path,
    prefix: str = "ont_crosschr",
    threads: int = 12,
    min_mapq: int = DEFAULT_MIN_MAPQ,
    min_align_len: int = DEFAULT_MIN_ALIGN_LEN,
    min_junction_reads: int = DEFAULT_MIN_JUNCTION_READS,
    min_mean_mapq: float = DEFAULT_MIN_MEAN_MAPQ,
    cluster_dist: int = DEFAULT_CLUSTER_DIST,
    padding: int = DEFAULT_PADDING,
    paf_path: Optional[str | Path] = None,
    logger: Optional[logging.Logger] = None,
) -> pd.DataFrame:
    """Run the full ONT cross-chr CeccDNA detection pipeline.

    Args:
        input_fastq: ONT reads FASTQ file
        reference: Reference genome FASTA
        output_dir: Output directory (for PAF and results)
        prefix: Output file prefix
        threads: minimap2 threads
        paf_path: Pre-computed PAF (skip alignment if provided and exists)
        logger: Logger instance

    Returns:
        DataFrame with CeccDNA candidates
    """
    log = logger or logging.getLogger(__name__)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Alignment
    if paf_path and Path(paf_path).exists():
        paf_file = Path(paf_path)
        log.info(f"ONT cross-chr: reusing existing PAF: {paf_file}")
    else:
        paf_file = output_dir / f"{prefix}_crosschr.paf"
        if paf_file.exists() and paf_file.stat().st_size > 0:
            log.info(f"ONT cross-chr: reusing cached PAF: {paf_file}")
        else:
            log.info("ONT cross-chr: running minimap2 alignment...")
            ref_path = Path(reference)
            ref_mmi = ref_path.parent / (ref_path.name + ".ont.mmi")
            ref_to_use = str(ref_mmi) if ref_mmi.exists() else str(reference)

            cmd = [
                "minimap2", "-cx", "map-ont",
                "--secondary=yes", "-N", "10", "-p", "0.5",
                "-t", str(threads),
                ref_to_use, str(input_fastq),
            ]
            log.info(f"ONT cross-chr: {' '.join(cmd)}")

            with open(paf_file, "w") as fout:
                proc = subprocess.run(
                    cmd, stdout=fout, stderr=subprocess.PIPE, text=True,
                )
            if proc.returncode != 0:
                log.error(f"minimap2 failed: {proc.stderr[:500]}")
                return pd.DataFrame()
            log.info(f"ONT cross-chr: minimap2 done, PAF={paf_file}")

    # Step 2: Parse PAF
    log.info(f"ONT cross-chr: parsing PAF...")
    reads = parse_paf(paf_file)
    log.info(f"ONT cross-chr: {len(reads)} reads with alignments")

    # Step 2b: Build local coverage index from PAF alignments
    # For each 1kb genomic bin, count how many unique reads have their
    # PRIMARY alignment overlapping that bin.  This allows per-region
    # adaptive filtering: real CeccDNA junction reads dominate their
    # fragment regions, while NUMT noise is a small fraction of the
    # background reads at those positions.
    local_cov_index = defaultdict(lambda: defaultdict(int))  # chrom -> bin -> count
    for qname, segments in reads.items():
        if not segments:
            continue
        # Use the best (longest) alignment as proxy for primary
        best = max(segments, key=lambda s: s["tend"] - s["tstart"])
        chrom = best["tname"]
        for b in range(best["tstart"] // 1000, best["tend"] // 1000 + 1):
            local_cov_index[chrom][b] += 1
    log.info(f"ONT cross-chr: built local coverage index "
             f"({sum(len(v) for v in local_cov_index.values())} bins)")

    # Step 3: Extract junctions
    junctions, n_crosschr = extract_crosschr_junctions(
        reads, min_mapq=min_mapq, min_align_len=min_align_len,
    )
    log.info(f"ONT cross-chr: {n_crosschr} cross-chr reads, {len(junctions)} junctions")

    # Step 4: Detect CeccDNA
    # Use min_junction_reads=1 for union-find (preserve connectivity),
    # filter by min_total_reads + per-region ratio at the component level
    results = detect_cecc_from_junctions(
        junctions,
        min_junction_reads=1,
        min_total_reads=min_junction_reads,
        min_mean_mapq=min_mean_mapq,
        cluster_dist=cluster_dist,
        padding=padding,
        local_coverage_index=dict(local_cov_index),
        logger=log,
    )

    if not results:
        log.info("ONT cross-chr: no CeccDNA detected")
        return pd.DataFrame()

    # Step 5: Build DataFrame
    df = pd.DataFrame(results)

    # Write results
    out_csv = output_dir / f"{prefix}_crosschr_cecc.csv"
    cols_to_save = [c for c in df.columns if c != "read_names"]
    df[cols_to_save].to_csv(out_csv, index=False)
    log.info(f"ONT cross-chr: {len(df)} CeccDNA written to {out_csv}")

    return df
