"""SplitReads-Core: Integrated split-read based eccDNA detection for HiFi data.

This module provides split-read based eccDNA detection, optimized for PacBio HiFi
sequencing data. It detects eccDNA by analyzing split-read alignment patterns
and breakpoint evidence using graph-based algorithms.

Algorithm inspired by CReSIL (https://github.com/Peppermint-Lab/CReSIL),
with significant optimizations for HiFi data:
- Uses map-hifi preset for accurate long-read alignment
- Adjusted mapq threshold for HiFi's high accuracy
- Removed ONT-specific consensus correction (not needed for HiFi)
- Streamlined for integration with CircleSeeker pipeline
"""

from __future__ import annotations

import datetime
import logging
import pathlib
import tempfile
from collections import Counter
from dataclasses import dataclass
from itertools import chain
from multiprocessing import Pool, cpu_count
from pathlib import Path
from typing import Any, Optional, cast

import intervaltree as tree
import mappy
import networkx as nx
import numpy as np
import pandas as pd
import pybedtools as bt
from Bio import SeqIO

from circleseeker.utils.logging import get_logger


# ============================================================================
# Configuration
# ============================================================================


@dataclass
class SplitReadsConfig:
    """Configuration for SplitReads-Core module."""

    # Alignment parameters
    preset: str = "map-hifi"  # Changed from map-ont
    mapq: int = 20  # Lowered from 30 for HiFi (higher accuracy allows lower threshold)
    exclude_chrs: str = ""

    # Trim parameters
    allow_gap: int = 10
    allow_overlap: int = 10

    # Identify parameters
    min_region_size: int = 200
    overlap_check_size: int = 50
    breakpoint_depth: int = 5  # Minimum breakpoint support
    average_depth: float = 5.0  # Minimum average coverage depth

    # Merge parameters
    ref_merge_distance: int = 1000  # Maximum reference distance (bp) for merging adjacent alignments

    # Processing
    threads: int = 0  # 0 = use all CPUs
    skip_variant: bool = True  # Always skip for HiFi (not needed)

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "SplitReadsConfig":
        """Create config from dictionary."""
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})


# ============================================================================
# Trim Module Functions
# ============================================================================


class Region:
    """Represents a genomic region."""

    def __init__(self, chrom: str, start: int, end: int):
        self.chrom = chrom
        self.start = start
        self.end = end

    def get_str_region(self) -> str:
        return f"{self.chrom}_{self.start}_{self.end}"


def is_overlapping(obj: Region, offset: int, chrom: str, start: int, end: int) -> bool:
    """Check if two regions overlap within offset tolerance."""
    if obj.chrom != chrom:
        return False
    return abs(obj.start - start) <= offset and abs(obj.end - end) <= offset


def any_overlapping_range(start1: int, end1: int, start2: int, end2: int) -> bool:
    """Check if two ranges overlap."""
    return start1 <= end2 and start2 <= end1


def combine_region_strand(value: pd.Series) -> str:
    """Combine region and strand into string."""
    return f"{value['ref']}_{value['r_start']}_{value['r_end']}_{value['strand']}"


def make_group_order(list_test: list[str], offset: int) -> list[str]:
    """Group regions by overlap and return ordered list."""
    list_group_order: list[str] = []
    list_group: list[Region] = []

    for region_strand in list_test:
        chrom, start, end, strand = region_strand.rsplit("_", 3)

        if len(list_group) == 0:
            obj_region = Region(chrom, int(start), int(end))
            list_group.append(obj_region)
            list_group_order.append(f"{chrom}_{start}_{end}_{strand}")
        else:
            found_flag = 0
            for obj in list_group:
                if is_overlapping(obj, offset, chrom, int(start), int(end)):
                    found_flag = 1
                    mod_region = f"{obj.get_str_region()}_{strand}"
                    list_group_order.append(mod_region)
                    break

            if found_flag == 0:
                obj_region = Region(chrom, int(start), int(end))
                list_group.append(obj_region)
                list_group_order.append(f"{chrom}_{start}_{end}_{strand}")

    return list_group_order


def split_list_step_window(list_check: list, step: int, window: int) -> list[list]:
    """Split list into windows with step."""
    return [list_check[i : i + window] for i in range(0, len(list_check), window - step)]


def get_sum_pattern(list_group_order: list[str]) -> list:
    """Get summary pattern from grouped order."""
    list_sum_pattern = []

    list_most_common = Counter(list_group_order).most_common()
    if list_most_common[0][1] > 1:
        if len(list_most_common) == 1:
            list_sum_pattern = [True, list_most_common[0][0], list_most_common[0][1]]
        else:
            idx_start = list_group_order.index(list_most_common[0][0])
            list_trim = list_group_order[idx_start:]

            list_trim_most_common = Counter(list_trim).most_common()
            if len(list_trim_most_common) == 1:
                list_sum_pattern = [True, list_trim_most_common[0][0], list_trim_most_common[0][1]]
            else:
                list_sum_pattern = [False, "", 0]

                for i in range(1, len(list_trim)):
                    list_split = split_list_step_window(list_trim, 0, i)
                    list_region = [",".join(x) for x in list_split]

                    list_rev_split_most_common = Counter(list_region).most_common()[::-1]

                    for tup_ in list_rev_split_most_common:
                        if tup_[1] > 1:
                            list_sum_pattern = [True, tup_[0], tup_[1]]
    else:
        list_sum_pattern = [False, "", 0]

    return list_sum_pattern


def get_idx_longest_pattern(list_group_order: list[str], list_pattern_region: list[str]) -> list[int]:
    """Get indices of longest pattern match."""
    len_list_group_order = len(list_group_order)
    len_list_pattern_region = len(list_pattern_region)

    list_idx_pattern: list[list[int]] = []
    idx = 0
    while True:
        if idx + len_list_pattern_region > len_list_group_order:
            break
        idx_end = idx + len_list_pattern_region
        if list_group_order[idx:idx_end] == list_pattern_region:
            list_idx_pattern.append(list(range(idx, idx_end)))
            idx = idx_end
        else:
            idx += 1

    list_merge_idx_pattern: list[list[int]] = []

    for list_ in list_idx_pattern:
        if len(list_merge_idx_pattern) > 0:
            if list_merge_idx_pattern[-1][-1] + 1 == list_[0]:
                list_merge_idx_pattern[-1] += list_
            else:
                list_merge_idx_pattern.append(list_)
        else:
            list_merge_idx_pattern.append(list_)

    return max(list_merge_idx_pattern, key=len) if list_merge_idx_pattern else []


class pdRegion:
    """Represents a mapped region from pandas row."""

    def __init__(self, pd_object: pd.Series):
        self.readid = pd_object["readid"]
        self.q_len = pd_object["q_len"]
        self.q_start = pd_object["q_start"]
        self.q_end = pd_object["q_end"]
        self.ref = pd_object["ref"]
        self.r_start = pd_object["r_start"]
        self.r_end = pd_object["r_end"]
        self.matchLen = pd_object["matchLen"]
        self.blockLen = pd_object["blockLen"]
        self.mapq = pd_object["mapq"]
        self.strand = pd_object["strand"]

    def cal_q_len(self) -> int:
        return int(self.q_end) - int(self.q_start)

    def cal_r_len(self) -> int:
        return int(self.r_end) - int(self.r_start)

    def get_list_attr(self) -> list:
        return [
            self.readid,
            self.q_len,
            self.q_start,
            self.q_end,
            self.ref,
            self.r_start,
            self.r_end,
            self.matchLen,
            self.blockLen,
            self.mapq,
            self.strand,
        ]


# Global variables for multiprocessing
_ref: Optional[mappy.Aligner] = None
_mapq: int = 0
_list_ex_chr: set[str] = set()
_allow_gap: int = 0
_allow_overlap: int = 0
_ref_merge_distance: int = 1000


def _init_trim_worker(
    fref: str,
    mapq: int,
    list_ex_chr: list[str],
    allow_gap: int,
    allow_overlap: int,
    preset: str,
    ref_merge_distance: int = 1000,
) -> None:
    """Initialize global variables in worker processes."""
    global _ref, _mapq, _list_ex_chr, _allow_gap, _allow_overlap, _ref_merge_distance
    MM_F_NO_LJOIN = 0x400
    _ref = mappy.Aligner(fref, preset=preset, extra_flags=MM_F_NO_LJOIN)
    _mapq = mapq
    _list_ex_chr = set(list_ex_chr)
    _allow_gap = allow_gap
    _allow_overlap = allow_overlap
    _ref_merge_distance = ref_merge_distance


def _check_read_pattern(name: str, seq: str) -> list:
    """Check for CTC (Circular Tandem Copy) pattern in read."""
    global _ref, _mapq, _list_ex_chr

    if _ref is None:
        raise RuntimeError("SplitReads-Core worker not initialized (reference aligner is None)")

    list_result = []
    len_seq = len(seq)

    list_hit = []
    for hit in _ref.map(seq):
        if hit.is_primary and hit.mapq >= _mapq and hit.ctg not in _list_ex_chr:
            list_local_hit = [
                name,
                hit.is_primary,
                hit.mapq,
                len_seq,
                hit.q_st,
                hit.q_en,
                hit.ctg,
                hit.r_st,
                hit.r_en,
                hit.mlen,
                hit.blen,
                hit.strand,
            ]
            list_hit.append(list_local_hit)

    if len(list_hit) > 1:
        header = [
            "readid",
            "is_primary",
            "mapq",
            "q_len",
            "q_start",
            "q_end",
            "ref",
            "r_start",
            "r_end",
            "matchLen",
            "blockLen",
            "strand",
        ]
        df_mapping = pd.DataFrame(list_hit, columns=header)

        df_mapping["mappedRegion"] = df_mapping.apply(combine_region_strand, axis=1)
        df_mapping = df_mapping.sort_values(by="q_start")
        df_mapping.reset_index(drop=True, inplace=True)

        list_check_pattern = df_mapping["mappedRegion"].tolist()
        list_group_order = make_group_order(list_check_pattern, 50)
        list_sum_pattern = get_sum_pattern(list_group_order)

        if list_sum_pattern[0]:
            list_pattern_region = list_sum_pattern[1].split(",")
            list_idx_pattern = get_idx_longest_pattern(list_group_order, list_pattern_region)

            if len(list_idx_pattern) > len(list_pattern_region):
                selected_cols = [
                    "readid",
                    "q_len",
                    "q_start",
                    "q_end",
                    "ref",
                    "r_start",
                    "r_end",
                    "matchLen",
                    "blockLen",
                    "mapq",
                    "strand",
                ]
                df_mod_mapping = df_mapping[df_mapping.index.isin(list_idx_pattern)][selected_cols].copy()
                df_mod_mapping["qlenTrimmed"] = (
                    df_mod_mapping["q_end"].tolist()[-1] - df_mod_mapping["q_start"].tolist()[0]
                )

                list_result = df_mod_mapping.values.tolist()

    return list_result


def _merge_pd_regions(prev: "pdRegion", curr: "pdRegion", columns: pd.Index) -> "pdRegion":
    """Merge two adjacent/overlapping pdRegion objects into one."""
    series_new = pd.Series(
        [
            curr.readid,
            curr.q_len,
            min(prev.q_start, curr.q_start),
            max(prev.q_end, curr.q_end),
            curr.ref,
            min(prev.r_start, curr.r_start),
            max(prev.r_end, curr.r_end),
            max(prev.matchLen, curr.matchLen),
            max(prev.blockLen, curr.blockLen),
            max(prev.mapq, curr.mapq),
            curr.strand,
        ],
        index=columns,
    )
    return pdRegion(series_new)


def _get_merge_all(name: str, seq: str) -> list:
    """Get all merged alignments for a read."""
    global _ref, _mapq, _list_ex_chr, _allow_gap, _allow_overlap, _ref_merge_distance

    if _ref is None:
        raise RuntimeError("SplitReads-Core worker not initialized (reference aligner is None)")

    list_merged_result = []
    len_seq = len(seq)

    list_hit = []
    for hit in _ref.map(seq):
        if hit.is_primary and hit.mapq >= _mapq and hit.ctg not in _list_ex_chr:
            list_local_hit = [
                name,
                len_seq,
                hit.q_st,
                hit.q_en,
                hit.ctg,
                hit.r_st,
                hit.r_en,
                hit.mlen,
                hit.blen,
                hit.mapq,
                hit.strand,
            ]
            list_hit.append(list_local_hit)

    if len(list_hit) > 0:
        header = [
            "readid",
            "q_len",
            "q_start",
            "q_end",
            "ref",
            "r_start",
            "r_end",
            "matchLen",
            "blockLen",
            "mapq",
            "strand",
        ]
        df_mapping = pd.DataFrame(list_hit, columns=header)

        df_mapping = df_mapping.sort_values(by="q_start")
        df_mapping.reset_index(drop=True, inplace=True)

        if len(df_mapping) > 1:
            # Make initial hit objects
            list_init_obj = [pdRegion(row) for _, row in df_mapping.iterrows()]

            # Merge overlapping/adjacent regions
            list_result_obj = [list_init_obj[0]]

            for obj_ in list_init_obj[1:]:
                prev = list_result_obj[-1]

                is_q_overlap = any_overlapping_range(
                    prev.q_start, prev.q_end, obj_.q_start, obj_.q_end
                )
                significant_overlap = (
                    is_q_overlap and abs(prev.q_end - obj_.q_start) >= _allow_overlap
                )
                same_locus = prev.ref == obj_.ref and prev.strand == obj_.strand
                close_on_ref = (
                    same_locus
                    and abs(prev.r_end - obj_.r_start) <= _ref_merge_distance
                )

                if significant_overlap:
                    # Deletion compensation: large query overlap
                    if close_on_ref:
                        list_result_obj[-1] = _merge_pd_regions(prev, obj_, df_mapping.columns)
                    elif same_locus:
                        list_result_obj.append(obj_)
                    elif obj_.cal_q_len() > prev.cal_q_len():
                        list_result_obj[-1] = obj_
                    # else: prev is bigger, skip obj_
                elif is_q_overlap:
                    # Minor overlap: keep as separate region
                    list_result_obj.append(obj_)
                else:
                    # Insertion compensation: check with gap tolerance
                    is_gap_overlap = any_overlapping_range(
                        prev.q_start, prev.q_end + _allow_gap,
                        obj_.q_start, obj_.q_end,
                    )
                    if is_gap_overlap and close_on_ref:
                        list_result_obj[-1] = _merge_pd_regions(prev, obj_, df_mapping.columns)
                    elif is_gap_overlap:
                        list_result_obj.append(obj_)
                    # else: no overlap at all, skip

            list_merged_result = [obj.get_list_attr() for obj in list_result_obj]
        else:
            list_merged_result = df_mapping.values.tolist()

    return list_merged_result


def _cal_pattern_trim(record: tuple[str, str]) -> list:
    """Calculate pattern and trim for a single record."""
    name, seq = record

    list_result = []

    list_chk_pattern = _check_read_pattern(name, seq)
    if len(list_chk_pattern) > 0:
        for idx, rTab in enumerate(list_chk_pattern, start=0):
            oTab = rTab + ["{:.2f}".format((rTab[3] - rTab[2]) / rTab[11]), idx, True]
            list_result.append(oTab)
    else:
        list_merged_result = _get_merge_all(name, seq)
        if len(list_merged_result) > 0:
            len_trim_read = list_merged_result[-1][3] - list_merged_result[0][2]

            for idx, rTab in enumerate(list_merged_result, start=0):
                oTab = rTab + [len_trim_read, "{:.2f}".format((rTab[3] - rTab[2]) / len_trim_read), idx, False]
                list_result.append(oTab)

    return list_result


# ============================================================================
# Identify Module Functions
# ============================================================================


def lenLoci(loci: str) -> int:
    """Calculate total length of loci string."""
    length = 0
    for region in loci.split(","):
        chrom, start, end, strand = region.rsplit("_", 3)
        length += int(end) - int(start)
    return length


def format_merge_region(merge_region: str) -> str:
    """Format merge region string."""
    list_formatted = []
    for region in merge_region.split(","):
        chrom, start, end, strand = region.rsplit("_", 3)
        list_formatted.append(f"{chrom}:{int(start) + 1}-{end}_{strand}")
    return ",".join(list_formatted)


def majority_strand(df_merged: pd.DataFrame) -> str:
    """Get majority strand from merged dataframe."""
    strand = "+"
    list_str_strand = df_merged["strand"].astype(str).tolist()
    if list_str_strand.count("-1") > list_str_strand.count("1"):
        strand = "-"
    return strand


def reverse_strand(strand: str) -> str:
    """Reverse strand characters."""
    return "".join(["+" if x == "-" else "-" if x == "+" else x for x in list(strand)])


def check_breakpoint_direction(df_check: pd.DataFrame) -> list[tuple]:
    """Check breakpoint direction for consecutive alignments."""
    list_pairs = []

    df_check = df_check.reset_index(drop=True)

    list_order_check = [(i, i + 1) for i in df_check.index.tolist() if i + 1 < len(df_check)]

    for tup_check in list_order_check:
        l1 = df_check.index.isin(tup_check)

        q_starts = df_check[l1]["q_start"].tolist()
        q_ends = df_check[l1]["q_end"].tolist()

        if any_overlapping_range(q_starts[0], q_ends[0] + 50, q_starts[1], q_ends[1]):
            list_mergeid = df_check[l1]["mergeid"].tolist()
            list_check = df_check[l1]["strand"].tolist()

            if list_check[0] == 1 and list_check[1] == 1:
                list_5 = df_check[l1]["ovl_5end"].tolist()
                list_3 = df_check[l1]["ovl_3end"].tolist()
                if list_3[0] == 1 and list_5[1] == 1:
                    pair = (list_mergeid[0], list_mergeid[1], "+_+", True)
                    list_pairs.append(pair)

            if list_check[0] == -1 and list_check[1] == -1:
                list_5 = df_check[l1]["ovl_5end"].tolist()
                list_3 = df_check[l1]["ovl_3end"].tolist()
                if list_3[1] == 1 and list_5[0] == 1:
                    pair = (list_mergeid[1], list_mergeid[0], "-_-", True)
                    list_pairs.append(pair)

            if list_check[0] == -1 and list_check[1] == 1:
                list_5 = df_check[l1]["ovl_5end"].tolist()
                list_3 = df_check[l1]["ovl_3end"].tolist()
                if list_5[0] == 1 and list_5[1] == 1:
                    pair = (list_mergeid[0], list_mergeid[1], "-_+", True)
                    list_pairs.append(pair)

            if list_check[0] == 1 and list_check[1] == -1:
                list_5 = df_check[l1]["ovl_5end"].tolist()
                list_3 = df_check[l1]["ovl_3end"].tolist()
                if list_3[0] == 1 and list_3[1] == 1:
                    pair = (list_mergeid[0], list_mergeid[1], "+_-", True)
                    list_pairs.append(pair)
        else:
            break

    return list_pairs


def chk_circular_subgraph(
    G: nx.MultiDiGraph, graph: nx.Graph, dict_pair_strand: dict
) -> tuple[str, int, bool, bool, bool]:
    """Check if subgraph represents a circular structure."""
    list_nodes = list(graph.nodes)
    num_nodes = len(list_nodes)
    print_nodes = ",".join(list_nodes)
    sl = len(list(nx.selfloop_edges(graph))) > 0
    solved = True
    cyclic = True

    if num_nodes == 2:
        rm_sl_subgraph = graph.copy()
        rm_sl_subgraph.remove_edges_from(nx.selfloop_edges(rm_sl_subgraph))
        solved = all(x <= 2 for x in [rm_sl_subgraph.degree[node] for node in list_nodes])
        key_forward = repr(list_nodes)
        key_reverse = repr(list_nodes[::-1])
        cyclic = False
        if all(x in dict_pair_strand.keys() for x in [key_forward, key_reverse]):
            for fstrands in dict_pair_strand[key_forward]:
                fl, fr = fstrands.split("_")
                key_check = f"{fr}{key_reverse}{fl}"
                for rstrands in dict_pair_strand[key_reverse]:
                    rl, rr = rstrands.split("_")
                    key_reverse_check = f"{rl}{key_reverse}{rr}"
                    if key_check == key_reverse_check:
                        cyclic = True
                        break
    elif num_nodes > 2:
        rm_sl_subgraph = graph.copy()
        rm_sl_subgraph.remove_edges_from(nx.selfloop_edges(rm_sl_subgraph))
        solved = all(x <= 2 for x in [rm_sl_subgraph.degree[node] for node in list_nodes])
        cyclic = len(nx.cycle_basis(nx.DiGraph(rm_sl_subgraph).to_undirected())) > 0
    else:
        cyclic = len(nx.cycle_basis(nx.DiGraph(graph).to_undirected())) > 0

    if not solved:
        cyclic = False

    test_graph = graph.copy()
    test_graph.remove_edges_from(nx.selfloop_edges(test_graph))
    list_traversal = nx.cycle_basis(nx.DiGraph(test_graph).to_undirected())
    if len(list_traversal) > 0:
        list_traversal_ini = list_traversal[0]
        if len(list_traversal_ini) == len(list_nodes):
            print_nodes = ",".join(list_traversal_ini)

    return (print_nodes, num_nodes, solved, sl, cyclic)


# ============================================================================
# Main SplitReads-Core Class
# ============================================================================


class SplitReadsCore:
    """Integrated split-read based eccDNA detection for HiFi data."""

    def __init__(
        self,
        reference_mmi: Path,
        config: Optional[SplitReadsConfig] = None,
        logger: Optional[logging.Logger] = None,
        threads: int = 0,
    ):
        """Initialize SplitReadsCore.

        Args:
            reference_mmi: Path to reference minimap2 index (.mmi)
            config: Configuration object
            logger: Logger instance
            threads: Number of threads (0 = auto)
        """
        self.reference_mmi = Path(reference_mmi)
        self.config = config or SplitReadsConfig()
        self.logger = logger or get_logger(self.__class__.__name__)
        self.threads = threads if threads > 0 else cpu_count()

        # Parse excluded chromosomes
        self.exclude_chrs: list[str] = []
        if self.config.exclude_chrs:
            self.exclude_chrs = [c.strip() for c in self.config.exclude_chrs.split(",") if c.strip()]

    def run(
        self,
        input_fasta: Path,
        chrom_sizes_file: Path,
        output_dir: Path,
    ) -> Path:
        """Run the full SplitReads-Core pipeline.

        Args:
            input_fasta: Input FASTA file
            chrom_sizes_file: Chromosome sizes file (.fai format)
            output_dir: Output directory

        Returns:
            Path to eccDNA_final.txt output file
        """
        input_fasta = Path(input_fasta)
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        self.logger.info("Starting SplitReads-Core pipeline (HiFi optimized)")
        self.logger.info(f"  Input: {input_fasta}")
        self.logger.info(f"  Reference: {self.reference_mmi}")
        self.logger.info(f"  Threads: {self.threads}")

        # Step 1: Trim
        self.logger.info("Phase 1/2: Trim (alignment and pattern detection)")
        trim_df = self._run_trim(input_fasta)
        self.logger.info(f"Trim complete: {len(trim_df)} aligned regions")

        if len(trim_df) == 0:
            self.logger.warning("No aligned regions found; creating empty output")
            return self._write_empty_output(output_dir)

        # Step 2: Identify
        self.logger.info("Phase 2/2: Identify (eccDNA detection)")
        eccDNA_file = self._run_identify(trim_df, chrom_sizes_file, output_dir)
        self.logger.info(f"SplitReads-Core complete. Output: {eccDNA_file}")

        return eccDNA_file

    def _run_trim(self, input_fasta: Path) -> pd.DataFrame:
        """Run trim phase."""
        # Read all sequences
        records: list[tuple[str, str]] = []
        for record in SeqIO.parse(str(input_fasta), "fasta"):
            records.append((record.id, str(record.seq)))

        if not records:
            return pd.DataFrame()

        self.logger.info(f"Processing {len(records)} reads")

        # Process in parallel
        batch_size = 10000
        all_results: list[list[Any]] = []

        for i in range(0, len(records), batch_size):
            batch = records[i : i + batch_size]

            with Pool(
                processes=self.threads,
                initializer=_init_trim_worker,
                initargs=(
                    str(self.reference_mmi),
                    self.config.mapq,
                    self.exclude_chrs,
                    self.config.allow_gap,
                    self.config.allow_overlap,
                    self.config.preset,
                    self.config.ref_merge_distance,
                ),
            ) as pool:
                batch_results = pool.map(_cal_pattern_trim, batch)

            all_results.extend(chain.from_iterable(batch_results))

            if (i + batch_size) % 10000 == 0:
                self.logger.debug(f"Processed {min(i + batch_size, len(records))} reads")

        if not all_results:
            return pd.DataFrame()

        # Create dataframe
        columns = [
            "readid",
            "q_len",
            "q_start",
            "q_end",
            "ref",
            "r_start",
            "r_end",
            "match",
            "mapBlock",
            "mapq",
            "strand",
            "qlenTrimmed",
            "freqCov",
            "order",
            "ctc_read",
        ]
        df = pd.DataFrame(all_results, columns=columns)

        return df

    def _run_identify(self, trim_df: pd.DataFrame, chrom_sizes_file: Path, output_dir: Path) -> Path:
        """Run identify phase."""
        tmp_dir = output_dir / "tmp"
        tmp_dir.mkdir(parents=True, exist_ok=True)

        # Prepare aligned reads
        ord_header = [
            "ref",
            "r_start",
            "r_end",
            "readid",
            "q_start",
            "q_end",
            "match",
            "mapBlock",
            "mapq",
            "strand",
            "qlenTrimmed",
            "freqCov",
            "order",
        ]
        read_trim = trim_df.loc[:, ord_header].copy()

        # Create BedTool from aligned reads
        aln_reads = bt.BedTool.from_dataframe(read_trim).sort()

        # Calculate coverage
        genome_cov = aln_reads.genome_coverage(bg=True, g=str(chrom_sizes_file))

        # Merge adjacent regions
        aln_reads_merge = aln_reads.merge()

        # Join merged regions with coverage
        merge_genomecov_path = tmp_dir / "merge_genomecov.txt"
        aln_reads_merge.intersect(genome_cov, output=str(merge_genomecov_path), wo=True)

        merge_genomecov_df = pd.read_csv(merge_genomecov_path, sep="\t", header=None, dtype=str)
        merge_genomecov_df.columns = [
            "m_chrom",
            "m_start",
            "m_end",
            "bg_chrom",
            "bg_start",
            "bg_end",
            "depth",
            "d_ovl",
        ]
        merge_genomecov_df[["depth", "bg_start", "d_ovl"]] = merge_genomecov_df[["depth", "bg_start", "d_ovl"]].apply(
            pd.to_numeric, errors="coerce"
        )

        # Generate merge ID
        merge_genomecov_df["mergeid"] = (
            merge_genomecov_df["m_chrom"] + "_" + merge_genomecov_df["m_start"] + "_" + merge_genomecov_df["m_end"]
        )

        # Filter by depth
        merge_genomecov_df_filt = merge_genomecov_df[merge_genomecov_df["depth"] >= self.config.average_depth]

        if merge_genomecov_df_filt.empty:
            self.logger.warning("No regions with sufficient depth")
            return self._write_empty_output(output_dir)

        # Group and aggregate
        merge_bg_filt = (
            merge_genomecov_df_filt.groupby(["mergeid"])
            .agg({"bg_chrom": "max", "bg_start": "min", "bg_end": "max", "depth": "mean"})
            .reset_index()
        )

        merge_bg_filt[["bg_start", "bg_end"]] = merge_bg_filt[["bg_start", "bg_end"]].apply(
            pd.to_numeric, errors="coerce"
        )
        merge_bg_filt["length"] = merge_bg_filt["bg_end"] - merge_bg_filt["bg_start"]
        merge_bg_filt["bg_start"] = merge_bg_filt["bg_start"].astype(str)
        merge_bg_filt["bg_end"] = merge_bg_filt["bg_end"].astype(str)
        merge_bg_filt["mergeid"] = merge_bg_filt["bg_chrom"] + "_" + merge_bg_filt["bg_start"] + "_" + merge_bg_filt["bg_end"]
        merge_bg_filt = merge_bg_filt.loc[:, ["bg_chrom", "bg_start", "bg_end", "depth", "length", "mergeid"]]

        # Filter by minimum size
        merge_bg_filt = merge_bg_filt[merge_bg_filt["length"] >= self.config.min_region_size]

        if merge_bg_filt.empty:
            self.logger.warning("No regions meeting size criteria")
            return self._write_empty_output(output_dir)

        self.logger.info(f"Found {len(merge_bg_filt)} potential eccDNA regions")

        # Create merge region BedTool
        aln_reads_merge = bt.BedTool.from_dataframe(merge_bg_filt)
        reads_merge_intersect_path = tmp_dir / "reads_merge_intersect.bed"
        aln_reads_merge.intersect(aln_reads, output=str(reads_merge_intersect_path), wo=True)

        if reads_merge_intersect_path.stat().st_size == 0:
            self.logger.warning("No overlapping reads found")
            return self._write_empty_output(output_dir)

        read_merged_ins_df_org = pd.read_csv(reads_merge_intersect_path, sep="\t", header=None, dtype=str)
        header = list(merge_bg_filt.columns) + list(read_trim.columns) + ["ovl"]
        read_merged_ins_df_org.columns = header
        header_select = [
            "ref",
            "r_start",
            "r_end",
            "mergeid",
            "depth",
            "length",
            "readid",
            "q_start",
            "q_end",
            "match",
            "mapBlock",
            "mapq",
            "strand",
            "qlenTrimmed",
            "freqCov",
            "order",
            "ovl",
        ]
        read_merged_ins_df_org = read_merged_ins_df_org[header_select]

        # Annotate 5' and 3' ends
        end5_path = tmp_dir / "end5_merge_region.bed"
        end3_path = tmp_dir / "end3_merge_region.bed"

        with open(end5_path, "w") as end5, open(end3_path, "w") as end3:
            for idx_, row in merge_bg_filt.iterrows():
                chrom = row["bg_chrom"]
                start = row["bg_start"]
                end = row["bg_end"]
                end_size = (
                    self.config.overlap_check_size
                    if self.config.min_region_size >= 200
                    else int(round(self.config.min_region_size * 0.3, 0))
                )
                end_size = max(end_size, 15)
                end5.write(f"{chrom}\t{start}\t{int(start) + end_size}\n")
                end3.write(f"{chrom}\t{int(end) - end_size}\t{end}\n")

        read_merged_ins_df_org_bt = bt.BedTool.from_dataframe(read_merged_ins_df_org)
        read_merged_final_bt = read_merged_ins_df_org_bt.annotate(files=[str(end5_path), str(end3_path)], counts=True)

        read_merged_ins_df = read_merged_final_bt.to_dataframe(names=header_select + ["ovl_5end", "ovl_3end"])
        read_merged_ins_df["sum_ends"] = read_merged_ins_df["ovl_5end"] + read_merged_ins_df["ovl_3end"]

        # Create strand dictionary
        dict_majority_strand = {}
        for group_name, group in read_merged_ins_df.groupby(by="mergeid"):
            dict_majority_strand[group_name] = majority_strand(group)

        # Build graph
        self.logger.info("Analyzing breakpoint graphs")
        dict_pair_strand: dict[str, set] = {}
        graph: dict[str, int] = {}

        for readid, group in read_merged_ins_df.groupby(by="readid"):
            df_check = group.sort_values(by="order").copy()
            df_check.reset_index(drop=True, inplace=True)

            if len(df_check) > 1:
                for tup_result in check_breakpoint_direction(df_check):
                    if tup_result[3]:
                        repr_mergeid = repr([tup_result[0], tup_result[1]])

                        if repr_mergeid in dict_pair_strand:
                            dict_pair_strand[repr_mergeid].add(tup_result[2])
                        else:
                            dict_pair_strand[repr_mergeid] = {tup_result[2]}

                        if repr_mergeid in graph:
                            graph[repr_mergeid] += 1
                        else:
                            graph[repr_mergeid] = 1

        # Filter by breakpoint depth
        graph_filt = {k: v for k, v in graph.items() if v >= self.config.breakpoint_depth}

        # Create NetworkX graph
        G = nx.MultiDiGraph()
        for k in graph_filt:
            nodes = eval(k)
            G.add_edge(nodes[0], nodes[1], weight=graph_filt[k])

        # Get connected components (subgraphs)
        subgraphs = list(nx.connected_components(G.to_undirected()))

        self.logger.info(f"Found {len(subgraphs)} potential circular subgraphs")

        if not subgraphs:
            return self._write_empty_output(output_dir)

        # Analyze each subgraph
        list_graph_summary = []

        for idx, comp_nodes in enumerate(subgraphs, start=1):
            gname = f"ec{idx}"
            subgraph = G.to_undirected().subgraph(comp_nodes)
            nodes = list(subgraph.nodes())

            regions, num_nodes, can_be_solved, contain_selfloop, is_cyclic = chk_circular_subgraph(
                G, subgraph, dict_pair_strand
            )

            # Determine region string with strands
            if len(nodes) == 1:
                select_repr = repr([nodes[0], nodes[0]])
                if select_repr in dict_pair_strand:
                    regions = f"{eval(select_repr)[0]}_{list(dict_pair_strand[select_repr])[0][0]}"
                else:
                    regions = f"{nodes[0]}_{dict_majority_strand.get(nodes[0], '+')}"

            elif len(nodes) == 2:
                select_repr = repr(nodes)
                list_region = eval(select_repr)

                if select_repr not in dict_pair_strand:
                    select_repr = repr(nodes[::-1])

                if select_repr in dict_pair_strand:
                    list_strand = list(dict_pair_strand[select_repr])[0].split("_")
                    regions = ",".join(f"{x[0]}_{x[1]}" for x in zip(list_region, list_strand))
                else:
                    regions = ",".join(f"{n}_{dict_majority_strand.get(n, '+')}" for n in nodes)

            else:
                test_graph = subgraph.copy()
                test_graph.remove_edges_from(nx.selfloop_edges(test_graph))

                list_traversal = nx.cycle_basis(nx.DiGraph(test_graph).to_undirected())
                if len(list_traversal) == 0:
                    regions = ",".join(f"{node}_{dict_majority_strand.get(node, '+')}" for node in nodes)
                else:
                    list_traversal = list_traversal[0]
                    list_order_check = [(i, i + 1) for i in range(len(list_traversal)) if i + 1 < len(list_traversal)]

                    list_temp_all: list[list[list[str]]] = []
                    for tup_check in list_order_check:
                        l_region = list_traversal[tup_check[0]]
                        r_region = list_traversal[tup_check[1]]
                        repr_mergeid = repr([l_region, r_region])

                        list_this_level: list[list[str]] = []

                        if repr_mergeid in dict_pair_strand:
                            for str_strand in list(dict_pair_strand[repr_mergeid]):
                                list_str_strand = [
                                    "_".join(x) for x in zip([l_region, r_region], str_strand.split("_"))
                                ]
                                list_this_level.append(list_str_strand)
                        else:
                            rev_repr_mergeid = repr([r_region, l_region])
                            if rev_repr_mergeid in dict_pair_strand:
                                for str_strand in [reverse_strand(x) for x in list(dict_pair_strand[rev_repr_mergeid])]:
                                    list_str_strand = [
                                        "_".join(x) for x in zip([l_region, r_region], str_strand.split("_"))
                                    ]
                                    list_this_level.append(list_str_strand)

                        if list_this_level:
                            list_temp_all.append(list_this_level)

                    if list_temp_all:
                        list_traverse: list[list[str]] = list_temp_all[0]
                        for list_level in list_temp_all[1:]:
                            for list_ in list_level:
                                for index, list_traverse_level in enumerate(list_traverse):
                                    if list_[0] == list_traverse_level[-1]:
                                        list_traverse[index].append(list_[1])

                        regions = ",".join(max(list_traverse, key=lambda xs: len(list(xs))))
                    else:
                        regions = ",".join(f"{node}_{dict_majority_strand.get(node, '+')}" for node in nodes)

            tup_graph_summary = (gname, regions, num_nodes, can_be_solved, contain_selfloop, is_cyclic)
            list_graph_summary.append(tup_graph_summary)

        if not list_graph_summary:
            return self._write_empty_output(output_dir)

        # Create summary dataframe
        header = ["id", "regions", "num_nodes", "can_be_solved", "contain_selfloop", "is_cyclic"]
        df_graph_summary = pd.DataFrame(list_graph_summary, columns=header)

        # Calculate statistics for each eccDNA
        list_identify_results = []

        for idx, value in df_graph_summary.iterrows():
            gname = value["id"]
            merge_region = value["regions"]
            num_region = value["num_nodes"]
            is_cyclic = value["is_cyclic"]

            length_loci = lenLoci(merge_region)
            ctc = "True" if is_cyclic else "False"

            # Collect reads
            total_base = 0
            list_readid = []

            for node in merge_region.split(","):
                node_no_strand = node.rsplit("_", 1)[0] if "_" in node else node
                l1 = read_merged_ins_df["mergeid"] == node_no_strand
                read_o = read_merged_ins_df[l1].copy()
                for _, read_o_row in read_o.iterrows():
                    read_o_readid = str(read_o_row["readid"])
                    read_o_start = int(read_o_row["q_start"])
                    read_o_end = int(read_o_row["q_end"])
                    total_base += read_o_end - read_o_start
                    list_readid.append(read_o_readid)

            expect_cov = f"{float(total_base) / length_loci:.2f}" if length_loci > 0 else "0.00"

            list_identify_results.append(
                (gname, merge_region, length_loci, num_region, ctc, len(set(list_readid)), total_base, expect_cov)
            )

        # Write output
        header = ["id", "merge_region", "merge_len", "num_region", "ctc", "numreads", "totalbase", "coverage"]
        df_identify_summary = pd.DataFrame(list_identify_results, columns=header)
        df_identify_summary["merge_region"] = df_identify_summary["merge_region"].apply(format_merge_region)

        eccdna_final_path = output_dir / "eccDNA_final.txt"
        df_identify_summary.to_csv(eccdna_final_path, sep="\t", index=None)

        self.logger.info(f"Identified {len(df_identify_summary)} eccDNA candidates")

        # Cleanup
        bt.cleanup()

        return eccdna_final_path

    def _write_empty_output(self, output_dir: Path) -> Path:
        """Write empty output file."""
        eccdna_final_path = output_dir / "eccDNA_final.txt"
        header = ["id", "merge_region", "merge_len", "num_region", "ctc", "numreads", "totalbase", "coverage"]
        df_empty = pd.DataFrame(columns=header)
        df_empty.to_csv(eccdna_final_path, sep="\t", index=None)
        return eccdna_final_path
