from pathlib import Path
import sys

import pandas as pd
import pytest
import networkx as nx

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.modules.splitreads_core import (
    SplitReadsConfig,
    Region,
    PdRegion,
    is_overlapping,
    any_overlapping_range,
    combine_region_strand,
    make_group_order,
    split_list_step_window,
    get_sum_pattern,
    get_idx_longest_pattern,
    len_loci,
    format_merge_region,
    majority_strand,
    reverse_strand,
    check_breakpoint_direction,
    chk_circular_subgraph,
    SplitReadsCore,
)

import logging


# ============================================================================
# TestSplitReadsConfig
# ============================================================================


class TestSplitReadsConfig:
    def test_from_dict_basic(self):
        cfg = SplitReadsConfig.from_dict({"preset": "map-hifi", "mapq": 30})
        assert cfg.preset == "map-hifi"
        assert cfg.mapq == 30

    def test_from_dict_defaults(self):
        cfg = SplitReadsConfig.from_dict({})
        assert cfg.preset == "map-hifi"
        assert cfg.mapq == 20
        assert cfg.exclude_chrs == ""
        assert cfg.allow_gap == 10
        assert cfg.allow_overlap == 10
        assert cfg.min_region_size == 200
        assert cfg.overlap_check_size == 50
        assert cfg.breakpoint_depth == 5
        assert cfg.average_depth == 5.0
        assert cfg.threads == 0
        assert cfg.skip_variant is True

    def test_from_dict_extra_keys_ignored(self):
        cfg = SplitReadsConfig.from_dict({"foo": "bar", "preset": "map-hifi"})
        assert cfg.preset == "map-hifi"

    def test_from_dict_partial(self):
        cfg = SplitReadsConfig.from_dict({"mapq": 50})
        assert cfg.mapq == 50
        assert cfg.preset == "map-hifi"
        assert cfg.threads == 0
        assert cfg.allow_gap == 10

    def test_default_threads_zero(self):
        cfg = SplitReadsConfig()
        assert cfg.threads == 0


# ============================================================================
# TestRegion
# ============================================================================


class TestRegion:
    def test_region_attributes(self):
        r = Region("chr1", 100, 200)
        assert r.chrom == "chr1"
        assert r.start == 100
        assert r.end == 200

    def test_get_str_region(self):
        r = Region("chr1", 100, 200)
        assert r.get_str_region() == "chr1_100_200"

    def test_region_different_chrom(self):
        r = Region("chrX", 0, 500)
        assert r.chrom == "chrX"
        assert r.start == 0
        assert r.end == 500


# ============================================================================
# TestPdRegion
# ============================================================================


class TestPdRegion:
    @pytest.fixture()
    def sample_series(self):
        return pd.Series({
            "readid": "read1", "q_len": 1000, "q_start": 100, "q_end": 300,
            "ref": "chr1", "r_start": 1000, "r_end": 1500,
            "matchLen": 190, "blockLen": 200, "mapq": 60, "strand": 1,
        })

    def test_calc_q_len(self, sample_series):
        pr = PdRegion(sample_series)
        assert pr.calc_q_len() == 200

    def test_cal_r_len(self, sample_series):
        pr = PdRegion(sample_series)
        assert pr.cal_r_len() == 500

    def test_get_list_attr(self, sample_series):
        pr = PdRegion(sample_series)
        attrs = pr.get_list_attr()
        assert len(attrs) == 11
        assert attrs == ["read1", 1000, 100, 300, "chr1", 1000, 1500, 190, 200, 60, 1]

    def test_pdregion_stores_all_attributes(self, sample_series):
        pr = PdRegion(sample_series)
        assert pr.readid == "read1"
        assert pr.q_len == 1000
        assert pr.q_start == 100
        assert pr.q_end == 300
        assert pr.ref == "chr1"
        assert pr.r_start == 1000
        assert pr.r_end == 1500
        assert pr.matchLen == 190
        assert pr.blockLen == 200
        assert pr.mapq == 60
        assert pr.strand == 1


# ============================================================================
# TestIsOverlapping
# ============================================================================


class TestIsOverlapping:
    def test_exact_match(self):
        r = Region("chr1", 100, 200)
        assert is_overlapping(r, 0, "chr1", 100, 200) is True

    def test_within_offset(self):
        r = Region("chr1", 100, 200)
        assert is_overlapping(r, 10, "chr1", 105, 195) is True

    def test_outside_offset(self):
        r = Region("chr1", 100, 200)
        assert is_overlapping(r, 5, "chr1", 200, 300) is False

    def test_different_chrom(self):
        r = Region("chr1", 100, 200)
        assert is_overlapping(r, 100, "chr2", 100, 200) is False

    def test_zero_offset_exact(self):
        r = Region("chr1", 100, 200)
        assert is_overlapping(r, 0, "chr1", 100, 200) is True

    def test_boundary_at_offset(self):
        r = Region("chr1", 100, 200)
        assert is_overlapping(r, 10, "chr1", 110, 210) is True


# ============================================================================
# TestAnyOverlappingRange
# ============================================================================


class TestAnyOverlappingRange:
    def test_complete_overlap(self):
        assert any_overlapping_range(100, 200, 100, 200) is True

    def test_partial_overlap(self):
        assert any_overlapping_range(100, 200, 150, 250) is True

    def test_adjacent_touch(self):
        assert any_overlapping_range(100, 200, 200, 300) is True

    def test_no_overlap(self):
        assert any_overlapping_range(100, 200, 300, 400) is False

    def test_contained(self):
        assert any_overlapping_range(100, 400, 200, 300) is True


# ============================================================================
# TestCombineRegionStrand
# ============================================================================


class TestCombineRegionStrand:
    def test_positive_strand(self):
        s = pd.Series({"ref": "chr1", "r_start": 100, "r_end": 200, "strand": "+"})
        assert combine_region_strand(s) == "chr1_100_200_+"

    def test_negative_strand(self):
        s = pd.Series({"ref": "chr1", "r_start": 100, "r_end": 200, "strand": "-"})
        assert combine_region_strand(s) == "chr1_100_200_-"


# ============================================================================
# TestMakeGroupOrder
# ============================================================================


class TestMakeGroupOrder:
    def test_single_region(self):
        result = make_group_order(["chr1_100_200_+"], 50)
        assert result == ["chr1_100_200_+"]

    def test_overlapping_regions(self):
        result = make_group_order(["chr1_100_200_+", "chr1_110_210_-"], 50)
        # Second region overlaps first within offset=50, so it maps to first region's coords + its own strand
        assert result[0] == "chr1_100_200_+"
        assert result[1] == "chr1_100_200_-"

    def test_non_overlapping(self):
        result = make_group_order(["chr1_100_200_+", "chr1_500_600_+"], 50)
        assert len(result) == 2
        assert result[0] == "chr1_100_200_+"
        assert result[1] == "chr1_500_600_+"

    def test_empty_list(self):
        result = make_group_order([], 50)
        assert result == []

    def test_cross_chromosome(self):
        result = make_group_order(["chr1_100_200_+", "chr2_100_200_+"], 50)
        assert len(result) == 2
        assert result[0] == "chr1_100_200_+"
        assert result[1] == "chr2_100_200_+"


# ============================================================================
# TestSplitListStepWindow
# ============================================================================


class TestSplitListStepWindow:
    def test_basic_split(self):
        result = split_list_step_window([1, 2, 3, 4, 5, 6], step=1, window=3)
        assert result == [[1, 2, 3], [3, 4, 5], [5, 6]]

    def test_no_overlap(self):
        result = split_list_step_window([1, 2, 3, 4], step=0, window=2)
        assert result == [[1, 2], [3, 4]]

    def test_single_element(self):
        result = split_list_step_window([1], step=0, window=1)
        assert result == [[1]]

    def test_empty_list(self):
        result = split_list_step_window([], step=0, window=2)
        assert result == []


# ============================================================================
# TestGetSumPattern
# ============================================================================


class TestGetSumPattern:
    def test_repeated_single_pattern(self):
        result = get_sum_pattern(["A", "A", "A"])
        assert result == [True, "A", 3]

    def test_no_repeat(self):
        result = get_sum_pattern(["A", "B", "C"])
        assert result == [False, "", 0]

    def test_ctc_pattern(self):
        result = get_sum_pattern(["A", "B", "A", "B"])
        assert result[0] is True

    def test_single_element(self):
        result = get_sum_pattern(["A"])
        assert result == [False, "", 0]

    def test_repeated_multi_pattern(self):
        result = get_sum_pattern(["A", "B", "A", "B", "A", "B"])
        assert result[0] is True


# ============================================================================
# TestGetIdxLongestPattern
# ============================================================================


class TestGetIdxLongestPattern:
    def test_simple_match(self):
        result = get_idx_longest_pattern(["A", "B", "A", "B"], ["A", "B"])
        assert len(result) > 0
        # Should match all four indices as two consecutive AB patterns merge
        assert result == [0, 1, 2, 3]

    def test_no_match(self):
        result = get_idx_longest_pattern(["A", "B", "C"], ["X"])
        assert result == []

    def test_intermittent_match(self):
        result = get_idx_longest_pattern(["A", "B", "C", "A", "B"], ["A", "B"])
        # Two matches at [0,1] and [3,4] but not contiguous, return longest
        assert len(result) == 2

    def test_consecutive_merge(self):
        result = get_idx_longest_pattern(["A", "B", "A", "B", "A", "B"], ["A", "B"])
        # Three consecutive AB matches -> all merged into one range
        assert result == [0, 1, 2, 3, 4, 5]


# ============================================================================
# TestLenLoci
# ============================================================================


class TestLenLoci:
    def test_single_region(self):
        assert len_loci("chr1_100_200_+") == 100

    def test_multiple_regions(self):
        assert len_loci("chr1_100_200_+,chr2_300_500_-") == 300

    def test_zero_length(self):
        assert len_loci("chr1_100_100_+") == 0


# ============================================================================
# TestFormatMergeRegion
# ============================================================================


class TestFormatMergeRegion:
    def test_single_region(self):
        assert format_merge_region("chr1_100_200_+") == "chr1:101-200_+"

    def test_multiple_regions(self):
        result = format_merge_region("chr1_100_200_+,chr2_300_500_-")
        assert result == "chr1:101-200_+,chr2:301-500_-"


# ============================================================================
# TestMajorityStrand
# ============================================================================


class TestMajorityStrand:
    def test_majority_positive(self):
        df = pd.DataFrame({"strand": [1, 1, -1]})
        assert majority_strand(df) == "+"

    def test_majority_negative(self):
        df = pd.DataFrame({"strand": [-1, -1, 1]})
        assert majority_strand(df) == "-"

    def test_all_same(self):
        df = pd.DataFrame({"strand": [1, 1]})
        assert majority_strand(df) == "+"


# ============================================================================
# TestReverseStrand
# ============================================================================


class TestReverseStrand:
    def test_plus_to_minus(self):
        assert reverse_strand("+") == "-"

    def test_minus_to_plus(self):
        assert reverse_strand("-") == "+"

    def test_mixed(self):
        assert reverse_strand("+-") == "-+"


# ============================================================================
# TestCheckBreakpointDirection
# ============================================================================


class TestCheckBreakpointDirection:
    def test_plus_plus_pattern(self):
        df = pd.DataFrame({
            "q_start": [100, 150],
            "q_end": [200, 250],
            "strand": [1, 1],
            "mergeid": ["regionA", "regionB"],
            "ovl_5end": [0, 1],
            "ovl_3end": [1, 0],
        })
        result = check_breakpoint_direction(df)
        assert len(result) == 1
        assert result[0][2] == "+_+"

    def test_minus_minus_pattern(self):
        df = pd.DataFrame({
            "q_start": [100, 150],
            "q_end": [200, 250],
            "strand": [-1, -1],
            "mergeid": ["regionA", "regionB"],
            "ovl_5end": [1, 0],
            "ovl_3end": [0, 1],
        })
        result = check_breakpoint_direction(df)
        assert len(result) == 1
        assert result[0][2] == "-_-"

    def test_minus_plus_pattern(self):
        df = pd.DataFrame({
            "q_start": [100, 150],
            "q_end": [200, 250],
            "strand": [-1, 1],
            "mergeid": ["regionA", "regionB"],
            "ovl_5end": [1, 1],
            "ovl_3end": [0, 0],
        })
        result = check_breakpoint_direction(df)
        assert len(result) == 1
        assert result[0][2] == "-_+"

    def test_plus_minus_pattern(self):
        df = pd.DataFrame({
            "q_start": [100, 150],
            "q_end": [200, 250],
            "strand": [1, -1],
            "mergeid": ["regionA", "regionB"],
            "ovl_5end": [0, 0],
            "ovl_3end": [1, 1],
        })
        result = check_breakpoint_direction(df)
        assert len(result) == 1
        assert result[0][2] == "+_-"

    def test_no_overlap_breaks(self):
        df = pd.DataFrame({
            "q_start": [100, 500],
            "q_end": [200, 600],
            "strand": [1, 1],
            "mergeid": ["regionA", "regionB"],
            "ovl_5end": [0, 1],
            "ovl_3end": [1, 0],
        })
        result = check_breakpoint_direction(df)
        assert result == []


# ============================================================================
# TestChkCircularSubgraph
# ============================================================================


class TestChkCircularSubgraph:
    def test_self_loop(self):
        G = nx.MultiDiGraph()
        G.add_edge("A", "A")
        subgraph = G.to_undirected().subgraph(["A"])
        dict_pair_strand = {}
        result = chk_circular_subgraph(G, subgraph, dict_pair_strand)
        # Single node with self-loop
        assert result[3] is True  # sl (contains self-loop)
        assert result[4] is True  # cyclic (cycle_basis on single node with self-loop)

    def test_two_node_cycle(self):
        G = nx.MultiDiGraph()
        G.add_edge("A", "B")
        G.add_edge("B", "A")
        subgraph = G.to_undirected().subgraph(["A", "B"])
        dict_pair_strand = {
            repr(["A", "B"]): {"+_+"},
            repr(["B", "A"]): {"+_+"},
        }
        result = chk_circular_subgraph(G, subgraph, dict_pair_strand)
        assert result[4] is True  # cyclic

    def test_non_cycle(self):
        G = nx.MultiDiGraph()
        G.add_edge("A", "B")
        subgraph = G.to_undirected().subgraph(["A", "B"])
        # Only one direction, strand pairs don't form a cycle match
        dict_pair_strand = {repr(["A", "B"]): {"+_-"}}
        result = chk_circular_subgraph(G, subgraph, dict_pair_strand)
        assert result[4] is False  # not cyclic

    def test_three_node_cycle(self):
        G = nx.MultiDiGraph()
        G.add_edge("A", "B")
        G.add_edge("B", "C")
        G.add_edge("C", "A")
        subgraph = G.to_undirected().subgraph(["A", "B", "C"])
        dict_pair_strand = {}
        result = chk_circular_subgraph(G, subgraph, dict_pair_strand)
        assert result[4] is True  # cyclic (triangle forms cycle)


# ============================================================================
# TestWriteEmptyOutput
# ============================================================================


class TestWriteEmptyOutput:
    def _make_core(self):
        core = SplitReadsCore.__new__(SplitReadsCore)
        core.logger = logging.getLogger("test")
        return core

    def test_creates_file(self, tmp_path):
        core = self._make_core()
        result_path = core._write_empty_output(tmp_path)
        assert result_path.exists()
        assert result_path.name == "eccDNA_final.txt"

    def test_correct_columns(self, tmp_path):
        core = self._make_core()
        result_path = core._write_empty_output(tmp_path)
        df = pd.read_csv(result_path, sep="\t")
        expected_columns = [
            "id", "merge_region", "merge_len", "num_region",
            "ctc", "numreads", "totalbase", "coverage",
        ]
        assert list(df.columns) == expected_columns
