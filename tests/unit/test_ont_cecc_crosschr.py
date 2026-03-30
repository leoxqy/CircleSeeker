"""Tests for ont_cecc_crosschr module — unit functions only."""

from pathlib import Path
import sys
import math
import pytest
import numpy as np
import pandas as pd
from unittest.mock import MagicMock

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.modules.ont_cecc_crosschr import (
    parse_paf,
    extract_crosschr_junctions,
    detect_cecc_from_junctions,
    _sigmoid,
    _CECC_CROSSCHR_THRESHOLD,
    _CECC_CROSSCHR_FEATURE_NAMES,
    _CECC_CROSSCHR_COEFS,
    _CECC_CROSSCHR_SCALER_MEAN,
    _CECC_CROSSCHR_SCALER_STD,
    _CECC_CROSSCHR_INTERCEPT,
    DEFAULT_MIN_MAPQ,
    DEFAULT_MIN_ALIGN_LEN,
    DEFAULT_CLUSTER_DIST,
)


# ════════════════════════════════════════════════════════════════════
# _sigmoid
# ════════════════════════════════════════════════════════════════════


class TestSigmoid:
    """Test the _sigmoid function (ONT module copy)."""

    def test_sigmoid_zero(self):
        assert _sigmoid(0.0) == pytest.approx(0.5)

    def test_sigmoid_large_positive(self):
        assert _sigmoid(100.0) == pytest.approx(1.0, abs=1e-10)

    def test_sigmoid_underflow_guard(self):
        assert _sigmoid(-600.0) == 0.0

    def test_sigmoid_minus500_boundary(self):
        assert _sigmoid(-500.0) == 0.0

    def test_sigmoid_positive_value(self):
        expected = 1.0 / (1.0 + math.exp(-2.0))
        assert _sigmoid(2.0) == pytest.approx(expected)


# ════════════════════════════════════════════════════════════════════
# parse_paf
# ════════════════════════════════════════════════════════════════════


class TestParsePaf:
    """Test PAF file parsing."""

    def _write_paf(self, tmp_path, lines):
        """Helper to write PAF lines to a temp file."""
        paf = tmp_path / "test.paf"
        paf.write_text("\n".join(lines) + "\n")
        return paf

    def test_parse_single_alignment(self, tmp_path):
        # PAF format: qname qlen qstart qend strand tname tlen tstart tend nmatch blocklen mapq
        line = "read1\t5000\t0\t1000\t+\tchr1\t100000\t500\t1500\t950\t1000\t60"
        paf = self._write_paf(tmp_path, [line])
        reads = parse_paf(paf)

        assert "read1" in reads
        assert len(reads["read1"]) == 1
        seg = reads["read1"][0]
        assert seg["qname"] == "read1"
        assert seg["qlen"] == 5000
        assert seg["qstart"] == 0
        assert seg["qend"] == 1000
        assert seg["strand"] == "+"
        assert seg["tname"] == "chr1"
        assert seg["tlen"] == 100000
        assert seg["tstart"] == 500
        assert seg["tend"] == 1500
        assert seg["n_match"] == 950
        assert seg["block_len"] == 1000
        assert seg["mapq"] == 60

    def test_parse_multiple_alignments_same_read(self, tmp_path):
        lines = [
            "read1\t5000\t0\t1000\t+\tchr1\t100000\t500\t1500\t950\t1000\t60",
            "read1\t5000\t1500\t2500\t-\tchr2\t80000\t2000\t3000\t900\t1000\t30",
        ]
        paf = self._write_paf(tmp_path, lines)
        reads = parse_paf(paf)

        assert len(reads["read1"]) == 2
        assert reads["read1"][0]["tname"] == "chr1"
        assert reads["read1"][1]["tname"] == "chr2"
        assert reads["read1"][1]["strand"] == "-"

    def test_parse_multiple_reads(self, tmp_path):
        lines = [
            "read1\t5000\t0\t1000\t+\tchr1\t100000\t500\t1500\t950\t1000\t60",
            "read2\t3000\t0\t800\t+\tchr2\t80000\t100\t900\t780\t800\t40",
        ]
        paf = self._write_paf(tmp_path, lines)
        reads = parse_paf(paf)

        assert len(reads) == 2
        assert "read1" in reads
        assert "read2" in reads

    def test_parse_skip_comments(self, tmp_path):
        lines = [
            "# This is a comment",
            "read1\t5000\t0\t1000\t+\tchr1\t100000\t500\t1500\t950\t1000\t60",
        ]
        paf = self._write_paf(tmp_path, lines)
        reads = parse_paf(paf)

        assert len(reads) == 1

    def test_parse_skip_short_lines(self, tmp_path):
        lines = [
            "read1\t5000\t0",  # too few columns
            "read2\t5000\t0\t1000\t+\tchr1\t100000\t500\t1500\t950\t1000\t60",
        ]
        paf = self._write_paf(tmp_path, lines)
        reads = parse_paf(paf)

        assert "read1" not in reads
        assert "read2" in reads

    def test_parse_empty_file(self, tmp_path):
        paf = tmp_path / "empty.paf"
        paf.write_text("")
        reads = parse_paf(paf)
        assert len(reads) == 0


# ════════════════════════════════════════════════════════════════════
# extract_crosschr_junctions
# ════════════════════════════════════════════════════════════════════


class TestExtractCrosschrJunctions:
    """Test cross-chromosome junction extraction from parsed PAF data."""

    def test_simple_crosschr_junction(self):
        """A read with alignments on 2 chromosomes should produce junctions."""
        reads = {
            "read1": [
                {"qname": "read1", "qlen": 5000, "qstart": 0, "qend": 1000,
                 "strand": "+", "tname": "chr1", "tlen": 100000,
                 "tstart": 500, "tend": 1500, "n_match": 950, "block_len": 1000, "mapq": 60},
                {"qname": "read1", "qlen": 5000, "qstart": 2000, "qend": 3000,
                 "strand": "-", "tname": "chr2", "tlen": 80000,
                 "tstart": 2000, "tend": 3000, "n_match": 900, "block_len": 1000, "mapq": 40},
            ]
        }
        junctions, n_crosschr = extract_crosschr_junctions(reads, min_mapq=1, min_align_len=200)

        assert n_crosschr == 1
        assert len(junctions) == 1
        j = junctions[0]
        assert j[0] == "chr1"  # chr_a (sorted)
        assert j[3] == "chr2"  # chr_b
        assert j[6] == "read1"  # read_name
        assert j[7] == "split"  # type

    def test_same_chromosome_no_junction(self):
        """Alignments on same chromosome should not produce cross-chr junctions."""
        reads = {
            "read1": [
                {"qname": "read1", "qlen": 5000, "qstart": 0, "qend": 1000,
                 "strand": "+", "tname": "chr1", "tlen": 100000,
                 "tstart": 500, "tend": 1500, "n_match": 950, "block_len": 1000, "mapq": 60},
                {"qname": "read1", "qlen": 5000, "qstart": 2000, "qend": 3000,
                 "strand": "+", "tname": "chr1", "tlen": 100000,
                 "tstart": 5000, "tend": 6000, "n_match": 900, "block_len": 1000, "mapq": 40},
            ]
        }
        junctions, n_crosschr = extract_crosschr_junctions(reads, min_mapq=1, min_align_len=200)

        assert n_crosschr == 0
        assert len(junctions) == 0

    def test_low_mapq_filtered(self):
        """Alignments below min_mapq should be filtered."""
        reads = {
            "read1": [
                {"qname": "read1", "qlen": 5000, "qstart": 0, "qend": 1000,
                 "strand": "+", "tname": "chr1", "tlen": 100000,
                 "tstart": 500, "tend": 1500, "n_match": 950, "block_len": 1000, "mapq": 60},
                {"qname": "read1", "qlen": 5000, "qstart": 2000, "qend": 3000,
                 "strand": "-", "tname": "chr2", "tlen": 80000,
                 "tstart": 2000, "tend": 3000, "n_match": 900, "block_len": 1000, "mapq": 0},
            ]
        }
        junctions, n_crosschr = extract_crosschr_junctions(reads, min_mapq=5, min_align_len=200)

        assert n_crosschr == 0
        assert len(junctions) == 0

    def test_short_alignment_filtered(self):
        """Alignments shorter than min_align_len should be filtered."""
        reads = {
            "read1": [
                {"qname": "read1", "qlen": 5000, "qstart": 0, "qend": 1000,
                 "strand": "+", "tname": "chr1", "tlen": 100000,
                 "tstart": 500, "tend": 1500, "n_match": 950, "block_len": 1000, "mapq": 60},
                {"qname": "read1", "qlen": 5000, "qstart": 2000, "qend": 2100,
                 "strand": "+", "tname": "chr2", "tlen": 80000,
                 "tstart": 2000, "tend": 2100, "n_match": 95, "block_len": 100, "mapq": 40},
            ]
        }
        # Second alignment is only 100bp (below 200bp threshold)
        junctions, n_crosschr = extract_crosschr_junctions(reads, min_mapq=1, min_align_len=200)

        assert n_crosschr == 0
        assert len(junctions) == 0

    def test_three_chromosomes_pairwise_junctions(self):
        """A read with alignments on 3 chromosomes should produce pairwise junctions."""
        reads = {
            "read1": [
                {"qname": "read1", "qlen": 10000, "qstart": 0, "qend": 1000,
                 "strand": "+", "tname": "chr1", "tlen": 100000,
                 "tstart": 500, "tend": 1500, "n_match": 950, "block_len": 1000, "mapq": 60},
                {"qname": "read1", "qlen": 10000, "qstart": 2000, "qend": 3000,
                 "strand": "+", "tname": "chr2", "tlen": 80000,
                 "tstart": 2000, "tend": 3000, "n_match": 900, "block_len": 1000, "mapq": 50},
                {"qname": "read1", "qlen": 10000, "qstart": 5000, "qend": 6000,
                 "strand": "-", "tname": "chr3", "tlen": 60000,
                 "tstart": 4000, "tend": 5000, "n_match": 900, "block_len": 1000, "mapq": 40},
            ]
        }
        junctions, n_crosschr = extract_crosschr_junctions(reads, min_mapq=1, min_align_len=200)

        assert n_crosschr == 1
        # 3 chromosomes -> C(3,2) = 3 pairs (chr1-chr2, chr1-chr3, chr2-chr3)
        assert len(junctions) == 3

    def test_single_alignment_no_junction(self):
        """A read with only one alignment should not produce junctions."""
        reads = {
            "read1": [
                {"qname": "read1", "qlen": 5000, "qstart": 0, "qend": 1000,
                 "strand": "+", "tname": "chr1", "tlen": 100000,
                 "tstart": 500, "tend": 1500, "n_match": 950, "block_len": 1000, "mapq": 60},
            ]
        }
        junctions, n_crosschr = extract_crosschr_junctions(reads, min_mapq=1, min_align_len=200)

        assert n_crosschr == 0
        assert len(junctions) == 0

    def test_empty_reads(self):
        junctions, n_crosschr = extract_crosschr_junctions({}, min_mapq=1, min_align_len=200)
        assert n_crosschr == 0
        assert len(junctions) == 0


# ════════════════════════════════════════════════════════════════════
# detect_cecc_from_junctions — union-find clustering
# ════════════════════════════════════════════════════════════════════


class TestDetectCeccFromJunctions:
    """Test the union-find based CeccDNA detection."""

    def test_empty_junctions(self):
        """No junctions -> no results."""
        results = detect_cecc_from_junctions([])
        assert results == []

    def test_single_junction_pair_sufficient_reads(self):
        """A junction cluster with enough reads should produce a CeccDNA.

        Requires strong enough signal to pass the LR quality filter:
        many reads, multiple junctions in the same union-find component,
        high MAPQ, strand closure, and non-zero local coverage.

        All chr1 positions share the same bin (bin 1), but chr2 positions
        are spread across bins 5, 6, 7. This creates 3 junction bins that
        share node (chr1, 1) and thus merge into one component via union-find.
        """
        junctions = []
        # All chr1 positions in bin 1 (1000-1999)
        # chr2 positions spread across bins 5, 6, 7
        for i in range(10):
            junctions.append(
                ("chr1", 1500 + i, "+", "chr2", 5500 + i, "+",
                 f"readA{i}", "split", 60)
            )
        for i in range(10):
            junctions.append(
                ("chr1", 1600 + i, "+", "chr2", 6500 + i, "+",
                 f"readB{i}", "split", 60)
            )
        for i in range(10):
            junctions.append(
                ("chr1", 1700 + i, "+", "chr2", 7500 + i, "+",
                 f"readC{i}", "split", 60)
            )
        # Local coverage index for non-zero local_reads
        local_cov = {
            "chr1": {b: 30 for b in range(0, 10)},
            "chr2": {b: 30 for b in range(4, 12)},
        }
        results = detect_cecc_from_junctions(
            junctions,
            min_junction_reads=1,
            min_total_reads=2,
            min_mean_mapq=10,
            cluster_dist=1000,
            local_coverage_index=local_cov,
        )
        assert len(results) >= 1
        result = results[0]
        assert result["eccdna_type"] == "Cecc"
        assert result["state"] == "Confirmed"
        assert result["n_reads"] >= 2
        assert result["n_chromosomes"] >= 2

    def test_insufficient_reads_filtered(self):
        """Junctions with too few reads should be filtered."""
        junctions = [
            ("chr1", 1000, "+", "chr2", 5000, "-", "read1", "split", 60),
        ]
        results = detect_cecc_from_junctions(
            junctions,
            min_junction_reads=1,
            min_total_reads=5,  # require 5 reads but only 1
            min_mean_mapq=10,
            cluster_dist=1000,
        )
        assert len(results) == 0

    def test_low_mapq_component_filtered(self):
        """Components with mean MAPQ below threshold should be filtered."""
        junctions = [
            ("chr1", 1000, "+", "chr2", 5000, "-", "read1", "split", 5),
            ("chr1", 1010, "+", "chr2", 5010, "-", "read2", "split", 3),
        ]
        results = detect_cecc_from_junctions(
            junctions,
            min_junction_reads=1,
            min_total_reads=1,
            min_mean_mapq=20,  # mean mapq is ~4, below 20
            cluster_dist=1000,
        )
        assert len(results) == 0

    def test_union_find_merges_related_junctions(self):
        """Junctions sharing chromosomal bins should be merged into one component.

        Uses 3 chromosomes with reads spread across junction bins.
        Provides local_coverage_index for strong LR signal.
        """
        junctions = []
        # chr1-chr2 junctions: spread across 3 bin pairs
        for i in range(15):
            bin_offset = (i % 3) * 1000
            junctions.append(
                ("chr1", 1500 + bin_offset + (i // 3), "+",
                 "chr2", 5500 + bin_offset + (i // 3), "+",
                 f"readA{i}", "split", 60)
            )
        # chr1-chr3 junctions: share chr1 bins with the above
        for i in range(15):
            bin_offset = (i % 3) * 1000
            junctions.append(
                ("chr1", 1500 + bin_offset + (i // 3), "+",
                 "chr3", 8500 + bin_offset + (i // 3), "+",
                 f"readB{i}", "split", 55)
            )
        local_cov = {
            "chr1": {b: 20 for b in range(0, 10)},
            "chr2": {b: 20 for b in range(4, 12)},
            "chr3": {b: 20 for b in range(7, 15)},
        }
        results = detect_cecc_from_junctions(
            junctions,
            min_junction_reads=1,
            min_total_reads=2,
            min_mean_mapq=10,
            cluster_dist=1000,
            local_coverage_index=local_cov,
        )
        # All share chr1 bins -> should form one component with 3 chromosomes
        assert len(results) >= 1
        # The single component should have fragments on chr1, chr2, chr3
        total_chroms = sum(r["n_chromosomes"] for r in results)
        assert total_chroms >= 2

    def test_strand_closure_flag(self):
        """Verify strand_closure is computed (even if it does not filter)."""
        junctions = [
            ("chr1", 1000, "+", "chr2", 5000, "+", "read1", "split", 60),
            ("chr1", 1010, "+", "chr2", 5010, "+", "read2", "split", 55),
        ]
        results = detect_cecc_from_junctions(
            junctions,
            min_junction_reads=1,
            min_total_reads=1,
            min_mean_mapq=10,
            cluster_dist=1000,
        )
        if results:
            # strand_closure should be an integer (0 or 1)
            assert results[0]["strand_closure"] in (0, 1)


# ════════════════════════════════════════════════════════════════════
# LR quality-control filter (sigmoid + threshold)
# ════════════════════════════════════════════════════════════════════


class TestCrossChrLRFilter:
    """Test that the LR scoring/filtering logic works correctly."""

    def test_lr_score_range(self):
        """LR score should be between 0 and 1."""
        # Build a mock candidate with typical features
        features = {f: 0.0 for f in _CECC_CROSSCHR_FEATURE_NAMES}
        features["n_reads"] = 10
        features["n_junctions"] = 4
        features["n_fragments"] = 2
        features["n_chromosomes"] = 2
        features["mean_junction_mapq"] = 50
        features["strand_closure"] = 1
        features["log_length"] = 3.5
        features["ratio"] = 0.1
        features["reads_per_junction"] = 2.5
        features["local_reads"] = 100

        feats = np.array([features.get(f, 0) for f in _CECC_CROSSCHR_FEATURE_NAMES], dtype=float)
        normed = (feats - _CECC_CROSSCHR_SCALER_MEAN) / np.maximum(_CECC_CROSSCHR_SCALER_STD, 1e-6)
        logit = float(np.dot(normed, _CECC_CROSSCHR_COEFS) + _CECC_CROSSCHR_INTERCEPT)
        score = _sigmoid(logit)

        assert 0.0 <= score <= 1.0

    def test_high_quality_candidate_passes(self):
        """A high-quality candidate should pass the LR threshold."""
        features = {
            "n_reads": 30,
            "n_junctions": 8,
            "n_fragments": 2,
            "n_chromosomes": 2,
            "mean_junction_mapq": 60,
            "strand_closure": 1,
            "log_length": 4.0,
            "ratio": 0.2,
            "reads_per_junction": 4.0,
            "local_reads": 200,
        }

        feats = np.array([features.get(f, 0) for f in _CECC_CROSSCHR_FEATURE_NAMES], dtype=float)
        normed = (feats - _CECC_CROSSCHR_SCALER_MEAN) / np.maximum(_CECC_CROSSCHR_SCALER_STD, 1e-6)
        logit = float(np.dot(normed, _CECC_CROSSCHR_COEFS) + _CECC_CROSSCHR_INTERCEPT)
        score = _sigmoid(logit)

        assert score >= _CECC_CROSSCHR_THRESHOLD

    def test_low_quality_candidate_fails(self):
        """A low-quality candidate should fail the LR threshold."""
        features = {
            "n_reads": 1,
            "n_junctions": 1,
            "n_fragments": 2,
            "n_chromosomes": 2,
            "mean_junction_mapq": 5,
            "strand_closure": 0,
            "log_length": 2.0,
            "ratio": 0.001,
            "reads_per_junction": 1.0,
            "local_reads": 5000,
        }

        feats = np.array([features.get(f, 0) for f in _CECC_CROSSCHR_FEATURE_NAMES], dtype=float)
        normed = (feats - _CECC_CROSSCHR_SCALER_MEAN) / np.maximum(_CECC_CROSSCHR_SCALER_STD, 1e-6)
        logit = float(np.dot(normed, _CECC_CROSSCHR_COEFS) + _CECC_CROSSCHR_INTERCEPT)
        score = _sigmoid(logit)

        assert score < _CECC_CROSSCHR_THRESHOLD
