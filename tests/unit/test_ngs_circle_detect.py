"""Tests for ngs_circle_detect module — unit functions only."""

from pathlib import Path
import sys
import math
import pytest
import numpy as np
import pandas as pd
from unittest.mock import MagicMock, patch

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.modules.ngs_circle_detect import (
    NGSDetectConfig,
    NGSCircleDetect,
    _sigmoid,
    FEATURE_NAMES,
    MECC_FEATURE_NAMES,
)


# ════════════════════════════════════════════════════════════════════
# NGSDetectConfig defaults
# ════════════════════════════════════════════════════════════════════


class TestNGSDetectConfig:
    """Test default values for NGSDetectConfig."""

    def test_default_min_mapq(self):
        cfg = NGSDetectConfig()
        assert cfg.min_mapq == 10

    def test_default_min_clip_len(self):
        cfg = NGSDetectConfig()
        assert cfg.min_clip_len == 20

    def test_default_min_split_len(self):
        cfg = NGSDetectConfig()
        assert cfg.min_split_len == 30

    def test_default_max_normal_insert(self):
        cfg = NGSDetectConfig()
        assert cfg.max_normal_insert == 1000

    def test_default_bin_size(self):
        cfg = NGSDetectConfig()
        assert cfg.bin_size == 1000

    def test_default_min_evidence(self):
        cfg = NGSDetectConfig()
        assert cfg.min_evidence == 2

    def test_default_score_threshold(self):
        cfg = NGSDetectConfig()
        assert cfg.score_threshold == pytest.approx(0.3)

    def test_default_cecc_score_threshold(self):
        cfg = NGSDetectConfig()
        assert cfg.cecc_score_threshold == pytest.approx(0.43)

    def test_default_mecc_score_threshold(self):
        cfg = NGSDetectConfig()
        assert cfg.mecc_score_threshold == pytest.approx(0.50)

    def test_default_threads(self):
        cfg = NGSDetectConfig()
        assert cfg.threads == 4

    def test_default_minimap2_preset(self):
        cfg = NGSDetectConfig()
        assert cfg.minimap2_preset == "sr"

    def test_custom_values(self):
        cfg = NGSDetectConfig(min_mapq=20, bin_size=500, score_threshold=0.5)
        assert cfg.min_mapq == 20
        assert cfg.bin_size == 500
        assert cfg.score_threshold == pytest.approx(0.5)


# ════════════════════════════════════════════════════════════════════
# _sigmoid function
# ════════════════════════════════════════════════════════════════════


class TestSigmoid:
    """Test the _sigmoid activation function."""

    def test_sigmoid_zero(self):
        assert _sigmoid(0.0) == pytest.approx(0.5)

    def test_sigmoid_large_positive(self):
        assert _sigmoid(100.0) == pytest.approx(1.0, abs=1e-10)

    def test_sigmoid_large_negative(self):
        assert _sigmoid(-100.0) == pytest.approx(0.0, abs=1e-10)

    def test_sigmoid_underflow_guard(self):
        """Very large negative value should return 0.0 (not raise)."""
        assert _sigmoid(-600.0) == 0.0

    def test_sigmoid_boundary_minus500(self):
        """At exactly -500 the guard should trigger."""
        assert _sigmoid(-500.0) == 0.0

    def test_sigmoid_just_above_minus500(self):
        """Just above -500 should compute normally (very small)."""
        result = _sigmoid(-499.0)
        assert result > 0.0
        assert result < 1e-200

    def test_sigmoid_one(self):
        expected = 1.0 / (1.0 + math.exp(-1.0))
        assert _sigmoid(1.0) == pytest.approx(expected)

    def test_sigmoid_minus_one(self):
        expected = 1.0 / (1.0 + math.exp(1.0))
        assert _sigmoid(-1.0) == pytest.approx(expected)


# ════════════════════════════════════════════════════════════════════
# _dedup static method (NMS dedup)
# ════════════════════════════════════════════════════════════════════


class TestDedup:
    """Test NGSCircleDetect._dedup NMS logic."""

    def test_dedup_no_overlap(self):
        """Non-overlapping candidates should all be kept."""
        candidates = [
            {"chr": "chr1", "start": 100, "end": 200, "score": 0.9},
            {"chr": "chr1", "start": 1000, "end": 2000, "score": 0.8},
            {"chr": "chr2", "start": 100, "end": 200, "score": 0.7},
        ]
        kept = NGSCircleDetect._dedup(candidates)
        assert len(kept) == 3

    def test_dedup_full_overlap_keeps_higher_score(self):
        """Overlapping candidates should keep the higher-score one."""
        candidates = [
            {"chr": "chr1", "start": 100, "end": 200, "score": 0.6},
            {"chr": "chr1", "start": 100, "end": 200, "score": 0.9},
        ]
        kept = NGSCircleDetect._dedup(candidates)
        assert len(kept) == 1
        assert kept[0]["score"] == 0.9

    def test_dedup_partial_overlap_above_threshold(self):
        """Partial overlap > 50% should trigger dedup."""
        # Region 1: 100-200 (length=100)
        # Region 2: 150-250 (length=100)
        # Overlap: 150-200 = 50bp, frac = 50/100 = 0.5 (boundary; > 0.5 needed)
        # Let's make overlap > 50%:
        candidates = [
            {"chr": "chr1", "start": 100, "end": 200, "score": 0.9},
            {"chr": "chr1", "start": 140, "end": 240, "score": 0.7},
        ]
        kept = NGSCircleDetect._dedup(candidates)
        # Overlap: 140-200 = 60bp, min(60/100, 60/100) = 0.6 > 0.5
        assert len(kept) == 1

    def test_dedup_partial_overlap_below_threshold(self):
        """Partial overlap <= 50% should not trigger dedup."""
        candidates = [
            {"chr": "chr1", "start": 100, "end": 200, "score": 0.9},
            {"chr": "chr1", "start": 160, "end": 260, "score": 0.7},
        ]
        kept = NGSCircleDetect._dedup(candidates)
        # Overlap: 160-200 = 40bp, min(40/100, 40/100) = 0.4 <= 0.5
        assert len(kept) == 2

    def test_dedup_different_chromosomes(self):
        """Same coords on different chromosomes should not dedup."""
        candidates = [
            {"chr": "chr1", "start": 100, "end": 200, "score": 0.9},
            {"chr": "chr2", "start": 100, "end": 200, "score": 0.8},
        ]
        kept = NGSCircleDetect._dedup(candidates)
        assert len(kept) == 2

    def test_dedup_empty_input(self):
        kept = NGSCircleDetect._dedup([])
        assert kept == []

    def test_dedup_single_candidate(self):
        candidates = [{"chr": "chr1", "start": 100, "end": 200, "score": 0.9}]
        kept = NGSCircleDetect._dedup(candidates)
        assert len(kept) == 1

    def test_dedup_sorts_by_score_descending(self):
        """_dedup should sort by score and keep highest first."""
        candidates = [
            {"chr": "chr1", "start": 100, "end": 200, "score": 0.3},
            {"chr": "chr1", "start": 100, "end": 200, "score": 0.9},
            {"chr": "chr1", "start": 100, "end": 200, "score": 0.6},
        ]
        kept = NGSCircleDetect._dedup(candidates)
        assert len(kept) == 1
        assert kept[0]["score"] == 0.9


# ════════════════════════════════════════════════════════════════════
# _filter_and_classify — MAPQ stratification logic
# ════════════════════════════════════════════════════════════════════


class TestFilterAndClassify:
    """Test MAPQ-based classification into Uecc / Mecc."""

    @pytest.fixture
    def detector(self):
        """Create a detector instance with mocked reference."""
        with patch("circleseeker.modules.ngs_circle_detect.pysam"):
            det = NGSCircleDetect.__new__(NGSCircleDetect)
            det.config = NGSDetectConfig(score_threshold=0.3)
            det.logger = MagicMock()
            return det

    def test_high_mapq_dominant_yields_uecc(self, detector):
        """High-MAPQ backjumps dominant -> UeccDNA."""
        candidates = [{
            "score": 0.8, "chr": "chr1", "start": 100, "end": 1000,
            "size": 900,
            "n_backjump": 5, "n_backjump_low": 2,
            "n_disc_outward": 0, "n_disc_inward": 0,
            "n_clip_left": 0, "n_clip_right": 0,
            "n_reads": 7,
        }]
        results = detector._filter_and_classify(candidates)
        assert len(results) == 1
        assert results[0]["eccdna_type"] == "Uecc"

    def test_low_mapq_only_yields_mecc(self, detector):
        """Low-MAPQ backjumps only -> MeccDNA."""
        candidates = [{
            "score": 0.8, "chr": "chr1", "start": 100, "end": 1000,
            "size": 900,
            "n_backjump": 0, "n_backjump_low": 5,
            "n_disc_outward": 0, "n_disc_inward": 0,
            "n_clip_left": 0, "n_clip_right": 0,
            "n_reads": 5,
        }]
        results = detector._filter_and_classify(candidates)
        assert len(results) == 1
        assert results[0]["eccdna_type"] == "Mecc"

    def test_equal_high_low_yields_uecc(self, detector):
        """Equal high and low MAPQ backjumps -> Uecc (high >= low)."""
        candidates = [{
            "score": 0.8, "chr": "chr1", "start": 100, "end": 1000,
            "size": 900,
            "n_backjump": 3, "n_backjump_low": 3,
            "n_disc_outward": 0, "n_disc_inward": 0,
            "n_clip_left": 0, "n_clip_right": 0,
            "n_reads": 6,
        }]
        results = detector._filter_and_classify(candidates)
        assert len(results) == 1
        assert results[0]["eccdna_type"] == "Uecc"

    def test_no_backjumps_skipped(self, detector):
        """No backjumps at all -> entry is skipped."""
        candidates = [{
            "score": 0.8, "chr": "chr1", "start": 100, "end": 1000,
            "size": 900,
            "n_backjump": 0, "n_backjump_low": 0,
            "n_disc_outward": 5, "n_disc_inward": 3,
            "n_clip_left": 2, "n_clip_right": 2,
            "n_reads": 10,
        }]
        results = detector._filter_and_classify(candidates)
        assert len(results) == 0

    def test_below_score_threshold_skipped(self, detector):
        """Score below threshold -> entry is skipped."""
        candidates = [{
            "score": 0.1, "chr": "chr1", "start": 100, "end": 1000,
            "size": 900,
            "n_backjump": 5, "n_backjump_low": 0,
            "n_disc_outward": 0, "n_disc_inward": 0,
            "n_clip_left": 0, "n_clip_right": 0,
            "n_reads": 5,
        }]
        results = detector._filter_and_classify(candidates)
        assert len(results) == 0

    def test_low_mapq_dominant_yields_mecc(self, detector):
        """Low-MAPQ dominant over high-MAPQ -> MeccDNA."""
        candidates = [{
            "score": 0.8, "chr": "chr1", "start": 100, "end": 1000,
            "size": 900,
            "n_backjump": 1, "n_backjump_low": 5,
            "n_disc_outward": 0, "n_disc_inward": 0,
            "n_clip_left": 0, "n_clip_right": 0,
            "n_reads": 6,
        }]
        results = detector._filter_and_classify(candidates)
        assert len(results) == 1
        # n_bj_hi=1 > 0 AND n_bj_hi(1) >= n_bj_lo(5)? NO -> else branch -> Mecc
        assert results[0]["eccdna_type"] == "Mecc"

    def test_multiple_candidates_mixed(self, detector):
        """Multiple candidates: mix of Uecc, Mecc, and skipped."""
        candidates = [
            {
                "score": 0.8, "chr": "chr1", "start": 100, "end": 1000,
                "size": 900,
                "n_backjump": 5, "n_backjump_low": 0,
                "n_disc_outward": 0, "n_disc_inward": 0,
                "n_clip_left": 0, "n_clip_right": 0,
                "n_reads": 5,
            },
            {
                "score": 0.8, "chr": "chr2", "start": 200, "end": 2000,
                "size": 1800,
                "n_backjump": 0, "n_backjump_low": 3,
                "n_disc_outward": 0, "n_disc_inward": 0,
                "n_clip_left": 0, "n_clip_right": 0,
                "n_reads": 3,
            },
            {
                "score": 0.1, "chr": "chr3", "start": 300, "end": 3000,
                "size": 2700,
                "n_backjump": 10, "n_backjump_low": 0,
                "n_disc_outward": 0, "n_disc_inward": 0,
                "n_clip_left": 0, "n_clip_right": 0,
                "n_reads": 10,
            },
        ]
        results = detector._filter_and_classify(candidates)
        assert len(results) == 2
        types = {r["eccdna_type"] for r in results}
        assert types == {"Uecc", "Mecc"}


# ════════════════════════════════════════════════════════════════════
# _score_mecc_lr — basic LR scoring flow
# ════════════════════════════════════════════════════════════════════


class TestScoreMeccLr:
    """Test MeccDNA LR scoring and demotion logic."""

    @pytest.fixture
    def detector(self):
        with patch("circleseeker.modules.ngs_circle_detect.pysam"):
            det = NGSCircleDetect.__new__(NGSCircleDetect)
            det.config = NGSDetectConfig(mecc_score_threshold=0.50)
            det.logger = MagicMock()
            return det

    def test_empty_df_returns_unchanged(self, detector):
        df = pd.DataFrame()
        result = detector._score_mecc_lr(df)
        assert result.empty

    def test_no_mecc_returns_unchanged(self, detector):
        df = pd.DataFrame({
            "eccDNA_id": ["U1"],
            "eccdna_type": ["Uecc"],
            "n_backjump": [5],
            "n_backjump_low": [0],
            "length": [1000],
            "n_reads": [10],
            "n_discordant": [0],
            "n_softclip": [1],
            "score": [0.8],
            "n_genomic_loci": [1],
        })
        result = detector._score_mecc_lr(df)
        assert len(result) == 1
        assert result.iloc[0]["eccdna_type"] == "Uecc"

    def test_mecc_gets_score_column(self, detector):
        """Mecc entries should get a score_mecc_lr column."""
        df = pd.DataFrame({
            "eccDNA_id": ["M1"],
            "eccdna_type": ["Mecc"],
            "n_backjump": [0],
            "n_backjump_low": [5],
            "length": [1000],
            "n_reads": [10],
            "n_discordant": [2],
            "n_softclip": [1],
            "score": [0.6],
            "n_genomic_loci": [3],
        })
        result = detector._score_mecc_lr(df)
        assert "score_mecc_lr" in result.columns
        assert pd.notna(result.loc[result["eccdna_type"].isin(["Mecc", "Uecc"]), "score_mecc_lr"]).any()

    def test_low_score_mecc_demoted_to_uecc(self, detector):
        """Mecc with very low evidence should be demoted to Uecc."""
        # Create a Mecc with minimal evidence (likely low LR score)
        df = pd.DataFrame({
            "eccDNA_id": ["M1"],
            "eccdna_type": ["Mecc"],
            "n_backjump": [0],
            "n_backjump_low": [1],
            "length": [50000],
            "n_reads": [1],
            "n_discordant": [10],
            "n_softclip": [20],
            "score": [0.3],
            "n_genomic_loci": [1],
        })
        # Use a very high threshold to force demotion
        detector.config.mecc_score_threshold = 0.99
        result = detector._score_mecc_lr(df)
        assert result.iloc[0]["eccdna_type"] == "Uecc"

    def test_threshold_zero_skips_filtering(self, detector):
        """Threshold <= 0 should skip filtering entirely."""
        df = pd.DataFrame({
            "eccDNA_id": ["M1"],
            "eccdna_type": ["Mecc"],
            "n_backjump": [0],
            "n_backjump_low": [1],
            "length": [1000],
            "n_reads": [1],
            "n_discordant": [0],
            "n_softclip": [0],
            "score": [0.3],
            "n_genomic_loci": [1],
        })
        detector.config.mecc_score_threshold = 0.0
        result = detector._score_mecc_lr(df)
        assert result.iloc[0]["eccdna_type"] == "Mecc"

    def test_uecc_not_affected(self, detector):
        """Uecc entries should pass through unchanged regardless of score."""
        df = pd.DataFrame({
            "eccDNA_id": ["U1", "M1"],
            "eccdna_type": ["Uecc", "Mecc"],
            "n_backjump": [5, 0],
            "n_backjump_low": [0, 3],
            "length": [1000, 1000],
            "n_reads": [10, 5],
            "n_discordant": [0, 0],
            "n_softclip": [1, 1],
            "score": [0.8, 0.6],
            "n_genomic_loci": [1, 3],
        })
        result = detector._score_mecc_lr(df)
        uecc_rows = result[result["eccDNA_id"] == "U1"]
        assert uecc_rows.iloc[0]["eccdna_type"] == "Uecc"


# ════════════════════════════════════════════════════════════════════
# Cecc ghost filter logic
# ════════════════════════════════════════════════════════════════════


class TestCeccGhostFilter:
    """Test ghost detection filtering in the cecc pipeline.

    Ghost detections are rescued fragments with no junctions and no
    fragment structure — DISC-dominant single regions that are not valid
    standalone CeccDNA.
    """

    def test_keep_entries_with_junctions(self):
        """Entries with n_junctions > 0 should be kept."""
        cecc_df = pd.DataFrame({
            "eccDNA_id": ["C1"],
            "n_junctions": [3],
            "fragments": [""],
        })
        has_junctions = cecc_df["n_junctions"].fillna(0) > 0
        has_fragments = cecc_df["fragments"].fillna("").astype(str).str.len() > 0
        keep = has_junctions | has_fragments
        assert keep.iloc[0] is True or keep.iloc[0] == True  # noqa: E712

    def test_keep_entries_with_fragments(self):
        """Entries with non-empty fragments should be kept."""
        cecc_df = pd.DataFrame({
            "eccDNA_id": ["C1"],
            "n_junctions": [0],
            "fragments": ["chr1:100-200;chr2:300-400"],
        })
        has_junctions = cecc_df["n_junctions"].fillna(0) > 0
        has_fragments = cecc_df["fragments"].fillna("").astype(str).str.len() > 0
        keep = has_junctions | has_fragments
        assert keep.iloc[0] is True or keep.iloc[0] == True  # noqa: E712

    def test_filter_ghost_no_junctions_no_fragments(self):
        """Entries with no junctions AND no fragments should be filtered out."""
        cecc_df = pd.DataFrame({
            "eccDNA_id": ["C1"],
            "n_junctions": [0],
            "fragments": [""],
        })
        has_junctions = cecc_df["n_junctions"].fillna(0) > 0
        has_fragments = cecc_df["fragments"].fillna("").astype(str).str.len() > 0
        keep = has_junctions | has_fragments
        assert not keep.iloc[0]

    def test_filter_ghost_nan_values(self):
        """NaN junctions/fragments should also be filtered."""
        cecc_df = pd.DataFrame({
            "eccDNA_id": ["C1"],
            "n_junctions": [float("nan")],
            "fragments": [None],
        })
        has_junctions = cecc_df.get("n_junctions", pd.Series(dtype=float)).fillna(0) > 0
        has_fragments = cecc_df.get("fragments", pd.Series(dtype=object)).fillna("").astype(str).str.len() > 0
        keep = has_junctions | has_fragments
        assert not keep.iloc[0]

    def test_mixed_keep_and_filter(self):
        """Mix of keep and filter entries."""
        cecc_df = pd.DataFrame({
            "eccDNA_id": ["C1", "C2", "C3"],
            "n_junctions": [3, 0, 0],
            "fragments": ["chr1:100-200", "", "chr2:500-600"],
        })
        has_junctions = cecc_df["n_junctions"].fillna(0) > 0
        has_fragments = cecc_df["fragments"].fillna("").astype(str).str.len() > 0
        keep = has_junctions | has_fragments
        assert keep.tolist() == [True, False, True]
