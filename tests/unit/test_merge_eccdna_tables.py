"""Unit tests for merge_eccdna_tables module."""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import sys

# Add src to path
ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.modules.merge_eccdna_tables import (
    parse_region,
    reciprocal_overlap_ok,
    build_chr_index,
    find_redundant_simple,
    find_redundant_chimeric,
    prepare_inferred_table,
    renumber_eccdna,
    merge_eccdna_tables
)


class TestHelperFunctions:
    """Test helper functions."""
    
    def test_parse_region_valid(self):
        """Test parsing valid region strings."""
        assert parse_region("chr1:1000-2000") == ("chr1", 1000, 2000)
        assert parse_region("chrX:0-1000") == ("chrX", 0, 1000)
        assert parse_region("2:500-600") == ("2", 500, 600)
        assert parse_region("chrMT:100-200") == ("chrMT", 100, 200)
    
    def test_parse_region_invalid(self):
        """Test parsing invalid region strings."""
        with pytest.raises(ValueError, match="Bad region"):
            parse_region("invalid")
        
        with pytest.raises(ValueError, match="Bad region"):
            parse_region("chr1:1000")
        
        with pytest.raises(ValueError, match="Bad region"):
            parse_region("chr1-1000-2000")
    
    def test_reciprocal_overlap_ok_perfect_match(self):
        """Test perfect overlap detection."""
        assert reciprocal_overlap_ok(1000, 2000, 1000, 2000) is True
        assert reciprocal_overlap_ok(0, 100, 0, 100) is True
    
    def test_reciprocal_overlap_ok_no_overlap(self):
        """Test non-overlapping regions."""
        assert reciprocal_overlap_ok(1000, 2000, 3000, 4000) is False
        assert reciprocal_overlap_ok(1000, 2000, 2001, 3000) is False
    
    def test_reciprocal_overlap_ok_partial_overlap(self):
        """Test partial overlap scenarios."""
        # 50% overlap - should fail with 0.99 threshold
        assert reciprocal_overlap_ok(1000, 2000, 1500, 2500, thr=0.99) is False
        
        # 50% overlap - should pass with 0.5 threshold
        assert reciprocal_overlap_ok(1000, 2000, 1500, 2500, thr=0.5) is True
        
        # Near-complete overlap
        assert reciprocal_overlap_ok(1000, 2000, 1001, 1999, thr=0.95) is True
    
    def test_reciprocal_overlap_ok_with_tolerance(self):
        """Test overlap with coordinate tolerance."""
        # Slightly offset but within tolerance
        assert reciprocal_overlap_ok(1000, 2000, 995, 2005, thr=0.99, tol=10) is True
        
        # Outside tolerance
        assert reciprocal_overlap_ok(1000, 2000, 980, 2020, thr=0.99, tol=10) is False
        
        # Edge case: touching regions with tolerance
        assert reciprocal_overlap_ok(1000, 2000, 2000, 3000, thr=0.99, tol=10) is False
    
    def test_reciprocal_overlap_ok_zero_length(self):
        """Test handling of zero-length regions."""
        # Zero-length regions should not crash
        assert reciprocal_overlap_ok(1000, 1000, 1000, 2000) is False
        assert reciprocal_overlap_ok(1000, 2000, 1500, 1500) is False


class TestBuildChrIndex:
    """Test chromosome index building."""
    
    def test_build_chr_index_basic(self):
        """Test basic index building."""
        regions = pd.Series([
            "chr1:1000-2000",
            "chr1:3000-4000",
            "chr2:5000-6000",
            "chrX:100-200"
        ])
        
        idx = build_chr_index(regions)
        
        assert set(idx.keys()) == {"chr1", "chr2", "chrX"}
        assert len(idx["chr1"]) == 2
        assert len(idx["chr2"]) == 1
        assert len(idx["chrX"]) == 1
        
        # Check sorting
        chr1_entries = idx["chr1"]
        assert chr1_entries[0][0] < chr1_entries[1][0]  # Sorted by start
    
    def test_build_chr_index_invalid_regions(self):
        """Test index building with invalid regions."""
        regions = pd.Series([
            "chr1:1000-2000",
            "invalid_region",
            "chr2:3000-4000",
            None,
            "chr3:missing-end"
        ])
        
        idx = build_chr_index(regions)
        
        # Only valid regions should be included
        assert set(idx.keys()) == {"chr1", "chr2"}
        assert len(idx["chr1"]) == 1
        assert len(idx["chr2"]) == 1
    
    def test_build_chr_index_empty(self):
        """Test index building with empty series."""
        regions = pd.Series(dtype=str)
        idx = build_chr_index(regions)
        assert idx == {}


class TestFindRedundantSimple:
    """Test finding redundant simple eccDNA."""
    
    @pytest.fixture
    def confirmed_u_core(self):
        """Create sample confirmed UeccDNA core table."""
        return pd.DataFrame({
            "eccDNA_id": ["UeccDNA1", "UeccDNA2", "UeccDNA3"],
            "Regions": ["chr1:1000-2000", "chr2:3000-4000", "chrX:5000-6000"],
            "Length": [1000, 1000, 1000],
            "eccDNA_type": ["UeccDNA", "UeccDNA", "UeccDNA"],
            "state": ["Confirmed", "Confirmed", "Confirmed"]
        })
    
    def test_find_redundant_simple_exact_match(self, confirmed_u_core):
        """Test finding exact matches."""
        inferred = pd.DataFrame({
            "eccDNA_id": ["IUeccDNA1", "IUeccDNA2"],
            "chr": ["chr1", "chr3"],
            "start0": [1000, 7000],
            "end0": [2000, 8000]
        })
        
        redundant = find_redundant_simple(inferred, confirmed_u_core)
        
        assert redundant == {"IUeccDNA1"}  # Only chr1:1000-2000 matches
    
    def test_find_redundant_simple_with_tolerance(self, confirmed_u_core):
        """Test finding matches within tolerance."""
        inferred = pd.DataFrame({
            "eccDNA_id": ["IUeccDNA1", "IUeccDNA2", "IUeccDNA3"],
            "eChr": ["chr1", "chr1", "chr2"],
            "eStart0": [995, 900, 3005],
            "eEnd0": [2005, 2100, 3995]
        })
        
        redundant = find_redundant_simple(inferred, confirmed_u_core, tol=10)
        
        assert redundant == {"IUeccDNA1", "IUeccDNA3"}  # Within tolerance
        assert "IUeccDNA2" not in redundant  # Outside tolerance
    
    def test_find_redundant_simple_no_matches(self, confirmed_u_core):
        """Test when no redundant entries exist."""
        inferred = pd.DataFrame({
            "eccDNA_id": ["IUeccDNA1", "IUeccDNA2"],
            "eChr": ["chr10", "chr11"],
            "eStart0": [1000, 2000],
            "eEnd0": [2000, 3000]
        })
        
        redundant = find_redundant_simple(inferred, confirmed_u_core)
        assert redundant == set()
    
    def test_find_redundant_simple_region_parsing(self, confirmed_u_core):
        """Test redundancy detection using regions column."""
        inferred = pd.DataFrame({
            "eccDNA_id": ["IUeccDNA1", "IUeccDNA2"],
            "Regions": ["chr1:1000-2000", "chr4:8000-9000"]
        })
        
        redundant = find_redundant_simple(inferred, confirmed_u_core)
        assert redundant == {"IUeccDNA1"}
    
    def test_find_redundant_simple_mixed_format(self, confirmed_u_core):
        """Test with mixed coordinate formats."""
        inferred = pd.DataFrame({
            "eccDNA_id": ["IUeccDNA1", "IUeccDNA2", "IUeccDNA3"],
            "eChr": ["chr1", None, "chr5"],
            "eStart0": [1000, None, 10000],
            "eEnd0": [2000, None, 11000],
            "Regions": [None, "chr2:3000-4000", None]
        })
        
        redundant = find_redundant_simple(inferred, confirmed_u_core)
        assert redundant == {"IUeccDNA1", "IUeccDNA2"}


class TestFindRedundantChimeric:
    """Test finding redundant chimeric eccDNA."""
    
    @pytest.fixture
    def confirmed_c_core(self):
        """Create sample confirmed CeccDNA core table."""
        return pd.DataFrame({
            "eccDNA_id": ["CeccDNA1", "CeccDNA2"],
            "Regions": [
                "chr1:1000-2000;chr2:3000-4000",
                "chr3:5000-6000;chr4:7000-8000;chr5:9000-10000"
            ],
            "Length": [3000, 6000],
            "eccDNA_type": ["CeccDNA", "CeccDNA"],
            "State": ["Confirmed", "Confirmed"]
        })
    
    def test_find_redundant_chimeric_exact_match(self, confirmed_c_core):
        """Test finding exact chimeric matches."""
        inferred = pd.DataFrame({
            "eccDNA_id": ["ICeccDNA1", "ICeccDNA2"],
            "Regions": [
                "chr1:1000-2000;chr2:3000-4000",  # Exact match
                "chr6:1000-2000;chr7:3000-4000"   # No match
            ]
        })
        
        redundant = find_redundant_chimeric(inferred, confirmed_c_core)
        assert redundant == {"ICeccDNA1"}
    
    def test_find_redundant_chimeric_no_coordinate_tolerance(self, confirmed_c_core):
        """Test that chimeric matching doesn't use coordinate tolerance."""
        inferred = pd.DataFrame({
            "eccDNA_id": ["ICeccDNA1"],
            "Regions": ["chr1:995-2005;chr2:2995-4005"]  # Slightly different
        })
        
        # Even with tolerance, should not match (exact string matching)
        redundant = find_redundant_chimeric(inferred, confirmed_c_core, tol=10)
        assert redundant == set()
    
    def test_find_redundant_chimeric_missing_regions(self, confirmed_c_core):
        """Test handling of missing regions."""
        inferred = pd.DataFrame({
            "eccDNA_id": ["ICeccDNA1", "ICeccDNA2"],
            "Regions": [None, pd.NA]
        })
        
        redundant = find_redundant_chimeric(inferred, confirmed_c_core)
        assert redundant == set()


class TestPrepareInferredTable:
    """Test inferred table preparation."""
    
    def test_prepare_inferred_simple(self):
        """Test preparing simple inferred table."""
        df = pd.DataFrame({
            "eccDNA_id": ["IUeccDNA1", "IUeccDNA2", "IUeccDNA3"],
            "eChr": ["chr1", "chr2", "chr3"],
            "eStart0": [1000, 2000, 3000],
            "eEnd0": [2000, 3000, 4000],
            "Length": [1000, 1000, 1000]
        })
        
        redundant = {"IUeccDNA2"}
        
        result = prepare_inferred_table(df, "UeccDNA", redundant)
        
        # Check filtering
        assert len(result) == 2
        assert "IUeccDNA2" not in result["eccDNA_id"].values
        
        # Check standard columns
        expected_cols = ["eccDNA_id", "Regions", "Length", "eccDNA_type", "State", "Seg_total", "Hit_count"]
        assert list(result.columns) == expected_cols
        
        # Check values
        assert all(result["eccDNA_type"] == "UeccDNA")
        assert all(result["State"] == "Inferred")
        assert all(result["Hit_count"] == 1)
        assert all(result["Seg_total"] == 1)
    
    def test_prepare_inferred_chimeric(self):
        """Test preparing chimeric inferred table."""
        df = pd.DataFrame({
            "eccDNA_id": ["ICeccDNA1", "ICeccDNA2"],
            "Regions": [
                "chr1:1000-2000;chr2:3000-4000",
                "chr3:5000-6000;chr4:7000-8000"
            ],
            "Length": [3000, 3000],
            "Seg_total": [2, 2]
        })
        
        result = prepare_inferred_table(df, "CeccDNA", set())
        
        assert all(result["eccDNA_type"] == "CeccDNA")
        assert all(result["Seg_total"] == 2)
    
    def test_prepare_inferred_build_regions(self):
        """Test building regions column when missing."""
        df = pd.DataFrame({
            "eccDNA_id": ["IUeccDNA1"],
            "eChr": ["chr1"],
            "eStart0": [1000],
            "eEnd0": [2000],
            "Length": [1000]
        })
        
        result = prepare_inferred_table(df, "UeccDNA", set())
        
        assert "Regions" in result.columns
        assert result.iloc[0]["Regions"] == "chr1:1000-2000"
    
    def test_prepare_inferred_all_redundant(self):
        """Test when all entries are redundant."""
        df = pd.DataFrame({
            "eccDNA_id": ["IUeccDNA1", "IUeccDNA2"],
            "eChr": ["chr1", "chr2"],
            "eStart0": [1000, 2000],
            "eEnd0": [2000, 3000],
            "Length": [1000, 1000]
        })
        
        redundant = {"IUeccDNA1", "IUeccDNA2"}
        
        result = prepare_inferred_table(df, "UeccDNA", redundant)
        assert len(result) == 0


class TestRenumberEccDNA:
    """Test eccDNA renumbering."""
    
    def test_renumber_basic(self):
        """Test basic renumbering functionality."""
        df = pd.DataFrame({
            "eccDNA_id": ["U3", "M2", "C1", "U1", "M1", "U2"],
            "Regions": ["r1", "r2", "r3", "r4", "r5", "r6"],
            "Length": [100, 200, 300, 400, 500, 600],
            "eccDNA_type": ["UeccDNA", "MeccDNA", "CeccDNA", "UeccDNA", "MeccDNA", "UeccDNA"]
        })
        
        result = renumber_eccdna(df)
        
        # Check new column exists
        assert "new_eccDNA_id" in result.columns
        
        # Check ordering: U, M, C
        new_ids = result["new_eccDNA_id"].tolist()
        assert new_ids == ["UeccDNA1", "UeccDNA2", "UeccDNA3", "MeccDNA1", "MeccDNA2", "CeccDNA1"]
        
        # Check that new_eccDNA_id is first column
        assert result.columns[0] == "new_eccDNA_id"
    
    def test_renumber_single_class(self):
        """Test renumbering with only one eccDNA class."""
        df = pd.DataFrame({
            "eccDNA_id": ["U5", "U2", "U9"],
            "Regions": ["r1", "r2", "r3"],
            "Length": [100, 200, 300],
            "eccDNA_type": ["UeccDNA", "UeccDNA", "UeccDNA"]
        })
        
        result = renumber_eccdna(df)
        
        assert list(result["new_eccDNA_id"]) == ["UeccDNA1", "UeccDNA2", "UeccDNA3"]
    
    def test_renumber_empty(self):
        """Test renumbering empty DataFrame."""
        df = pd.DataFrame(columns=["eccDNA_id", "Regions", "Length", "eccDNA_type"])
        result = renumber_eccdna(df)
        assert len(result) == 0
        assert "new_eccDNA_id" in result.columns
    
    def test_renumber_preserves_data(self):
        """Test that renumbering preserves all data."""
        df = pd.DataFrame({
            "eccDNA_id": ["U1", "M1"],
            "Regions": ["chr1:100-200", "chr2:300-400"],
            "Length": [100, 100],
            "eccDNA_type": ["UeccDNA", "MeccDNA"],
            "extra_col": ["data1", "data2"]
        })
        
        result = renumber_eccdna(df)
        
        # All columns should be preserved
        assert "extra_col" in result.columns
        assert list(result["extra_col"]) == ["data1", "data2"]


class TestMergeEccDNATables:
    """Test the main merge function."""
    
    @pytest.fixture
    def confirmed_cores(self):
        """Create sample confirmed core tables."""
        u_core = pd.DataFrame({
            "eccDNA_id": ["UeccDNA1", "UeccDNA2"],
            "Regions": ["chr1:1000-2000", "chr2:3000-4000"],
            "Length": [1000, 1000],
            "eccDNA_type": ["UeccDNA", "UeccDNA"],
            "State": ["Confirmed", "Confirmed"],
            "Seg_total": [1, 1],
            "Hit_count": [1, 1]
        })
        
        m_core = pd.DataFrame({
            "eccDNA_id": ["MeccDNA1"],
            "Regions": ["chr3:5000-6000|chr3:7000-8000"],
            "Length": [3000],
            "eccDNA_type": ["MeccDNA"],
            "State": ["Confirmed"],
            "Seg_total": [1],
            "Hit_count": [2]
        })
        
        c_core = pd.DataFrame({
            "eccDNA_id": ["CeccDNA1"],
            "Regions": ["chr4:9000-10000;chr5:11000-12000"],
            "Length": [2000],
            "eccDNA_type": ["CeccDNA"],
            "State": ["Confirmed"],
            "Seg_total": [2],
            "Hit_count": [1]
        })
        
        return {"U": u_core, "M": m_core, "C": c_core}
    
    def test_merge_confirmed_only(self, confirmed_cores):
        """Test merging only confirmed tables."""
        result = merge_eccdna_tables(confirmed_cores)
        
        # Check total count
        assert len(result) == 4  # 2 U + 1 M + 1 C
        
        # Check renumbering
        assert list(result["new_eccDNA_id"]) == ["UeccDNA1", "UeccDNA2", "MeccDNA1", "CeccDNA1"]
        
        # Check all are confirmed
        assert all(result["State"] == "Confirmed")
    
    def test_merge_with_inferred_simple(self, confirmed_cores):
        """Test merging with inferred simple eccDNA."""
        inferred_simple = pd.DataFrame({
            "eccDNA_id": ["IUeccDNA1", "IUeccDNA2", "IUeccDNA3"],
            "eChr": ["chr1", "chr6", "chr7"],
            "eStart0": [1000, 13000, 15000],
            "eEnd0": [2000, 14000, 16000],
            "Length": [1000, 1000, 1000]
        })
        
        result = merge_eccdna_tables(
            confirmed_cores,
            inferred_simple=inferred_simple
        )
        
        # IUeccDNA1 should be redundant (matches UeccDNA1)
        # So we should have 4 confirmed + 2 inferred = 6 total
        assert len(result) == 6
        
        # Check states
        confirmed_count = (result["State"] == "Confirmed").sum()
        inferred_count = (result["State"] == "Inferred").sum()
        assert confirmed_count == 4
        assert inferred_count == 2
    
    def test_merge_with_inferred_chimeric(self, confirmed_cores):
        """Test merging with inferred chimeric eccDNA."""
        inferred_chimeric = pd.DataFrame({
            "eccDNA_id": ["ICeccDNA1", "ICeccDNA2"],
            "Regions": [
                "chr4:9000-10000;chr5:11000-12000",  # Matches CeccDNA1
                "chr8:17000-18000;chr9:19000-20000"  # New
            ],
            "Length": [2000, 2000],
            "Seg_total": [2, 2]
        })
        
        result = merge_eccdna_tables(
            confirmed_cores,
            inferred_chimeric=inferred_chimeric
        )
        
        # ICeccDNA1 should be redundant
        assert len(result) == 5  # 4 confirmed + 1 inferred
    
    def test_merge_all_types(self, confirmed_cores):
        """Test merging all types together."""
        inferred_simple = pd.DataFrame({
            "eccDNA_id": ["IUeccDNA1"],
            "eChr": ["chr10"],
            "eStart0": [20000],
            "eEnd0": [21000],
            "Length": [1000]
        })
        
        inferred_chimeric = pd.DataFrame({
            "eccDNA_id": ["ICeccDNA1"],
            "Regions": ["chr11:22000-23000;chr12:24000-25000"],
            "Length": [2000],
            "Seg_total": [2]
        })
        
        result = merge_eccdna_tables(
            confirmed_cores,
            inferred_simple=inferred_simple,
            inferred_chimeric=inferred_chimeric,
            overlap_threshold=0.95,
            tolerance=5
        )
        
        assert len(result) == 6  # 4 confirmed + 2 inferred
        
        # Check final ordering
        classes = result["eccDNA_type"].tolist()
        assert classes == ["UeccDNA", "UeccDNA", "UeccDNA", "MeccDNA", "CeccDNA", "CeccDNA"]
    
    def test_merge_no_renumber(self, confirmed_cores):
        """Test merging without renumbering."""
        result = merge_eccdna_tables(confirmed_cores, renumber=False)
        
        assert "new_eccDNA_id" not in result.columns
        assert len(result) == 4
        
        # Original IDs should be preserved
        assert "UeccDNA1" in result["eccDNA_id"].values
        assert "MeccDNA1" in result["eccDNA_id"].values
    
    def test_merge_empty_confirmed(self):
        """Test error handling with empty confirmed tables."""
        with pytest.raises(ValueError, match="No data to merge"):
            merge_eccdna_tables({})
    
    def test_merge_partial_confirmed(self):
        """Test with only some confirmed types."""
        cores = {
            "U": pd.DataFrame({
                "eccDNA_id": ["UeccDNA1"],
                "Regions": ["chr1:100-200"],
                "Length": [100],
                "eccDNA_type": ["UeccDNA"],
                "State": ["Confirmed"],
                "Seg_total": [1],
                "Hit_count": [1]
            })
        }
        
        result = merge_eccdna_tables(cores)
        assert len(result) == 1
        assert result.iloc[0]["new_eccDNA_id"] == "UeccDNA1"
    
    def test_merge_custom_thresholds(self, confirmed_cores):
        """Test with custom overlap thresholds."""
        inferred_simple = pd.DataFrame({
            "eccDNA_id": ["IUeccDNA1"],
            "eChr": ["chr1"],
            "eStart0": [990],  # Slightly offset
            "eEnd0": [2010],
            "Length": [1020]
        })
        
        # With strict threshold, should not be redundant
        result_strict = merge_eccdna_tables(
            confirmed_cores,
            inferred_simple=inferred_simple,
            overlap_threshold=0.99,
            tolerance=5
        )
        
        # With relaxed threshold, should be redundant
        result_relaxed = merge_eccdna_tables(
            confirmed_cores,
            inferred_simple=inferred_simple,
            overlap_threshold=0.90,
            tolerance=15
        )
        
        assert len(result_strict) > len(result_relaxed)


class TestIntegration:
    """Integration tests for the complete workflow."""
    
    def test_complete_workflow(self):
        """Test a complete eccDNA merging workflow."""
        # Create confirmed cores
        confirmed = {
            "U": pd.DataFrame({
                "eccDNA_id": ["UeccDNA1", "UeccDNA2", "UeccDNA3"],
                "Regions": ["chr1:1000-2000", "chr2:3000-4000", "chr3:5000-6000"],
                "Length": [1000, 1000, 1000],
                "eccDNA_type": ["UeccDNA", "UeccDNA", "UeccDNA"],
                "state": ["Confirmed", "Confirmed", "Confirmed"],
                "Seg_total": [1, 1, 1],
                "Hit_count": [1, 2, 3]
            }),
            "M": pd.DataFrame({
                "eccDNA_id": ["MeccDNA1", "MeccDNA2"],
                "Regions": [
                    "chr4:7000-8000|chr4:9000-10000",
                    "chr5:11000-12000|chr6:13000-14000"
                ],
                "Length": [3000, 3000],
                "eccDNA_type": ["MeccDNA", "MeccDNA"],
                "State": ["Confirmed", "Confirmed"],
                "Seg_total": [1, 1],
                "Hit_count": [1, 1]
            }),
            "C": pd.DataFrame({
                "eccDNA_id": ["CeccDNA1"],
                "Regions": ["chr7:15000-16000;chr8:17000-18000;chr9:19000-20000"],
                "Length": [3000],
                "eccDNA_type": ["CeccDNA"],
                "State": ["Confirmed"],
                "Seg_total": [3],
                "Hit_count": [1]
            })
        }
        
        # Create inferred tables with some redundancy
        inferred_simple = pd.DataFrame({
            "eccDNA_id": ["IU1", "IU2", "IU3", "IU4"],
            "eChr": ["chr1", "chr10", "chr11", "chr2"],
            "eStart0": [995, 21000, 23000, 3005],  
            "eEnd0": [2005, 22000, 24000, 3995],
            "Length": [1010, 1000, 1000, 990]
        })
        
        inferred_chimeric = pd.DataFrame({
            "eccDNA_id": ["IC1", "IC2"],
            "Regions": [
                "chr7:15000-16000;chr8:17000-18000;chr9:19000-20000",  # Exact match
                "chr12:25000-26000;chr13:27000-28000"  # New
            ],
            "Length": [3000, 2000],
            "Seg_total": [3, 2]
        })
        
        # Merge all
        result = merge_eccdna_tables(
            confirmed,
            inferred_simple=inferred_simple,
            inferred_chimeric=inferred_chimeric,
            overlap_threshold=0.99,
            tolerance=10
        )
        
        # Verify counts
        # Confirmed: 3 U + 2 M + 1 C = 6
        # Inferred redundant: IU1 (matches U1), IU4 (matches U2), IC1 (matches C1)
        # Inferred kept: IU2, IU3, IC2 = 3
        # Total: 6 + 3 = 9
        assert len(result) == 9
        
        # Verify ordering and renumbering
        u_entries = result[result["eccDNA_type"] == "UeccDNA"]
        m_entries = result[result["eccDNA_type"] == "MeccDNA"]
        c_entries = result[result["eccDNA_type"] == "CeccDNA"]
        
        assert len(u_entries) == 5  # 3 confirmed + 2 inferred
        assert len(m_entries) == 2  # 2 confirmed
        assert len(c_entries) == 2  # 1 confirmed + 1 inferred
        
        # Check new IDs are sequential
        assert list(u_entries["new_eccDNA_id"]) == ["UeccDNA1", "UeccDNA2", "UeccDNA3", "UeccDNA4", "UeccDNA5"]
        assert list(m_entries["new_eccDNA_id"]) == ["MeccDNA1", "MeccDNA2"]
        assert list(c_entries["new_eccDNA_id"]) == ["CeccDNA1", "CeccDNA2"]
        
        # Verify states
        confirmed_count = (result["State"] == "Confirmed").sum()
        inferred_count = (result["State"] == "Inferred").sum()
        assert confirmed_count == 6
        assert inferred_count == 3
    
    def test_edge_cases(self):
        """Test various edge cases."""
        # Single confirmed entry
        confirmed = {
            "U": pd.DataFrame({
                "eccDNA_id": ["U1"],
                "Regions": ["chr1:100-200"],
                "Length": [100],
                "eccDNA_type": ["UeccDNA"],
                "State": ["Confirmed"],
                "Seg_total": [1],
                "Hit_count": [1]
            })
        }
        
        # All inferred entries are redundant
        inferred_simple = pd.DataFrame({
            "eccDNA_id": ["IU1", "IU2"],
            "eChr": ["chr1", "chr1"],
            "eStart0": [100, 95],
            "eEnd0": [200, 205],
            "Length": [100, 110]
        })
        
        result = merge_eccdna_tables(
            confirmed,
            inferred_simple=inferred_simple,
            tolerance=10
        )
        
        # Should only have the confirmed entry
        assert len(result) == 1
        assert result.iloc[0]["State"] == "Confirmed"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])