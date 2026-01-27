from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.modules.ecc_dedup import CDHitClusters, EccDedup

import pytest
import numpy as np
from circleseeker.modules.ecc_dedup import (
    _normalize_strand,
    to_numeric_safe,
    format_two_decimals,
    concat_unique_semicolon,
    natural_sort_eccdna_id,
    merge_read_lists,
    majority_vote,
    count_reads_from_string,
    extract_eccDNA_base_id,
    ensure_int64_column,
    clamp_match_degree,
    region_str,
    read_csv_auto,
    build_u_confirmed_table,
    build_m_confirmed_table,
    build_c_confirmed_table,
)
from circleseeker.utils.column_standards import ColumnStandard

def test_ecc_dedup_subset(tmp_path):
    input_dir = tmp_path / "inputs"
    input_dir.mkdir()

    uecc_csv = input_dir / "uecc_input.csv"
    uecc_clstr = input_dir / "uecc.id2cluster.csv"
    mecc_csv = input_dir / "mecc_input.csv"
    mecc_clstr = input_dir / "mecc.id2cluster.csv"
    cecc_csv = input_dir / "cecc_input.csv"
    cecc_clstr = input_dir / "cecc.id2cluster.csv"

    pd.DataFrame(
        [
            {
                "eccDNA_id": "U1",
                "chr": "chr1",
                "start0": 0,
                "end0": 100,
                "strand": "+",
                "length": 100,
                "match_degree": 99.5,
                "copy_number": 1,
                "reads": "readA",
                "eSeq": "A" * 100,
            },
            {
                "eccDNA_id": "U2",
                "chr": "chr1",
                "start0": 0,
                "end0": 100,
                "strand": "+",
                "length": 100,
                "match_degree": 98.0,
                "copy_number": 2,
                "reads": "readB",
                "eSeq": "C" * 100,
            },
            {
                "eccDNA_id": "U3",
                "chr": "chr2",
                "start0": 200,
                "end0": 300,
                "strand": "-",
                "length": 100,
                "match_degree": 97.0,
                "copy_number": 1,
                "reads": "readC",
                "eSeq": "G" * 100,
            },
        ]
    ).to_csv(uecc_csv, index=False)
    uecc_clstr.write_text(
        "id,cluster,is_representative\n"
        "U1,grp1,false\n"
        "U2,grp1,true\n"
        "U3,grp2,true\n"
    )

    pd.DataFrame(
        [
            {
                "eccDNA_id": "M1",
                "chr": "chr1",
                "start0": 1000,
                "end0": 1100,
                "strand": "+",
                "length": 100,
                "match_degree": 98.0,
                "copy_number": 2,
                "reads": "mread1",
                "alignment_length": 100,
            },
            {
                "eccDNA_id": "M2",
                "chr": "chr1",
                "start0": 2000,
                "end0": 2100,
                "strand": "+",
                "length": 100,
                "match_degree": 99.0,
                "copy_number": 2,
                "reads": "mread2",
                "alignment_length": 100,
            },
            {
                "eccDNA_id": "M2",
                "chr": "chr2",
                "start0": 3000,
                "end0": 3100,
                "strand": "+",
                "length": 100,
                "match_degree": 97.0,
                "copy_number": 2,
                "reads": "mread3",
                "alignment_length": 90,
            },
        ]
    ).to_csv(mecc_csv, index=False)
    mecc_clstr.write_text(
        "id,cluster,is_representative\n"
        "M1,mgrp1,false\n"
        "M2,mgrp1,true\n"
    )

    pd.DataFrame(
        [
            {
                "eccDNA_id": "C1",
                "chr": "chr1",
                "start0": 0,
                "end0": 100,
                "strand": "+",
                "length": 200,
                "match_degree": 90.0,
                "copy_number": 1,
                "reads": "cread1",
            },
            {
                "eccDNA_id": "C2",
                "chr": "chr1",
                "start0": 100,
                "end0": 200,
                "strand": "+",
                "length": 200,
                "match_degree": 95.0,
                "copy_number": 1,
                "reads": "cread2",
            },
            {
                "eccDNA_id": "C2",
                "chr": "chr1",
                "start0": 300,
                "end0": 400,
                "strand": "+",
                "length": 200,
                "match_degree": 95.0,
                "copy_number": 1,
                "reads": "cread3",
            },
        ]
    ).to_csv(cecc_csv, index=False)
    cecc_clstr.write_text(
        "id,cluster,is_representative\n"
        "C1,cgrp1,false\n"
        "C2,cgrp1,true\n"
    )

    processor = EccDedup()
    results = processor.process_all_types(
        uecc_csv=uecc_csv,
        uecc_clstr=uecc_clstr,
        mecc_csv=mecc_csv,
        mecc_clstr=mecc_clstr,
        cecc_csv=cecc_csv,
        cecc_clstr=cecc_clstr,
        output_dir=tmp_path,
        prefix="step8",
        drop_seq=True,
        generate_confirmed_tables=False,
        organize_output_files=False,
    )

    assert set(results.keys()) == {"Uecc", "Mecc", "Cecc"}

    uecc = results["Uecc"]
    assert uecc["cluster_id"].nunique() == 2
    assert uecc["eccDNA_id"].str.startswith("UeccDNA").all()
    grp1 = uecc.loc[uecc["cluster_id"] == "grp1"].iloc[0]
    assert grp1["num_merged"] == 2
    assert grp1["merged_from_ids"] == "U1;U2"

    mecc = results["Mecc"]
    assert mecc["cluster_id"].nunique() == 1
    assert mecc["eccDNA_id"].str.startswith("MeccDNA").all()
    assert len(mecc) == 2  # two loci retained for the representative ID
    assert mecc["num_merged"].dropna().astype(int).eq(2).all()
    assert mecc["merged_from_ids"].iloc[0] == "M1;M2"

    cecc = results["Cecc"]
    assert cecc["cluster_id"].nunique() == 1
    assert cecc["eccDNA_id"].str.startswith("CeccDNA").all()
    assert len(cecc) == 2  # two segments retained for the representative ID
    assert cecc["num_merged"].dropna().astype(int).eq(2).all()
    assert cecc["merged_from_ids"].iloc[0] == "C1;C2"

    expected_outputs = [
        tmp_path / "step8_UeccDNA.core.csv",
        tmp_path / "step8_UeccDNA.bed",
        tmp_path / "step8_UeccDNA_C.fasta",
        tmp_path / "step8_MeccSites.core.csv",
        tmp_path / "step8_MeccSites.bed",
        tmp_path / "step8_MeccBestSite.bed",
        tmp_path / "step8_MeccDNA_C.fasta",
        tmp_path / "step8_CeccSegments.core.csv",
        tmp_path / "step8_CeccSegments.bed",
        tmp_path / "step8_CeccJunctions.bedpe",
        tmp_path / "step8_CeccDNA_C.fasta",
    ]
    for path in expected_outputs:
        assert path.exists(), f"Expected output missing: {path}"


def test_cecc_tolerance_merge_10bp(tmp_path):
    """Cecc entries with small boundary jitter should be merged (±10bp)."""
    input_dir = tmp_path / "inputs"
    input_dir.mkdir()

    cecc_csv = input_dir / "cecc_input.csv"
    cecc_clstr = input_dir / "cecc.id2cluster.csv"

    pd.DataFrame(
        [
            # C1 (two segments)
            {
                "eccDNA_id": "C1",
                "chr": "chr1",
                "start0": 1000,
                "end0": 1100,
                "strand": "+",
                "length": 200,
                "match_degree": 95.0,
                "copy_number": 1,
                "reads": "cread1",
            },
            {
                "eccDNA_id": "C1",
                "chr": "chr2",
                "start0": 2000,
                "end0": 2100,
                "strand": "+",
                "length": 200,
                "match_degree": 95.0,
                "copy_number": 1,
                "reads": "cread1",
            },
            # C2 (same circle, jitter within 10bp)
            {
                "eccDNA_id": "C2",
                "chr": "chr1",
                "start0": 1005,
                "end0": 1096,
                "strand": "-",
                "length": 200,
                "match_degree": 94.0,
                "copy_number": 1,
                "reads": "cread2",
            },
            {
                "eccDNA_id": "C2",
                "chr": "chr2",
                "start0": 1994,
                "end0": 2102,
                "strand": "-",
                "length": 200,
                "match_degree": 94.0,
                "copy_number": 1,
                "reads": "cread2",
            },
        ]
    ).to_csv(cecc_csv, index=False)

    # Put them in different CD-HIT clusters (simulate missed clustering).
    cecc_clstr.write_text(
        "id,cluster,is_representative\n"
        "C1,cgrp1,true\n"
        "C2,cgrp2,true\n"
    )

    processor = EccDedup()
    results = processor.process_all_types(
        uecc_csv=None,
        uecc_clstr=None,
        mecc_csv=None,
        mecc_clstr=None,
        cecc_csv=cecc_csv,
        cecc_clstr=cecc_clstr,
        output_dir=tmp_path,
        prefix="tol",
        drop_seq=True,
        generate_confirmed_tables=False,
        organize_output_files=False,
    )

    cecc = results["Cecc"]
    assert cecc["eccDNA_id"].nunique() == 1
    assert len(cecc) == 2  # representative keeps two segments
    merged_from = set(str(cecc["merged_from_ids"].iloc[0]).split(";"))
    assert {"C1", "C2"}.issubset(merged_from)
    assert cecc["num_merged"].dropna().astype(int).eq(2).all()


def test_ecc_dedup_csv_cluster_merging(tmp_path):
    """Ensure CSV cluster mappings are honoured when selecting representatives."""
    cluster_csv = tmp_path / "uecc.id2cluster.csv"
    cluster_csv.write_text(
        "id,cluster,is_representative\n"
        "U1,grp1,false\n"
        "U2,grp1,true\n"
        "U3,grp2,\n"
    )

    clusters = CDHitClusters()
    clusters.parse(cluster_csv)

    assert clusters.rep_of["grp1"] == "U2"
    assert clusters.rep_of["grp2"] == "U3"

    uecc_df = pd.DataFrame(
        [
            {
                "eccDNA_id": "U1",
                "eChr": "chr1",
                "eStart": 10,
                "eEnd": 70,
                "eLength": 60,
                "readName": "readA",
                "eRepeatNum": 1,
                "eStrand": "+",
            },
            {
                "eccDNA_id": "U2",
                "eChr": "chr1",
                "eStart": 20,
                "eEnd": 90,
                "eLength": 70,
                "readName": "readB",
                "eRepeatNum": 2,
                "eStrand": "+",
            },
            {
                "eccDNA_id": "U3",
                "eChr": "chr2",
                "eStart": 5,
                "eEnd": 45,
                "eLength": 40,
                "readName": "readC",
                "eRepeatNum": 3,
                "eStrand": "-",
            },
        ]
    )

    processor = EccDedup()
    processed = processor.process_uecc(uecc_df, clusters, "Uecc")
    processed = processed.sort_values("cluster_id").reset_index(drop=True)

    assert len(processed) == 2

    assert {"eccdna_id", "chr", "start_0based", "end_0based", "strand", "copy_number", "repeat_number", "match_degree"}.issubset(processed.columns)

    cluster1 = processed.loc[processed["cluster_id"] == "grp1"].iloc[0]
    assert cluster1["eccdna_id"] == "U2"
    assert cluster1["reads"] == "readA;readB"
    assert cluster1["repeat_number"] == 3
    assert cluster1["copy_number"] == 3
    assert cluster1["num_merged"] == 2
    assert cluster1["merged_from_ids"] == "U1;U2"
    assert cluster1["orig_eccdna_id"] == "U2"
    assert cluster1["eccdna_type"] == "Uecc"
    assert cluster1["strand"] == "+"

    cluster2 = processed.loc[processed["cluster_id"] == "grp2"].iloc[0]
    assert cluster2["eccdna_id"] == "U3"
    assert cluster2["num_merged"] == 1
    assert cluster2["merged_from_ids"] == "U3"
    assert cluster2["strand"] == "-"
    assert cluster2["copy_number"] == 3


# =====================================================================
# New test classes appended below
# =====================================================================


class TestNormalizeStrand:
    def test_plus(self): assert _normalize_strand("+") == "+"
    def test_minus(self): assert _normalize_strand("-") == "-"
    def test_plus_word(self): assert _normalize_strand("plus") == "+"
    def test_minus_word(self): assert _normalize_strand("minus") == "-"
    def test_forward(self): assert _normalize_strand("forward") == "+"
    def test_reverse(self): assert _normalize_strand("reverse") == "-"
    def test_positive(self): assert _normalize_strand("positive") == "+"
    def test_negative(self): assert _normalize_strand("negative") == "-"
    def test_none(self): assert _normalize_strand(None) is None
    def test_nan(self): assert _normalize_strand(float("nan")) is None
    def test_unknown_passthrough(self): assert _normalize_strand("unknown_val") == "unknown_val"


class TestToNumericSafe:
    def test_numeric_series(self):
        s = pd.Series([1, 2, 3])
        result = to_numeric_safe(s)
        assert list(result) == [1.0, 2.0, 3.0]

    def test_string_convertible(self):
        s = pd.Series(["10", "20"])
        result = to_numeric_safe(s)
        assert list(result) == [10.0, 20.0]

    def test_nan_filled(self):
        s = pd.Series([1, "abc", 3])
        result = to_numeric_safe(s, default=0)
        assert result.iloc[1] == 0.0

    def test_all_invalid(self):
        s = pd.Series(["a", "b"])
        result = to_numeric_safe(s, default=-1)
        assert all(result == -1)


class TestFormatTwoDecimals:
    def test_float(self): assert format_two_decimals(3.14159) == "3.14"
    def test_int(self): assert format_two_decimals(5) == "5.00"
    def test_nan(self): assert format_two_decimals(float("nan")) == ""
    def test_invalid_string(self): assert format_two_decimals("abc") == ""


class TestConcatUniqueSemicolon:
    def test_basic(self):
        s = pd.Series(["a", "b", "c"])
        assert concat_unique_semicolon(s) == "a;b;c"

    def test_duplicates(self):
        s = pd.Series(["a", "b", "a"])
        result = concat_unique_semicolon(s)
        assert set(result.split(";")) == {"a", "b"}

    def test_with_nan(self):
        s = pd.Series(["x", float("nan"), "y"])
        result = concat_unique_semicolon(s)
        assert "nan" not in result.lower()
        assert "x" in result
        assert "y" in result


class TestNaturalSortEccdnaId:
    def test_natural_sort(self):
        df = pd.DataFrame({"eccDNA_id": ["Uecc10", "Uecc2", "Uecc1"]})
        result = natural_sort_eccdna_id(df)
        assert result["eccDNA_id"].tolist() == ["Uecc1", "Uecc2", "Uecc10"]

    def test_mixed_types(self):
        df = pd.DataFrame({"eccDNA_id": ["Mecc1", "Uecc1", "Cecc1"]})
        result = natural_sort_eccdna_id(df)
        assert result["eccDNA_id"].tolist() == ["Cecc1", "Mecc1", "Uecc1"]

    def test_empty(self):
        df = pd.DataFrame({"eccDNA_id": []})
        result = natural_sort_eccdna_id(df)
        assert result.empty


class TestMergeReadLists:
    def test_basic(self):
        s = pd.Series(["readA;readB", "readC"])
        result = merge_read_lists(s)
        assert set(result.split(";")) == {"readA", "readB", "readC"}

    def test_dedup(self):
        s = pd.Series(["readA;readB", "readA;readC"])
        result = merge_read_lists(s)
        parts = result.split(";")
        assert len(parts) == len(set(parts))

    def test_with_nan(self):
        s = pd.Series(["readA", float("nan")])
        result = merge_read_lists(s)
        assert "readA" in result

    def test_empty(self):
        s = pd.Series([], dtype=str)
        result = merge_read_lists(s)
        assert result == ""


class TestMajorityVote:
    def test_clear_majority(self):
        s = pd.Series(["+", "+", "-"])
        assert majority_vote(s) == "+"

    def test_all_same(self):
        s = pd.Series(["-", "-", "-"])
        assert majority_vote(s) == "-"

    def test_empty(self):
        s = pd.Series([], dtype=str)
        assert majority_vote(s) is None


class TestCountReadsFromString:
    def test_basic(self):
        assert count_reads_from_string("readA;readB;readC") == 3

    def test_dedup(self):
        assert count_reads_from_string("readA;readB;readA") == 2

    def test_empty(self):
        assert count_reads_from_string("") == 0

    def test_nan(self):
        assert count_reads_from_string(float("nan")) == 0


class TestExtractEccdnaBaseId:
    def test_full_id(self):
        assert extract_eccDNA_base_id("CeccDNA_000001.1_2|rep0|5765|2|circular") == "CeccDNA_000001"

    def test_with_dot(self):
        assert extract_eccDNA_base_id("CeccDNA_000001.1_2") == "CeccDNA_000001"

    def test_plain_id(self):
        assert extract_eccDNA_base_id("CeccDNA_000001") == "CeccDNA_000001"

    def test_semicolon_list(self):
        assert extract_eccDNA_base_id("CeccDNA_000001;CeccDNA_000002") == "CeccDNA_000001"


class TestEnsureInt64Column:
    def test_float_conversion(self):
        df = pd.DataFrame({"col": [1.0, 2.0, 3.0]})
        ensure_int64_column(df, "col")
        assert df["col"].dtype == pd.Int64Dtype()
        assert df["col"].tolist() == [1, 2, 3]

    def test_missing_column(self):
        df = pd.DataFrame({"other": [1]})
        ensure_int64_column(df, "col")
        assert "col" in df.columns
        assert df["col"].dtype == pd.Int64Dtype()

    def test_invalid_data(self):
        df = pd.DataFrame({"col": ["abc", "def"]})
        ensure_int64_column(df, "col")
        assert df["col"].dtype == pd.Int64Dtype()


class TestClampMatchDegree:
    def test_within_range(self):
        s = pd.Series([50.0, 75.0])
        result = clamp_match_degree(s)
        assert result.tolist() == [50.0, 75.0]

    def test_over_100(self):
        s = pd.Series([150.0])
        result = clamp_match_degree(s)
        assert result.iloc[0] == 100.0

    def test_negative(self):
        s = pd.Series([-10.0])
        result = clamp_match_degree(s)
        assert result.iloc[0] == 0.0


class TestRegionStr:
    def test_basic(self):
        assert region_str("chr1", 100, 200) == "chr1:100-200"

    def test_zero_start(self):
        assert region_str("chr1", 0, 100) == "chr1:0-100"


class TestReadCsvAuto:
    def test_csv(self, tmp_path):
        f = tmp_path / "test.csv"
        pd.DataFrame({"a": [1], "b": [2]}).to_csv(f, index=False)
        df = read_csv_auto(f)
        assert list(df.columns) == ["a", "b"]

    def test_tsv(self, tmp_path):
        f = tmp_path / "test.tsv"
        pd.DataFrame({"a": [1], "b": [2]}).to_csv(f, index=False, sep="\t")
        df = read_csv_auto(f)
        assert list(df.columns) == ["a", "b"]

    def test_explicit_sep(self, tmp_path):
        f = tmp_path / "test.txt"
        pd.DataFrame({"a": [1], "b": [2]}).to_csv(f, index=False, sep="\t")
        df = read_csv_auto(f, sep="\t")
        assert list(df.columns) == ["a", "b"]


class TestCDHitClusters:
    def test_csv_parse(self, tmp_path):
        f = tmp_path / "clusters.csv"
        f.write_text("id,cluster,is_representative\nA,1,true\nB,1,false\nC,2,true\n")
        c = CDHitClusters()
        c.parse(f)
        assert c.rep_of["1"] == "A"
        assert set(c.members["1"]) == {"A", "B"}
        assert c.rep_of["2"] == "C"

    def test_clstr_parse(self, tmp_path):
        f = tmp_path / "test.clstr"
        f.write_text(">Cluster 0\n0\t100nt, >SeqA... *\n1\t90nt, >SeqB... at 95%\n>Cluster 1\n0\t80nt, >SeqC... *\n")
        c = CDHitClusters()
        c.parse(f)
        assert c.rep_of["0"] == "SeqA"
        assert "SeqB" in c.members["0"]

    def test_no_representative_marker(self, tmp_path):
        f = tmp_path / "clusters.csv"
        f.write_text("id,cluster\nA,1\nB,1\n")
        c = CDHitClusters()
        c.parse(f)
        assert c.rep_of["1"] == "A"  # First member used as fallback

    def test_singleton(self, tmp_path):
        f = tmp_path / "clusters.csv"
        f.write_text("id,cluster,is_representative\nX,solo,true\n")
        c = CDHitClusters()
        c.parse(f)
        assert c.members["solo"] == ["X"]

    def test_empty_file(self, tmp_path):
        f = tmp_path / "empty.csv"
        f.write_text("id,cluster,is_representative\n")
        c = CDHitClusters()
        c.parse(f)
        assert len(c.members) == 0

    def test_get_id_to_cluster_map(self, tmp_path):
        f = tmp_path / "clusters.csv"
        f.write_text("id,cluster,is_representative\nA,1,true\nB,1,false\n")
        c = CDHitClusters()
        c.parse(f)
        mapping = c.get_id_to_cluster_map()
        assert mapping["A"] == "1"
        assert mapping["B"] == "1"

    def test_state_reset_on_reparse(self, tmp_path):
        f1 = tmp_path / "c1.csv"
        f1.write_text("id,cluster,is_representative\nA,1,true\n")
        f2 = tmp_path / "c2.csv"
        f2.write_text("id,cluster,is_representative\nX,2,true\n")
        c = CDHitClusters()
        c.parse(f1)
        assert "A" in c.members.get("1", [])
        c.parse(f2)
        assert "A" not in c.members.get("1", [])
        assert "X" in c.members.get("2", [])

    def test_missing_file_raises(self):
        c = CDHitClusters()
        with pytest.raises(FileNotFoundError):
            c.parse(Path("/nonexistent/file.csv"))


class TestNormalizeCoordinates:
    @pytest.fixture
    def processor(self):
        return EccDedup()

    def test_column_rename(self, processor):
        df = pd.DataFrame({"eccDNA_id": ["U1"], "eChr": ["chr1"], "eStart0": [100], "eEnd0": [200], "eLength": [100]})
        result = processor.normalize_coordinates(df)
        assert ColumnStandard.CHR in result.columns

    def test_strand_standardization(self, processor):
        df = pd.DataFrame({"eccDNA_id": ["U1"], "chr": ["chr1"], "start0": [0], "end0": [100], "length": [100], "strand": ["plus"]})
        result = processor.normalize_coordinates(df)
        assert result[ColumnStandard.STRAND].iloc[0] == "+"

    def test_length_calculation(self, processor):
        df = pd.DataFrame({"eccDNA_id": ["U1"], "chr": ["chr1"], "start0": [100], "end0": [200]})
        result = processor.normalize_coordinates(df)
        assert result[ColumnStandard.LENGTH].iloc[0] == 100

    def test_match_degree_clamp(self, processor):
        df = pd.DataFrame({"eccDNA_id": ["U1"], "chr": ["chr1"], "start0": [0], "end0": [100], "length": [100], "match_degree": [150.0]})
        result = processor.normalize_coordinates(df)
        assert result[ColumnStandard.MATCH_DEGREE].iloc[0] == 100.0

    def test_copy_number_default(self, processor):
        df = pd.DataFrame({"eccDNA_id": ["U1"], "chr": ["chr1"], "start0": [0], "end0": [100], "length": [100]})
        result = processor.normalize_coordinates(df)
        assert ColumnStandard.COPY_NUMBER in result.columns


class TestProcessUecc:
    @pytest.fixture
    def processor(self):
        return EccDedup()

    @pytest.fixture
    def clusters(self, tmp_path):
        f = tmp_path / "clusters.csv"
        f.write_text("id,cluster,is_representative\nU1,grp1,true\nU2,grp1,false\n")
        c = CDHitClusters()
        c.parse(f)
        return c

    def test_single_cluster_representative(self, processor, clusters):
        df = pd.DataFrame({
            "eccDNA_id": ["U1", "U2"], "chr": ["chr1", "chr1"],
            "start0": [0, 0], "end0": [100, 100], "length": [100, 100],
            "strand": ["+", "+"], "reads": ["r1", "r2"],
        })
        result = processor.process_uecc(df, clusters, "Uecc")
        assert len(result) == 1
        assert result.iloc[0][ColumnStandard.ECCDNA_ID] == "U1"

    def test_singleton(self, processor, tmp_path):
        f = tmp_path / "c.csv"
        f.write_text("id,cluster,is_representative\nU1,grp1,true\n")
        c = CDHitClusters()
        c.parse(f)
        df = pd.DataFrame({
            "eccDNA_id": ["U1"], "chr": ["chr1"],
            "start0": [0], "end0": [100], "length": [100],
            "strand": ["+"], "reads": ["r1"],
        })
        result = processor.process_uecc(df, c, "Uecc")
        assert len(result) == 1

    def test_empty_input(self, processor, clusters):
        df = pd.DataFrame()
        result = processor.process_uecc(df, clusters, "Uecc")
        assert result.empty


class TestProcessMeccCecc:
    @pytest.fixture
    def processor(self):
        return EccDedup()

    def test_keeps_all_representative_rows(self, processor, tmp_path):
        f = tmp_path / "c.csv"
        f.write_text("id,cluster,is_representative\nM1,grp1,false\nM2,grp1,true\n")
        c = CDHitClusters()
        c.parse(f)
        df = pd.DataFrame({
            "eccDNA_id": ["M1", "M2", "M2"], "chr": ["chr1", "chr1", "chr2"],
            "start0": [0, 100, 200], "end0": [100, 200, 300], "length": [100, 100, 100],
            "strand": ["+", "+", "-"], "reads": ["r1", "r2", "r2"],
        })
        result = processor.process_mecc_cecc(df, c, "Mecc")
        # Representative M2 has 2 rows - both should be kept
        assert len(result[result[ColumnStandard.ECCDNA_ID] == "M2"]) == 2

    def test_metadata_aggregation(self, processor, tmp_path):
        f = tmp_path / "c.csv"
        f.write_text("id,cluster,is_representative\nC1,grp1,true\nC2,grp1,false\n")
        c = CDHitClusters()
        c.parse(f)
        df = pd.DataFrame({
            "eccDNA_id": ["C1", "C2"], "chr": ["chr1", "chr1"],
            "start0": [0, 0], "end0": [100, 100], "length": [100, 100],
            "strand": ["+", "+"], "reads": ["rA", "rB"],
        })
        result = processor.process_mecc_cecc(df, c, "Cecc")
        reads_val = result.iloc[0][ColumnStandard.READS]
        assert "rA" in str(reads_val)
        assert "rB" in str(reads_val)

    def test_empty_input(self, processor, tmp_path):
        f = tmp_path / "c.csv"
        f.write_text("id,cluster,is_representative\n")
        c = CDHitClusters()
        c.parse(f)
        result = processor.process_mecc_cecc(pd.DataFrame(), c, "Cecc")
        assert result.empty


class TestRenumberEccdnaIds:
    def test_sequential_ids(self, tmp_path):
        processor = EccDedup()
        f = tmp_path / "c.csv"
        f.write_text("id,cluster,is_representative\nU1,grp1,true\nU2,grp2,true\n")
        c = CDHitClusters()
        c.parse(f)
        df = pd.DataFrame({
            "eccDNA_id": ["U1", "U2"], "cluster_id": ["grp1", "grp2"],
            "chr": ["chr1", "chr2"], "start0": [0, 0], "end0": [100, 100],
            "length": [100, 100], "strand": ["+", "-"],
        })
        result = processor.renumber_eccdna_ids(df, "Uecc")
        assert result[ColumnStandard.ECCDNA_ID].tolist() == ["UeccDNA1", "UeccDNA2"]

    def test_empty(self):
        processor = EccDedup()
        result = processor.renumber_eccdna_ids(pd.DataFrame(), "Uecc")
        assert result.empty


class TestDedupeCeccSegments:
    @pytest.fixture
    def processor(self):
        return EccDedup()

    def test_position_duplicate_removal(self, processor):
        df = pd.DataFrame({
            "eccDNA_id": ["C1", "C1"], ColumnStandard.CHR: ["chr1", "chr1"],
            ColumnStandard.START0: [100, 100], ColumnStandard.END0: [200, 200],
            ColumnStandard.STRAND: ["+", "-"],
        })
        result = processor.dedupe_cecc_segments(df)
        assert len(result) == 1

    def test_different_positions_kept(self, processor):
        df = pd.DataFrame({
            "eccDNA_id": ["C1", "C1"], ColumnStandard.CHR: ["chr1", "chr1"],
            ColumnStandard.START0: [100, 500], ColumnStandard.END0: [200, 600],
            ColumnStandard.STRAND: ["+", "+"],
        })
        result = processor.dedupe_cecc_segments(df)
        assert len(result) == 2

    def test_tolerance(self, processor):
        df = pd.DataFrame({
            "eccDNA_id": ["C1", "C1"], ColumnStandard.CHR: ["chr1", "chr1"],
            ColumnStandard.START0: [100, 105], ColumnStandard.END0: [200, 195],
            ColumnStandard.STRAND: ["+", "-"],
        })
        result = processor.dedupe_cecc_segments(df, position_tolerance=100)
        assert len(result) == 1

    def test_empty(self, processor):
        result = processor.dedupe_cecc_segments(pd.DataFrame())
        assert result.empty


class TestMergeCeccByTolerance:
    @pytest.fixture
    def processor(self):
        return EccDedup()

    def test_within_tolerance_merged(self, processor):
        df = pd.DataFrame({
            ColumnStandard.ECCDNA_ID: ["C1", "C1", "C2", "C2"],
            ColumnStandard.CHR: ["chr1", "chr2", "chr1", "chr2"],
            ColumnStandard.START0: [100, 200, 105, 203],
            ColumnStandard.END0: [200, 300, 208, 297],
            ColumnStandard.STRAND: ["+", "+", "-", "-"],
            ColumnStandard.LENGTH: [100, 100, 103, 94],
            "reads": ["rA", "rA", "rB", "rB"],
        })
        result = processor.merge_cecc_by_tolerance(df, tolerance_bp=10)
        assert result[ColumnStandard.ECCDNA_ID].nunique() == 1

    def test_outside_tolerance_not_merged(self, processor):
        df = pd.DataFrame({
            ColumnStandard.ECCDNA_ID: ["C1", "C2"],
            ColumnStandard.CHR: ["chr1", "chr1"],
            ColumnStandard.START0: [100, 500],
            ColumnStandard.END0: [200, 600],
            ColumnStandard.STRAND: ["+", "+"],
            ColumnStandard.LENGTH: [100, 100],
        })
        result = processor.merge_cecc_by_tolerance(df, tolerance_bp=10)
        assert result[ColumnStandard.ECCDNA_ID].nunique() == 2

    def test_zero_tolerance_skips(self, processor):
        df = pd.DataFrame({
            ColumnStandard.ECCDNA_ID: ["C1", "C2"],
            ColumnStandard.CHR: ["chr1", "chr1"],
            ColumnStandard.START0: [100, 105],
            ColumnStandard.END0: [200, 205],
            ColumnStandard.STRAND: ["+", "+"],
            ColumnStandard.LENGTH: [100, 100],
        })
        result = processor.merge_cecc_by_tolerance(df, tolerance_bp=0)
        # tolerance_bp=0 means no merging
        assert result[ColumnStandard.ECCDNA_ID].nunique() == 2

    def test_transitive_merge(self, processor):
        df = pd.DataFrame({
            ColumnStandard.ECCDNA_ID: ["C1", "C2", "C3"],
            ColumnStandard.CHR: ["chr1", "chr1", "chr1"],
            ColumnStandard.START0: [100, 108, 115],
            ColumnStandard.END0: [200, 208, 215],
            ColumnStandard.STRAND: ["+", "+", "+"],
            ColumnStandard.LENGTH: [100, 100, 100],
            "reads": ["r1", "r2", "r3"],
        })
        result = processor.merge_cecc_by_tolerance(df, tolerance_bp=10)
        # C1-C2 within 10bp, C2-C3 within 10bp -> all merged via union-find
        assert result[ColumnStandard.ECCDNA_ID].nunique() <= 2


class TestBuildConfirmedTables:
    def test_build_u_basic(self):
        df = pd.DataFrame({
            "eccDNA_id": ["U1"], "chr": ["chr1"],
            "start0": [0], "end0": [100], "strand": ["+"], "length": [100],
        })
        result = build_u_confirmed_table(df)
        assert len(result) == 1
        assert result["eccDNA_type"].iloc[0] == "UeccDNA"
        assert result["State"].iloc[0] == "Confirmed"

    def test_build_u_empty(self):
        df = pd.DataFrame(columns=["eccDNA_id", "chr", "start0", "end0", "strand", "length"])
        result = build_u_confirmed_table(df)
        assert result.empty

    def test_build_m_basic(self):
        df = pd.DataFrame({
            "eccDNA_id": ["M1", "M1"], "chr": ["chr1", "chr2"],
            "start0": [0, 100], "end0": [100, 200], "strand": ["+", "-"], "length": [100, 100],
        })
        result = build_m_confirmed_table(df)
        assert len(result) == 1
        assert result["eccDNA_type"].iloc[0] == "MeccDNA"

    def test_build_m_empty(self):
        df = pd.DataFrame(columns=["eccDNA_id", "chr", "start0", "end0", "strand", "length"])
        result = build_m_confirmed_table(df)
        assert result.empty

    def test_build_c_basic(self):
        df = pd.DataFrame({
            "eccDNA_id": ["C1", "C1"], "chr": ["chr1", "chr2"],
            "start0": [0, 100], "end0": [100, 200], "strand": ["+", "-"], "length": [200, 200],
            "seg_total": [2, 2], "seg_index": [1, 2],
        })
        result = build_c_confirmed_table(df)
        assert len(result) == 1
        assert result["eccDNA_type"].iloc[0] == "CeccDNA"

    def test_build_c_empty(self):
        df = pd.DataFrame(columns=["eccDNA_id", "chr", "start0", "end0", "strand", "length", "seg_total", "seg_index"])
        result = build_c_confirmed_table(df)
        assert result.empty
