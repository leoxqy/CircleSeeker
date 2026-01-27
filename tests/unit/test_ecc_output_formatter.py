from pathlib import Path
import sys

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.modules.ecc_output_formatter import (
    format_eccdna_id,
    renumber_with_padding,
    format_location,
    format_mecc_location,
    format_cecc_location,
    generate_summary_table,
    generate_regions_table,
    generate_reads_table,
    generate_uecc_bed,
    generate_mecc_bed,
    generate_cecc_bed,
    generate_cecc_bedpe,
    generate_fasta_files,
)


class TestFormatEccdnaId:
    def test_uecc_type(self):
        assert format_eccdna_id("Uecc", 1) == "UeccDNA0001"

    def test_mecc_type(self):
        assert format_eccdna_id("Mecc", 5) == "MeccDNA0005"

    def test_cecc_type(self):
        assert format_eccdna_id("Cecc", 10) == "CeccDNA0010"

    def test_custom_width(self):
        assert format_eccdna_id("Uecc", 1, width=6) == "UeccDNA000001"


class TestRenumberWithPadding:
    def test_single_type(self):
        df = pd.DataFrame({
            "eccDNA_id": ["old1", "old2"],
            "eccDNA_type": ["UeccDNA", "UeccDNA"],
            "State": ["Confirmed", "Confirmed"],
        })
        result = renumber_with_padding(df, width=4)
        assert result["eccDNA_id"].tolist() == ["UeccDNA0001", "UeccDNA0002"]

    def test_confirmed_before_inferred(self):
        df = pd.DataFrame({
            "eccDNA_id": ["inf1", "conf1"],
            "eccDNA_type": ["UeccDNA", "UeccDNA"],
            "State": ["Inferred", "Confirmed"],
        })
        result = renumber_with_padding(df, width=4)
        # Confirmed should come first
        assert result.iloc[0]["State"] == "Confirmed"
        assert result.iloc[0]["eccDNA_id"] == "UeccDNA0001"

    def test_type_sorting(self):
        df = pd.DataFrame({
            "eccDNA_id": ["c1", "u1", "m1"],
            "eccDNA_type": ["CeccDNA", "UeccDNA", "MeccDNA"],
            "State": ["Confirmed", "Confirmed", "Confirmed"],
        })
        result = renumber_with_padding(df, width=4)
        assert result.iloc[0]["eccDNA_type"] == "UeccDNA"
        assert result.iloc[1]["eccDNA_type"] == "MeccDNA"
        assert result.iloc[2]["eccDNA_type"] == "CeccDNA"

    def test_preserves_old_id(self):
        df = pd.DataFrame({
            "eccDNA_id": ["original1"],
            "eccDNA_type": ["UeccDNA"],
            "State": ["Confirmed"],
        })
        result = renumber_with_padding(df)
        assert result["_old_eccDNA_id"].iloc[0] == "original1"

    def test_custom_width(self):
        df = pd.DataFrame({
            "eccDNA_id": ["x"],
            "eccDNA_type": ["UeccDNA"],
            "State": ["Confirmed"],
        })
        result = renumber_with_padding(df, width=6)
        assert result["eccDNA_id"].iloc[0] == "UeccDNA000001"


class TestFormatLocation:
    def test_basic(self):
        assert format_location("chr1", 100, 200, "+") == "chr1:100-200(+)"

    def test_negative_strand(self):
        assert format_location("chr2", 500, 600, "-") == "chr2:500-600(-)"


class TestFormatMeccLocation:
    def test_single_site(self):
        df = pd.DataFrame({
            "chr": ["chr1"], "start": [100], "end": [200], "strand": ["+"],
        })
        assert format_mecc_location(df) == "chr1:100-200(+)"

    def test_multiple_sites(self):
        df = pd.DataFrame({
            "chr": ["chr1", "chr2"], "start": [100, 300], "end": [200, 400], "strand": ["+", "-"],
        })
        result = format_mecc_location(df)
        assert "|" in result
        assert "chr1:100-200(+)" in result
        assert "chr2:300-400(-)" in result


class TestFormatCeccLocation:
    def test_two_segments(self):
        df = pd.DataFrame({
            "chr": ["chr1", "chr2"], "start": [100, 300], "end": [200, 400], "strand": ["+", "-"],
        })
        result = format_cecc_location(df)
        assert ";" in result
        assert "chr1:100-200(+)" in result

    def test_with_region_idx(self):
        df = pd.DataFrame({
            "chr": ["chr2", "chr1"], "start": [300, 100], "end": [400, 200],
            "strand": ["-", "+"], "region_idx": [2, 1],
        })
        result = format_cecc_location(df)
        # Should be sorted by region_idx
        parts = result.split(";")
        assert "chr1" in parts[0]

    def test_single_segment(self):
        df = pd.DataFrame({
            "chr": ["chr1"], "start": [100], "end": [200], "strand": ["+"],
        })
        result = format_cecc_location(df)
        assert ";" not in result


class TestGenerateSummaryTable:
    def test_uecc_entry(self):
        unified = pd.DataFrame({
            "eccDNA_id": ["U1"], "eccDNA_type": ["UeccDNA"],
            "State": ["Confirmed"], "Length": [100],
        })
        regions = pd.DataFrame({
            "eccDNA_id": ["U1"], "chr": ["chr1"], "start": [100], "end": [200],
            "strand": ["+"], "role": ["source"],
        })
        result = generate_summary_table(unified, regions)
        assert len(result) == 1
        assert result.iloc[0]["type"] == "Uecc"

    def test_mecc_entry(self):
        unified = pd.DataFrame({
            "eccDNA_id": ["M1"], "eccDNA_type": ["MeccDNA"],
            "State": ["Confirmed"], "Length": [100],
        })
        regions = pd.DataFrame({
            "eccDNA_id": ["M1", "M1"],
            "chr": ["chr1", "chr2"], "start": [100, 300], "end": [200, 400],
            "strand": ["+", "-"], "role": ["primary", "candidate"],
        })
        result = generate_summary_table(unified, regions)
        assert "|" in result.iloc[0]["location"]

    def test_cecc_entry(self):
        unified = pd.DataFrame({
            "eccDNA_id": ["C1"], "eccDNA_type": ["CeccDNA"],
            "State": ["Confirmed"], "Length": [200],
        })
        regions = pd.DataFrame({
            "eccDNA_id": ["C1", "C1"],
            "chr": ["chr1", "chr2"], "start": [100, 300], "end": [200, 400],
            "strand": ["+", "-"], "role": ["head", "tail"],
        })
        result = generate_summary_table(unified, regions)
        assert result.iloc[0]["segment_count"] == 2

    def test_empty_regions(self):
        unified = pd.DataFrame({
            "eccDNA_id": ["U1"], "eccDNA_type": ["UeccDNA"],
            "State": ["Confirmed"], "Length": [100],
        })
        regions = pd.DataFrame(columns=["eccDNA_id", "chr", "start", "end", "strand", "role"])
        result = generate_summary_table(unified, regions)
        assert result.iloc[0]["location"] == ""

    def test_mixed_types(self):
        unified = pd.DataFrame({
            "eccDNA_id": ["U1", "M1"], "eccDNA_type": ["UeccDNA", "MeccDNA"],
            "State": ["Confirmed", "Confirmed"], "Length": [100, 200],
        })
        regions = pd.DataFrame({
            "eccDNA_id": ["U1", "M1"],
            "chr": ["chr1", "chr2"], "start": [100, 300], "end": [200, 400],
            "strand": ["+", "+"], "role": ["source", "primary"],
        })
        result = generate_summary_table(unified, regions)
        assert len(result) == 2


class TestGenerateRegionsTable:
    def test_uecc_only(self):
        uecc = pd.DataFrame({
            "eccDNA_id": ["U1"], "chr": ["chr1"], "start0": [100], "end0": [200],
            "strand": ["+"], "length": [100],
        })
        result = generate_regions_table(uecc_df=uecc)
        assert len(result) == 1
        assert result.iloc[0]["role"] == "source"

    def test_mecc_multiple_sites(self):
        mecc = pd.DataFrame({
            "eccDNA_id": ["M1", "M1"], "chr": ["chr1", "chr2"],
            "start0": [100, 300], "end0": [200, 400],
            "strand": ["+", "-"], "length": [100, 100],
        })
        result = generate_regions_table(mecc_df=mecc)
        assert len(result) == 2
        assert "primary" in result["role"].values

    def test_cecc_head_middle_tail(self):
        cecc = pd.DataFrame({
            "eccDNA_id": ["C1", "C1", "C1"],
            "chr": ["chr1", "chr2", "chr3"],
            "start0": [100, 200, 300], "end0": [200, 300, 400],
            "strand": ["+", "+", "+"], "length": [100, 100, 100],
        })
        result = generate_regions_table(cecc_df=cecc)
        roles = result["role"].tolist()
        assert "head" in roles
        assert "tail" in roles
        assert "middle" in roles

    def test_all_none(self):
        result = generate_regions_table()
        assert result.empty

    def test_inferred_simple(self):
        inferred = pd.DataFrame({
            "eccDNA_id": ["IU1"], "chr": ["chr1"], "start0": [100], "end0": [200],
            "strand": ["+"], "length": [100],
        })
        result = generate_regions_table(inferred_simple_df=inferred)
        assert len(result) == 1
        assert result.iloc[0]["role"] == "source"


class TestGenerateReadsTable:
    def test_single_read(self):
        uecc = pd.DataFrame({
            "eccDNA_id": ["U1"], "read_name": ["readA"],
            "copy_number": [1], "match_degree": [99.0],
        })
        result = generate_reads_table(uecc_df=uecc)
        assert len(result) == 1
        assert result.iloc[0]["read_name"] == "readA"

    def test_semicolon_separated(self):
        uecc = pd.DataFrame({
            "eccDNA_id": ["U1"], "read_name": ["readA;readB"],
            "copy_number": [1], "match_degree": [99.0],
        })
        result = generate_reads_table(uecc_df=uecc)
        assert len(result) == 2

    def test_empty_input(self):
        result = generate_reads_table()
        assert result.empty


class TestGenerateUeccBed:
    def test_basic(self, tmp_path):
        regions = pd.DataFrame({
            "eccDNA_id": ["U1"], "chr": ["chr1"], "start": [100], "end": [200],
            "strand": ["+"], "role": ["source"],
        })
        out = tmp_path / "uecc.bed"
        generate_uecc_bed(regions, out)
        assert out.exists()
        content = out.read_text()
        assert "chr1" in content

    def test_empty(self, tmp_path):
        regions = pd.DataFrame(columns=["eccDNA_id", "chr", "start", "end", "strand", "role"])
        out = tmp_path / "uecc.bed"
        generate_uecc_bed(regions, out)
        assert not out.exists()


class TestGenerateMeccBed:
    def test_with_roles(self, tmp_path):
        regions = pd.DataFrame({
            "eccDNA_id": ["M1", "M1"],
            "chr": ["chr1", "chr2"], "start": [100, 300], "end": [200, 400],
            "strand": ["+", "-"], "role": ["primary", "candidate"],
        })
        out = tmp_path / "mecc.bed"
        generate_mecc_bed(regions, out)
        assert out.exists()

    def test_empty(self, tmp_path):
        regions = pd.DataFrame(columns=["eccDNA_id", "chr", "start", "end", "strand", "role"])
        out = tmp_path / "mecc.bed"
        generate_mecc_bed(regions, out)
        assert not out.exists()


class TestGenerateCeccBed:
    def test_multiple_segments(self, tmp_path):
        regions = pd.DataFrame({
            "eccDNA_id": ["C1", "C1"],
            "region_idx": [1, 2],
            "chr": ["chr1", "chr2"], "start": [100, 300], "end": [200, 400],
            "strand": ["+", "-"], "role": ["head", "tail"],
        })
        out = tmp_path / "cecc.bed"
        generate_cecc_bed(regions, out)
        assert out.exists()
        content = out.read_text()
        assert "seg1/2" in content

    def test_empty(self, tmp_path):
        regions = pd.DataFrame(columns=["eccDNA_id", "region_idx", "chr", "start", "end", "strand", "role"])
        out = tmp_path / "cecc.bed"
        generate_cecc_bed(regions, out)
        assert not out.exists()


class TestGenerateCeccBedpe:
    def test_two_junctions(self, tmp_path):
        regions = pd.DataFrame({
            "eccDNA_id": ["C1", "C1", "C1"],
            "region_idx": [1, 2, 3],
            "chr": ["chr1", "chr2", "chr3"],
            "start": [100, 300, 500], "end": [200, 400, 600],
            "strand": ["+", "+", "-"], "role": ["head", "middle", "tail"],
        })
        out = tmp_path / "cecc.bedpe"
        generate_cecc_bedpe(regions, out)
        assert out.exists()
        lines = out.read_text().strip().split("\n")
        assert len(lines) == 2  # 2 junctions for 3 segments

    def test_single_segment(self, tmp_path):
        regions = pd.DataFrame({
            "eccDNA_id": ["C1"], "region_idx": [1],
            "chr": ["chr1"], "start": [100], "end": [200],
            "strand": ["+"], "role": ["head"],
        })
        out = tmp_path / "cecc.bedpe"
        generate_cecc_bedpe(regions, out)
        # Single segment -> no junctions, file is created but empty
        content = out.read_text().strip()
        assert content == ""


class TestGenerateFastaFiles:
    def test_all_types(self, tmp_path):
        sequences = {"U1": "ATCG", "M1": "GGCC", "C1": "TTAA"}
        summary = pd.DataFrame({
            "eccDNA_id": ["U1", "M1", "C1"],
            "type": ["Uecc", "Mecc", "Cecc"],
            "location": ["chr1:1-4(+)", "chr2:1-4(+)", "chr3:1-4(+)"],
            "length": [4, 4, 4],
            "state": ["Confirmed", "Confirmed", "Confirmed"],
        })
        generate_fasta_files(sequences, tmp_path, summary)
        assert (tmp_path / "eccDNA_all.fasta").exists()
        assert (tmp_path / "Uecc" / "uecc.fasta").exists()
        assert (tmp_path / "Mecc" / "mecc.fasta").exists()
        assert (tmp_path / "Cecc" / "cecc.fasta").exists()

    def test_header_format(self, tmp_path):
        sequences = {"U1": "ATCG"}
        summary = pd.DataFrame({
            "eccDNA_id": ["U1"], "type": ["Uecc"],
            "location": ["chr1:1-4(+)"], "length": [4], "state": ["Confirmed"],
        })
        generate_fasta_files(sequences, tmp_path, summary)
        content = (tmp_path / "eccDNA_all.fasta").read_text()
        assert ">U1|" in content
        assert "length=4" in content

    def test_empty_sequences(self, tmp_path):
        sequences = {}
        summary = pd.DataFrame(columns=["eccDNA_id", "type", "location", "length", "state"])
        generate_fasta_files(sequences, tmp_path, summary)
        content = (tmp_path / "eccDNA_all.fasta").read_text()
        assert content == ""
