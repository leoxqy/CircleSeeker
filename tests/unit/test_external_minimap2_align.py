"""Tests for minimap2 PAF to alignment TSV conversion."""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.external.minimap2_align import paf_to_alignment_tsv


def test_paf_to_alignment_tsv_plus_strand_with_cs(tmp_path):
    paf_path = tmp_path / "alignments.paf"
    out_path = tmp_path / "alignments.tsv"

    paf_line = "\t".join(
        [
            "query1",
            "200",
            "0",
            "155",
            "+",
            "ref1",
            "1000",
            "100",
            "255",
            "150",
            "155",
            "60",
            "cs:Z::100*ac+tt-gg:50",
        ]
    )
    paf_path.write_text(f"{paf_line}\n")

    written = paf_to_alignment_tsv(paf_path, out_path)
    assert written == 1

    fields = out_path.read_text().strip().split("\t")
    assert len(fields) == 14
    assert fields[0] == "query1"
    assert fields[1] == "ref1"
    assert fields[2] == "96.77"
    assert fields[3] == "155"
    assert fields[4] == "1"
    assert fields[5] == "2"
    assert fields[6] == "1"
    assert fields[7] == "155"
    assert fields[8] == "101"
    assert fields[9] == "255"
    assert fields[10] == "0"
    assert fields[11] == "0"
    assert fields[12] == "plus"
    assert fields[13] == "60"


def test_paf_to_alignment_tsv_minus_strand_without_cs(tmp_path):
    paf_path = tmp_path / "alignments.paf"
    out_path = tmp_path / "alignments.tsv"

    paf_line = "\t".join(
        [
            "query2",
            "500",
            "10",
            "60",
            "-",
            "ref2",
            "2000",
            "500",
            "550",
            "40",
            "50",
            "40",
        ]
    )
    paf_path.write_text(f"{paf_line}\n")

    written = paf_to_alignment_tsv(paf_path, out_path)
    assert written == 1

    fields = out_path.read_text().strip().split("\t")
    assert len(fields) == 14
    assert fields[0] == "query2"
    assert fields[1] == "ref2"
    assert fields[2] == "80.00"
    assert fields[3] == "50"
    assert fields[4] == "10"
    assert fields[5] == "0"
    assert fields[6] == "11"
    assert fields[7] == "60"
    assert fields[8] == "550"
    assert fields[9] == "501"
    assert fields[10] == "0"
    assert fields[11] == "0"
    assert fields[12] == "minus"
    assert fields[13] == "40"
