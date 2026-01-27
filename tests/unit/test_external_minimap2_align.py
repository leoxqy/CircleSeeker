"""Tests for minimap2 PAF to alignment TSV conversion."""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.external.minimap2_align import paf_to_alignment_tsv, _parse_cs_tag, _extract_cs_tag


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


# ---------------------------------------------------------------------------
# TestParseCsTag
# ---------------------------------------------------------------------------


class TestParseCsTag:
    """Tests for _parse_cs_tag()."""

    def test_match_only(self):
        """Match-only ':100' produces zero mismatches/gaps."""
        assert _parse_cs_tag(":100") == (0, 0, 0)

    def test_single_mismatch(self):
        """'*ac' counts as one mismatch, no gaps."""
        assert _parse_cs_tag("*ac") == (1, 0, 0)

    def test_insertion(self):
        """'+acgt' is one gap-open with 4 gap bases."""
        assert _parse_cs_tag("+acgt") == (0, 1, 4)

    def test_deletion(self):
        """'-acgt' is one gap-open with 4 gap bases."""
        assert _parse_cs_tag("-acgt") == (0, 1, 4)

    def test_splice(self):
        """'~gt1000ag' is one gap-open with 1000 gap bases."""
        assert _parse_cs_tag("~gt1000ag") == (0, 1, 1000)

    def test_mixed_operations(self):
        """':50*ac+tt-gg:30' has 1 mismatch, 2 gap-opens, 4 gap bases."""
        assert _parse_cs_tag(":50*ac+tt-gg:30") == (1, 2, 4)

    def test_empty_string(self):
        """Empty cs string returns zeros."""
        assert _parse_cs_tag("") == (0, 0, 0)


# ---------------------------------------------------------------------------
# TestExtractCsTag
# ---------------------------------------------------------------------------


class TestExtractCsTag:
    """Tests for _extract_cs_tag()."""

    def test_found(self):
        """cs:Z: prefix is stripped and value returned."""
        assert _extract_cs_tag(["cs:Z::100"]) == ":100"

    def test_not_found(self):
        """Returns None when no cs tag present."""
        assert _extract_cs_tag(["NM:i:5", "ms:i:200"]) is None

    def test_multiple_tags_cs_present(self):
        """Finds cs tag among several optional tags."""
        tags = ["NM:i:5", "ms:i:200", "cs:Z::50*ac+tt-gg:30", "de:f:0.01"]
        assert _extract_cs_tag(tags) == ":50*ac+tt-gg:30"


# ---------------------------------------------------------------------------
# TestPafToAlignmentTsv
# ---------------------------------------------------------------------------


def _make_paf_line(
    qname="q1",
    qlen="200",
    qstart="0",
    qend="155",
    strand="+",
    tname="ref1",
    tlen="1000",
    tstart="100",
    tend="255",
    nmatch="150",
    alen="155",
    mapq="60",
    extra_tags=None,
):
    """Helper to build a PAF line from keyword args."""
    fields = [
        qname, qlen, qstart, qend, strand,
        tname, tlen, tstart, tend, nmatch, alen, mapq,
    ]
    if extra_tags:
        fields.extend(extra_tags)
    return "\t".join(str(f) for f in fields)


class TestPafToAlignmentTsv:
    """Tests for paf_to_alignment_tsv()."""

    def test_empty_paf(self, tmp_path):
        """Empty PAF file produces 0 records."""
        paf = tmp_path / "empty.paf"
        out = tmp_path / "empty.tsv"
        paf.write_text("")
        assert paf_to_alignment_tsv(paf, out) == 0

    def test_malformed_line_insufficient_fields(self, tmp_path):
        """Line with fewer than 12 fields is skipped."""
        paf = tmp_path / "bad.paf"
        out = tmp_path / "bad.tsv"
        paf.write_text("only\tfive\tfields\there\tnow\n")
        assert paf_to_alignment_tsv(paf, out) == 0

    def test_multiple_alignments(self, tmp_path):
        """Multiple valid lines each produce a record."""
        paf = tmp_path / "multi.paf"
        out = tmp_path / "multi.tsv"
        lines = [
            _make_paf_line(qname="q1"),
            _make_paf_line(qname="q2"),
            _make_paf_line(qname="q3"),
        ]
        paf.write_text("\n".join(lines) + "\n")
        assert paf_to_alignment_tsv(paf, out) == 3

    def test_identity_filter_rejects_low(self, tmp_path):
        """Alignment with identity 96.77% filtered when min_identity=99.0."""
        paf = tmp_path / "id.paf"
        out = tmp_path / "id.tsv"
        # nmatch=150, alen=155 → 96.77%
        paf.write_text(_make_paf_line() + "\n")
        assert paf_to_alignment_tsv(paf, out, min_identity=99.0) == 0

    def test_identity_filter_with_decay(self, tmp_path):
        """Decay lowers the threshold for longer sequences."""
        paf = tmp_path / "decay.paf"
        out = tmp_path / "decay.tsv"
        # qlen=20000 → decay = (20/10)*0.5 = 1.0 → threshold = 99.0-1.0 = 98.0
        # identity = 150/155*100 = 96.77% — still below 98.0 → filtered
        paf.write_text(
            _make_paf_line(qlen="20000") + "\n"
        )
        assert paf_to_alignment_tsv(
            paf, out, min_identity=99.0, identity_decay_per_10kb=0.5
        ) == 0

    def test_identity_filter_decay_passes(self, tmp_path):
        """With enough decay a borderline alignment passes."""
        paf = tmp_path / "decay_pass.paf"
        out = tmp_path / "decay_pass.tsv"
        # nmatch=154, alen=155 → 99.35%
        # qlen=40000 → decay = (40/10)*2.0 = 8.0 → threshold = max(99.0-8.0, 97.0) = 97.0
        # 99.35% > 97.0 → passes
        paf.write_text(
            _make_paf_line(qlen="40000", nmatch="154", alen="155") + "\n"
        )
        assert paf_to_alignment_tsv(
            paf, out, min_identity=99.0, identity_decay_per_10kb=2.0
        ) == 1

    def test_identity_floor(self, tmp_path):
        """Threshold never drops below min_identity_floor."""
        paf = tmp_path / "floor.paf"
        out = tmp_path / "floor.tsv"
        # qlen=1000000 → decay = (1000/10)*10 = 1000 → threshold clamped to floor 98.5
        # nmatch=153, alen=155 → 98.71% > 98.5 → passes
        paf.write_text(
            _make_paf_line(qlen="1000000", nmatch="153", alen="155") + "\n"
        )
        assert paf_to_alignment_tsv(
            paf, out, min_identity=99.0, identity_decay_per_10kb=10.0, min_identity_floor=98.5
        ) == 1

    def test_de_f_tag_preferred(self, tmp_path):
        """de:f tag identity is preferred over nmatch/alen."""
        paf = tmp_path / "de.paf"
        out = tmp_path / "de.tsv"
        # de:f:0.005 → identity = 99.5% — passes 99.0 threshold
        # nmatch/alen alone would give 96.77%
        line = _make_paf_line(extra_tags=["de:f:0.005"])
        paf.write_text(line + "\n")
        assert paf_to_alignment_tsv(paf, out, min_identity=99.0) == 1
        fields = out.read_text().strip().split("\t")
        assert fields[2] == "99.50"

    def test_de_f_absent_fallback(self, tmp_path):
        """Without de:f tag, identity comes from nmatch/alen."""
        paf = tmp_path / "noDe.paf"
        out = tmp_path / "noDe.tsv"
        paf.write_text(_make_paf_line() + "\n")
        assert paf_to_alignment_tsv(paf, out) == 1
        fields = out.read_text().strip().split("\t")
        # 150/155 = 96.774...
        assert fields[2] == "96.77"

    def test_output_directory_created(self, tmp_path):
        """Non-existent output directory is created automatically."""
        paf = tmp_path / "dir.paf"
        out = tmp_path / "new_dir" / "sub" / "out.tsv"
        paf.write_text(_make_paf_line() + "\n")
        assert paf_to_alignment_tsv(paf, out) == 1
        assert out.exists()

    def test_blank_lines_skipped(self, tmp_path):
        """Blank lines in the PAF are silently ignored."""
        paf = tmp_path / "blank.paf"
        out = tmp_path / "blank.tsv"
        content = "\n\n" + _make_paf_line() + "\n\n"
        paf.write_text(content)
        assert paf_to_alignment_tsv(paf, out) == 1

    def test_plus_strand_coordinates(self, tmp_path):
        """Plus strand: s_start=tstart+1, s_end=tend, sstrand='plus'."""
        paf = tmp_path / "plus.paf"
        out = tmp_path / "plus.tsv"
        paf.write_text(_make_paf_line(strand="+", tstart="100", tend="255") + "\n")
        assert paf_to_alignment_tsv(paf, out) == 1
        fields = out.read_text().strip().split("\t")
        assert fields[8] == "101"  # s_start = tstart+1
        assert fields[9] == "255"  # s_end = tend
        assert fields[12] == "plus"

    def test_minus_strand_coordinates(self, tmp_path):
        """Minus strand: s_start=tend, s_end=tstart+1, sstrand='minus'."""
        paf = tmp_path / "minus.paf"
        out = tmp_path / "minus.tsv"
        paf.write_text(
            _make_paf_line(strand="-", tstart="500", tend="550") + "\n"
        )
        assert paf_to_alignment_tsv(paf, out) == 1
        fields = out.read_text().strip().split("\t")
        assert fields[8] == "550"  # s_start = tend
        assert fields[9] == "501"  # s_end = tstart+1
        assert fields[12] == "minus"

    def test_mapq_passed_through(self, tmp_path):
        """MAPQ value appears as the 14th field."""
        paf = tmp_path / "mq.paf"
        out = tmp_path / "mq.tsv"
        paf.write_text(_make_paf_line(mapq="42") + "\n")
        assert paf_to_alignment_tsv(paf, out) == 1
        fields = out.read_text().strip().split("\t")
        assert fields[13] == "42"

    def test_evalue_and_bitscore_always_zero(self, tmp_path):
        """evalue (field 11) and bit_score (field 12) are always '0'."""
        paf = tmp_path / "ev.paf"
        out = tmp_path / "ev.tsv"
        paf.write_text(_make_paf_line() + "\n")
        assert paf_to_alignment_tsv(paf, out) == 1
        fields = out.read_text().strip().split("\t")
        assert fields[10] == "0"
        assert fields[11] == "0"

    def test_mixed_pass_fail_identity(self, tmp_path):
        """Lines with mixed identities: only those passing threshold are kept."""
        paf = tmp_path / "mix.paf"
        out = tmp_path / "mix.tsv"
        lines = [
            # 150/155 = 96.77% — fail at 99%
            _make_paf_line(qname="low"),
            # 155/155 = 100% — pass
            _make_paf_line(qname="high", nmatch="155"),
            # 154/155 = 99.35% — pass
            _make_paf_line(qname="mid", nmatch="154"),
        ]
        paf.write_text("\n".join(lines) + "\n")
        assert paf_to_alignment_tsv(paf, out, min_identity=99.0) == 2
        text = out.read_text()
        assert "high" in text
        assert "mid" in text
        assert "low" not in text
