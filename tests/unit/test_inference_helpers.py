"""Unit tests for standalone helper functions in circleseeker.core.steps.inference."""

from pathlib import Path
import sys
import logging

import pytest

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.core.steps.inference import (
    _get_xecc_source_read_names,
    _append_reads_from_fasta,
)


@pytest.fixture()
def logger():
    return logging.getLogger("test")


# ---------------------------------------------------------------------------
# _get_xecc_source_read_names
# ---------------------------------------------------------------------------


class TestGetXeccSourceReadNames:
    """Tests for _get_xecc_source_read_names."""

    def test_basic_two_reads(self, tmp_path, logger):
        """Two distinct reads yield a set of two names."""
        fasta = tmp_path / "xecc.fasta"
        fasta.write_text(
            ">readA|rep1|100|2|circular\nACGT\n"
            ">readB|rep2|200|3|circular\nTTTT\n"
        )
        result = _get_xecc_source_read_names(fasta, logger)
        assert result == {"readA", "readB"}

    def test_duplicate_read_names_are_deduplicated(self, tmp_path, logger):
        """Same read name appearing twice yields one entry."""
        fasta = tmp_path / "xecc.fasta"
        fasta.write_text(
            ">readA|rep1|100|2|circular\nACGT\n"
            ">readA|rep2|200|3|circular\nTTTT\n"
        )
        result = _get_xecc_source_read_names(fasta, logger)
        assert result == {"readA"}

    def test_header_with_suffix(self, tmp_path, logger):
        """Suffix like __X2 after pipe-delimited fields is handled correctly."""
        fasta = tmp_path / "xecc.fasta"
        fasta.write_text(">readName|rep1|100|1|circular__X2\nACGT\n")
        result = _get_xecc_source_read_names(fasta, logger)
        assert result == {"readName"}

    def test_empty_file(self, tmp_path, logger):
        """An empty file returns an empty set."""
        fasta = tmp_path / "xecc.fasta"
        fasta.write_text("")
        result = _get_xecc_source_read_names(fasta, logger)
        assert result == set()

    def test_nonexistent_file(self, tmp_path, logger):
        """A non-existent path returns an empty set."""
        fasta = tmp_path / "does_not_exist.fasta"
        result = _get_xecc_source_read_names(fasta, logger)
        assert result == set()

    def test_malformed_header_no_pipe(self, tmp_path, logger):
        """Header without pipe uses the full header word as read name."""
        fasta = tmp_path / "xecc.fasta"
        fasta.write_text(">simpleReadName\nACGT\n")
        result = _get_xecc_source_read_names(fasta, logger)
        assert result == {"simpleReadName"}

    def test_header_with_spaces(self, tmp_path, logger):
        """Only the first whitespace-delimited token is used for parsing."""
        fasta = tmp_path / "xecc.fasta"
        fasta.write_text(">readName|rep1|100|1|circular extra_info comment\nACGT\n")
        result = _get_xecc_source_read_names(fasta, logger)
        assert result == {"readName"}

    def test_zero_size_file(self, tmp_path, logger):
        """A file that exists but has zero bytes returns an empty set."""
        fasta = tmp_path / "xecc.fasta"
        fasta.touch()
        assert fasta.stat().st_size == 0
        result = _get_xecc_source_read_names(fasta, logger)
        assert result == set()


# ---------------------------------------------------------------------------
# _append_reads_from_fasta
# ---------------------------------------------------------------------------


class TestAppendReadsFromFasta:
    """Tests for _append_reads_from_fasta."""

    def _make_source(self, tmp_path, content):
        src = tmp_path / "source.fasta"
        src.write_text(content)
        return src

    def _make_target(self, tmp_path, content=""):
        tgt = tmp_path / "target.fasta"
        if content:
            tgt.write_text(content)
        else:
            tgt.touch()
        return tgt

    def test_basic_append_one_of_three(self, tmp_path, logger):
        """Append 1 matching read out of 3; returns 1."""
        source = self._make_source(
            tmp_path,
            ">readA\nAAAA\n>readB\nBBBB\n>readC\nCCCC\n",
        )
        target = self._make_target(tmp_path)
        count = _append_reads_from_fasta(source, target, {"readA"}, logger)
        assert count == 1
        text = target.read_text()
        assert ">readA" in text
        assert ">readB" not in text

    def test_append_two_of_three(self, tmp_path, logger):
        """Append 2 matching reads out of 3; returns 2."""
        source = self._make_source(
            tmp_path,
            ">readA\nAAAA\n>readB\nBBBB\n>readC\nCCCC\n",
        )
        target = self._make_target(tmp_path)
        count = _append_reads_from_fasta(source, target, {"readA", "readC"}, logger)
        assert count == 2

    def test_no_matching_reads(self, tmp_path, logger):
        """No matching reads returns 0, target is unchanged."""
        source = self._make_source(tmp_path, ">readA\nAAAA\n")
        target = self._make_target(tmp_path, ">existing\nTTTT\n")
        original = target.read_bytes()
        count = _append_reads_from_fasta(source, target, {"noMatch"}, logger)
        assert count == 0
        assert target.read_bytes() == original

    def test_empty_read_names(self, tmp_path, logger):
        """Empty read_names set returns 0 immediately."""
        source = self._make_source(tmp_path, ">readA\nAAAA\n")
        target = self._make_target(tmp_path)
        count = _append_reads_from_fasta(source, target, set(), logger)
        assert count == 0

    def test_nonexistent_source(self, tmp_path, logger):
        """Non-existent source FASTA returns 0."""
        source = tmp_path / "missing.fasta"
        target = self._make_target(tmp_path)
        count = _append_reads_from_fasta(source, target, {"readA"}, logger)
        assert count == 0

    def test_target_already_has_content(self, tmp_path, logger):
        """Appending does not overwrite existing target content."""
        source = self._make_source(tmp_path, ">readA\nAAAA\n")
        target = self._make_target(tmp_path, ">existing\nTTTT\n")
        count = _append_reads_from_fasta(source, target, {"readA"}, logger)
        assert count == 1
        text = target.read_text()
        assert ">existing" in text
        assert ">readA" in text

    def test_target_not_ending_with_newline(self, tmp_path, logger):
        """A newline is inserted before appending when target lacks trailing newline."""
        source = self._make_source(tmp_path, ">readA\nAAAA\n")
        target = self._make_target(tmp_path, ">existing\nTTTT")
        count = _append_reads_from_fasta(source, target, {"readA"}, logger)
        assert count == 1
        raw = target.read_bytes()
        # After "TTTT" there should be a newline before ">readA"
        idx = raw.find(b">readA")
        assert idx > 0
        assert raw[idx - 1 : idx] == b"\n"

    def test_empty_target_file(self, tmp_path, logger):
        """Appending to a zero-byte target works normally."""
        source = self._make_source(tmp_path, ">readA\nAAAA\n")
        target = self._make_target(tmp_path)
        assert target.stat().st_size == 0
        count = _append_reads_from_fasta(source, target, {"readA"}, logger)
        assert count == 1
        assert target.read_text().startswith(">readA")

    def test_multiline_sequences_preserved(self, tmp_path, logger):
        """Multi-line FASTA sequences are fully preserved."""
        source = self._make_source(
            tmp_path,
            ">readA\nAAAA\nGGGG\nCCCC\n>readB\nTTTT\n",
        )
        target = self._make_target(tmp_path)
        count = _append_reads_from_fasta(source, target, {"readA"}, logger)
        assert count == 1
        text = target.read_text()
        assert "AAAA" in text
        assert "GGGG" in text
        assert "CCCC" in text
        assert ">readB" not in text

    def test_exact_match_no_partial(self, tmp_path, logger):
        """Read names are matched exactly, not as substrings."""
        source = self._make_source(
            tmp_path,
            ">read\nAAAA\n>readA\nBBBB\n>readABC\nCCCC\n",
        )
        target = self._make_target(tmp_path)
        count = _append_reads_from_fasta(source, target, {"readA"}, logger)
        assert count == 1
        text = target.read_text()
        assert ">readA\n" in text
        assert ">read\n" not in text
        assert ">readABC" not in text
