from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import pytest

from circleseeker.utils.display import (
    supports_emoji,
    EmojiHandler,
    ConsoleFormatter,
    ProgressDisplay,
    E,
    format_stats,
)


# ---------------------------------------------------------------------------
# TestSupportsEmoji
# ---------------------------------------------------------------------------
class TestSupportsEmoji:
    """Tests for the supports_emoji() helper."""

    def test_returns_bool(self):
        result = supports_emoji()
        assert isinstance(result, bool)

    def test_ci_env_returns_false(self, monkeypatch):
        monkeypatch.setenv("CI", "true")
        assert supports_emoji() is False

    def test_term_dumb_returns_false(self, monkeypatch):
        monkeypatch.delenv("CI", raising=False)
        monkeypatch.setenv("TERM", "dumb")
        assert supports_emoji() is False


# ---------------------------------------------------------------------------
# TestEmojiHandler
# ---------------------------------------------------------------------------
class TestEmojiHandler:
    """Tests for EmojiHandler."""

    def test_get_returns_fallback_when_no_support(self, monkeypatch):
        # Force no emoji support via CI env
        monkeypatch.setenv("CI", "true")
        handler = EmojiHandler()
        assert handler.get("✅", "[OK]") == "[OK]"

    def test_emoji_map_contains_known_mappings(self):
        assert "✅" in EmojiHandler.EMOJI_MAP
        assert "❌" in EmojiHandler.EMOJI_MAP
        assert "📄" in EmojiHandler.EMOJI_MAP

    def test_unknown_emoji_returns_fallback(self, monkeypatch):
        monkeypatch.setenv("CI", "true")
        handler = EmojiHandler()
        # Use an emoji that is NOT in EMOJI_MAP
        result = handler.get("\U0001f9ca", "ICE")  # ice cube emoji, unlikely in map
        assert result == "ICE"


# ---------------------------------------------------------------------------
# TestConsoleFormatter
# ---------------------------------------------------------------------------
class TestConsoleFormatter:
    """Tests for ConsoleFormatter."""

    def test_header_contains_circleseeker(self):
        fmt = ConsoleFormatter()
        header = fmt.header()
        assert "CircleSeeker" in header

    def test_separator_returns_string_of_given_char(self):
        fmt = ConsoleFormatter(width=40)
        sep = fmt.separator("=")
        assert sep == "=" * 40

    def test_separator_default_dash(self):
        fmt = ConsoleFormatter(width=20)
        sep = fmt.separator()
        assert sep == "-" * 20

    def test_format_line_formats_with_bullet(self):
        fmt = ConsoleFormatter()
        line = fmt.format_line("", "Label", "Value")
        # Should contain the middle-dot bullet and label/value
        assert "Label" in line
        assert "Value" in line
        assert "\u00b7" in line  # middle dot

    def test_format_config_formats_dict(self):
        fmt = ConsoleFormatter()
        config = {
            "input_file": "sample.fasta",
            "reference": "ref.fa",
            "output_dir": "/tmp/output",
            "threads": 4,
        }
        result = fmt.format_config(config)
        assert "sample.fasta" in result
        assert "ref.fa" in result
        assert "4" in result


# ---------------------------------------------------------------------------
# TestProgressDisplay
# ---------------------------------------------------------------------------
class TestProgressDisplay:
    """Tests for ProgressDisplay."""

    def test_update_increments_and_returns_string(self):
        progress = ProgressDisplay(total=10, desc="Test")
        result = progress.update(1)
        assert progress.current == 1
        assert isinstance(result, str)
        assert "10%" in result or "Test" in result

    def test_format_step_returns_formatted_line(self):
        progress = ProgressDisplay(total=5)
        line = progress.format_step("tidehunter", status="running")
        assert "tidehunter" in line
        assert isinstance(line, str)
