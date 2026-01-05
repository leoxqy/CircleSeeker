"""Enhanced display utilities for CircleSeeker with emoji support and beautiful formatting."""

from __future__ import annotations

import os
import sys
import locale
from datetime import datetime
from pathlib import Path
from typing import Optional, Dict, Any
import time

# Import version carefully to avoid circular imports
try:
    from circleseeker.__version__ import __version__
except ImportError:
    __version__ = "unknown"  # Fallback version


def supports_emoji() -> bool:
    """Check if the current terminal supports emoji display."""
    # Check if we're in a CI environment
    if os.environ.get("CI"):
        return False

    # Check terminal encoding
    try:
        encoding = locale.getpreferredencoding() or "ascii"
        if "utf" not in encoding.lower():
            return False
    except (LookupError, TypeError, ValueError):
        return False

    # Check terminal type
    term = os.environ.get("TERM", "").lower()
    if term in ["dumb", "unknown"]:
        return False

    # Check if stdout is a TTY
    if not sys.stdout.isatty():
        return False

    # Platform-specific checks
    if sys.platform == "win32":
        # Windows Terminal and modern Windows 10/11 support emoji
        if os.environ.get("WT_SESSION"):  # Windows Terminal
            return True
        # Check Windows version
        try:
            import platform

            version = tuple(map(int, platform.version().split(".")))
            return version >= (10, 0, 14393)  # Windows 10 Anniversary Update
        except (ValueError, AttributeError, TypeError):
            return False

    return True  # Most Unix-like systems support emoji


class EmojiHandler:
    """Smart emoji handler with ASCII fallbacks."""

    # Emoji to ASCII mapping
    EMOJI_MAP = {
        "📄": "[F]",  # File
        "🧬": "[R]",  # Reference
        "📁": "[D]",  # Directory
        "🔧": "[#]",  # Threads
        "☕": ">",  # Processing
        "✅": "[OK]",  # Success
        "❌": "[X]",  # Error
        "⌛": "[T]",  # Time (hourglass)
        "🧹": "[C]",  # Cleanup
        "🔥": "[>]",  # Start (fire)
        "💾": "[S]",  # Save
        "📊": "[#]",  # Stats
        "🔍": "[?]",  # Search/Analysis
        "⚠️": "[!]",  # Warning
        "ℹ️": "[i]",  # Info
        "🎯": "[*]",  # Target
        "📝": "[W]",  # Write
        "🔄": "[~]",  # Process/Convert
        "✨": "[+]",  # New/Create
        "🏁": "[=]",  # Finish
    }

    def __init__(self):
        self._supports_emoji = supports_emoji()

    def get(self, emoji: str, fallback: str = "") -> str:
        """Get emoji or its ASCII fallback."""
        if self._supports_emoji:
            return emoji
        return self.EMOJI_MAP.get(emoji, fallback)


# Global emoji handler instance
emoji_handler = EmojiHandler()


def E(emoji: str, fallback: str = "") -> str:
    """Convenience function for emoji handling."""
    return emoji_handler.get(emoji, fallback)


class ConsoleFormatter:
    """Beautiful console output formatter for CircleSeeker."""

    def __init__(self, width: int = 60, use_emoji: bool = False):
        self.width = width
        self.start_time = None
        self.use_emoji = use_emoji  # Allow disabling emoji even if supported
        self._emoji_supported = supports_emoji() if use_emoji else False

    def header(self, sample_name: Optional[str] = None) -> str:
        """Generate a beautiful header."""
        lines = []
        lines.append("=" * self.width)

        # Title line with version
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        title = f"CircleSeeker Analysis  •  {timestamp}  •  v{__version__}"
        if len(title) < self.width:
            padding = (self.width - len(title)) // 2
            title = " " * padding + title
        lines.append(title)

        lines.append("-" * self.width)  # Use single dash line instead of equals
        return "\n".join(lines)

    def separator(self, char: str = "-") -> str:
        """Generate a separator line."""
        return char * self.width

    def format_line(self, emoji: str, label: str, value: str, label_width: int = 12) -> str:
        """Format a single information line with bullet point."""
        # Use middle dot instead of emoji for cleaner look
        bullet = "·"  # Middle dot character
        return f"{bullet} {label:<{label_width}} : {value}"

    def _is_emoji(self, s: str) -> bool:
        """Check if string contains actual emoji (not ASCII fallback)."""
        # If emoji is not supported globally, it's definitely ASCII
        if not self._emoji_supported:
            return False
        # Check if it's an ASCII fallback pattern
        return not s.startswith("[") and s not in [">", "-", "*", "~", "=", "+"]

    def _format_path(self, path_str: str) -> str:
        """Format path for user-friendly display."""
        from pathlib import Path

        path = Path(path_str)
        home = Path.home()

        try:
            # Try to make path relative to home directory
            if path.is_absolute():
                try:
                    rel_to_home = path.relative_to(home)
                    return f"~/{rel_to_home}"
                except ValueError:
                    # Path is not under home directory
                    pass

            # Try to make path relative to current directory
            cwd = Path.cwd()
            try:
                rel_path = path.relative_to(cwd)
                # If it's in current dir or subdirs, show as relative
                if str(rel_path) != str(path):
                    return str(rel_path)
            except ValueError:
                pass

        except Exception as e:
            # Path resolution failed (e.g., permission denied, invalid path)
            # Log at debug level and fall back to original string
            import logging
            logging.getLogger("display").debug(f"Path formatting failed for '{path_str}': {e}")

        # Fallback to original path
        return str(path_str)

    def format_config(self, config: Dict[str, Any]) -> str:
        """Format configuration information."""
        lines = []

        # Input/Output info - no emoji, just clean formatting
        if "input_file" in config:
            lines.append(self.format_line("", "Input sample", str(config["input_file"])))
        if "reference" in config:
            lines.append(self.format_line("", "Reference", str(config["reference"])))
        if "output_dir" in config:
            # Convert to user-friendly path
            output_path = self._format_path(config["output_dir"])
            lines.append(self.format_line("", "Output dir", output_path))
        if "threads" in config:
            lines.append(self.format_line("", "Threads", str(config["threads"])))
        if "inference" in config:
            lines.append(self.format_line("", "Inference", str(config["inference"])))

        return "\n".join(lines)

    def start_message(self) -> str:
        """Generate start message."""
        self.start_time = time.time()
        return self.format_line("", "Status", "Starting analysis...", label_width=12)

    def success_message(self) -> str:
        """Generate success message with timing."""
        lines = []
        lines.append(self._format_status_line("", "Pipeline completed successfully!"))

        if self.start_time:
            elapsed = time.time() - self.start_time
            minutes = int(elapsed // 60)
            seconds = int(elapsed % 60)
            if minutes > 0:
                time_str = f"{minutes}m {seconds}s"
            else:
                time_str = f"{seconds}s"
            lines.append(self._format_status_line("", f"Time elapsed: {time_str}"))

        return "\n".join(lines)

    def _format_status_line(self, emoji: str, message: str) -> str:
        """Format a status line with bullet point."""
        bullet = "·"  # Middle dot character
        return f"{bullet} {message}"

    def error_message(self, error: str) -> str:
        """Generate error message."""
        return self._format_status_line("", f"Pipeline failed: {error}")

    def cleanup_message(self) -> str:
        """Generate cleanup message."""
        return self._format_status_line("", "Temporary workdir cleaned up")

    def final_output_message(self, output_dir: Path) -> str:
        """Generate final output location message."""
        friendly_path = self._format_path(output_dir)
        return self._format_status_line("", f"Finalized outputs saved in: {friendly_path}")


class ProgressDisplay:
    """Progress bar display with emoji support."""

    def __init__(self, total: int, desc: str = "Pipeline"):
        self.total = total
        self.desc = desc
        self.current = 0
        self.bar_width = 30

    def update(self, n: int = 1) -> str:
        """Update progress and return formatted string."""
        self.current += n
        progress = self.current / self.total

        # Build progress bar
        filled = int(self.bar_width * progress)
        bar = "█" * filled + "░" * (self.bar_width - filled)

        # Format with emoji
        emoji = E("☕") if progress < 1.0 else E("✅")
        percent = int(progress * 100)

        return f"{emoji} {self.desc}: {percent:3d}%|{bar}| {self.current}/{self.total}"

    def format_step(self, step_name: str, status: str = "running") -> str:
        """Format a step status line."""
        status_emoji = {
            "running": E("🔄"),
            "complete": E("✅"),
            "failed": E("❌"),
            "skipped": E("⏭️", "[>]"),
        }.get(status, E("🔄"))

        return f"  {status_emoji} {step_name}"


def format_stats(stats: Dict[str, int]) -> str:
    """Format statistics output."""
    lines = []
    lines.append(f"\n{E('📊')} Final eccDNA counts:")

    for ecc_type, count in stats.items():
        lines.append(f"  {ecc_type:<10} : {count:>6} sequences")

    return "\n".join(lines)


def print_formatted(message: str, file=sys.stdout) -> None:
    """Print formatted message to console."""
    print(message, file=file, flush=True)


# Example usage function
def demo():
    """Demonstrate the display utilities."""
    formatter = ConsoleFormatter()

    # Print header
    print(formatter.header())

    # Print config
    config = {
        "input_file": "HeLa_rep1_1k.fasta",
        "reference": "chm13v2.0.fa",
        "output_dir": "/abs/path/to/10kkkk",
        "threads": 8,
    }
    print(formatter.format_config(config))
    print(formatter.separator())

    # Start message
    print(formatter.start_message())

    # Progress simulation
    progress = ProgressDisplay(16)
    for i in range(16):
        print(f"\r{progress.update()}", end="")
        time.sleep(0.1)
    print()  # New line after progress

    print(formatter.separator())

    # Success message
    print(formatter.success_message())
    print(formatter.cleanup_message())
    print(formatter.separator("="))


if __name__ == "__main__":
    demo()
