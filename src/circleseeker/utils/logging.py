"""Centralized logging utilities for CircleSeeker.

Provides a single place to configure logging and fetch namespaced loggers.
"""

from __future__ import annotations

import logging
from logging.handlers import RotatingFileHandler
from pathlib import Path
from typing import Optional


DEFAULT_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
CONSOLE_FORMAT = "%(levelname)s: %(message)s"

# Log rotation settings
DEFAULT_MAX_BYTES = 10 * 1024 * 1024  # 10 MB
DEFAULT_BACKUP_COUNT = 5


def setup_logging(
    level: int = logging.INFO,
    log_file: Optional[Path] = None,
    max_bytes: int = DEFAULT_MAX_BYTES,
    backup_count: int = DEFAULT_BACKUP_COUNT,
) -> None:
    """Configure logging for the 'circleseeker' namespace.

    Args:
        level: Logging level for the application logger
        log_file: Optional path for log file output
        max_bytes: Maximum size per log file before rotation (default: 10MB)
        backup_count: Number of backup files to keep (default: 5)

    Notes:
        - Root logger kept at WARNING to suppress third-party noise
        - 'circleseeker' logger uses the requested level
        - Console handler uses concise format; file handler (if any) is detailed at DEBUG
        - File handler uses rotation to prevent unbounded log growth
    """
    # Keep root quiet
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.WARNING)

    # Configure namespaced logger
    app_logger = logging.getLogger("circleseeker")
    app_logger.setLevel(level)
    # Avoid duplicate logs if called multiple times
    if app_logger.handlers:
        for h in list(app_logger.handlers):
            app_logger.removeHandler(h)

    console = logging.StreamHandler()
    console.setLevel(level)
    console.setFormatter(logging.Formatter(CONSOLE_FORMAT))
    app_logger.addHandler(console)

    if log_file:
        try:
            log_file.parent.mkdir(parents=True, exist_ok=True)
            # Use rotating file handler to prevent unbounded log growth
            file_handler = RotatingFileHandler(
                log_file,
                maxBytes=max_bytes,
                backupCount=backup_count,
            )
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(logging.Formatter(DEFAULT_FORMAT))
            app_logger.addHandler(file_handler)
        except (OSError, IOError) as e:
            import warnings
            warnings.warn(f"Failed to create log file {log_file}: {e}")

    # Do not propagate to root to avoid double-printing
    app_logger.propagate = False


def get_logger(name: str) -> logging.Logger:
    """Return a namespaced logger under 'circleseeker' root."""
    base = logging.getLogger("circleseeker")
    return base.getChild(name)


class LogTemplates:
    """Standard log message templates for consistent logging across modules.

    Use these templates to ensure uniform log output throughout the application.

    Example usage:
        logger.info(LogTemplates.STEP_START.format(
            step_name="tidehunter",
            step_number=2,
            total=16
        ))
    """

    # Step lifecycle messages
    STEP_START = "Starting step: {step_name} (#{step_number}/{total})"
    STEP_SUCCESS = "Completed step: {step_name} in {duration:.1f}s"
    STEP_FAILURE = "Failed at step: {step_name} - {error}"
    STEP_SKIPPED = "Skipping step: {step_name} - {reason}"

    # File operations
    FILE_CREATED = "Created output file: {path} ({size:,} bytes)"
    FILE_LOADED = "Loaded {count:,} records from {path}"
    FILE_NOT_FOUND = "File not found: {path}"

    # Processing statistics
    PROCESSING_STATS = "Processed {input_count:,} items → {output_count:,} results"
    FILTERING_STATS = "Filtered: {kept:,} kept, {removed:,} removed ({percent:.1f}% pass rate)"

    # External tool execution
    TOOL_START = "Running {tool_name}: {description}"
    TOOL_SUCCESS = "{tool_name} completed in {duration:.1f}s"
    TOOL_FAILURE = "{tool_name} failed with exit code {exit_code}"
