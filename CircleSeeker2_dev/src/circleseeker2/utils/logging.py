# src/circleseeker2/utils/logging.py
"""Logging setup for CircleSeeker2."""
import logging
import sys
from pathlib import Path
from typing import Optional

def setup_logging(level: int = logging.INFO, log_file: Optional[Path] = None):
    """
    Configure logging for the application.

    Args:
        level: The logging level (e.g., logging.INFO).
        log_file: Optional path to a file to write logs to.
    """
    root_logger = logging.getLogger()
    root_logger.setLevel(level)

    # Remove existing handlers to avoid duplication
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)

    # Create formatter
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(level)
    console_handler.setFormatter(formatter)
    root_logger.addHandler(console_handler)

    # File handler
    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file, mode='a')
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        root_logger.addHandler(file_handler)

    logging.info(f"Logging configured with level {logging.getLevelName(level)}")
