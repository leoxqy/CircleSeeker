"""Centralized logging utilities for CircleSeeker.

Provides a single place to configure logging and fetch namespaced loggers.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional


DEFAULT_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
CONSOLE_FORMAT = "%(levelname)s: %(message)s"


def setup_logging(level: int = logging.INFO, log_file: Optional[Path] = None) -> None:
    """Configure logging for the 'circleseeker' namespace.

    - Root logger kept at WARNING to suppress third-party noise
    - 'circleseeker' logger uses the requested level
    - Console handler uses concise format; file handler (if any) is detailed at DEBUG
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
        log_file.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(logging.Formatter(DEFAULT_FORMAT))
        app_logger.addHandler(file_handler)

    # Do not propagate to root to avoid double-printing
    app_logger.propagate = False


def get_logger(name: str) -> logging.Logger:
    """Return a namespaced logger under 'circleseeker' root."""
    base = logging.getLogger("circleseeker")
    return base.getChild(name)

