"""Pytest configuration for CircleSeeker tests."""

import logging
import sys
from pathlib import Path

import pytest

# Add src directory to Python path
ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


@pytest.fixture(autouse=True)
def reset_logging_after_test():
    """Reset circleseeker logger state after each test.

    This prevents test pollution from tests that call setup_logging(),
    which sets propagate=False and breaks caplog in subsequent tests.
    """
    yield
    # Restore logger to clean state after each test
    app_logger = logging.getLogger("circleseeker")
    for handler in list(app_logger.handlers):
        app_logger.removeHandler(handler)
    app_logger.setLevel(logging.NOTSET)
    app_logger.propagate = True