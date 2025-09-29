"""Tests for utils logging module."""

from pathlib import Path
import sys
import pytest
import logging
import tempfile
from unittest.mock import patch, MagicMock

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.utils.logging import (
    setup_logging, get_logger, DEFAULT_FORMAT, CONSOLE_FORMAT
)


class TestLoggingUtilities:
    """Test cases for logging utilities."""

    def test_default_format_constants(self):
        """Test that format constants are defined."""
        assert DEFAULT_FORMAT is not None
        assert CONSOLE_FORMAT is not None
        assert "%(asctime)s" in DEFAULT_FORMAT
        assert "%(name)s" in DEFAULT_FORMAT
        assert "%(levelname)s" in DEFAULT_FORMAT
        assert "%(message)s" in DEFAULT_FORMAT
        assert "%(levelname)s" in CONSOLE_FORMAT
        assert "%(message)s" in CONSOLE_FORMAT

    def test_get_logger_basic(self):
        """Test get_logger returns a logger with correct name."""
        logger = get_logger("test_module")
        assert logger is not None
        assert isinstance(logger, logging.Logger)
        assert logger.name == "circleseeker.test_module"

    def test_get_logger_different_names(self):
        """Test get_logger with different module names."""
        logger1 = get_logger("module1")
        logger2 = get_logger("module2")
        logger3 = get_logger("external.blast")

        assert logger1.name == "circleseeker.module1"
        assert logger2.name == "circleseeker.module2"
        assert logger3.name == "circleseeker.external.blast"

        # Should return different logger instances
        assert logger1 != logger2
        assert logger1 != logger3

    def test_get_logger_same_name_returns_same_instance(self):
        """Test that get_logger returns the same instance for the same name."""
        logger1 = get_logger("same_module")
        logger2 = get_logger("same_module")

        assert logger1 is logger2

    def test_setup_logging_basic(self):
        """Test basic setup_logging functionality."""
        # Clear any existing handlers first
        app_logger = logging.getLogger("circleseeker")
        for handler in list(app_logger.handlers):
            app_logger.removeHandler(handler)

        setup_logging(level=logging.INFO)

        # Check that circleseeker logger is configured
        app_logger = logging.getLogger("circleseeker")
        assert app_logger.level == logging.INFO
        assert len(app_logger.handlers) == 1
        assert isinstance(app_logger.handlers[0], logging.StreamHandler)
        assert app_logger.propagate is False

        # Check root logger level
        root_logger = logging.getLogger()
        assert root_logger.level == logging.WARNING

    def test_setup_logging_with_file(self):
        """Test setup_logging with file handler."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            log_file = Path(tmp_dir) / "test.log"

            # Clear any existing handlers first
            app_logger = logging.getLogger("circleseeker")
            for handler in list(app_logger.handlers):
                app_logger.removeHandler(handler)

            setup_logging(level=logging.DEBUG, log_file=log_file)

            # Check that both console and file handlers are present
            app_logger = logging.getLogger("circleseeker")
            assert len(app_logger.handlers) == 2

            # Find console and file handlers
            console_handler = None
            file_handler = None
            for handler in app_logger.handlers:
                if isinstance(handler, logging.StreamHandler) and not isinstance(handler, logging.FileHandler):
                    console_handler = handler
                elif isinstance(handler, logging.FileHandler):
                    file_handler = handler

            assert console_handler is not None
            assert file_handler is not None

            # Check levels
            assert console_handler.level == logging.DEBUG
            assert file_handler.level == logging.DEBUG

            # Check that log file was created
            assert log_file.exists()

    def test_setup_logging_multiple_calls(self):
        """Test that multiple calls to setup_logging don't create duplicate handlers."""
        # Clear any existing handlers first
        app_logger = logging.getLogger("circleseeker")
        for handler in list(app_logger.handlers):
            app_logger.removeHandler(handler)

        # Call setup_logging multiple times
        setup_logging(level=logging.INFO)
        setup_logging(level=logging.DEBUG)
        setup_logging(level=logging.WARNING)

        # Should only have one handler
        app_logger = logging.getLogger("circleseeker")
        assert len(app_logger.handlers) == 1
        assert app_logger.level == logging.WARNING  # Last level set

    def test_setup_logging_file_directory_creation(self):
        """Test that setup_logging creates directory for log file."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            log_file = Path(tmp_dir) / "subdir" / "nested" / "test.log"

            # Directory should not exist initially
            assert not log_file.parent.exists()

            # Clear any existing handlers first
            app_logger = logging.getLogger("circleseeker")
            for handler in list(app_logger.handlers):
                app_logger.removeHandler(handler)

            setup_logging(level=logging.INFO, log_file=log_file)

            # Directory should be created
            assert log_file.parent.exists()
            assert log_file.exists()

    def test_logger_hierarchy(self):
        """Test that logger hierarchy works correctly."""
        # Setup logging first
        app_logger = logging.getLogger("circleseeker")
        for handler in list(app_logger.handlers):
            app_logger.removeHandler(handler)
        setup_logging(level=logging.DEBUG)

        # Get child loggers
        parent_logger = get_logger("external")
        child_logger = get_logger("external.blast")

        # Check names
        assert parent_logger.name == "circleseeker.external"
        assert child_logger.name == "circleseeker.external.blast"

        # Child logger should inherit from parent
        assert child_logger.parent == parent_logger

    def test_formatter_configuration(self):
        """Test that formatters are configured correctly."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            log_file = Path(tmp_dir) / "test.log"

            # Clear any existing handlers first
            app_logger = logging.getLogger("circleseeker")
            for handler in list(app_logger.handlers):
                app_logger.removeHandler(handler)

            setup_logging(level=logging.INFO, log_file=log_file)

            app_logger = logging.getLogger("circleseeker")

            # Find handlers
            console_handler = None
            file_handler = None
            for handler in app_logger.handlers:
                if isinstance(handler, logging.StreamHandler) and not isinstance(handler, logging.FileHandler):
                    console_handler = handler
                elif isinstance(handler, logging.FileHandler):
                    file_handler = handler

            # Check formatter formats
            console_format = console_handler.formatter._fmt
            file_format = file_handler.formatter._fmt

            assert console_format == CONSOLE_FORMAT
            assert file_format == DEFAULT_FORMAT

    @patch('logging.StreamHandler')
    @patch('logging.FileHandler')
    def test_setup_logging_handler_creation(self, mock_file_handler, mock_stream_handler):
        """Test that handlers are created correctly."""
        mock_console = MagicMock()
        mock_file = MagicMock()
        mock_stream_handler.return_value = mock_console
        mock_file_handler.return_value = mock_file

        with tempfile.TemporaryDirectory() as tmp_dir:
            log_file = Path(tmp_dir) / "test.log"

            # Clear any existing handlers first
            app_logger = logging.getLogger("circleseeker")
            for handler in list(app_logger.handlers):
                app_logger.removeHandler(handler)

            setup_logging(level=logging.INFO, log_file=log_file)

            # Verify handlers were created
            mock_stream_handler.assert_called_once()
            mock_file_handler.assert_called_once_with(log_file)

            # Verify setLevel was called on handlers
            mock_console.setLevel.assert_called_with(logging.INFO)
            mock_file.setLevel.assert_called_with(logging.DEBUG)

    def test_logger_usage_example(self):
        """Test typical logger usage pattern."""
        # Setup logging
        app_logger = logging.getLogger("circleseeker")
        for handler in list(app_logger.handlers):
            app_logger.removeHandler(handler)
        setup_logging(level=logging.INFO)

        # Get logger and use it
        logger = get_logger("pipeline")
        assert logger.name == "circleseeker.pipeline"

        # Logger should be properly configured (this would normally log to console)
        # We can't easily test the actual output without capturing it
        assert logger.isEnabledFor(logging.INFO)
        assert logger.getEffectiveLevel() <= logging.INFO

    def test_different_log_levels(self):
        """Test setup_logging with different log levels."""
        levels = [logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR, logging.CRITICAL]

        for level in levels:
            # Clear any existing handlers first
            app_logger = logging.getLogger("circleseeker")
            for handler in list(app_logger.handlers):
                app_logger.removeHandler(handler)

            setup_logging(level=level)

            app_logger = logging.getLogger("circleseeker")
            assert app_logger.level == level

            # Console handler should also have the same level
            console_handler = app_logger.handlers[0]
            assert console_handler.level == level


@pytest.fixture
def clean_logging():
    """Fixture to clean up logging configuration after each test."""
    yield
    # Clean up after test
    app_logger = logging.getLogger("circleseeker")
    for handler in list(app_logger.handlers):
        app_logger.removeHandler(handler)
    app_logger.setLevel(logging.NOTSET)


class TestLoggingIntegration:
    """Integration tests for logging functionality."""

    def test_complete_logging_workflow(self, clean_logging):
        """Test complete logging workflow from setup to usage."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            log_file = Path(tmp_dir) / "workflow.log"

            # Step 1: Setup logging
            setup_logging(level=logging.DEBUG, log_file=log_file)

            # Step 2: Get loggers for different modules
            pipeline_logger = get_logger("pipeline")
            blast_logger = get_logger("external.blast")
            utils_logger = get_logger("utils.progress")

            # Step 3: Verify logger names
            assert pipeline_logger.name == "circleseeker.pipeline"
            assert blast_logger.name == "circleseeker.external.blast"
            assert utils_logger.name == "circleseeker.utils.progress"

            # Step 4: Verify all loggers are properly configured
            for logger in [pipeline_logger, blast_logger, utils_logger]:
                assert logger.isEnabledFor(logging.DEBUG)
                # They should inherit handlers from parent
                effective_handlers = logger.handlers or logger.parent.handlers
                assert len(effective_handlers) >= 1

            # Step 5: Verify log file exists
            assert log_file.exists()