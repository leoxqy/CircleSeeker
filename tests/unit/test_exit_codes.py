"""Tests for exit_codes module."""

import pytest

from circleseeker.cli.exit_codes import (
    EXIT_SUCCESS,
    EXIT_ERROR,
    EXIT_USAGE,
    EXIT_SIGINT,
    EXIT_SIGTERM,
)


class TestExitCodes:
    """Test exit code constants."""

    def test_exit_success_is_zero(self):
        """EXIT_SUCCESS should be 0."""
        assert EXIT_SUCCESS == 0

    def test_exit_error_is_one(self):
        """EXIT_ERROR should be 1."""
        assert EXIT_ERROR == 1

    def test_exit_usage_is_two(self):
        """EXIT_USAGE should be 2 (shell convention for usage errors)."""
        assert EXIT_USAGE == 2

    def test_exit_sigint_is_130(self):
        """EXIT_SIGINT should be 130 (128 + SIGINT(2))."""
        assert EXIT_SIGINT == 130
        # Verify the formula
        assert EXIT_SIGINT == 128 + 2

    def test_exit_sigterm_is_143(self):
        """EXIT_SIGTERM should be 143 (128 + SIGTERM(15))."""
        assert EXIT_SIGTERM == 143
        # Verify the formula
        assert EXIT_SIGTERM == 128 + 15

    def test_all_exit_codes_are_integers(self):
        """All exit codes should be integers."""
        codes = [EXIT_SUCCESS, EXIT_ERROR, EXIT_USAGE, EXIT_SIGINT, EXIT_SIGTERM]
        for code in codes:
            assert isinstance(code, int)

    def test_all_exit_codes_are_non_negative(self):
        """All exit codes should be non-negative."""
        codes = [EXIT_SUCCESS, EXIT_ERROR, EXIT_USAGE, EXIT_SIGINT, EXIT_SIGTERM]
        for code in codes:
            assert code >= 0

    def test_all_exit_codes_fit_in_byte(self):
        """All exit codes should fit in a byte (0-255)."""
        codes = [EXIT_SUCCESS, EXIT_ERROR, EXIT_USAGE, EXIT_SIGINT, EXIT_SIGTERM]
        for code in codes:
            assert 0 <= code <= 255

    def test_exit_codes_are_unique(self):
        """All exit codes should be unique."""
        codes = [EXIT_SUCCESS, EXIT_ERROR, EXIT_USAGE, EXIT_SIGINT, EXIT_SIGTERM]
        assert len(codes) == len(set(codes))
