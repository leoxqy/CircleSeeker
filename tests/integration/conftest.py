"""Pytest configuration for integration tests."""

import pytest


def pytest_configure(config):
    """Register custom markers."""
    config.addinivalue_line(
        "markers", "integration: mark test as integration test (requires external data)"
    )
    config.addinivalue_line(
        "markers", "requires_data: mark test as requiring external data files"
    )
    config.addinivalue_line(
        "markers", "requires_reference: mark test as requiring reference genome"
    )
