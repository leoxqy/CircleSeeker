# tests/conftest.py
"""Shared fixtures for pytest."""

import pytest
from pathlib import Path
from circleseeker2.config import Config

@pytest.fixture
def temp_dir(tmp_path: Path) -> Path:
    """Create a temporary directory for testing."""
    return tmp_path

@pytest.fixture
def basic_config(temp_dir: Path) -> Config:
    """Return a basic, valid Config object."""
    input_file = temp_dir / "input.fasta"
    ref_file = temp_dir / "reference.fasta"
    input_file.touch()
    ref_file.touch()
    
    return Config(
        input_file=input_file,
        reference=ref_file,
        output_dir=temp_dir / "output"
    )
