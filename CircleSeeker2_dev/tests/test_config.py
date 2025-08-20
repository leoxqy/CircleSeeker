# tests/test_config.py
"""Tests for the configuration module."""

import pytest
import yaml
from pathlib import Path
from circleseeker2.config import Config, def_config, save_config, ToolConfig
from circleseeker2.exceptions import ConfigurationError

def test_config_defaults():
    """Test that Config initializes with correct default values."""
    cfg = Config()
    assert cfg.output_dir == Path("circleseeker2_output")
    assert cfg.prefix == "sample"
    assert cfg.performance.threads == 8
    assert isinstance(cfg.tools, ToolConfig)
    assert cfg.tools.blast["evalue"] == 1e-50

def test_config_validation(temp_dir: Path):
    """Test the configuration validation logic."""
    # Missing input file
    with pytest.raises(ConfigurationError, match="Input file is required"):
        Config(reference=temp_dir / "ref.fa").validate()

    # Non-existent reference file
    input_file = temp_dir / "input.fa"
    input_file.touch()
    with pytest.raises(ConfigurationError, match="Reference file not found"):
        Config(input_file=input_file, reference=temp_dir / "ref.fa").validate()

    # Invalid threads
    ref_file = temp_dir / "ref.fa"
    ref_file.touch()
    cfg = Config(input_file=input_file, reference=ref_file)
    cfg.performance.threads = 0
    with pytest.raises(ConfigurationError, match="Threads must be >= 1"):
        cfg.validate()

def test_save_and_load_config(temp_dir: Path):
    """Test saving a config to YAML and loading it back."""
    input_file = temp_dir / "input.fq"
    ref_file = temp_dir / "ref.fa"
    input_file.touch()
    ref_file.touch()

    cfg = Config(
        input_file=input_file,
        reference=ref_file,
        output_dir=temp_dir / "out",
        prefix="my_sample",
        enable_xecc=True
    )
    cfg.performance.threads = 16
    cfg.quality.min_coverage = 20
    cfg.tools.tidehunter["k"] = 24

    config_path = temp_dir / "config.yaml"
    save_config(cfg, config_path)

    assert config_path.exists()

    # Check content
    with open(config_path, 'r') as f:
        data = yaml.safe_load(f)
    assert data["prefix"] == "my_sample"
    assert data["performance"]["threads"] == 16
    assert data["tools"]["tidehunter"]["k"] == 24

    # Load it back
    loaded_cfg = def_config(config_path)
    assert loaded_cfg.prefix == "my_sample"
    assert loaded_cfg.performance.threads == 16
    assert loaded_cfg.tools.tidehunter["k"] == 24
    assert loaded_cfg.input_file == input_file
    assert loaded_cfg.reference == ref_file
