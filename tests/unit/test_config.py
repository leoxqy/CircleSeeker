"""Tests for config module."""

from pathlib import Path
import sys
import pytest
import tempfile
import yaml
from dataclasses import asdict

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.config import (
    RuntimeConfig, PerformanceConfig, ToolConfig, Config, load_config
)
from circleseeker.exceptions import ConfigurationError


class TestRuntimeConfig:
    """Test cases for RuntimeConfig."""

    def test_runtime_config_defaults(self):
        """Test RuntimeConfig default values."""
        config = RuntimeConfig()
        assert config.log_level == "WARNING"
        assert config.log_file is None
        assert config.tmp_dir == Path(".tmp_work")
        assert config.keep_tmp is False
        assert config.checkpoint_policy == "continue"
        assert config.enable_progress is True

    def test_runtime_config_custom_values(self):
        """Test RuntimeConfig with custom values."""
        log_file = Path("/tmp/test.log")
        tmp_dir = Path("/tmp/custom")

        config = RuntimeConfig(
            log_level="DEBUG",
            log_file=log_file,
            tmp_dir=tmp_dir,
            keep_tmp=True,
            checkpoint_policy="reset",
            enable_progress=False
        )

        assert config.log_level == "DEBUG"
        assert config.log_file == log_file
        assert config.tmp_dir == tmp_dir
        assert config.keep_tmp is True
        assert config.checkpoint_policy == "reset"
        assert config.enable_progress is False


class TestPerformanceConfig:
    """Test cases for PerformanceConfig."""

    def test_performance_config_defaults(self):
        """Test PerformanceConfig default values."""
        config = PerformanceConfig()
        assert config.threads == 8

    def test_performance_config_custom_values(self):
        """Test PerformanceConfig with custom values."""
        config = PerformanceConfig(threads=16)
        assert config.threads == 16


class TestToolConfig:
    """Test cases for ToolConfig."""

    def test_tool_config_defaults(self):
        """Test ToolConfig default values."""
        config = ToolConfig()

        # Check tidehunter defaults
        assert config.tidehunter["k"] == 16
        assert config.tidehunter["w"] == 1
        assert config.tidehunter["p"] == 100
        assert config.tidehunter["P"] == 2000000
        assert config.tidehunter["e"] == 0.1
        assert config.tidehunter["f"] == 2

        # Check blast defaults
        assert config.blast["word_size"] == 100
        assert config.blast["evalue"] == "1e-50"
        assert config.blast["perc_identity"] == 99.0

        # Check minimap2 defaults
        assert config.minimap2["preset"] == "map-hifi"
        assert config.minimap2["additional_args"] == ""

        # Check samtools defaults (empty dict)
        assert config.samtools == {}

    def test_tool_config_custom_values(self):
        """Test ToolConfig with custom values."""
        config = ToolConfig(
            tidehunter={"k": 20, "w": 2},
            blast={"word_size": 50, "evalue": "1e-10"},
            minimap2={"preset": "map-ont", "additional_args": "-k 15"},
            samtools={"sort_memory": "4G"}
        )

        assert config.tidehunter["k"] == 20
        assert config.tidehunter["w"] == 2
        assert config.blast["word_size"] == 50
        assert config.blast["evalue"] == "1e-10"
        assert config.minimap2["preset"] == "map-ont"
        assert config.minimap2["additional_args"] == "-k 15"
        assert config.samtools["sort_memory"] == "4G"


class TestConfig:
    """Test cases for main Config class."""

    def test_config_defaults(self):
        """Test Config default values."""
        config = Config()
        assert config.input_file is None
        assert config.reference is None
        assert config.output_dir == Path("circleseeker_output")
        assert config.prefix == "sample"
        assert config.enable_xecc is True

        # Check skip flags
        assert config.skip_make_db is False
        assert config.skip_tidehunter is False
        assert config.skip_carousel is False
        assert config.skip_blast is False
        assert config.skip_gatekeeper is False
        assert config.skip_organize is False

        # Check sub-configurations exist
        assert isinstance(config.runtime, RuntimeConfig)
        assert isinstance(config.performance, PerformanceConfig)
        assert isinstance(config.tools, ToolConfig)

    def test_config_custom_values(self):
        """Test Config with custom values."""
        input_file = Path("test_input.fasta")
        reference = Path("test_ref.fa")
        output_dir = Path("test_output")

        config = Config(
            input_file=input_file,
            reference=reference,
            output_dir=output_dir,
            prefix="test_sample",
            enable_xecc=False,
            skip_make_db=True,
            skip_tidehunter=True
        )

        assert config.input_file == input_file
        assert config.reference == reference
        assert config.output_dir == output_dir
        assert config.prefix == "test_sample"
        assert config.enable_xecc is False
        assert config.skip_make_db is True
        assert config.skip_tidehunter is True

    def test_config_threads_property(self):
        """Test threads property getter and setter."""
        config = Config()

        # Default threads
        assert config.threads == 8

        # Set threads via property
        config.threads = 16
        assert config.threads == 16
        assert config.performance.threads == 16

    def test_config_keep_tmp_property(self):
        """Test keep_tmp property getter and setter."""
        config = Config()

        # Default keep_tmp
        assert config.keep_tmp is False

        # Set keep_tmp via property
        config.keep_tmp = True
        assert config.keep_tmp is True
        assert config.runtime.keep_tmp is True

    def test_config_canonical_skip_properties(self):
        """Test canonical skip flag properties."""
        config = Config()

        # Test skip_make_blastdb property
        assert config.skip_make_blastdb is False
        config.skip_make_blastdb = True
        assert config.skip_make_blastdb is True
        assert config.skip_make_db is True

        # Test skip_tandem_to_ring property
        assert config.skip_tandem_to_ring is False
        config.skip_tandem_to_ring = True
        assert config.skip_tandem_to_ring is True
        assert config.skip_carousel is True

        # Test skip_run_blast property
        assert config.skip_run_blast is False
        config.skip_run_blast = True
        assert config.skip_run_blast is True
        assert config.skip_blast is True

        # Test skip_um_classify property
        assert config.skip_um_classify is False
        config.skip_um_classify = True
        assert config.skip_um_classify is True
        assert config.skip_gatekeeper is True

    def test_config_validate_missing_required(self):
        """Test config validation with missing required fields."""
        config = Config()

        # Should fail validation - missing required fields
        with pytest.raises(ConfigurationError) as exc_info:
            config.validate()
        assert "Input file is required" in str(exc_info.value)

    def test_config_validate_missing_reference(self):
        """Test config validation with missing reference."""
        config = Config()
        config.input_file = Path("dummy.fasta")

        with pytest.raises(ConfigurationError) as exc_info:
            config.validate()
        assert "Reference genome is required" in str(exc_info.value)

    def test_config_validate_nonexistent_files(self):
        """Test config validation with nonexistent files."""
        config = Config()
        config.input_file = Path("nonexistent_input.fasta")
        config.reference = Path("nonexistent_ref.fa")

        with pytest.raises(ConfigurationError) as exc_info:
            config.validate()
        assert "Input file not found" in str(exc_info.value)

    def test_config_validate_invalid_threads(self):
        """Test config validation with invalid thread count."""
        with tempfile.NamedTemporaryFile(suffix=".fasta") as input_f, \
             tempfile.NamedTemporaryFile(suffix=".fa") as ref_f:

            config = Config()
            config.input_file = Path(input_f.name)
            config.reference = Path(ref_f.name)
            config.performance.threads = 0

            with pytest.raises(ConfigurationError) as exc_info:
                config.validate()
            assert "Threads must be >= 1" in str(exc_info.value)

    def test_config_validate_success(self):
        """Test successful config validation."""
        with tempfile.NamedTemporaryFile(suffix=".fasta") as input_f, \
             tempfile.NamedTemporaryFile(suffix=".fa") as ref_f:

            config = Config()
            config.input_file = Path(input_f.name)
            config.reference = Path(ref_f.name)

            # Should not raise any exception
            config.validate()

    def test_config_to_dict(self):
        """Test config to_dict conversion."""
        config = Config()
        config.input_file = Path("test.fasta")
        config.reference = Path("ref.fa")
        config.threads = 4

        config_dict = config.to_dict()

        assert isinstance(config_dict, dict)
        assert config_dict["input_file"] == "test.fasta"  # Path converted to string
        assert config_dict["reference"] == "ref.fa"
        assert config_dict["performance"]["threads"] == 4
        assert "runtime" in config_dict
        assert "tools" in config_dict

    def test_config_to_dict_nested_paths(self):
        """Test config to_dict with nested Path objects."""
        config = Config()
        config.runtime.log_file = Path("/tmp/test.log")
        config.runtime.tmp_dir = Path("/tmp/custom")

        config_dict = config.to_dict()

        assert config_dict["runtime"]["log_file"] == "/tmp/test.log"
        assert config_dict["runtime"]["tmp_dir"] == "/tmp/custom"


class TestLoadConfig:
    """Test cases for load_config function."""

    def test_load_config_basic(self):
        """Test basic config loading from YAML."""
        config_data = {
            "input_file": "test_input.fasta",
            "reference": "test_ref.fa",
            "output_dir": "test_output",
            "prefix": "test_sample"
        }

        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml.dump(config_data, f)
            config_path = Path(f.name)

        try:
            config = load_config(config_path)

            assert config.input_file == Path("test_input.fasta")
            assert config.reference == Path("test_ref.fa")
            assert config.output_dir == Path("test_output")
            assert config.prefix == "test_sample"
        finally:
            config_path.unlink()

    def test_load_config_empty_file(self):
        """Test loading config from empty YAML file."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write("")  # Empty file
            config_path = Path(f.name)

        try:
            config = load_config(config_path)

            # Should create config with defaults
            assert isinstance(config, Config)
            assert config.prefix == "sample"  # Default value
        finally:
            config_path.unlink()

    def test_load_config_with_nested_sections(self):
        """Test loading config with nested sections."""
        config_data = {
            "input_file": "test.fasta",
            "reference": "ref.fa",
            "performance": {
                "threads": 16
            },
            "tools": {
                "blast": {
                    "evalue": "1e-10",
                    "word_size": 50
                }
            }
        }

        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml.dump(config_data, f)
            config_path = Path(f.name)

        try:
            config = load_config(config_path)

            assert config.input_file == Path("test.fasta")
            assert config.performance.threads == 16
            assert config.tools.blast["evalue"] == "1e-10"
            assert config.tools.blast["word_size"] == 50
        finally:
            config_path.unlink()

    def test_load_config_skip_flags(self):
        """Test loading config with skip flags."""
        config_data = {
            "input_file": "test.fasta",
            "reference": "ref.fa",
            "skip_make_db": True,
            "skip_tidehunter": True,
            "skip_blast": True
        }

        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml.dump(config_data, f)
            config_path = Path(f.name)

        try:
            config = load_config(config_path)

            assert config.skip_make_db is True
            assert config.skip_tidehunter is True
            assert config.skip_blast is True
            # Test canonical properties
            assert config.skip_make_blastdb is True
            assert config.skip_run_blast is True
        finally:
            config_path.unlink()

    def test_load_config_nonexistent_file(self):
        """Test loading config from nonexistent file."""
        nonexistent_path = Path("nonexistent_config.yaml")

        with pytest.raises(FileNotFoundError):
            load_config(nonexistent_path)


class TestConfigIntegration:
    """Integration tests for config functionality."""

    def test_complete_config_workflow(self):
        """Test complete config workflow from creation to validation."""
        # Create config
        config = Config()

        # Set required fields
        with tempfile.NamedTemporaryFile(suffix=".fasta") as input_f, \
             tempfile.NamedTemporaryFile(suffix=".fa") as ref_f:

            config.input_file = Path(input_f.name)
            config.reference = Path(ref_f.name)
            config.prefix = "integration_test"
            config.threads = 8
            config.keep_tmp = True

            # Validate
            config.validate()

            # Convert to dict
            config_dict = config.to_dict()
            assert config_dict["prefix"] == "integration_test"
            assert config_dict["performance"]["threads"] == 8
            assert config_dict["runtime"]["keep_tmp"] is True

    def test_config_yaml_roundtrip(self):
        """Test config save and load roundtrip."""
        # Create config with custom values
        original_config = Config()
        original_config.prefix = "roundtrip_test"
        original_config.threads = 12
        original_config.skip_make_blastdb = True

        # Save to YAML
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml.dump(original_config.to_dict(), f)
            config_path = Path(f.name)

        try:
            # Load from YAML
            loaded_config = load_config(config_path)

            # Compare important values
            assert loaded_config.prefix == "roundtrip_test"
            assert loaded_config.performance.threads == 12
            assert loaded_config.skip_make_db is True  # Note: uses legacy field name
            assert loaded_config.skip_make_blastdb is True  # But canonical property works
        finally:
            config_path.unlink()
