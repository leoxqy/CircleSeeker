"""Tests for exceptions module."""

from pathlib import Path
import sys
import pytest

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.exceptions import (
    CircleSeekerError, ConfigurationError, ExternalToolError,
    PipelineError, ValidationError, FileFormatError
)


class TestCircleSeekerError:
    """Test cases for CircleSeekerError base exception."""

    def test_circleseeker_error_creation(self):
        """Test CircleSeekerError basic creation."""
        error = CircleSeekerError("Test error message")
        assert str(error) == "Test error message"
        assert isinstance(error, Exception)

    def test_circleseeker_error_inheritance(self):
        """Test that CircleSeekerError inherits from Exception."""
        error = CircleSeekerError("Test")
        assert isinstance(error, Exception)
        assert isinstance(error, CircleSeekerError)

    def test_circleseeker_error_raise(self):
        """Test raising CircleSeekerError."""
        with pytest.raises(CircleSeekerError) as exc_info:
            raise CircleSeekerError("Test error")

        assert str(exc_info.value) == "Test error"

    def test_circleseeker_error_empty_message(self):
        """Test CircleSeekerError with empty message."""
        error = CircleSeekerError("")
        assert str(error) == ""

    def test_circleseeker_error_no_message(self):
        """Test CircleSeekerError with no message."""
        error = CircleSeekerError()
        assert str(error) == ""


class TestConfigurationError:
    """Test cases for ConfigurationError."""

    def test_configuration_error_creation(self):
        """Test ConfigurationError creation."""
        error = ConfigurationError("Invalid configuration")
        assert str(error) == "Invalid configuration"
        assert isinstance(error, ConfigurationError)
        assert isinstance(error, CircleSeekerError)

    def test_configuration_error_inheritance(self):
        """Test ConfigurationError inheritance chain."""
        error = ConfigurationError("Test")
        assert isinstance(error, Exception)
        assert isinstance(error, CircleSeekerError)
        assert isinstance(error, ConfigurationError)

    def test_configuration_error_raise(self):
        """Test raising ConfigurationError."""
        with pytest.raises(ConfigurationError) as exc_info:
            raise ConfigurationError("Config validation failed")

        assert "Config validation failed" in str(exc_info.value)

    def test_configuration_error_catch_as_base(self):
        """Test catching ConfigurationError as base CircleSeekerError."""
        with pytest.raises(CircleSeekerError):
            raise ConfigurationError("Config error")


class TestExternalToolError:
    """Test cases for ExternalToolError."""

    def test_external_tool_error_basic(self):
        """Test ExternalToolError basic creation."""
        error = ExternalToolError("Tool failed")
        assert str(error) == "Tool failed"
        assert isinstance(error, ExternalToolError)
        assert isinstance(error, CircleSeekerError)

    def test_external_tool_error_with_command(self):
        """Test ExternalToolError with command information."""
        command = ["blast", "-query", "input.fasta"]
        error = ExternalToolError("BLAST failed", command=command)

        assert str(error) == "BLAST failed"
        assert error.command == command
        assert error.returncode is None
        assert error.stderr is None

    def test_external_tool_error_with_all_params(self):
        """Test ExternalToolError with all parameters."""
        command = ["minimap2", "ref.fa", "reads.fq"]
        returncode = 1
        stderr = "Error: file not found"

        error = ExternalToolError(
            "Minimap2 execution failed",
            command=command,
            returncode=returncode,
            stderr=stderr
        )

        assert str(error) == "Minimap2 execution failed"
        assert error.command == command
        assert error.returncode == returncode
        assert error.stderr == stderr

    def test_external_tool_error_inheritance(self):
        """Test ExternalToolError inheritance."""
        error = ExternalToolError("Tool error")
        assert isinstance(error, Exception)
        assert isinstance(error, CircleSeekerError)
        assert isinstance(error, ExternalToolError)

    def test_external_tool_error_raise(self):
        """Test raising ExternalToolError."""
        command = ["samtools", "sort", "input.bam"]

        with pytest.raises(ExternalToolError) as exc_info:
            raise ExternalToolError("Samtools failed", command=command, returncode=2)

        error = exc_info.value
        assert str(error) == "Samtools failed"
        assert error.command == command
        assert error.returncode == 2


class TestPipelineError:
    """Test cases for PipelineError."""

    def test_pipeline_error_creation(self):
        """Test PipelineError creation."""
        error = PipelineError("Pipeline step failed")
        assert str(error) == "Pipeline step failed"
        assert isinstance(error, PipelineError)
        assert isinstance(error, CircleSeekerError)

    def test_pipeline_error_inheritance(self):
        """Test PipelineError inheritance."""
        error = PipelineError("Test")
        assert isinstance(error, Exception)
        assert isinstance(error, CircleSeekerError)
        assert isinstance(error, PipelineError)

    def test_pipeline_error_raise(self):
        """Test raising PipelineError."""
        with pytest.raises(PipelineError) as exc_info:
            raise PipelineError("Step 3 failed: missing input file")

        assert "Step 3 failed" in str(exc_info.value)


class TestValidationError:
    """Test cases for ValidationError."""

    def test_validation_error_creation(self):
        """Test ValidationError creation."""
        error = ValidationError("Data validation failed")
        assert str(error) == "Data validation failed"
        assert isinstance(error, ValidationError)
        assert isinstance(error, CircleSeekerError)

    def test_validation_error_inheritance(self):
        """Test ValidationError inheritance."""
        error = ValidationError("Test")
        assert isinstance(error, Exception)
        assert isinstance(error, CircleSeekerError)
        assert isinstance(error, ValidationError)

    def test_validation_error_raise(self):
        """Test raising ValidationError."""
        with pytest.raises(ValidationError) as exc_info:
            raise ValidationError("Invalid data format in column 'chr'")

        assert "Invalid data format" in str(exc_info.value)


class TestFileFormatError:
    """Test cases for FileFormatError."""

    def test_file_format_error_creation(self):
        """Test FileFormatError creation."""
        error = FileFormatError("Invalid FASTA format")
        assert str(error) == "Invalid FASTA format"
        assert isinstance(error, FileFormatError)
        assert isinstance(error, CircleSeekerError)

    def test_file_format_error_inheritance(self):
        """Test FileFormatError inheritance."""
        error = FileFormatError("Test")
        assert isinstance(error, Exception)
        assert isinstance(error, CircleSeekerError)
        assert isinstance(error, FileFormatError)

    def test_file_format_error_raise(self):
        """Test raising FileFormatError."""
        with pytest.raises(FileFormatError) as exc_info:
            raise FileFormatError("Cannot parse BED file: missing columns")

        assert "Cannot parse BED file" in str(exc_info.value)


class TestExceptionHierarchy:
    """Test cases for exception hierarchy and interactions."""

    def test_catch_all_with_base_exception(self):
        """Test catching all CircleSeeker exceptions with base class."""
        exceptions_to_test = [
            ConfigurationError("Config error"),
            ExternalToolError("Tool error"),
            PipelineError("Pipeline error"),
            ValidationError("Validation error"),
            FileFormatError("Format error")
        ]

        for exc in exceptions_to_test:
            with pytest.raises(CircleSeekerError):
                raise exc

    def test_specific_exception_catching(self):
        """Test catching specific exception types."""
        # Test that we can catch specific exceptions
        with pytest.raises(ConfigurationError):
            raise ConfigurationError("Specific config error")

        with pytest.raises(ExternalToolError):
            raise ExternalToolError("Specific tool error")

        with pytest.raises(PipelineError):
            raise PipelineError("Specific pipeline error")

    def test_exception_type_checking(self):
        """Test exception type checking."""
        config_error = ConfigurationError("Config")
        tool_error = ExternalToolError("Tool")

        # Same type
        assert isinstance(config_error, ConfigurationError)
        assert isinstance(tool_error, ExternalToolError)

        # Different types
        assert not isinstance(config_error, ExternalToolError)
        assert not isinstance(tool_error, ConfigurationError)

        # Base type
        assert isinstance(config_error, CircleSeekerError)
        assert isinstance(tool_error, CircleSeekerError)

    def test_exception_message_preservation(self):
        """Test that exception messages are preserved through inheritance."""
        test_cases = [
            (ConfigurationError, "Configuration failed"),
            (ExternalToolError, "External tool failed"),
            (PipelineError, "Pipeline failed"),
            (ValidationError, "Validation failed"),
            (FileFormatError, "File format failed")
        ]

        for exception_class, message in test_cases:
            try:
                raise exception_class(message)
            except CircleSeekerError as e:
                assert str(e) == message
                assert isinstance(e, exception_class)


class TestExternalToolErrorSpecialFeatures:
    """Test special features of ExternalToolError."""

    def test_external_tool_error_attributes_default(self):
        """Test ExternalToolError default attribute values."""
        error = ExternalToolError("Simple error")
        assert error.command is None
        assert error.returncode is None
        assert error.stderr is None

    def test_external_tool_error_attributes_access(self):
        """Test ExternalToolError attribute access."""
        command = ["tool", "--arg", "value"]
        returncode = 127
        stderr = "command not found"

        error = ExternalToolError(
            "Command failed",
            command=command,
            returncode=returncode,
            stderr=stderr
        )

        # Test attribute access
        assert error.command == command
        assert error.returncode == returncode
        assert error.stderr == stderr

        # Test that we can modify attributes
        error.returncode = 1
        assert error.returncode == 1

    def test_external_tool_error_with_complex_command(self):
        """Test ExternalToolError with complex command structure."""
        complex_command = [
            "blast", "-query", "input.fasta", "-db", "reference.fa",
            "-out", "results.txt", "-evalue", "1e-5", "-num_threads", "8"
        ]

        error = ExternalToolError("BLAST execution failed", command=complex_command)
        assert error.command == complex_command
        assert len(error.command) == 10

    def test_external_tool_error_string_representation(self):
        """Test ExternalToolError string representation."""
        error1 = ExternalToolError("Simple error")
        assert str(error1) == "Simple error"

        error2 = ExternalToolError("Complex error", command=["cmd"], returncode=1)
        assert str(error2) == "Complex error"  # Additional info not in str representation