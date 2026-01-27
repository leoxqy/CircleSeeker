"""Custom exceptions for CircleSeeker."""


class CircleSeekerError(Exception):
    """Base exception for all CircleSeeker errors."""

    pass


class ConfigurationError(CircleSeekerError):
    """Raised when configuration is invalid or missing."""

    pass


class ExternalToolError(CircleSeekerError):
    """Raised when an external tool execution fails."""

    def __init__(self, message="", command=None, returncode=None, stderr=None):
        """Initialize ExternalToolError with optional command details.

        Args:
            message: Error message
            command: Command that was executed (list of strings)
            returncode: Exit code from the command
            stderr: Standard error output from the command
        """
        super().__init__(message)
        self.command = command
        self.returncode = returncode
        self.stderr = stderr


class PipelineError(CircleSeekerError):
    """Raised when a pipeline step fails."""

    pass


class ValidationError(CircleSeekerError):
    """Raised when data validation fails."""

    pass


class FileFormatError(CircleSeekerError):
    """Raised when file format is invalid or unsupported."""

    pass


class DependencyError(CircleSeekerError):
    """Raised when required external dependencies are missing or incompatible."""

    pass
