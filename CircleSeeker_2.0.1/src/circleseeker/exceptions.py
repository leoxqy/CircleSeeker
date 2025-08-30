"""Custom exceptions for CircleSeeker."""


class CircleSeekerError(Exception):
    """Base exception for CircleSeeker."""
    pass


class ConfigurationError(CircleSeekerError):
    """Configuration validation error."""
    pass


class ExternalToolError(CircleSeekerError):
    """External tool execution error."""
    
    def __init__(self, message, command=None, returncode=None, stderr=None):
        super().__init__(message)
        self.command = command
        self.returncode = returncode
        self.stderr = stderr


class PipelineError(CircleSeekerError):
    """Pipeline execution error."""
    pass


class ValidationError(CircleSeekerError):
    """Data validation error."""
    pass


class FileFormatError(CircleSeekerError):
    """File format or parsing error."""
    pass
