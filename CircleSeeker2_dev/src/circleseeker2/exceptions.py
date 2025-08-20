"""Custom exceptions for CircleSeeker2."""

class CircleSeekerError(Exception):
    """Base exception for all CircleSeeker2 errors."""
    pass

class ConfigurationError(CircleSeekerError):
    """Error related to configuration."""
    pass

class ExternalToolError(CircleSeekerError):
    """Error related to an external tool execution."""
    def __init__(self, message, command=None, returncode=None, stderr=None):
        super().__init__(message)
        self.command = command
        self.returncode = returncode
        self.stderr = stderr

    def __str__(self):
        msg = super().__str__()
        if self.command:
            msg += f"\n  Command: {' '.join(map(str, self.command))}"
        if self.returncode is not None:
            msg += f"\n  Return Code: {self.returncode}"
        if self.stderr:
            msg += f"\n  Stderr: {self.stderr[:1000]}" # Truncate for display
        return msg
