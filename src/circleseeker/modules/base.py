"""
Base module interface for CircleSeeker pipeline modules.

This module defines the base classes and interfaces that all pipeline modules
should inherit from to ensure consistency and maintainability.
"""

from __future__ import annotations

import logging
import time
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional, Union

from circleseeker.utils.logging import get_logger


@dataclass
class ModuleResult:
    """Standard result container for all modules."""

    success: bool
    module_name: str
    output_files: dict[str, Path] = field(default_factory=dict)
    metrics: dict[str, Any] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    error_message: Optional[str] = None
    execution_time: float = 0.0

    def add_output(self, key: str, path: Union[str, Path]) -> None:
        """Add an output file to the result."""
        self.output_files[key] = Path(path)

    def add_metric(self, key: str, value: Any) -> None:
        """Add a metric to the result."""
        self.metrics[key] = value

    def add_warning(self, message: str) -> None:
        """Add a warning message."""
        self.warnings.append(message)


class ModuleBase(ABC):
    """Base class for all CircleSeeker pipeline modules."""

    def __init__(
        self,
        name: Optional[str] = None,
        logger: Optional[logging.Logger] = None,
        debug: bool = False,
    ):
        """
        Initialize the module.

        Args:
            name: Module name (defaults to class name)
            logger: Logger instance (creates new if None)
            debug: Enable debug mode
        """
        self.name = name or self.__class__.__name__
        self.logger = logger or get_logger(self.name)
        self.debug = debug
        self._start_time: Optional[float] = None

    def validate_input_file(self, file_path: Union[str, Path], file_type: str = "input") -> Path:
        """
        Validate that an input file exists and is readable.

        Args:
            file_path: Path to the file
            file_type: Description of file type for error messages

        Returns:
            Path object for the validated file

        Raises:
            FileNotFoundError: If file doesn't exist
            PermissionError: If file isn't readable
        """
        path = Path(file_path)

        if not path.exists():
            raise FileNotFoundError(f"{file_type} file not found: {path}")

        if not path.is_file():
            raise ValueError(f"{file_type} is not a file: {path}")

        if not path.stat().st_size > 0:
            self.logger.warning(f"{file_type} file is empty: {path}")

        return path

    def validate_output_dir(self, output_dir: Union[str, Path]) -> Path:
        """
        Validate/create output directory.

        Args:
            output_dir: Path to output directory

        Returns:
            Path object for the validated directory
        """
        path = Path(output_dir)
        path.mkdir(parents=True, exist_ok=True)
        return path

    @abstractmethod
    def validate_inputs(self, **kwargs: Any) -> bool:
        """
        Validate all required inputs for the module.

        Returns:
            True if validation passes

        Raises:
            ValueError: If validation fails
        """
        pass

    @abstractmethod
    def execute(self, **kwargs: Any) -> ModuleResult:
        """
        Execute the module's main logic.

        Returns:
            ModuleResult object with outputs and metrics
        """
        pass

    def run(self, **kwargs: Any) -> ModuleResult:
        """
        Main entry point for running the module.

        This method handles:
        1. Input validation
        2. Execution timing
        3. Error handling
        4. Result packaging
        """
        self._start_time = time.time()
        result = ModuleResult(success=False, module_name=self.name)

        try:
            # Validate inputs
            self.logger.info(f"Starting {self.name}")
            self.validate_inputs(**kwargs)

            # Execute main logic
            result = self.execute(**kwargs)
            result.module_name = self.name

            # Calculate execution time
            result.execution_time = time.time() - self._start_time

            if result.success:
                self.logger.info(
                    f"{self.name} completed successfully in " f"{result.execution_time:.2f} seconds"
                )
            else:
                self.logger.error(f"{self.name} failed: {result.error_message}")

            # Log warnings if any
            for warning in result.warnings:
                self.logger.warning(warning)

        except Exception as e:
            result.success = False
            result.error_message = str(e)
            result.execution_time = time.time() - self._start_time
            self.logger.error(f"{self.name} failed with error: {e}", exc_info=True)

        return result

    def cleanup(self, temp_files: Optional[list[Path]] = None) -> None:
        """
        Clean up temporary files.

        Args:
            temp_files: List of temporary files to remove
        """
        if temp_files:
            for file_path in temp_files:
                try:
                    if file_path.exists():
                        file_path.unlink()
                        self.logger.debug(f"Removed temporary file: {file_path}")
                except OSError as e:
                    self.logger.warning(f"Failed to remove {file_path}: {e}")


class ExternalToolModule(ModuleBase):
    """Base class for modules that wrap external tools."""

    def __init__(self, tool_name: str, required_version: Optional[str] = None, **kwargs: Any):
        """
        Initialize external tool module.

        Args:
            tool_name: Name of the external tool
            required_version: Minimum required version
        """
        super().__init__(**kwargs)
        self.tool_name = tool_name
        self.required_version = required_version

    @abstractmethod
    def check_tool_availability(self) -> bool:
        """Check if the external tool is available."""
        pass

    @abstractmethod
    def get_tool_version(self) -> str:
        """Get the version of the external tool."""
        pass

    def validate_tool(self) -> None:
        """Validate that the tool is available and meets version requirements."""
        if not self.check_tool_availability():
            raise RuntimeError(f"{self.tool_name} is not available in PATH")

        if self.required_version:
            version = self.get_tool_version()
            self.logger.debug(f"{self.tool_name} version: {version}")
