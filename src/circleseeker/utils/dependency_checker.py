"""Dependency checker for CircleSeeker.

Performs comprehensive pre-flight checks for all required and optional tools.
"""

import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from typing import Optional

from circleseeker.utils.logging import get_logger
from circleseeker.exceptions import DependencyError


@dataclass
class Tool:
    """Tool dependency definition."""

    name: str
    required: bool
    purpose: str
    install_hint: str
    min_version: Optional[str] = None
    check_command: Optional[str] = None
    alt_names: Optional[list[str]] = None  # Alternative executable names


def find_tool(name: str, alt_names: Optional[list[str]] = None) -> Optional[str]:
    """Find a tool by name, checking alternative names if provided.

    Args:
        name: Primary tool name
        alt_names: Alternative names to check

    Returns:
        The name that was found, or None if not found
    """
    if shutil.which(name) is not None:
        return name
    if alt_names:
        for alt in alt_names:
            if shutil.which(alt) is not None:
                return alt
    return None


def get_tool_version(tool_name: str, version_arg: str = "--version") -> Optional[str]:
    """Get version string from a tool.

    Args:
        tool_name: Name of the tool executable
        version_arg: Argument to get version (default: --version)

    Returns:
        Version string if found, None otherwise
    """
    try:
        result = subprocess.run(
            [tool_name, version_arg],
            capture_output=True,
            text=True,
            timeout=5,
        )
        output = result.stdout + result.stderr
        # Common version patterns
        match = re.search(r"(\d+\.\d+(?:\.\d+)?)", output)
        if match:
            return match.group(1)
    except (OSError, subprocess.SubprocessError, subprocess.TimeoutExpired):
        pass
    return None


def _basic_version_compare(current: str, minimum: str) -> bool:
    """Basic version comparison fallback.

    Extracts numeric parts and compares them as tuples.

    Args:
        current: Current version string
        minimum: Minimum required version string

    Returns:
        True if current >= minimum, False if parsing fails (conservative)
    """
    def parse(v: str) -> tuple[int, ...]:
        # Extract first 3 numeric components
        return tuple(int(x) for x in re.findall(r"\d+", v)[:3])

    try:
        return parse(current) >= parse(minimum)
    except (ValueError, TypeError):
        return False  # Conservative: parsing failure means requirement not met


def compare_versions(current: str, minimum: str) -> bool:
    """Compare version strings with robust handling.

    Uses packaging.version if available for proper semantic versioning,
    falls back to basic comparison for pre-release versions or when
    packaging is not available.

    Args:
        current: Current version string (e.g., "2.24.0", "1.17-rc1")
        minimum: Minimum required version string

    Returns:
        True if current >= minimum, False otherwise
    """
    try:
        from packaging import version

        try:
            return version.parse(current) >= version.parse(minimum)
        except version.InvalidVersion:
            # Fall back to basic comparison for unusual version formats
            return _basic_version_compare(current, minimum)
    except ImportError:
        # packaging not available, use basic comparison
        return _basic_version_compare(current, minimum)


# Tool dependency definitions
TOOLS = [
    # Core required tools
    Tool(
        name="TideHunter",
        required=True,
        purpose="Tandem repeat detection",
        install_hint="conda install -c bioconda tidehunter",
        min_version="1.5",
        alt_names=["tidehunter"],  # bioconda may install as lowercase
    ),
    Tool(
        name="minimap2",
        required=True,
        purpose="Sequence alignment",
        install_hint="conda install -c bioconda minimap2",
        min_version="2.24",
    ),
    Tool(
        name="samtools",
        required=True,
        purpose="BAM file processing",
        install_hint="conda install -c bioconda samtools",
        min_version="1.17",
    ),
    Tool(
        name="cd-hit-est",
        required=True,
        purpose="Sequence clustering",
        install_hint="conda install -c bioconda cd-hit",
    ),
    Tool(
        name="lastal",
        required=True,
        purpose="High-accuracy CeccDNA detection (LAST aligner)",
        install_hint="conda install -c bioconda last",
        # Note: cecc_build module checks lastal/lastdb/last-split availability
    ),
    Tool(
        name="bedtools",
        required=True,
        purpose="Genomic interval operations (required by pybedtools)",
        install_hint="conda install -c bioconda bedtools",
    ),
    # Note: SplitReads-Core is built-in and doesn't require external tools
]


class DependencyChecker:
    """Check and report on tool dependencies."""

    def __init__(self, logger=None):
        self.logger = logger or get_logger("dependency_checker")
        self.missing_required: list[Tool] = []
        self.missing_optional: list[Tool] = []
        self.found_tools: list[str] = []
        self.version_warnings: list[str] = []

    def check_all(self, skip_tools: Optional[set[str]] = None) -> bool:
        """Check all dependencies.

        Returns:
            True if all required tools are available
        """
        self.logger.info("Checking dependencies...")

        skip_set = {name.lower() for name in (skip_tools or set())}

        def is_skipped(tool: Tool) -> bool:
            if not skip_set:
                return False
            names = [tool.name]
            if tool.alt_names:
                names.extend(tool.alt_names)
            return any(name.lower() in skip_set for name in names)

        for tool in TOOLS:
            if is_skipped(tool):
                self.logger.debug("Skipping dependency check for %s", tool.name)
                continue
            is_required = tool.required
            found_name = find_tool(tool.name, tool.alt_names)
            if found_name is not None:
                self.found_tools.append(found_name)
                # Check version if min_version is specified
                if tool.min_version:
                    current_version = get_tool_version(found_name)
                    if current_version:
                        if not compare_versions(current_version, tool.min_version):
                            warning = (
                                f"{tool.name}: version {current_version} < "
                                f"recommended {tool.min_version}"
                            )
                            self.version_warnings.append(warning)
                            self.logger.warning(f"⚠ {warning}")
                        else:
                            self.logger.debug(f"✓ {tool.name} v{current_version}")
                    else:
                        self.logger.debug(f"✓ {tool.name} found (version unknown)")
                elif found_name != tool.name:
                    self.logger.debug(f"✓ {tool.name} found as '{found_name}'")
                else:
                    self.logger.debug(f"✓ {tool.name} found")
            else:
                if is_required:
                    self.missing_required.append(tool)
                    self.logger.error(f"✗ {tool.name} not found (REQUIRED)")
                else:
                    self.missing_optional.append(tool)
                    self.logger.warning(f"⚠ {tool.name} not found (optional)")

        # SplitReads-Core is built-in; no external inference tools required.
        return len(self.missing_required) == 0

    def print_report(self) -> None:
        """Print a detailed dependency report."""
        print("\n" + "=" * 70)
        print("CircleSeeker Dependency Check")
        print("=" * 70)

        if self.found_tools:
            print("\n✓ Found tools:")
            for tool_name in sorted(self.found_tools):
                print(f"  - {tool_name}")

        if self.version_warnings:
            print("\n⚠ Version warnings:")
            for warning in self.version_warnings:
                print(f"  - {warning}")

        if self.missing_optional:
            print("\n⚠ Missing optional tools:")
            for tool in self.missing_optional:
                print(f"  - {tool.name}")
                print(f"    Purpose: {tool.purpose}")
                print(f"    Install: {tool.install_hint}")

        if self.missing_required:
            print("\n✗ Missing REQUIRED tools:")
            for tool in self.missing_required:
                print(f"  - {tool.name}")
                print(f"    Purpose: {tool.purpose}")
                print(f"    Install: {tool.install_hint}")

            print("\n" + "=" * 70)
            print("ERROR: Cannot proceed without required dependencies.")
            print("Please install missing tools and try again.")
            print("=" * 70 + "\n")
        else:
            print("\n" + "=" * 70)
            print("✓ All required dependencies satisfied!")
            print("=" * 70 + "\n")

    def raise_if_missing_required(self) -> None:
        """Raise error if required dependencies are missing."""
        if self.missing_required:
            self.print_report()
            raise DependencyError("Missing required dependencies")


def check_dependencies(logger=None, skip_tools: Optional[set[str]] = None) -> bool:
    """Convenience function to check all dependencies.

    Args:
        logger: Optional logger instance
        skip_tools: Optional set of tool names to skip

    Returns:
        True if all required tools are available

    Raises:
        DependencyError: If required dependencies are missing
    """
    checker = DependencyChecker(logger=logger)
    all_ok = checker.check_all(skip_tools=skip_tools)

    if not all_ok:
        checker.print_report()
        raise DependencyError("Missing required dependencies")

    # Print warnings for missing optional tools
    if checker.missing_optional:
        checker.print_report()

    return True


if __name__ == "__main__":
    # Can run as standalone tool
    try:
        check_dependencies()
    except DependencyError:
        sys.exit(1)
