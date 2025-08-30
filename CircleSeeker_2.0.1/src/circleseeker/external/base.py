"""Base class for external tool execution."""

from __future__ import annotations

import subprocess
import logging
import shutil
import re
from pathlib import Path
from typing import Sequence, Optional, Tuple, Dict, Any
from packaging import version
from circleseeker.exceptions import ExternalToolError


class ExternalTool:
    """Base class for external tool wrappers."""
    
    tool_name: str = ""
    required_version: Optional[str] = None
    version_command: Optional[str] = "--version"
    version_regex: Optional[str] = r"(\d+\.\d+(?:\.\d+)*)"
    
    def __init__(self, logger: Optional[logging.Logger] = None, threads: int = 1):
        self.threads = threads
        self.logger = logger or logging.getLogger(f"{__name__}.{self.tool_name}")
        self._check_installation()
    
    def check_tool_availability(self, tool_name: str) -> bool:
        """Check if a tool is available in PATH."""
        return shutil.which(tool_name) is not None
    
    def run_command(self, cmd: Sequence[str], **kwargs) -> subprocess.CompletedProcess:
        """Run a command and return the completed process."""
        return subprocess.run(cmd, **kwargs)
    
    def _check_installation(self) -> None:
        """Check if the tool is installed and meets version requirements."""
        if not self.check_tool_availability(self.tool_name):
            raise ExternalToolError(
                f"{self.tool_name} not found in PATH. "
                f"Please install it via: conda install {self.tool_name}"
            )
        
        if self.required_version:
            current_version = self.get_tool_version(self.tool_name)
            if current_version and not self.check_minimum_version(current_version, self.required_version):
                raise ExternalToolError(
                    f"{self.tool_name} version {current_version} is below required version {self.required_version}"
                )
            self.logger.debug(f"{self.tool_name} version: {current_version}")
        
        # Check additional requirements if defined
        additional_tools = getattr(self, '_get_required_tools', lambda: [])()
        for tool in additional_tools:
            if not self.check_tool_availability(tool):
                raise ExternalToolError(
                    f"Required dependency '{tool}' not found for {self.tool_name}"
                )
    
    def get_tool_version(self, tool_name: str) -> Optional[str]:
        """Get tool version string."""
        if not self.version_command:
            return None
            
        try:
            # Try different version command formats
            version_commands = [
                [tool_name, self.version_command],
                [tool_name, "-v"],
                [tool_name, "version"],
                [tool_name, "--help"]  # Some tools show version in help
            ]
            
            for cmd in version_commands:
                try:
                    result = subprocess.run(
                        cmd,
                        capture_output=True,
                        text=True,
                        check=False,
                        timeout=10
                    )
                    
                    output = result.stdout + result.stderr
                    if self.version_regex:
                        match = re.search(self.version_regex, output)
                        if match:
                            return match.group(1)
                    
                    # Fallback: return first line that looks like version
                    for line in output.split('\n'):
                        line = line.strip()
                        if re.search(r'\d+\.\d+', line):
                            return line
                            
                except (subprocess.TimeoutExpired, OSError):
                    continue
                    
        except Exception as e:
            self.logger.debug(f"Could not get version for {tool_name}: {e}")
            
        return None
    
    def check_minimum_version(self, current_version: str, required_version: str) -> bool:
        """Check if current version meets minimum requirement."""
        try:
            # Extract version numbers using regex
            current_match = re.search(r'(\d+\.\d+(?:\.\d+)*)', current_version)
            required_match = re.search(r'(\d+\.\d+(?:\.\d+)*)', required_version)
            
            if not current_match or not required_match:
                self.logger.warning(f"Could not parse versions: {current_version} vs {required_version}")
                return True  # Assume OK if we can't parse
            
            current_ver = version.parse(current_match.group(1))
            required_ver = version.parse(required_match.group(1))
            
            return current_ver >= required_ver
            
        except Exception as e:
            self.logger.debug(f"Version comparison failed: {e}")
            return True  # Assume OK on error
    
    def get_tool_info(self) -> Dict[str, Any]:
        """Get comprehensive tool information."""
        return {
            'name': self.tool_name,
            'available': self.check_tool_availability(self.tool_name),
            'version': self.get_tool_version(self.tool_name),
            'required_version': self.required_version,
            'path': shutil.which(self.tool_name)
        }
    
    def run(
        self,
        cmd: Sequence[str],
        cwd: Optional[Path] = None,
        check: bool = True,
        capture_output: bool = True,
        timeout: Optional[int] = None
    ) -> Tuple[str, str]:
        """Execute command with enhanced error handling."""
        cmd_str = " ".join(str(c) for c in cmd)
        self.logger.info(f"Running: {cmd_str}")
        
        try:
            result = subprocess.run(
                cmd,
                cwd=cwd,
                capture_output=capture_output,
                text=True,
                check=check,
                timeout=timeout
            )
            
            if result.stderr and not result.returncode:
                self.logger.debug(f"Command stderr: {result.stderr[:500]}")
            
            if capture_output:
                return result.stdout, result.stderr
            return "", ""
            
        except subprocess.TimeoutExpired as e:
            self.logger.error(f"Command timed out after {timeout}s: {cmd_str}")
            raise ExternalToolError(
                f"{self.tool_name} timed out",
                command=cmd,
                returncode=-1,
                stderr=f"Process timed out after {timeout} seconds"
            )
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Command failed: {cmd_str}")
            self.logger.error(f"Return code: {e.returncode}")
            self.logger.error(f"Error: {e.stderr[:1000] if e.stderr else 'No error output'}")
            raise ExternalToolError(
                f"{self.tool_name} failed",
                command=cmd,
                returncode=e.returncode,
                stderr=e.stderr
            )
        except OSError as e:
            self.logger.error(f"OS error running command: {cmd_str}")
            self.logger.error(f"Error: {e}")
            raise ExternalToolError(
                f"Failed to execute {self.tool_name}",
                command=cmd,
                returncode=-1,
                stderr=str(e)
            )
