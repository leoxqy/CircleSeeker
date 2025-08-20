# src/circleseeker2/external/base.py
"""Base class for external tool wrappers."""

import logging
import shutil
import subprocess
from abc import ABC, abstractmethod
from typing import Sequence, Tuple, Optional

from circleseeker2.config import Config
from circleseeker2.exceptions import ExternalToolError

logger = logging.getLogger(__name__)


class ExternalTool(ABC):
    """Abstract base class for wrapping external command-line tools."""

    def __init__(self, tool_name: str, cfg: Config):
        self.tool_name = tool_name
        self.cfg = cfg
        self.threads = cfg.performance.threads
        self._check_installation()

    @abstractmethod
    def get_version(self) -> Optional[str]:
        """Get the version of the external tool."""
        pass

    def _check_installation(self) -> None:
        """Check if the tool is available in the system's PATH."""
        if not shutil.which(self.tool_name):
            raise ExternalToolError(
                f"Tool '{self.tool_name}' not found in PATH. "
                f"Please ensure it is installed and accessible."
            )
        logger.debug(f"Found external tool: {self.tool_name}")

    def run(
        self, 
        cmd: Sequence[str],
        check: bool = True, 
        capture_output: bool = True
    ) -> Tuple[str, str]:
        """
        Run a command using subprocess.

        Args:
            cmd: The command to execute as a sequence of strings.
            check: If True, raise an exception on non-zero exit codes.
            capture_output: If True, capture and return stdout and stderr.

        Returns:
            A tuple containing stdout and stderr as strings.
        """
        logger.info(f"Running command: {' '.join(cmd)}")
        try:
            process = subprocess.run(
                cmd,
                check=check,
                capture_output=capture_output,
                text=True,
                encoding='utf-8'
            )
            stdout = process.stdout or ""
            stderr = process.stderr or ""
            
            if stdout:
                logger.debug(f"[{self.tool_name} stdout]\n{stdout.strip()}")
            if stderr:
                logger.debug(f"[{self.tool_name} stderr]\n{stderr.strip()}")
                
            return stdout, stderr
        except FileNotFoundError:
            raise ExternalToolError(f"Command '{cmd[0]}' not found.")
        except subprocess.CalledProcessError as e:
            raise ExternalToolError(
                f"Error executing {self.tool_name}",
                command=cmd,
                returncode=e.returncode,
                stderr=e.stderr
            ) from e
