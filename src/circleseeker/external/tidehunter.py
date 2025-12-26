"""TideHunter wrapper."""

import shutil
from pathlib import Path
from circleseeker.external.base import ExternalTool


def _find_tidehunter_executable() -> str:
    """Find the TideHunter executable, checking both name variants.

    Returns:
        The executable name found in PATH

    Raises:
        FileNotFoundError: If TideHunter is not found
    """
    for name in ("TideHunter", "tidehunter"):
        if shutil.which(name) is not None:
            return name
    raise FileNotFoundError(
        "TideHunter not found. Install via: conda install -c bioconda tidehunter"
    )


class TideHunter(ExternalTool):
    """TideHunter tandem repeat finder."""

    # Detect the actual binary name (case-sensitive on some systems)
    tool_name = "TideHunter"  # Default, will be updated in __init__

    def __init__(self, threads: int = 1, **kwargs):
        """Initialize TideHunter wrapper.

        Args:
            threads: Number of threads to use (default: 1)
            **kwargs: Additional arguments for ExternalTool
        """
        # Detect actual executable name before super().__init__
        try:
            self.tool_name = _find_tidehunter_executable()
        except FileNotFoundError:
            pass  # Will fail in _check_installation
        super().__init__(threads=threads, **kwargs)

    def run_analysis(
        self,
        input_file: Path,
        output_file: Path,
        k: int = 16,
        w: int = 1,
        p: int = 100,
        P: int = 2000000,
        e: float = 0.1,
        f: int = 2,
    ) -> None:
        """Run TideHunter analysis."""
        cmd = [
            self.tool_name,
            "-f",
            str(f),
            "-t",
            str(self.threads),
            "-k",
            str(k),
            "-w",
            str(w),
            "-p",
            str(p),
            "-P",
            str(P),
            "-e",
            str(e),
            str(input_file),
        ]

        # Create output directory
        output_file.parent.mkdir(parents=True, exist_ok=True)

        # Run TideHunter with output redirection
        with open(output_file, "w") as out_handle:
            import subprocess

            try:
                subprocess.run(
                    cmd, stdout=out_handle, stderr=subprocess.PIPE, text=True, check=True
                )

                self.logger.info("TideHunter completed successfully")
                self.logger.info(f"Output saved to: {output_file}")

                # Log output file size
                output_size = output_file.stat().st_size
                self.logger.debug(f"Output file size: {output_size} bytes")

            except subprocess.CalledProcessError as e:
                self.logger.error(f"TideHunter failed with exit code {e.returncode}")
                self.logger.error(f"Error message: {e.stderr}")
                from circleseeker.exceptions import ExternalToolError

                raise ExternalToolError(
                    "TideHunter failed", command=cmd, returncode=e.returncode, stderr=e.stderr
                )


class TideHunterRunner(TideHunter):
    """Backward compatibility wrapper for existing code."""

    def __init__(self, num_threads=8):
        super().__init__(threads=num_threads)
        self.num_threads = num_threads

    def run(self, input_fasta, output_path):
        """Run TideHunter with legacy interface."""
        return self.run_analysis(input_file=Path(input_fasta), output_file=Path(output_path))
