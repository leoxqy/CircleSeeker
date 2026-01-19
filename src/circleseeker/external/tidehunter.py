"""TideHunter wrapper."""

import shutil
from pathlib import Path
from typing import Any, Optional

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

    def __init__(self, threads: int = 1, **kwargs: Any) -> None:
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
        c: int = 2,
        out_fmt: int = 2,
        # Legacy parameter name (maps to out_fmt)
        f: Optional[int] = None,
    ) -> None:
        """Run TideHunter analysis.

        Args:
            k: k-mer length (max 16)
            w: window size for minimizer seeding
            p: minimum period size
            P: maximum period size
            e: maximum divergence rate
            c: minimum copy number (TideHunter -c)
            out_fmt: output format (1=FASTA, 2=Tabular, 3=FASTQ, 4=Tabular+qual)
            f: deprecated alias for out_fmt (for backward compatibility)
        """
        import subprocess

        # Handle legacy 'f' parameter (was incorrectly used for output format)
        output_format = f if f is not None else out_fmt

        cmd = [
            self.tool_name,
            "-f",
            str(output_format),
            "-c",
            str(c),
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

        log_file = output_file.parent / "tidehunter.log"

        def _read_tail(path: Path, max_bytes: int = 32_000) -> str:
            try:
                with open(path, "rb") as handle:
                    try:
                        handle.seek(0, 2)
                        size = handle.tell()
                        offset = max(0, size - max_bytes)
                        handle.seek(offset)
                    except OSError:
                        pass
                    data = handle.read()
                return data.decode(errors="replace")
            except OSError:
                return ""

        # Run TideHunter with stdout redirected to output file and stderr to a persistent log.
        try:
            with open(output_file, "w") as out_handle, open(log_file, "w") as log_handle:
                subprocess.run(cmd, stdout=out_handle, stderr=log_handle, text=True, check=True)

            self.logger.info("TideHunter completed successfully")
            self.logger.info(f"Output saved to: {output_file}")

            # Log output file size
            output_size = output_file.stat().st_size
            self.logger.debug(f"Output file size: {output_size} bytes")

        except subprocess.CalledProcessError as e:
            stderr_tail = _read_tail(log_file)
            self.logger.error(f"TideHunter failed with exit code {e.returncode}")
            if stderr_tail:
                self.logger.error("Last TideHunter log output:\n%s", stderr_tail[-2000:])
            from circleseeker.exceptions import ExternalToolError

            raise ExternalToolError(
                "TideHunter failed",
                command=cmd,
                returncode=e.returncode,
                stderr=stderr_tail or None,
            ) from e


class TideHunterRunner(TideHunter):
    """Backward compatibility wrapper for existing code."""

    def __init__(self, num_threads: int = 8) -> None:
        super().__init__(threads=num_threads)
        self.num_threads = num_threads

    def run(self, input_fasta: Path, output_path: Path) -> None:  # type: ignore[override]
        """Run TideHunter with legacy interface."""
        return self.run_analysis(input_file=Path(input_fasta), output_file=Path(output_path))
