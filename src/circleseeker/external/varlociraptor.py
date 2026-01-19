"""Varlociraptor external tool wrapper.

Covers the subcommands used in CircleSeeker's cyrcular_calling step.
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Optional

from circleseeker.external.base import ExternalTool
from circleseeker.exceptions import PipelineError


class Varlociraptor(ExternalTool):
    """Wrapper for the `varlociraptor` CLI."""

    tool_name = "varlociraptor"

    def estimate_alignment_properties(self, reference: Path, bam: Path, output_json: Path) -> None:
        """Estimate alignment properties and write JSON to file."""
        output_json.parent.mkdir(parents=True, exist_ok=True)
        cmd = [
            self.tool_name,
            "estimate",
            "alignment-properties",
            str(reference),
            "--bams",
            str(bam),
        ]
        # Capture JSON text and write to file
        stdout, _ = self.run(cmd, capture_output=True)
        output_json.write_text(stdout)

    def preprocess_variants(
        self,
        reference: Path,
        candidates_bcf_sorted: Path,
        alignprops_json: Path,
        bam: Path,
        output_obs_bcf: Path,
        *,
        max_depth: int = 200,
    ) -> None:
        """Preprocess variants to observation BCF."""
        output_obs_bcf.parent.mkdir(parents=True, exist_ok=True)
        cmd = [
            self.tool_name,
            "preprocess",
            "variants",
            str(reference),
            "--candidates",
            str(candidates_bcf_sorted),
            "--alignment-properties",
            str(alignprops_json),
            "--max-depth",
            str(max_depth),
            "--bam",
            str(bam),
            "--output",
            str(output_obs_bcf),
        ]
        stdout, stderr = self.run(cmd, capture_output=True)
        if stderr:
            self.logger.debug(f"varlociraptor preprocess stderr: {stderr[:500]}")
        self.logger.info(f"Preprocessed variants saved to: {output_obs_bcf}")

    def call_variants_generic(
        self,
        obs_bcf_sorted: Path,
        sample_name: str,
        scenario_yaml: Path,
        output_calls_bcf: Path,
    ) -> None:
        """Call variants (generic) and write BCF to file safely (no shell redirects)."""
        output_calls_bcf.parent.mkdir(parents=True, exist_ok=True)
        cmd = [
            self.tool_name,
            "call",
            "variants",
            "generic",
            "--obs",
            f"{sample_name}={obs_bcf_sorted}",
            "--scenario",
            str(scenario_yaml),
        ]
        # Write binary BCF directly via file handle
        with open(output_calls_bcf, "wb") as fout:
            try:
                subprocess.run(cmd, stdout=fout, stderr=subprocess.PIPE, check=True)
            except subprocess.CalledProcessError as e:
                err_text = e.stderr.decode(errors="ignore") if e.stderr else ""
                raise PipelineError(f"varlociraptor call variants failed: {err_text}")

    def filter_calls_fdr_local_smart(
        self,
        input_calls_bcf: Path,
        output_calls_fdr_bcf: Path,
        *,
        fdr: float = 0.2,
        memory_limit: str = "4G",
        timeout: Optional[int] = None,
    ) -> None:
        """Filter calls with control-FDR, decode phred, and sort to BCF.

        Implements the pipeline:
          varlociraptor filter-calls control-fdr ... input | \
          varlociraptor decode-phred | \
          bcftools sort -m MEM -Ob -o OUTPUT -

        Args:
            input_calls_bcf: Input BCF file with variant calls
            output_calls_fdr_bcf: Output BCF file path
            fdr: False discovery rate threshold (default: 0.2)
            memory_limit: Memory limit for bcftools sort (default: "4G")
            timeout: Maximum time in seconds for the pipeline (default: None = no timeout)
        """
        from circleseeker.exceptions import ExternalToolError

        output_calls_fdr_bcf.parent.mkdir(parents=True, exist_ok=True)

        p1 = None
        p2 = None
        p3 = None

        def _cleanup_processes() -> None:
            """Terminate and clean up all pipeline processes gracefully."""
            for proc in [p1, p2, p3]:
                if proc is None:
                    continue
                try:  # type: ignore[unreachable]
                    # First try graceful termination
                    proc.terminate()
                    proc.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    # Force kill if graceful termination fails
                    try:
                        proc.kill()
                        proc.wait(timeout=2)
                    except OSError:
                        pass
                except OSError:
                    pass  # Process already terminated
                # Close any remaining pipes
                for pipe in [proc.stdout, proc.stderr, proc.stdin]:
                    if pipe is not None:
                        try:
                            pipe.close()
                        except OSError:
                            pass

        try:
            p1 = subprocess.Popen(
                [
                    self.tool_name,
                    "filter-calls",
                    "control-fdr",
                    "--mode",
                    "local-smart",
                    "--events",
                    "PRESENT",
                    "--var",
                    "BND",
                    "--fdr",
                    str(fdr),
                    str(input_calls_bcf),
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )

            p2 = subprocess.Popen(
                [self.tool_name, "decode-phred"],
                stdin=p1.stdout,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )

            p3 = subprocess.Popen(
                [
                    "bcftools",
                    "sort",
                    "-m",
                    memory_limit,
                    "-O",
                    "b",
                    "-o",
                    str(output_calls_fdr_bcf),
                    "-",
                ],
                stdin=p2.stdout,
                stderr=subprocess.PIPE,
            )

            # Allow p1 to receive SIGPIPE if p2 exits
            if p1.stdout:
                p1.stdout.close()
            # Allow p2 to receive SIGPIPE if p3 exits
            if p2.stdout:
                p2.stdout.close()

            # Wait for processes in reverse order to avoid deadlock
            _, stderr3 = p3.communicate(timeout=timeout)
            _, stderr2 = p2.communicate(timeout=timeout)
            _, stderr1 = p1.communicate(timeout=timeout)

            # Check return codes
            if p1.returncode != 0:
                err_text = stderr1.decode(errors="ignore") if stderr1 else ""
                raise ExternalToolError(
                    f"varlociraptor filter-calls failed: {err_text}",
                    command=[self.tool_name, "filter-calls"],
                    returncode=p1.returncode,
                    stderr=err_text,
                )
            if p2.returncode != 0:
                err_text = stderr2.decode(errors="ignore") if stderr2 else ""
                raise ExternalToolError(
                    f"varlociraptor decode-phred failed: {err_text}",
                    command=[self.tool_name, "decode-phred"],
                    returncode=p2.returncode,
                    stderr=err_text,
                )
            if p3.returncode != 0:
                err_text = stderr3.decode(errors="ignore") if stderr3 else ""
                raise ExternalToolError(
                    f"bcftools sort failed: {err_text}",
                    command=["bcftools", "sort"],
                    returncode=p3.returncode,
                    stderr=err_text,
                )

        except subprocess.TimeoutExpired as e:
            _cleanup_processes()
            raise PipelineError(
                f"Pipeline timed out after {timeout} seconds: {e}"
            ) from e
        except ExternalToolError:
            raise
        except Exception as e:
            _cleanup_processes()
            raise PipelineError(f"Pipeline execution failed: {e}") from e
        finally:
            # Ensure stderr pipes are closed even on success
            for proc in [p1, p2, p3]:
                if proc is not None and proc.stderr is not None:
                    try:
                        proc.stderr.close()
                    except OSError:
                        pass
