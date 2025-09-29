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
                raise PipelineError(f"varlociraptor call variants failed: {e.stderr}")

    def filter_calls_fdr_local_smart(
        self,
        input_calls_bcf: Path,
        output_calls_fdr_bcf: Path,
        *,
        fdr: float = 0.05,
        memory_limit: str = "4G",
    ) -> None:
        """Filter calls with control-FDR, decode phred, and sort to BCF.

        Implements the pipeline:
          varlociraptor filter-calls control-fdr ... input | \
          varlociraptor decode-phred | \
          bcftools sort -m MEM -Ob -o OUTPUT -
        """
        output_calls_fdr_bcf.parent.mkdir(parents=True, exist_ok=True)

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

        # Ensure upstream pipes close
        if p1.stdout:
            p1.stdout.close()
        if p2.stdout:
            p2.stdout.close()

        rc3 = p3.wait()
        rc2 = p2.wait()
        rc1 = p1.wait()

        if rc1 != 0:
            _, err = p1.communicate()
            raise PipelineError(f"varlociraptor filter-calls failed: {err.decode(errors='ignore')}")
        if rc2 != 0:
            _, err = p2.communicate()
            raise PipelineError(f"varlociraptor decode-phred failed: {err.decode(errors='ignore')}")
        if rc3 != 0:
            err = p3.stderr.read().decode(errors="ignore") if p3.stderr else ""
            raise PipelineError(f"bcftools sort failed: {err}")

