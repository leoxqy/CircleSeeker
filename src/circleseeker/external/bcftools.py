"""Bcftools external tool wrapper."""

from __future__ import annotations

from pathlib import Path
from circleseeker.external.base import ExternalTool


class Bcftools(ExternalTool):
    """Wrapper for `bcftools`."""

    tool_name = "bcftools"

    def sort(
        self,
        input_vcf_bcf: Path,
        output_bcf: Path,
        *,
        memory_limit: str = "4G",
    ) -> None:
        """Sort VCF/BCF and write BCF output."""
        output_bcf.parent.mkdir(parents=True, exist_ok=True)
        cmd = [
            self.tool_name,
            "sort",
            "-m",
            memory_limit,
            "-O",
            "b",
            "-o",
            str(output_bcf),
            str(input_vcf_bcf),
        ]
        stdout, stderr = self.run(cmd, capture_output=True)
        if stdout:
            self.logger.debug(f"bcftools sort output: {stdout[:500]}")
        self.logger.info(f"Sorted BCF saved to: {output_bcf}")

    def index(self, bcf_file: Path, *, force: bool = False) -> None:
        cmd = [self.tool_name, "index"]
        if force:
            cmd.append("-f")
        cmd.append(str(bcf_file))
        stdout, stderr = self.run(cmd, capture_output=True)
        if stdout:
            self.logger.debug(f"bcftools index output: {stdout[:500]}")
        self.logger.info(f"BCF index created for: {bcf_file}")
