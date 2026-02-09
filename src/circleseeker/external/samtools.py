"""Samtools wrapper."""

from pathlib import Path
from circleseeker.external.base import ExternalTool


class Samtools(ExternalTool):
    """Samtools BAM/SAM manipulation."""

    tool_name = "samtools"

    def _write_log(self, log_file: Path, stdout: str, stderr: str) -> None:
        """Append stdout/stderr to a log file."""
        with open(log_file, "a") as f:
            if stdout:
                f.write(stdout)
            if stderr:
                f.write(stderr)

    def sort_bam(self, input_sam: Path, output_bam: Path, log_prefix: str = "") -> None:
        """Sort SAM/BAM file."""
        cmd = [
            self.tool_name,
            "sort",
            "-@",
            str(self.threads),
            "-o",
            str(output_bam),
            str(input_sam),
        ]

        output_bam.parent.mkdir(parents=True, exist_ok=True)
        stdout, stderr = self.run(cmd, capture_output=True)
        log_name = f"{log_prefix}_samtools.log" if log_prefix else "samtools.log"
        self._write_log(output_bam.parent / log_name, stdout, stderr)
        if stdout:
            self.logger.debug(f"samtools sort output: {stdout[:500]}")
        self.logger.debug(f"Sorted BAM saved to: {output_bam}")

    def index_bam(self, bam_file: Path, log_prefix: str = "") -> None:
        """Index BAM file."""
        cmd = [self.tool_name, "index", str(bam_file)]

        stdout, stderr = self.run(cmd, capture_output=True)
        log_name = f"{log_prefix}_samtools.log" if log_prefix else "samtools.log"
        self._write_log(bam_file.parent / log_name, stdout, stderr)
        if stdout:
            self.logger.debug(f"samtools index output: {stdout[:500]}")
        self.logger.debug(f"BAM index created: {bam_file}.bai")

    def faidx(self, reference_fasta: Path, log_prefix: str = "") -> None:
        """Create FASTA index (.fai) for a reference."""
        cmd = [
            self.tool_name,
            "faidx",
            str(reference_fasta),
        ]
        stdout, stderr = self.run(cmd, capture_output=True)
        log_name = f"{log_prefix}_samtools.log" if log_prefix else "samtools.log"
        self._write_log(reference_fasta.parent / log_name, stdout, stderr)
        if stdout:
            self.logger.debug(f"samtools faidx output: {stdout[:500]}")
        self.logger.debug(f"Reference index created: {reference_fasta}.fai")
