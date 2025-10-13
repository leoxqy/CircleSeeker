"""Samtools wrapper."""

from pathlib import Path
from typing import Optional
from circleseeker.external.base import ExternalTool


class Samtools(ExternalTool):
    """Samtools BAM/SAM manipulation."""
    
    tool_name = "samtools"
    
    def sort_bam(self, input_sam: Path, output_bam: Path) -> None:
        """Sort SAM/BAM file."""
        cmd = [
            self.tool_name, "sort",
            "-@", str(self.threads),
            "-o", str(output_bam),
            str(input_sam)
        ]
        
        output_bam.parent.mkdir(parents=True, exist_ok=True)
        stdout, stderr = self.run(cmd, capture_output=True)
        if stdout:
            self.logger.debug(f"samtools sort output: {stdout[:500]}")
        self.logger.info(f"Sorted BAM saved to: {output_bam}")
    
    def index_bam(self, bam_file: Path) -> None:
        """Index BAM file."""
        cmd = [
            self.tool_name, "index",
            str(bam_file)
        ]

        stdout, stderr = self.run(cmd, capture_output=True)
        if stdout:
            self.logger.debug(f"samtools index output: {stdout[:500]}")
        self.logger.info(f"BAM index created: {bam_file}.bai")

    def faidx(self, reference_fasta: Path) -> None:
        """Create FASTA index (.fai) for a reference."""
        cmd = [
            self.tool_name,
            "faidx",
            str(reference_fasta),
        ]
        stdout, stderr = self.run(cmd, capture_output=True)
        if stdout:
            self.logger.debug(f"samtools faidx output: {stdout[:500]}")
        self.logger.info(f"Reference index created: {reference_fasta}.fai")
