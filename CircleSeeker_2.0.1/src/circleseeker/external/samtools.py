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
        self.run(cmd, capture_output=False)
        self.logger.info(f"Sorted BAM saved to: {output_bam}")
    
    def index_bam(self, bam_file: Path) -> None:
        """Index BAM file."""
        cmd = [
            self.tool_name, "index",
            str(bam_file)
        ]
        
        self.run(cmd, capture_output=False)
        self.logger.info(f"BAM index created: {bam_file}.bai")
