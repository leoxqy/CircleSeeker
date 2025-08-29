"""
Mosdepth - Depth Analysis Tool Wrapper

This module provides a wrapper for the mosdepth tool used for calculating
depth coverage across genomic regions. Mosdepth is efficient for large BAM
files and provides windowed depth calculations.

Key features:
- Fast depth calculation using mosdepth
- Window-based depth analysis
- Automatic decompression of output files
- BAM indexing integration with samtools
- Comprehensive error handling and logging

Migrated from step11_run_mosdepth.py to the new CircleSeeker2 architecture.
"""

from __future__ import annotations

import gzip
import logging
import subprocess
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from circleseeker2.exceptions import PipelineError
from circleseeker2.external.base import ExternalTool


@dataclass
class MosdepthConfig:
   """Configuration for Mosdepth depth analysis."""
   # Analysis parameters
   window_size: int = 50
   threads: int = 8
   
   # Output options
   keep_intermediate: bool = False
   keep_csi: bool = False
   
   # Quality filters (applied during analysis)
   min_mapping_quality: int = 0
   include_flag: Optional[int] = None
   exclude_flag: Optional[int] = None


class Mosdepth(ExternalTool):
   """Wrapper for mosdepth depth calculation tool."""
   
   tool_name = "mosdepth"  # Set as class attribute for base class
   
   def __init__(self, config: Optional[MosdepthConfig] = None,
                logger: Optional[logging.Logger] = None):
       """Initialize Mosdepth wrapper."""
       super().__init__(logger)
       self.config = config or MosdepthConfig()
       self._check_samtools()
   
   def _check_samtools(self) -> None:
       """Check if samtools is available for BAM indexing."""
       self.has_samtools = shutil.which("samtools") is not None
       if not self.has_samtools:
           self.logger.warning("samtools not found - BAM indexing may fail")
   
   def _check_bam_index(self, bam_file: Path) -> None:
       """Check and create BAM index if needed."""
       bai_file = Path(f"{bam_file}.bai")
       
       if not bai_file.exists():
           if not self.has_samtools:
               self.logger.warning(f"BAM index not found and samtools not available: {bai_file}")
               return
               
           self.logger.info("Creating BAM index...")
           try:
               subprocess.run(
                   ["samtools", "index", str(bam_file)],
                   check=True,
                   stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL
               )
               self.logger.debug(f"Created BAM index: {bai_file}")
           except subprocess.CalledProcessError as e:
               raise PipelineError(f"Failed to create BAM index: {e}")
   
   def _build_command(self, bam_file: Path, output_prefix: str) -> list[str]:
       """Build mosdepth command."""
       command = [
           "mosdepth",
           "-t", str(self.config.threads),
           "--by", str(self.config.window_size)
       ]
       
       # Add quality filters if specified
       if self.config.min_mapping_quality > 0:
           command.extend(["-Q", str(self.config.min_mapping_quality)])
       
       if self.config.include_flag is not None:
           command.extend(["--include-flag", str(self.config.include_flag)])
       
       if self.config.exclude_flag is not None:
           command.extend(["--exclude-flag", str(self.config.exclude_flag)])
       
       # Add output prefix and input BAM
       command.extend([output_prefix, str(bam_file)])
       
       return command
   
   def _decompress_per_base(self, output_prefix: str) -> Path:
       """Decompress per-base.bed.gz file."""
       per_base_gz = Path(f"{output_prefix}.per-base.bed.gz")
       per_base_bed = Path(f"{output_prefix}.per-base.bed")
       
       if per_base_gz.exists():
           self.logger.debug("Decompressing per-base.bed.gz...")
           try:
               with gzip.open(per_base_gz, 'rt') as gz_file:
                   with open(per_base_bed, 'w') as out_file:
                       for line in gz_file:
                           out_file.write(line)
               self.logger.debug(f"Decompressed to: {per_base_bed}")
           except Exception as e:
               raise PipelineError(f"Failed to decompress per-base file: {e}")
       elif not per_base_bed.exists():
           raise PipelineError("Mosdepth per-base output not found")
       
       return per_base_bed
   
   def _cleanup_files(self, output_prefix: str) -> None:
       """Clean up intermediate mosdepth files."""
       if self.config.keep_intermediate:
           return
       
       # Files to potentially remove
       files_to_remove = [
           f"{output_prefix}.per-base.bed.gz",
           f"{output_prefix}.regions.bed.gz",
           f"{output_prefix}.mosdepth.global.dist.txt",
           f"{output_prefix}.mosdepth.region.dist.txt",
           f"{output_prefix}.mosdepth.summary.txt"
       ]
       
       # Include CSI files unless requested to keep
       if not self.config.keep_csi:
           files_to_remove.extend([
               f"{output_prefix}.per-base.bed.gz.csi",
               f"{output_prefix}.regions.bed.gz.csi"
           ])
       
       removed_count = 0
       for file_path in files_to_remove:
           path = Path(file_path)
           if path.exists():
               try:
                   path.unlink()
                   removed_count += 1
                   self.logger.debug(f"Removed: {file_path}")
               except Exception as e:
                   self.logger.warning(f"Failed to remove {file_path}: {e}")
       
       if removed_count > 0:
           self.logger.debug(f"Cleaned up {removed_count} intermediate files")
   
   def calculate_depth(self, bam_file: Path, output_prefix: str,
                      regions_bed: Optional[Path] = None) -> Path:
       """
       Calculate depth coverage using mosdepth.
       
       Args:
           bam_file: Input BAM file
           output_prefix: Output prefix for all files
           regions_bed: Optional BED file to restrict analysis to specific regions
           
       Returns:
           Path to the decompressed per-base BED file
       """
       bam_file = Path(bam_file)
       if not bam_file.exists():
           raise PipelineError(f"BAM file not found: {bam_file}")
       
       # Check and create BAM index
       self._check_bam_index(bam_file)
       
       # Build full output prefix with window size
       full_prefix = f"{output_prefix}.win{self.config.window_size}"
       
       # Build command
       command = self._build_command(bam_file, full_prefix)
       
       # Add regions if specified
       if regions_bed:
           if not regions_bed.exists():
               raise PipelineError(f"Regions BED file not found: {regions_bed}")
           command.extend(["--by", str(regions_bed)])
       
       self.logger.info(f"Running mosdepth: {bam_file.name} -> {full_prefix}")
       self.logger.debug(f"Command: {' '.join(command)}")
       
       # Run mosdepth
       try:
           result = subprocess.run(
               command,
               capture_output=True,
               text=True,
               check=True
           )
           self.logger.debug("Mosdepth completed successfully")
           
           # Log any stderr output (mosdepth info)
           if result.stderr:
               for line in result.stderr.strip().split('\n'):
                   if line:
                       self.logger.debug(f"mosdepth: {line}")
                       
       except subprocess.CalledProcessError as e:
           error_msg = f"Mosdepth failed with exit code {e.returncode}"
           if e.stderr:
               error_msg += f": {e.stderr}"
           raise PipelineError(error_msg)
       
       # Decompress per-base output
       per_base_bed = self._decompress_per_base(full_prefix)
       
       # Clean up intermediate files
       self._cleanup_files(full_prefix)
       
       self.logger.info(f"Depth analysis complete: {per_base_bed}")
       return per_base_bed
   
   def run(self, input_bam: Path, output_prefix: str,
           regions_bed: Optional[Path] = None) -> Path:
       """
       Run depth calculation (convenience method).
       
       This is an alias for calculate_depth to maintain compatibility.
       """
       return self.calculate_depth(input_bam, output_prefix, regions_bed)


def run_mosdepth(input_bam: Path, output_prefix: str, 
               window_size: int = 50, threads: int = 8,
               keep_intermediate: bool = False) -> Path:
   """
   Convenience function for running mosdepth analysis.
   
   This function maintains compatibility with the original step11 interface
   while using the new architecture.
   """
   config = MosdepthConfig(
       window_size=window_size,
       threads=threads,
       keep_intermediate=keep_intermediate
   )
   
   mosdepth = Mosdepth(config)
   return mosdepth.run(input_bam, output_prefix)