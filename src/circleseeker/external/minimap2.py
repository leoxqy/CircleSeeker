"""
Minimap2 wrapper module for CircleSeeker

Provides comprehensive minimap2 alignment functionality including:
- Reference genome indexing
- Sequence alignment with sorting
- BAM index generation
- Alignment statistics

Fixed version: avoids pipe data loss by using temporary files
"""

from __future__ import annotations

import shlex
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional
import logging

from circleseeker.external.base import ExternalTool
from circleseeker.exceptions import ExternalToolError, PipelineError


@dataclass
class MiniMapConfig:
    """Configuration for Minimap2 wrapper."""

    preset: str = "map-hifi"
    threads: int = 8
    output_format: str = "bam"  # "sam" or "bam"
    allow_secondary: bool = True
    build_index: bool = True
    sort_bam: bool = True
    index_bam: bool = True
    additional_args: str = ""


class Minimap2(ExternalTool):
    """Minimap2 wrapper with full BAM processing pipeline."""

    tool_name = "minimap2"

    def __init__(
        self,
        config: Optional[MiniMapConfig] = None,
        logger: Optional[logging.Logger] = None,
        threads: Optional[int] = None,
    ) -> None:
        if config is None:
            cfg_threads = threads if threads is not None else 8
            config = MiniMapConfig(threads=cfg_threads)
        # Set config BEFORE super().__init__ because _check_installation() may access it
        self.config = config
        self.preset = config.preset
        self.allow_secondary = config.allow_secondary
        super().__init__(logger=logger, threads=config.threads)
        # Note: super().__init__ calls _check_installation() via method resolution

    def _check_installation(self) -> None:
        """Check if required tools are available."""
        required_tools = ["minimap2"]
        missing_tools = []

        for tool in required_tools:
            if not shutil.which(tool):
                missing_tools.append(tool)
                self.logger.error(f"{tool} not found in PATH")

        if missing_tools:
            raise PipelineError(f"Required tools not installed: {', '.join(missing_tools)}")

        # Check samtools separately (optional but recommended)
        self.has_samtools = shutil.which("samtools") is not None
        if not self.has_samtools:
            if self.config.sort_bam or self.config.index_bam:
                self.logger.warning("samtools not found - BAM sorting/indexing will be disabled")
                self.config.sort_bam = False
                self.config.index_bam = False

        self.logger.debug("All required tools found in PATH")

    def build_index(self, reference_fasta: Path, index_path: Optional[Path] = None) -> Path:
        """
        Build minimap2 index for reference genome.

        Args:
            reference_fasta: Reference genome FASTA file
            index_path: Output index path (optional)

        Returns:
            Path to index file
        """
        reference_fasta = Path(reference_fasta)

        if index_path is None:
            # Auto-generate index filename (check .fasta first to avoid partial replacement)
            ref_str = str(reference_fasta)
            if ref_str.endswith(".fasta"):
                index_path = Path(ref_str[:-6] + ".mmi")
            elif ref_str.endswith(".fa"):
                index_path = Path(ref_str[:-3] + ".mmi")
            elif ref_str.endswith(".fna"):
                index_path = Path(ref_str[:-4] + ".mmi")
            else:
                index_path = Path(ref_str + ".mmi")
        else:
            index_path = Path(index_path)

        # Skip if index already exists
        if index_path.exists():
            self.logger.info(f"Index already exists: {index_path}")
            return index_path

        command = [
            "minimap2",
            "-t",
            str(self.config.threads),
            "-d",
            str(index_path),
            str(reference_fasta),
        ]

        self.logger.info(f"Building minimap2 index for {reference_fasta}")
        self.logger.debug(f"Command: {' '.join(command)}")

        try:
            subprocess.run(command, stderr=subprocess.PIPE, text=True, check=True)

            self.logger.info(f"Index built successfully: {index_path}")
            return index_path

        except subprocess.CalledProcessError as e:
            self.logger.error(f"minimap2 index building failed with exit code {e.returncode}")
            self.logger.error(f"Error message: {e.stderr}")
            raise ExternalToolError(
                "Failed to build minimap2 index",
                command=e.cmd,
                returncode=e.returncode,
                stderr=e.stderr,
            )

    def align(
        self,
        reference: Path,
        query: Path,
        output_file: Path,
        force_format: Optional[str] = None,
        build_index_first: Optional[bool] = None,
    ) -> Path:
        """
        Run minimap2 alignment with optional BAM processing.

        Args:
            reference: Reference genome (FASTA or .mmi index)
            query: Query reads file
            output_file: Output file path
            force_format: Force output format ("sam" or "bam")
            build_index_first: Whether to build index first (overrides config)

        Returns:
            Path to output file
        """
        reference = Path(reference)
        query = Path(query)
        output_file = Path(output_file)

        # Check input files
        if not query.exists():
            raise FileNotFoundError(f"Query file not found: {query}")

        # Determine if we should build index
        should_build_index = (
            build_index_first if build_index_first is not None else self.config.build_index
        )

        # Determine reference (index or FASTA)
        if should_build_index and reference.suffix in [".fa", ".fasta", ".fna"]:
            ref_index = self.build_index(reference)
        else:
            ref_index = reference
            if not ref_index.exists():
                raise FileNotFoundError(f"Reference file not found: {ref_index}")

        # Create output directory
        output_file.parent.mkdir(parents=True, exist_ok=True)

        # Determine output format
        output_format = force_format or self.config.output_format

        # Create temporary directory for intermediate files
        temp_dir = output_file.parent / "tmp_sort"
        temp_dir.mkdir(parents=True, exist_ok=True)

        try:
            if output_format.lower() == "bam":
                if not self.has_samtools:
                    raise PipelineError("samtools is required to produce BAM output")
                if self.config.sort_bam:
                    # Use temporary file approach to avoid pipe data loss
                    self._run_alignment_safe(ref_index, query, output_file, temp_dir)
                else:
                    self.logger.warning(
                        "BAM output requested without sorting; producing unsorted BAM"
                    )
                    self._run_alignment_unsorted(ref_index, query, output_file, temp_dir)
            else:
                self._run_alignment_simple(ref_index, query, output_file)

            # Create BAM index if requested and the BAM is sorted
            if output_format.lower() == "bam" and self.config.index_bam and self.has_samtools:
                if self.config.sort_bam:
                    self._index_bam(output_file)
                else:
                    self.logger.warning(
                        "Skipping BAM index creation because the BAM output is unsorted"
                    )

            self.logger.info("Alignment completed successfully")
            self.logger.info(f"Output file: {output_file}")

            if output_format.lower() == "bam" and self.config.index_bam:
                self.logger.info(f"Output index: {output_file}.bai")

            # Output statistics
            if output_format.lower() == "bam" and self.has_samtools:
                self._print_stats(output_file)

            return output_file

        finally:
            # Clean up temporary directory with error logging
            if temp_dir.exists():
                try:
                    shutil.rmtree(temp_dir)
                except OSError as e:
                    self.logger.warning(f"Failed to clean up temp directory {temp_dir}: {e}")

    def _run_alignment_safe(
        self, reference: Path, reads: Path, output_bam: Path, temp_dir: Path
    ) -> None:
        """
        Run alignment using temporary files to avoid pipe data loss.

        Args:
            reference: Reference sequence (index or FASTA)
            reads: Query reads file
            output_bam: Output BAM file
            temp_dir: Temporary directory for intermediate files
        """
        # Temporary SAM file
        temp_sam = temp_dir / "temp_alignment.sam"

        # Step 1: Run minimap2 to generate SAM file
        minimap2_cmd = ["minimap2", "-t", str(self.config.threads), "-ax", self.preset]

        # Add secondary alignment parameter
        if self.allow_secondary:
            minimap2_cmd.extend(["--secondary", "yes"])  # Use space instead of =
        else:
            minimap2_cmd.extend(["--secondary", "no"])

        # Add additional arguments if provided (use shlex for proper handling of quoted args)
        if self.config.additional_args:
            minimap2_cmd.extend(shlex.split(self.config.additional_args))

        minimap2_cmd.extend([str(reference), str(reads), "-o", str(temp_sam)])

        self.logger.info("Running minimap2 alignment")
        self.logger.debug(f"Command: {' '.join(minimap2_cmd)}")

        try:
            # Run minimap2
            result = subprocess.run(minimap2_cmd, stderr=subprocess.PIPE, text=True, check=True)

            # Log minimap2 stderr (contains progress info)
            if result.stderr:
                for line in result.stderr.strip().split("\n"):
                    if line:
                        self.logger.debug(f"minimap2: {line}")

            # Verify SAM file was created
            if not temp_sam.exists():
                raise PipelineError("minimap2 did not create output file")

            # Check SAM file size
            sam_size = temp_sam.stat().st_size
            if sam_size == 0:
                raise PipelineError("minimap2 created empty output file")

            self.logger.debug(f"SAM file created: {sam_size:,} bytes")

        except subprocess.CalledProcessError as e:
            self.logger.error(f"minimap2 failed with exit code {e.returncode}")
            self.logger.error(f"Error: {e.stderr}")
            raise ExternalToolError(
                "minimap2 alignment failed",
                command=e.cmd,
                returncode=e.returncode,
                stderr=e.stderr,
            )

        # Step 2: Sort SAM to BAM
        sort_threads = max(1, self.config.threads // 2)
        samtools_cmd = [
            "samtools",
            "sort",
            "-@",
            str(sort_threads),
            "-T",
            str(temp_dir / "sort_tmp"),
            "-o",
            str(output_bam),
            str(temp_sam),
        ]

        self.logger.info("Sorting alignment to BAM")
        self.logger.debug(f"Command: {' '.join(samtools_cmd)}")

        try:
            result = subprocess.run(samtools_cmd, stderr=subprocess.PIPE, text=True, check=True)

            if result.stderr:
                for line in result.stderr.strip().split("\n"):
                    if line:
                        self.logger.debug(f"samtools: {line}")

            # Verify BAM file was created
            if not output_bam.exists():
                raise PipelineError("samtools sort did not create output file")

            bam_size = output_bam.stat().st_size
            if bam_size == 0:
                raise PipelineError("samtools sort created empty output file")

            self.logger.debug(f"BAM file created: {bam_size:,} bytes")

        except subprocess.CalledProcessError as e:
            self.logger.error(f"samtools sort failed with exit code {e.returncode}")
            self.logger.error(f"Error: {e.stderr}")
            raise ExternalToolError(
                "BAM sorting failed",
                command=e.cmd,
                returncode=e.returncode,
                stderr=e.stderr,
            )

        finally:
            # Clean up temporary SAM file
            if temp_sam.exists():
                temp_sam.unlink()

    def _run_alignment_simple(self, reference: Path, reads: Path, output_file: Path) -> None:
        """
        Run simple alignment to SAM/BAM without sorting.

        Args:
            reference: Reference sequence
            reads: Query reads file
            output_file: Output file
        """
        command = ["minimap2", "-t", str(self.config.threads), "-ax", self.preset]

        if self.allow_secondary:
            command.extend(["--secondary", "yes"])
        else:
            command.extend(["--secondary", "no"])

        if self.config.additional_args:
            command.extend(shlex.split(self.config.additional_args))

        command.extend([str(reference), str(reads), "-o", str(output_file)])

        self.logger.info("Running minimap2 alignment")
        self.logger.debug(f"Command: {' '.join(command)}")

        try:
            subprocess.run(command, stderr=subprocess.PIPE, text=True, check=True)
        except subprocess.CalledProcessError as e:
            self.logger.error(f"minimap2 failed: {e.stderr}")
            raise PipelineError(f"Alignment failed: {e}")

    def _run_alignment_unsorted(
        self, reference: Path, reads: Path, output_bam: Path, temp_dir: Path
    ) -> None:
        """Run alignment and convert to BAM without sorting."""
        temp_sam = temp_dir / "temp_alignment.sam"

        minimap2_cmd = ["minimap2", "-t", str(self.config.threads), "-ax", self.preset]
        if self.allow_secondary:
            minimap2_cmd.extend(["--secondary", "yes"])
        else:
            minimap2_cmd.extend(["--secondary", "no"])

        if self.config.additional_args:
            minimap2_cmd.extend(shlex.split(self.config.additional_args))

        minimap2_cmd.extend([str(reference), str(reads), "-o", str(temp_sam)])

        self.logger.info("Running minimap2 alignment")
        self.logger.debug(f"Command: {' '.join(minimap2_cmd)}")

        try:
            subprocess.run(minimap2_cmd, stderr=subprocess.PIPE, text=True, check=True)
        except subprocess.CalledProcessError as e:
            self.logger.error(f"minimap2 failed: {e.stderr}")
            raise PipelineError(f"Alignment failed: {e}")

        samtools_cmd = [
            "samtools",
            "view",
            "-b",
            "-o",
            str(output_bam),
            str(temp_sam),
        ]

        self.logger.info("Converting SAM to BAM (unsorted)")
        self.logger.debug(f"Command: {' '.join(samtools_cmd)}")

        try:
            subprocess.run(samtools_cmd, stderr=subprocess.PIPE, text=True, check=True)
            if not output_bam.exists() or output_bam.stat().st_size == 0:
                raise PipelineError("samtools view created empty BAM output")
        except subprocess.CalledProcessError as e:
            self.logger.error(f"samtools view failed: {e.stderr}")
            raise PipelineError(f"BAM conversion failed: {e}")
        finally:
            if temp_sam.exists():
                temp_sam.unlink()

    def _index_bam(self, bam_file: Path) -> None:
        """
        Create BAM index.

        Args:
            bam_file: BAM file path
        """
        command = ["samtools", "index", str(bam_file)]

        self.logger.info("Creating BAM index")
        self.logger.debug(f"Command: {' '.join(command)}")

        try:
            subprocess.run(command, stderr=subprocess.PIPE, text=True, check=True)
            self.logger.info("BAM index created successfully")

        except subprocess.CalledProcessError as e:
            self.logger.error(f"samtools index failed with exit code {e.returncode}")
            self.logger.error(f"Error message: {e.stderr}")
            raise PipelineError(f"Failed to create BAM index: {e}")

    def _print_stats(self, bam_file: Path) -> None:
        """
        Output BAM file statistics.

        Args:
            bam_file: BAM file path
        """
        try:
            # Get alignment statistics
            result = subprocess.run(
                ["samtools", "flagstat", str(bam_file)],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True,
            )

            self.logger.info("Alignment statistics:")
            for line in result.stdout.strip().split("\n"):
                if line:
                    self.logger.info(f"  {line}")

            # Get file size
            bam_size = bam_file.stat().st_size / (1024 * 1024)  # MB
            self.logger.info(f"Output BAM size: {bam_size:.2f} MB")

        except Exception as e:
            self.logger.warning(f"Could not get BAM statistics: {e}")

    # Backward compatibility method
    def map_reads(self, reference: Path, reads: Path, output_sam: Path) -> None:
        """
        Legacy method for pipeline compatibility.
        Maps reads and outputs SAM format.
        """
        self.align(reference, reads, output_sam, force_format="sam")
