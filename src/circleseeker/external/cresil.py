"""Cresil external tool wrapper.

Provides typed methods for cresil subcommands used by CircleSeeker.
Cresil is a circular DNA detection tool that works directly with FASTA files.
"""

from __future__ import annotations

from pathlib import Path

from circleseeker.external.base import ExternalTool
from circleseeker.external.samtools import Samtools
from circleseeker.exceptions import ExternalToolError


class Cresil(ExternalTool):
    """Wrapper for the `cresil` CLI."""

    tool_name = "cresil"

    def trim(
        self,
        fasta_query: Path,
        reference_mmi: Path,
        output_dir: Path,
        *,
        threads: int = 1,
    ) -> Path:
        """Run Cresil trim step.

        Args:
            fasta_query: Input FASTA file (query sequences)
            reference_mmi: Reference genome minimap2 index (.mmi file)
            output_dir: Output directory for trim results
            threads: Number of threads to use

        Returns:
            Path to trim.txt output file
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        output_dir = output_dir.resolve()

        fasta_query = Path(fasta_query).resolve()
        reference_mmi = Path(reference_mmi).resolve()

        cmd = [
            self.tool_name,
            "trim",
            "-t",
            str(threads),
            "-fq",
            str(fasta_query),
            "-r",
            str(reference_mmi),
            "-o",
            str(output_dir),
        ]

        stdout, stderr = self.run(cmd, capture_output=True)

        if stderr and ("ERROR" in stderr or "error" in stderr):
            self.logger.warning(f"Cresil trim stderr: {stderr[:1000]}")
        else:
            self.logger.debug(f"Cresil trim progress: {stderr[:500] if stderr else 'no output'}")

        trim_output = output_dir / "trim.txt"
        if trim_output.exists():
            self.logger.info(f"Cresil trim completed: {trim_output}")
        else:
            raise FileNotFoundError(f"Expected trim output not found: {trim_output}")

        return trim_output

    def identify(
        self,
        fasta_query: Path,
        reference_fasta: Path,
        reference_fai: Path,
        trim_file: Path,
        output_dir: Path,
        *,
        threads: int = 1,
        split_reads: bool = True,
    ) -> Path:
        """Run Cresil identify step.

        Args:
            fasta_query: Input FASTA file (query sequences)
            reference_fasta: Reference genome FASTA file
            reference_fai: Reference genome FASTA index (.fai file)
            trim_file: trim.txt file from cresil trim step
            output_dir: Output directory for identify results
            threads: Number of threads to use
            split_reads: Enable split reads mode (-s flag)

        Returns:
            Path to eccDNA_final.txt output file
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        output_dir = output_dir.resolve()

        fasta_query = Path(fasta_query).resolve()
        reference_fasta = Path(reference_fasta).resolve()
        reference_fai = Path(reference_fai).resolve()
        trim_file = Path(trim_file).resolve()

        cmd = [
            self.tool_name,
            "identify",
            "-t",
            str(threads),
            "-fa",
            str(reference_fasta),
            "-fai",
            str(reference_fai),
            "-fq",
            str(fasta_query),
            "-trim",
            str(trim_file),
        ]

        if split_reads:
            cmd.append("-s")

        # Cresil identify writes output to current directory by default
        # We'll run it from output_dir to control output location
        stdout, stderr = self.run(cmd, capture_output=True, cwd=output_dir)

        if stderr and ("ERROR" in stderr or "error" in stderr):
            self.logger.warning(f"Cresil identify stderr: {stderr[:1000]}")
        else:
            self.logger.debug(
                f"Cresil identify progress: {stderr[:500] if stderr else 'no output'}"
            )

        # Expected output file
        eccDNA_output = output_dir / "eccDNA_final.txt"
        if eccDNA_output.exists():
            self.logger.info(f"Cresil identify completed: {eccDNA_output}")
        else:
            raise FileNotFoundError(f"Expected eccDNA output not found: {eccDNA_output}")

        return eccDNA_output

    def run_full_pipeline(
        self,
        fasta_query: Path,
        reference_fasta: Path,
        reference_mmi: Path,
        output_dir: Path,
        *,
        threads: int = 1,
        split_reads: bool = True,
    ) -> Path:
        """Run complete Cresil pipeline (trim + identify).

        Args:
            fasta_query: Input FASTA file
            reference_fasta: Reference genome FASTA
            reference_mmi: Reference genome minimap2 index
            output_dir: Output directory
            threads: Number of threads
            split_reads: Enable split reads mode

        Returns:
            Path to eccDNA_final.txt
        """
        self.logger.info("Running Cresil full pipeline")

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        output_dir = output_dir.resolve()

        # Track symlinks created for cleanup
        symlinks_to_cleanup: list[Path] = []

        # Cresil substep 1: trim
        self.logger.info("  → Cresil substep 1/2: trim")
        trim_file = self.trim(
            fasta_query=fasta_query,
            reference_mmi=reference_mmi,
            output_dir=output_dir,
            threads=threads,
        )

        # Cresil substep 2: identify
        self.logger.info("  → Cresil substep 2/2: identify")

        # Check for .fai file
        reference_fasta = Path(reference_fasta).resolve()
        reference_fai = reference_fasta.with_suffix(reference_fasta.suffix + ".fai")
        reference_fasta_for_identify = reference_fasta
        if not reference_fai.exists():
            self.logger.info(f"Reference FASTA index not found: {reference_fai}")
            self.logger.info("Attempting to create FASTA index with samtools faidx")
            try:
                samtools = Samtools(logger=self.logger.getChild("samtools"), threads=max(1, threads))
                samtools.faidx(reference_fasta_for_identify)
            except ExternalToolError as exc:
                # Fallback: if reference is not writable, create a symlink in output_dir so samtools
                # can write the .fai alongside it without touching the original reference directory.
                self.logger.warning(
                    "samtools faidx failed for reference FASTA (%s). "
                    "Trying to create index in output_dir via symlink.",
                    exc,
                )

                link_candidate = output_dir / reference_fasta.name
                try:
                    if link_candidate.exists() and link_candidate.resolve() != reference_fasta.resolve():
                        link_candidate = output_dir / f"{reference_fasta.stem}.circleseeker{reference_fasta.suffix}"
                except Exception:
                    link_candidate = output_dir / f"{reference_fasta.stem}.circleseeker{reference_fasta.suffix}"

                if not link_candidate.exists():
                    try:
                        link_candidate.symlink_to(reference_fasta)
                        symlinks_to_cleanup.append(link_candidate)
                    except OSError as link_exc:
                        raise FileNotFoundError(
                            "Reference FASTA index missing and automatic creation failed; "
                            "could not create symlink for fallback faidx in output_dir"
                        ) from link_exc

                reference_fasta_for_identify = link_candidate
                reference_fai = reference_fasta_for_identify.with_suffix(
                    reference_fasta_for_identify.suffix + ".fai"
                )
                try:
                    samtools.faidx(reference_fasta_for_identify)
                except ExternalToolError as exc2:
                    raise FileNotFoundError(
                        f"Reference FASTA index missing and automatic creation failed: {reference_fai}"
                    ) from exc2

            if not reference_fai.exists():
                raise FileNotFoundError(
                    f"Reference FASTA index still missing after samtools faidx: {reference_fai}"
                )

        eccDNA_output = self.identify(
            fasta_query=fasta_query,
            reference_fasta=reference_fasta_for_identify,
            reference_fai=reference_fai,
            trim_file=trim_file,
            output_dir=output_dir,
            threads=threads,
            split_reads=split_reads,
        )

        self.logger.info(f"Cresil pipeline completed. Output: {eccDNA_output}")

        # Cleanup any symlinks created during the process
        for symlink in symlinks_to_cleanup:
            try:
                if symlink.is_symlink():
                    symlink.unlink()
                    self.logger.debug(f"Cleaned up temporary symlink: {symlink}")
                    # Also cleanup any .fai file created alongside the symlink
                    fai_file = symlink.with_suffix(symlink.suffix + ".fai")
                    if fai_file.exists():
                        fai_file.unlink()
                        self.logger.debug(f"Cleaned up temporary index: {fai_file}")
            except OSError as e:
                self.logger.warning(f"Failed to cleanup symlink {symlink}: {e}")

        return eccDNA_output
