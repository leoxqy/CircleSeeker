"""
Wrapper modules for external tools in the CircleSeeker pipeline.
These integrate external bioinformatics tools into the new module architecture.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Optional

from circleseeker.modules.base import ExternalToolModule, ModuleResult


class TideHunterModule(ExternalToolModule):
    """Wrapper for TideHunter tandem repeat detection tool."""

    def __init__(self, **kwargs):
        super().__init__(tool_name="TideHunter", **kwargs)
        self.period_range = "50,10000"
        self.consensus_threshold = 0.8
        self.format_type = 2

    def check_tool_availability(self) -> bool:
        """Check if TideHunter is available."""
        return shutil.which("TideHunter") is not None

    def get_tool_version(self) -> str:
        """Get TideHunter version."""
        try:
            result = subprocess.run(
                ["TideHunter", "--version"],
                capture_output=True,
                text=True
            )
            return result.stdout.strip()
        except:
            return "unknown"

    def validate_inputs(self, **kwargs) -> bool:
        """Validate inputs for TideHunter."""
        if 'input_file' not in kwargs:
            raise ValueError("input_file is required")

        input_file = Path(kwargs['input_file'])
        if not input_file.exists():
            raise FileNotFoundError(f"Input file not found: {input_file}")

        return True

    def execute(self, **kwargs) -> ModuleResult:
        """Run TideHunter on input FASTA file."""
        result = ModuleResult(
            success=False,
            module_name=self.name
        )

        # Check tool availability
        if not self.check_tool_availability():
            result.error_message = "TideHunter not found in PATH"
            result.add_warning("Please install TideHunter or add it to PATH")
            return result

        input_file = Path(kwargs['input_file'])
        output_file = Path(kwargs.get('output_file', 'tidehunter_output.txt'))
        threads = kwargs.get('threads', 4)

        # Update parameters from kwargs if provided
        self.period_range = kwargs.get('period_range', self.period_range)
        self.consensus_threshold = kwargs.get('consensus_threshold', self.consensus_threshold)

        # Build command
        cmd = [
            "TideHunter",
            "-f", str(self.format_type),  # Output format
            "-p", self.period_range,  # Period range
            "-c", str(self.consensus_threshold),  # Consensus threshold
            "-t", str(threads),
            str(input_file),
            "-o", str(output_file)
        ]

        self.logger.info(f"Running TideHunter: {' '.join(cmd)}")

        try:
            # Run TideHunter
            process = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )

            result.success = True
            result.add_output('ecc_candidates', output_file)

            # Count lines in output (rough estimate of candidates)
            if output_file.exists():
                with open(output_file, 'r') as f:
                    num_candidates = sum(1 for line in f if line.strip() and not line.startswith('#'))
                result.add_metric('num_candidates', num_candidates)

            self.logger.info(f"TideHunter found {num_candidates} tandem repeat candidates")

        except subprocess.CalledProcessError as e:
            result.error_message = f"TideHunter failed: {e.stderr}"
            self.logger.error(result.error_message)

        return result


class BlastModule(ExternalToolModule):
    """Wrapper for BLAST sequence similarity search."""

    def __init__(self, **kwargs):
        super().__init__(tool_name="blastn", **kwargs)
        self.evalue = 1e-5
        self.max_target_seqs = 100
        self.outfmt = 6  # Tabular format

    def check_tool_availability(self) -> bool:
        """Check if BLAST is available."""
        return shutil.which("blastn") is not None and shutil.which("makeblastdb") is not None

    def get_tool_version(self) -> str:
        """Get BLAST version."""
        try:
            result = subprocess.run(
                ["blastn", "-version"],
                capture_output=True,
                text=True
            )
            return result.stdout.split('\n')[0]
        except:
            return "unknown"

    def validate_inputs(self, **kwargs) -> bool:
        """Validate inputs for BLAST."""
        required = ['query_file', 'db_file']
        for req in required:
            if req not in kwargs:
                raise ValueError(f"{req} is required")

            file_path = Path(kwargs[req])
            if not file_path.exists():
                raise FileNotFoundError(f"File not found: {file_path}")

        return True

    def execute(self, **kwargs) -> ModuleResult:
        """Run BLAST search."""
        result = ModuleResult(
            success=False,
            module_name=self.name
        )

        # Check tool availability
        if not self.check_tool_availability():
            result.error_message = "BLAST tools not found in PATH"
            result.add_warning("Please install BLAST+ or add it to PATH")
            return result

        query_file = Path(kwargs['query_file'])
        db_file = Path(kwargs['db_file'])
        output_file = Path(kwargs.get('output_file', 'blast_results.tsv'))
        threads = kwargs.get('threads', 4)

        # Get parameters from kwargs
        self.evalue = kwargs.get('evalue', self.evalue)
        self.max_target_seqs = kwargs.get('max_target_seqs', self.max_target_seqs)

        # Create BLAST database
        db_path = output_file.parent / "blast_db"
        db_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            # Make BLAST database
            self.logger.info("Creating BLAST database...")
            cmd_makedb = [
                "makeblastdb",
                "-in", str(db_file),
                "-dbtype", "nucl",
                "-out", str(db_path)
            ]

            subprocess.run(cmd_makedb, check=True, capture_output=True)

            # Run BLAST
            self.logger.info("Running BLAST search...")
            cmd_blast = [
                "blastn",
                "-query", str(query_file),
                "-db", str(db_path),
                "-out", str(output_file),
                "-outfmt", str(self.outfmt),
                "-evalue", str(self.evalue),
                "-max_target_seqs", str(self.max_target_seqs),
                "-num_threads", str(threads)
            ]

            process = subprocess.run(
                cmd_blast,
                capture_output=True,
                text=True,
                check=True
            )

            result.success = True
            result.add_output('blast_results', output_file)

            # Count hits
            if output_file.exists():
                with open(output_file, 'r') as f:
                    num_hits = sum(1 for line in f if line.strip())
                result.add_metric('num_hits', num_hits)

            self.logger.info(f"BLAST found {num_hits} hits")

        except subprocess.CalledProcessError as e:
            result.error_message = f"BLAST failed: {e.stderr}"
            self.logger.error(result.error_message)

        return result


class CDHitModule(ExternalToolModule):
    """Wrapper for CD-HIT sequence clustering tool."""

    def __init__(self, **kwargs):
        super().__init__(tool_name="cd-hit-est", **kwargs)
        self.similarity = 0.99  # 99% similarity threshold
        self.memory = 0  # No memory limit

    def check_tool_availability(self) -> bool:
        """Check if CD-HIT is available."""
        return shutil.which("cd-hit-est") is not None

    def get_tool_version(self) -> str:
        """Get CD-HIT version."""
        try:
            result = subprocess.run(
                ["cd-hit-est"],
                capture_output=True,
                text=True
            )
            # CD-HIT prints version in stderr
            for line in result.stderr.split('\n'):
                if 'CD-HIT' in line and 'version' in line.lower():
                    return line.strip()
            return "unknown"
        except:
            return "unknown"

    def validate_inputs(self, **kwargs) -> bool:
        """Validate inputs for CD-HIT."""
        if 'input_file' not in kwargs:
            raise ValueError("input_file is required")

        input_file = Path(kwargs['input_file'])
        if not input_file.exists():
            raise FileNotFoundError(f"Input file not found: {input_file}")

        return True

    def execute(self, **kwargs) -> ModuleResult:
        """Run CD-HIT clustering."""
        result = ModuleResult(
            success=False,
            module_name=self.name
        )

        # Check tool availability
        if not self.check_tool_availability():
            result.error_message = "cd-hit-est not found in PATH"
            result.add_warning("Please install CD-HIT or add it to PATH")
            return result

        input_file = Path(kwargs['input_file'])
        output_prefix = kwargs.get('output_prefix', input_file.stem)
        output_fasta = Path(f"{output_prefix}.fasta")
        output_clstr = Path(f"{output_prefix}.clstr")
        threads = kwargs.get('threads', 4)

        # Get parameters
        self.similarity = kwargs.get('similarity', self.similarity)

        # Build command
        cmd = [
            "cd-hit-est",
            "-i", str(input_file),
            "-o", str(output_fasta),
            "-c", str(self.similarity),
            "-T", str(threads),
            "-d", "0",  # Use full sequence names
            "-M", str(self.memory)
        ]

        self.logger.info(f"Running CD-HIT: {' '.join(cmd)}")

        try:
            # Run CD-HIT
            process = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )

            result.success = True
            result.add_output('clustered_fasta', output_fasta)
            result.add_output('cluster_file', output_clstr)

            # Count clusters
            if output_clstr.exists():
                num_clusters = 0
                with open(output_clstr, 'r') as f:
                    for line in f:
                        if line.startswith('>Cluster'):
                            num_clusters += 1
                result.add_metric('num_clusters', num_clusters)

            self.logger.info(f"CD-HIT created {num_clusters} clusters")

        except subprocess.CalledProcessError as e:
            result.error_message = f"CD-HIT failed: {e.stderr}"
            self.logger.error(result.error_message)

        return result


class Minimap2Module(ExternalToolModule):
    """Wrapper for Minimap2 aligner."""

    def __init__(self, **kwargs):
        super().__init__(tool_name="minimap2", **kwargs)

    def check_tool_availability(self) -> bool:
        """Check if minimap2 and samtools are available."""
        return (shutil.which("minimap2") is not None and
                shutil.which("samtools") is not None)

    def get_tool_version(self) -> str:
        """Get minimap2 version."""
        try:
            result = subprocess.run(
                ["minimap2", "--version"],
                capture_output=True,
                text=True
            )
            return result.stdout.strip()
        except:
            return "unknown"

    def validate_inputs(self, **kwargs) -> bool:
        """Validate inputs for minimap2."""
        required = ['reads_file', 'reference_file']
        for req in required:
            if req not in kwargs:
                raise ValueError(f"{req} is required")

            file_path = Path(kwargs[req])
            if not file_path.exists():
                raise FileNotFoundError(f"File not found: {file_path}")

        return True

    def execute(self, **kwargs) -> ModuleResult:
        """Run minimap2 alignment."""
        result = ModuleResult(
            success=False,
            module_name=self.name
        )

        # Check tool availability
        if not self.check_tool_availability():
            result.error_message = "minimap2 or samtools not found in PATH"
            result.add_warning("Please install minimap2 and samtools")
            return result

        reads_file = Path(kwargs['reads_file'])
        reference_file = Path(kwargs['reference_file'])
        output_bam = Path(kwargs.get('output_bam', 'aligned.bam'))
        threads = kwargs.get('threads', 4)

        try:
            # Run minimap2 and pipe to samtools
            self.logger.info("Running minimap2 alignment...")

            cmd = (
                f"minimap2 -ax map-hifi -t {threads} "
                f"{reference_file} {reads_file} | "
                f"samtools sort -@ {threads} -o {output_bam} -"
            )

            process = subprocess.run(
                cmd,
                shell=True,
                capture_output=True,
                text=True,
                check=True
            )

            # Index the BAM file
            self.logger.info("Indexing BAM file...")
            subprocess.run(
                ["samtools", "index", str(output_bam)],
                check=True
            )

            result.success = True
            result.add_output('bam_file', output_bam)
            result.add_output('bam_index', Path(f"{output_bam}.bai"))

            # Get alignment stats
            stats_output = subprocess.run(
                ["samtools", "flagstat", str(output_bam)],
                capture_output=True,
                text=True
            ).stdout

            # Parse basic stats
            for line in stats_output.split('\n'):
                if 'mapped' in line and 'mate' not in line:
                    parts = line.split()
                    if parts:
                        mapped_reads = int(parts[0])
                        result.add_metric('mapped_reads', mapped_reads)
                        break

            self.logger.info(f"Alignment completed: {output_bam}")

        except subprocess.CalledProcessError as e:
            result.error_message = f"Minimap2 failed: {e.stderr}"
            self.logger.error(result.error_message)

        return result