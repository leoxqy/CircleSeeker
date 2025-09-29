#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CD-HIT-EST integrated wrapper:
- Run cd-hit-est with default parameters: c=0.99, aS=0.95, aL=0.95
- Parse .clstr and export intuitive CSV (id, cluster, is_representative, cluster_size)

"""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List, Dict, Any, Tuple
import csv
import re

# ExternalTool from your project (unchanged)
from circleseeker.external.base import ExternalTool


@dataclass
class CDHitConfig:
    """Configuration for CD-HIT-EST."""
    similarity_threshold: float = 0.99     # -c
    threads: int = 8                       # -T
    memory_limit: int = 8000               # -M
    word_length: int = 10                  # -n (8/9/10/11)
    length_diff_cutoff: float = 0.99       # -s
    alignment_coverage_short: float = 0.95 # -aS
    alignment_coverage_long: float = 0.95  # -aL
    include_gaps: int = 1                  # -G
    desc_len: int = 0                      # -d


class CDHitEst(ExternalTool):
    """CD-HIT-EST nucleotide sequence clustering tool."""

    tool_name = "cd-hit-est"

    def __init__(
        self,
        config: Optional[CDHitConfig] = None,
        logger=None,
        *,
        threads: int = 24,
        c: float = 0.99,
        n: int = 10,
        s: float = 0.99,
        aS: float = 0.95,
        aL: float = 0.95,
        G: int = 1,
        d: int = 0,
        M: int = 0,
        extra: Optional[List[str]] = None,
    ):
        """
        Initialize CD-HIT-EST with clustering parameters.

        Args:
            config: CDHitConfig object with parameters
            logger: Logger instance
            threads: Number of threads to use (fallback) -> -T
            c: Sequence identity threshold (0.0-1.0) -> -c
            n: Word length, choose from 8, 9, 10, 11 -> -n
            s: Length difference cutoff (0.0-1.0) -> -s
            aS: Alignment coverage for shorter sequence (0.0-1.0) -> -aS
            aL: Alignment coverage for longer sequence (0.0-1.0) -> -aL
            G: Include/exclude gaps (1/0) -> -G
            d: Description length (0 for unlimited) -> -d
            M: Memory limit in MB (0 for unlimited) -> -M
            extra: Additional command line arguments
        """
        if config is not None:
            super().__init__(threads=config.threads, logger=logger)
            self.c = float(config.similarity_threshold)
            self.n = int(config.word_length)
            self.s = float(config.length_diff_cutoff)
            self.aS = float(config.alignment_coverage_short)
            self.aL = float(config.alignment_coverage_long)
            self.G = int(config.include_gaps)
            self.d = int(config.desc_len)
            self.M = int(config.memory_limit)
            self.extra = []
        else:
            super().__init__(threads=threads, logger=logger)
            self.c = float(c)
            self.n = int(n)
            self.s = float(s)
            self.aS = float(aS)
            self.aL = float(aL)
            self.G = int(G)
            self.d = int(d)
            self.M = int(M)
            self.extra = extra or []

    def _build_command(self, input_fasta: Path, output_prefix: Path) -> List[str]:
        """Build cd-hit-est command line."""
        cmd = [
            self.tool_name,
            "-i", str(input_fasta),
            "-o", str(output_prefix),
            "-c", str(self.c),
            "-n", str(self.n),
            "-s", str(self.s),
            "-aS", str(self.aS),
            "-aL", str(self.aL),
            "-G", str(self.G),
            "-d", str(self.d),
            "-M", str(self.M),
            "-T", str(self.threads),
        ]

        if self.extra:
            cmd.extend(self.extra)

        return cmd

    def cluster_sequences(self, input_fasta: Path, output_prefix: Path) -> Path:
        """
        Run CD-HIT-EST clustering.

        Args:
            input_fasta: Input FASTA file to cluster
            output_prefix: Output prefix for clustered sequences

        Produces:
            - {output_prefix}: Representative sequences FASTA
            - {output_prefix}.clstr: Cluster information file

        Returns:
            clstr_path: Path to the generated .clstr file
        """
        command = self._build_command(input_fasta, output_prefix)

        self.logger.info(f"Running CD-HIT-EST on: {input_fasta}")
        self.logger.debug(f"Command: {' '.join(command)}")

        # Ensure output directory exists
        output_prefix.parent.mkdir(parents=True, exist_ok=True)

        try:
            # Run CD-HIT-EST
            stdout, stderr = self.run(command)

            # Log output for debugging
            if stdout:
                self.logger.debug("CD-HIT-EST stdout:\n" + stdout)
            if stderr:
                self.logger.debug("CD-HIT-EST stderr:\n" + stderr)

            self.logger.info("CD-HIT-EST completed successfully")
            self.logger.info(f"Representative sequences: {output_prefix}")

            # Check for cluster file
            clstr_path = Path(str(output_prefix) + ".clstr")
            if clstr_path.exists():
                self.logger.info(f"Cluster file: {clstr_path}")

                # Log file sizes for debugging
                try:
                    rep_size = output_prefix.stat().st_size
                    clstr_size = clstr_path.stat().st_size
                    self.logger.debug(f"Representative file size: {rep_size} bytes")
                    self.logger.debug(f"Cluster file size: {clstr_size} bytes")
                except Exception as e:
                    self.logger.debug(f"Could not get file sizes: {e}")
            else:
                self.logger.warning(f"Cluster file not found: {clstr_path}")
            return clstr_path

        except Exception as e:
            self.logger.error(f"CD-HIT-EST failed: {e}")
            raise

    @staticmethod
    def _parse_clstr(clstr_file: Path) -> Tuple[Dict[str, List[str]], Dict[str, str]]:
        """
        Parse CD-HIT .clstr file
        Returns:
            members: cluster_id -> [member_ids]
            representatives: cluster_id -> representative_id
        """
        p = Path(clstr_file)
        if not p.exists():
            raise FileNotFoundError(f"Cluster file not found: {clstr_file}")

        members: Dict[str, List[str]] = {}
        representatives: Dict[str, str] = {}

        current_cluster: Optional[str] = None
        current_members: List[str] = []
        current_rep: Optional[str] = None

        with p.open("r") as f:
            for raw in f:
                line = raw.strip()
                if not line:
                    continue

                # Cluster header: ">Cluster 0"
                if line.startswith(">Cluster"):
                    # Finish previous cluster
                    if current_cluster is not None and current_members:
                        members[current_cluster] = current_members.copy()
                        if current_rep:
                            representatives[current_cluster] = current_rep

                    # Parse cluster id
                    try:
                        current_cluster = str(int(line.split()[1]))
                    except Exception:
                        # Fallback: use suffix as ID
                        current_cluster = line.replace(">Cluster", "").strip()
                    current_members = []
                    current_rep = None
                    continue

                # Cluster member line: "... >seqid... *"
                # Common format: "0   12634nt, >U23805... *"
                m = re.search(r'>\s*([^\.>\s][^\.]*?)\.\.\.', line) or re.search(r'>\s*([^\s,]+)', line)
                if m:
                    sid = m.group(1)
                    current_members.append(sid)
                    if line.endswith("*"):
                        current_rep = sid

        # Finish last cluster
        if current_cluster is not None and current_members:
            members[current_cluster] = current_members.copy()
            if current_rep:
                representatives[current_cluster] = current_rep

        return members, representatives

    @staticmethod
    def _write_id2cluster_csv(
        members: Dict[str, List[str]],
        representatives: Dict[str, str],
        output_csv: Path,
    ) -> Path:
        """
        Write CSV with four columns:
        - id                : sequence ID
        - cluster           : cluster ID (string format)
        - is_representative : whether representative sequence (true/false)
        - cluster_size      : number of members in this cluster
        """
        output_csv = Path(output_csv)
        output_csv.parent.mkdir(parents=True, exist_ok=True)

        with output_csv.open("w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["id", "cluster", "is_representative", "cluster_size"])
            for cid, ids in members.items():
                rep = representatives.get(cid)
                size = len(ids)
                for sid in ids:
                    writer.writerow([sid, cid, str(sid == rep).lower(), size])

        return output_csv

    def export_id2cluster(self, cluster_file: Path, output_csv: Optional[Path] = None) -> Path:
        """
        Convert .clstr to intuitive CSV (id, cluster, is_representative, cluster_size)

        Args:
            cluster_file: .clstr file path
            output_csv: output CSV path; defaults to same prefix as .clstr with .id2cluster.csv extension
        Returns:
            output_csv path
        """
        cluster_file = Path(cluster_file)
        if output_csv is None:
            # {output_prefix}.clstr -> {output_prefix}.id2cluster.csv
            output_csv = cluster_file.with_suffix("")  # Remove .clstr
            output_csv = output_csv.with_suffix(".id2cluster.csv")

        self.logger.info(f"Parsing clusters from: {cluster_file}")
        members, representatives = self._parse_clstr(cluster_file)

        csv_path = self._write_id2cluster_csv(members, representatives, output_csv)
        self.logger.info(f"ID→Cluster CSV written: {csv_path}")
        return csv_path


    def get_cluster_stats(self, cluster_file: Path) -> Dict[str, Any]:
        """
        Parse cluster file and return statistics.

        Args:
            cluster_file: Path to .clstr file

        Returns:
            Dictionary with cluster statistics
        """
        cluster_file = Path(cluster_file)
        if not cluster_file.exists():
            return {}

        try:
            members, _ = self._parse_clstr(cluster_file)

            total_clusters = len(members)
            total_sequences = sum(len(v) for v in members.values())

            stats: Dict[str, Any] = {
                'total_clusters': total_clusters,
                'total_sequences': total_sequences,
                'singleton_clusters': sum(1 for v in members.values() if len(v) == 1),
                'multi_member_clusters': sum(1 for v in members.values() if len(v) > 1),
            }

            if total_clusters > 0:
                stats['avg_cluster_size'] = (total_sequences / total_clusters)
            if total_sequences > 0:
                stats['reduction_rate'] = (1 - total_clusters / total_sequences) * 100

            self.logger.info(f"Cluster statistics: {stats}")
            return stats

        except Exception as e:
            self.logger.error(f"Error parsing cluster file: {e}")
            return {}


class CDHitRunner(CDHitEst):
    """Backward compatibility wrapper for existing code."""

    def __init__(
        self,
        threads: int = 24,
        c: float = 0.99,
        n: int = 10,
        s: float = 0.99,
        aS: float = 0.95,
        aL: float = 0.95,
        G: int = 1,
        d: int = 0,
        M: int = 0,
        extra: Optional[List[str]] = None,
    ):
        """Initialize with legacy parameter names."""
        super().__init__(
            threads=threads, c=c, n=n, s=s, aS=aS, aL=aL, G=G, d=d, M=M, extra=extra
        )

    def run(self, input_fasta: Path, output_prefix: Path):
        """Run CD-HIT-EST with legacy interface."""
        return self.cluster_sequences(input_fasta, output_prefix)


def _parse_args():
    """Parse CLI arguments for direct script execution"""
    import argparse
    parser = argparse.ArgumentParser(
        description="CD-HIT-EST - Sequence clustering tool"
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input FASTA file"
    )
    parser.add_argument(
        "-o", "--output-prefix",
        required=True,
        help="Output prefix (no extension)"
    )
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
        help="Log level (default: INFO)"
    )
    return parser.parse_args()


def main():
    """Main function for CLI execution"""
    from pathlib import Path
    import logging as _logging
    from circleseeker.utils.logging import get_logger
    
    args = _parse_args()
    
    input_path = Path(args.input)
    output_prefix = Path(args.output_prefix)
    
    # Set root log level
    _logging.getLogger().setLevel(getattr(_logging, args.log_level))
    
    # Create logger
    logger = get_logger("CDHitEst")
    
    # Create and run CD-HIT-EST
    cdhit = CDHitEst(logger=logger)
    
    try:
        # Run clustering
        cluster_file = cdhit.cluster_sequences(input_path, output_prefix)
        
        # Export ID to cluster CSV
        csv_file = cdhit.export_id2cluster(cluster_file)
        
        # Get and display stats
        stats = cdhit.get_cluster_stats(cluster_file)
        logger.info(f"Clustering done! Output: {output_prefix} and {csv_file}")
        
    except Exception as e:
        logger.error(f"CD-HIT-EST failed: {e}")
        raise


if __name__ == "__main__":
    main()