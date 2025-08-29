"""CD-HIT wrapper for sequence clustering."""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List
from circleseeker2.external.base import ExternalTool


@dataclass
class CDHitConfig:
    """Configuration for CD-HIT-EST."""
    similarity_threshold: float = 0.99
    threads: int = 8
    memory_limit: int = 8000
    word_length: int = 10
    length_diff_cutoff: float = 0.99
    alignment_coverage: float = 0.99


class CDHitEst(ExternalTool):
    """CD-HIT-EST nucleotide sequence clustering tool."""
    
    tool_name = "cd-hit-est"
    
    def __init__(self, config: Optional[CDHitConfig] = None, logger=None, 
                 threads: int = 24, c: float = 0.99, n: int = 10,
                 s: float = 0.99, aS: float = 0.99, G: int = 1,
                 d: int = 0, M: int = 0, extra: Optional[List[str]] = None):
        """
        Initialize CD-HIT-EST with clustering parameters.
        
        Args:
            config: CDHitConfig object with parameters
            logger: Logger instance
            threads: Number of threads to use (fallback)
            c: Sequence identity threshold (0.0-1.0)
            n: Word_length, choose from 8, 9, 10, 11  
            s: Length difference cutoff (0.0-1.0)
            aS: Alignment coverage for shorter sequence (0.0-1.0)
            G: Include/exclude gaps (1/0)
            d: Description length (0 for unlimited)
            M: Memory limit in MB (0 for unlimited)
            extra: Additional command line arguments
        """
        if config is not None:
            super().__init__(threads=config.threads, logger=logger)
            self.c = config.similarity_threshold
            self.n = config.word_length
            self.s = config.length_diff_cutoff
            self.aS = config.alignment_coverage
            self.G = 1
            self.d = 0
            self.M = config.memory_limit
            self.extra = []
        else:
            super().__init__(threads=threads, logger=logger)
            self.c = float(c)
            self.n = int(n)
            self.s = float(s)
            self.aS = float(aS)
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
                self.logger.debug("CD-HIT-EST stdout:\\n" + stdout)
            if stderr:
                self.logger.debug("CD-HIT-EST stderr:\\n" + stderr)
            
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
        
        # Return cluster file path even if not found (caller may check existence)
        return Path(str(output_prefix) + ".clstr")
    
    def get_cluster_stats(self, cluster_file: Path) -> dict:
        """
        Parse cluster file and return statistics.
        
        Args:
            cluster_file: Path to .clstr file
            
        Returns:
            Dictionary with cluster statistics
        """
        if not cluster_file.exists():
            return {}
        
        try:
            clusters = {}
            current_cluster = None
            
            with open(cluster_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('>Cluster'):
                        current_cluster = line
                        clusters[current_cluster] = []
                    elif current_cluster and line:
                        clusters[current_cluster].append(line)
            
            stats = {
                'total_clusters': len(clusters),
                'total_sequences': sum(len(seqs) for seqs in clusters.values()),
                'singleton_clusters': sum(1 for seqs in clusters.values() if len(seqs) == 1),
                'multi_member_clusters': sum(1 for seqs in clusters.values() if len(seqs) > 1),
            }
            
            if stats['total_clusters'] > 0:
                stats['avg_cluster_size'] = stats['total_sequences'] / stats['total_clusters']
                stats['reduction_rate'] = (1 - stats['total_clusters'] / stats['total_sequences']) * 100
            
            self.logger.info(f"Cluster statistics: {stats}")
            return stats
            
        except Exception as e:
            self.logger.error(f"Error parsing cluster file: {e}")
            return {}


class CDHitRunner(CDHitEst):
    """Backward compatibility wrapper for existing code."""
    
    def __init__(self, threads=24, c=0.99, n=10, s=0.99, aS=0.99, G=1, d=0, M=0, extra=None):
        """Initialize with legacy parameter names."""
        super().__init__(
            threads=threads, c=c, n=n, s=s, aS=aS, G=G, d=d, M=M, extra=extra
        )
    
    def run(self, input_fasta: Path, output_prefix: Path):
        """Run CD-HIT-EST with legacy interface."""
        return self.cluster_sequences(input_fasta, output_prefix)
