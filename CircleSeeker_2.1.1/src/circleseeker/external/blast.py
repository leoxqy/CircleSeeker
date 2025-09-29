"""BLAST+ wrapper."""

from pathlib import Path
from typing import Optional
from circleseeker.external.base import ExternalTool


class MakeBlastDB(ExternalTool):
    """makeblastdb wrapper for database creation."""
    
    tool_name = "makeblastdb"
    
    def build_database(
        self,
        input_file: Path,
        output_db: Path,
        dbtype: str = "nucl",
        title: Optional[str] = None,
        parse_seqids: bool = False,
        taxid: Optional[int] = None,
        max_file_size: str = "1GB"
    ) -> None:
        """Build BLAST database."""
        cmd = [
            self.tool_name,
            "-in", str(input_file),
            "-dbtype", dbtype,
            "-out", str(output_db),
            "-max_file_sz", max_file_size
        ]
        
        if title:
            cmd.extend(["-title", title])
        if parse_seqids:
            cmd.append("-parse_seqids")
        if taxid:
            cmd.extend(["-taxid", str(taxid)])
        
        stdout, stderr = self.run(cmd, capture_output=True)
        if stdout:
            self.logger.debug(f"makeblastdb output: {stdout[:500]}")
        self.logger.info(f"BLAST database created: {output_db}")
    
    def verify_database(self, db_prefix: Path) -> bool:
        """Verify BLAST database exists and is valid."""
        # Check for required files
        required_exts = ['.nhr', '.nin', '.nsq']  # For nucleotide DB
        return all((db_prefix.with_suffix(ext)).exists() for ext in required_exts)


class BlastN(ExternalTool):
    """blastn wrapper for nucleotide BLAST."""
    
    tool_name = "blastn"
    
    def run_blast(
        self,
        query_file: Path,
        database: Path,
        output_file: Path,
        evalue: str = "1e-50",
        word_size: int = 100,
        perc_identity: float = 99.0,
        outfmt: str = "6 std sstrand",
        max_target_seqs: int = 1000
    ) -> None:
        """Run BLAST search."""
        cmd = [
            self.tool_name,
            "-query", str(query_file),
            "-db", str(database),
            "-out", str(output_file),
            "-evalue", str(evalue),
            "-word_size", str(word_size),
            "-perc_identity", str(perc_identity),
            "-outfmt", outfmt,
            "-max_target_seqs", str(max_target_seqs),
            "-num_threads", str(self.threads)
        ]
        
        output_file.parent.mkdir(parents=True, exist_ok=True)
        stdout, stderr = self.run(cmd, capture_output=True)
        if stdout:
            self.logger.debug(f"blastn output: {stdout[:500]}")
        self.logger.info(f"BLAST results saved to: {output_file}")


class BlastRunner:
    """High-level BLAST runner combining makeblastdb and blastn."""
    
    def __init__(
        self,
        num_threads: int = 8,
        word_size: int = 100,
        evalue: str = "1e-50",
        perc_identity: float = 99.0,
        outfmt: str = "6 std sstrand"
    ):
        self.makeblastdb = MakeBlastDB(threads=num_threads)
        self.blastn = BlastN(threads=num_threads)
        self.word_size = word_size
        self.evalue = evalue
        self.perc_identity = perc_identity
        self.outfmt = outfmt
    
    def run(
        self,
        database: Path,
        query_file: Path,
        output_file: Path,
        use_time_cmd: bool = False
    ) -> None:
        """Run BLAST search using existing database."""
        self.blastn.run_blast(
            query_file=query_file,
            database=database,
            output_file=output_file,
            evalue=self.evalue,
            word_size=self.word_size,
            perc_identity=self.perc_identity,
            outfmt=self.outfmt
        )
