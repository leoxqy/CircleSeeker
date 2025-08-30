"""Test utilities for CircleSeeker tests."""

import subprocess
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple
import pandas as pd
from Bio import SeqIO


class TestComparator:
    """Utility class for comparing test outputs."""
    
    @staticmethod
    def compare_csv_row_counts(file1: Path, file2: Path, tolerance: float = 0.1) -> Tuple[bool, str]:
        """Compare row counts between two CSV files."""
        try:
            df1 = pd.read_csv(file1)
            df2 = pd.read_csv(file2)
            
            count1 = len(df1)
            count2 = len(df2)
            
            if count1 == 0 and count2 == 0:
                return True, "Both files are empty"
            
            diff_ratio = abs(count1 - count2) / max(count1, count2, 1)
            
            if diff_ratio <= tolerance:
                return True, f"Row counts match within tolerance: {count1} vs {count2}"
            else:
                return False, f"Row count mismatch: {count1} vs {count2} (diff: {diff_ratio:.2%})"
                
        except Exception as e:
            return False, f"Error comparing CSV files: {str(e)}"
    
    @staticmethod
    def compare_csv_columns(file1: Path, file2: Path) -> Tuple[bool, str]:
        """Compare column names between two CSV files."""
        try:
            df1 = pd.read_csv(file1)
            df2 = pd.read_csv(file2)
            
            cols1 = set(df1.columns)
            cols2 = set(df2.columns)
            
            if cols1 == cols2:
                return True, "Column names match"
            
            missing_in_2 = cols1 - cols2
            missing_in_1 = cols2 - cols1
            
            msg = "Column mismatch:"
            if missing_in_2:
                msg += f"\nMissing in file2: {missing_in_2}"
            if missing_in_1:
                msg += f"\nMissing in file1: {missing_in_1}"
            
            return False, msg
            
        except Exception as e:
            return False, f"Error comparing columns: {str(e)}"
    
    @staticmethod
    def compare_fasta_records(file1: Path, file2: Path, tolerance: float = 0.1) -> Tuple[bool, str]:
        """Compare FASTA files by record count."""
        try:
            records1 = list(SeqIO.parse(file1, "fasta"))
            records2 = list(SeqIO.parse(file2, "fasta"))
            
            count1 = len(records1)
            count2 = len(records2)
            
            if count1 == 0 and count2 == 0:
                return True, "Both FASTA files are empty"
            
            diff_ratio = abs(count1 - count2) / max(count1, count2, 1)
            
            if diff_ratio <= tolerance:
                return True, f"Record counts match within tolerance: {count1} vs {count2}"
            else:
                return False, f"Record count mismatch: {count1} vs {count2} (diff: {diff_ratio:.2%})"
                
        except Exception as e:
            return False, f"Error comparing FASTA files: {str(e)}"
    
    @staticmethod
    def compare_bed_files(file1: Path, file2: Path, tolerance: float = 0.1) -> Tuple[bool, str]:
        """Compare BED files."""
        try:
            # Read BED files
            df1 = pd.read_csv(file1, sep='\t', header=None)
            df2 = pd.read_csv(file2, sep='\t', header=None)
            
            # Compare row counts
            count1 = len(df1)
            count2 = len(df2)
            
            if count1 == 0 and count2 == 0:
                return True, "Both BED files are empty"
            
            diff_ratio = abs(count1 - count2) / max(count1, count2, 1)
            
            if diff_ratio <= tolerance:
                return True, f"BED entries match within tolerance: {count1} vs {count2}"
            else:
                return False, f"BED entry count mismatch: {count1} vs {count2}"
                
        except Exception as e:
            return False, f"Error comparing BED files: {str(e)}"


class ExternalToolChecker:
    """Check availability of external tools."""
    
    @staticmethod
    def check_tool(tool_name: str, version_cmd: Optional[str] = None) -> bool:
        """Check if a tool is available in PATH."""
        try:
            cmd = version_cmd or f"{tool_name} --version"
            result = subprocess.run(
                cmd.split(), 
                capture_output=True, 
                text=True, 
                timeout=5
            )
            return result.returncode == 0
        except (subprocess.TimeoutExpired, FileNotFoundError):
            return False
    
    @staticmethod
    def get_available_tools() -> Dict[str, bool]:
        """Check all required external tools."""
        tools = {
            "tidehunter": "tidehunter -h",
            "blastn": "blastn -version",
            "makeblastdb": "makeblastdb -version",
            "cd-hit-est": "cd-hit-est -h",
            "minimap2": "minimap2 --version",
            "samtools": "samtools --version",
            "mosdepth": "mosdepth --version"
        }
        
        return {
            tool: ExternalToolChecker.check_tool(tool, cmd)
            for tool, cmd in tools.items()
        }


def create_test_fasta(path: Path, records: List[Tuple[str, str]]) -> None:
    """Create a test FASTA file."""
    with open(path, 'w') as f:
        for seq_id, seq in records:
            f.write(f">{seq_id}\n{seq}\n")


def create_test_bed(path: Path, regions: List[Tuple[str, int, int]]) -> None:
    """Create a test BED file."""
    with open(path, 'w') as f:
        for chrom, start, end in regions:
            f.write(f"{chrom}\t{start}\t{end}\n")


def run_circleseeker_module(
    module_name: str,
    args: List[str],
    cwd: Optional[Path] = None
) -> subprocess.CompletedProcess:
    """Run a CircleSeeker module and return the result."""
    cmd = ["python", "-m", f"circleseeker.modules.{module_name}"] + args
    
    result = subprocess.run(
        cmd,
        cwd=cwd,
        capture_output=True,
        text=True
    )
    
    return result