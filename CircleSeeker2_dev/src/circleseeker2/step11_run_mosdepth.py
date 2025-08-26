#!/usr/bin/env python3
"""
run_mosdepth_depth.py - Mosdepth depth calculation script
"""

import subprocess
import shutil
import sys
import logging
import argparse
import gzip
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class MosdepthRunner:
    """Mosdepth depth calculator"""
    
    def __init__(self, num_threads=8, window_size=50):
        self.num_threads = num_threads
        self.window_size = window_size
        self._check_installation()
    
    def _check_installation(self):
        """Check required tools"""
        if not shutil.which("mosdepth"):
            logger.error("mosdepth not found in PATH")
            logger.error("Please install: conda install -c bioconda mosdepth")
            raise RuntimeError("mosdepth not installed")
        if not shutil.which("samtools"):
            logger.error("samtools not found in PATH")
            raise RuntimeError("samtools not installed")
        logger.debug("Required tools found")
    
    def run(self, input_bam, output_prefix):
        """
        Calculate depth using mosdepth
        
        Args:
            input_bam: Input BAM file path
            output_prefix: Output prefix for files
        
        Returns:
            Path to decompressed per-base BED file
        """
        input_bam = Path(input_bam)
        if not input_bam.exists():
            raise FileNotFoundError(f"BAM file not found: {input_bam}")
        
        # Check BAM index
        bai_file = Path(f"{input_bam}.bai")
        if not bai_file.exists():
            logger.info("Creating BAM index...")
            subprocess.run(["samtools", "index", str(input_bam)], check=True)
        
        # Build mosdepth output prefix with window size
        mosdepth_prefix = f"{output_prefix}.win{self.window_size}"
        
        # Run mosdepth
        command = [
            "mosdepth",
            "-t", str(self.num_threads),
            "--by", str(self.window_size),
            mosdepth_prefix,
            str(input_bam)
        ]
        
        logger.info(f"Running mosdepth: {input_bam} -> {mosdepth_prefix}")
        logger.debug(f"Command: {' '.join(command)}")
        
        try:
            result = subprocess.run(
                command,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
            logger.info("Mosdepth completed")
        except subprocess.CalledProcessError as e:
            logger.error(f"mosdepth failed with exit code {e.returncode}")
            logger.error(f"Error message: {e.stderr}")
            raise
        
        # Decompress per-base.bed.gz
        per_base_gz = Path(f"{mosdepth_prefix}.per-base.bed.gz")
        per_base_bed = Path(f"{mosdepth_prefix}.per-base.bed")
        
        if per_base_gz.exists():
            logger.info("Decompressing per-base.bed.gz...")
            with gzip.open(per_base_gz, 'rt') as gz_file:
                with open(per_base_bed, 'w') as out_file:
                    for line in gz_file:
                        out_file.write(line)
            logger.info(f"Output: {per_base_bed}")
        elif not per_base_bed.exists():
            raise FileNotFoundError("mosdepth per-base output not found")
        
        # Clean up other files (keep .csi files for potential future use)
        self._cleanup_files(mosdepth_prefix)
        
        return per_base_bed
    
    def _cleanup_files(self, prefix):
        """Remove unnecessary mosdepth output files but keep .csi index files"""
        files_to_remove = [
            f"{prefix}.per-base.bed.gz",  # Remove gz after decompression
            f"{prefix}.regions.bed.gz",
            f"{prefix}.mosdepth.global.dist.txt",
            f"{prefix}.mosdepth.region.dist.txt",
            f"{prefix}.mosdepth.summary.txt"
        ]
        
        # Also remove .csi files if you don't need them
        files_to_remove.extend([
            f"{prefix}.per-base.bed.gz.csi",
            f"{prefix}.regions.bed.gz.csi"
        ])
        
        for file_path in files_to_remove:
            path = Path(file_path)
            if path.exists():
                path.unlink()
                logger.debug(f"Removed: {file_path}")

def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description="Run mosdepth for depth calculation"
    )
    parser.add_argument(
        "-i", "--input",
        type=Path,
        required=True,
        help="Input sorted BAM file"
    )
    parser.add_argument(
        "-o", "--output-prefix",
        type=str,
        required=True,
        help="Output prefix (will create: prefix.winN.per-base.bed)"
    )
    parser.add_argument(
        "-w", "--window-size",
        type=int,
        default=50,
        help="Window size for depth calculation (default: 50)"
    )
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=8,
        help="Number of threads (default: 8)"
    )
    parser.add_argument(
        "--keep-csi",
        action="store_true",
        help="Keep .csi index files"
    )
    parser.add_argument(
        "--keep-all",
        action="store_true",
        help="Keep all mosdepth output files"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose logging"
    )
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    if not args.input.exists():
        logger.error(f"Input BAM file not found: {args.input}")
        sys.exit(1)
    
    try:
        runner = MosdepthRunner(
            num_threads=args.threads,
            window_size=args.window_size
        )
        
        # Modify cleanup behavior if requested
        if args.keep_all:
            runner._cleanup_files = lambda x: None
        elif args.keep_csi:
            original_cleanup = runner._cleanup_files
            def cleanup_no_csi(prefix):
                files_to_remove = [
                    f"{prefix}.per-base.bed.gz",
                    f"{prefix}.regions.bed.gz",
                    f"{prefix}.mosdepth.global.dist.txt",
                    f"{prefix}.mosdepth.region.dist.txt",
                    f"{prefix}.mosdepth.summary.txt"
                ]
                for file_path in files_to_remove:
                    path = Path(file_path)
                    if path.exists():
                        path.unlink()
            runner._cleanup_files = cleanup_no_csi
        
        per_base_bed = runner.run(args.input, args.output_prefix)
        sys.exit(0)
    except Exception as e:
        logger.error(f"Failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
