#!/usr/bin/env python3
"""
Sieve: Optimized FASTA Read Filtering Module

This module filters FASTA sequences based on a CSV classification file,
specifically removing sequences classified as CtcR variants.

Key features:
- Reads CSV file with readName and readClass columns
- Filters out sequences with specific CtcR classifications
- Pure Python implementation with optional seqkit support
- Memory-efficient processing for large files
- Comprehensive logging and statistics

Version: 2.0.0
License: GNU General Public License v2
Copyright (c) 2024 CircleSeeker Team
"""

import argparse
import subprocess
import csv
import sys
import logging
from typing import Set, List, Optional
from dataclasses import dataclass
from pathlib import Path


@dataclass
class FilterStats:
    """Statistics for filtering operation"""
    total_reads: int = 0
    filtered_reads: int = 0
    retained_reads: int = 0
    csv_total_reads: int = 0
    csv_ctcr_reads: int = 0
    
    @property
    def filtered_percentage(self) -> float:
        if self.total_reads == 0:
            return 0.0
        return (self.filtered_reads / self.total_reads) * 100
    
    @property
    def retained_percentage(self) -> float:
        if self.total_reads == 0:
            return 0.0
        return (self.retained_reads / self.total_reads) * 100


class Sieve:
    """FASTA sequence filter based on CSV classification file"""
    
    # Default CtcR classes to filter out
    DEFAULT_CTCR_CLASSES = {
        "CtcR-perfect",
        "CtcR-inversion", 
        "CtcR-hybrid"
    }
    
    def __init__(self, 
                 input_fasta: str, 
                 read_list_file: str,
                 output_fasta: str,
                 delimiter: str = ',',
                 ctcr_classes: Optional[Set[str]] = None,
                 use_seqkit: bool = True):
        """
        Initialize Sieve filter
        
        Args:
            input_fasta: Path to input FASTA file
            read_list_file: Path to CSV file with readName and readClass columns
            output_fasta: Path to output FASTA file
            delimiter: CSV delimiter (default: comma)
            ctcr_classes: Set of classes to filter out (default: CtcR variants)
            use_seqkit: Whether to use seqkit for statistics (if available)
        """
        # Setup logging
        self._setup_logging()
        
        # File paths
        self.input_fasta = Path(input_fasta)
        self.read_list_file = Path(read_list_file)
        self.output_fasta = Path(output_fasta)
        
        # Validate input files
        if not self.input_fasta.exists():
            raise FileNotFoundError(f"Input FASTA file not found: {input_fasta}")
        if not self.read_list_file.exists():
            raise FileNotFoundError(f"Read list file not found: {read_list_file}")
        
        # Settings
        self.delimiter = delimiter
        self.ctcr_classes = ctcr_classes or self.DEFAULT_CTCR_CLASSES
        self.use_seqkit = use_seqkit and self._check_seqkit_available()
        
        # Data structures
        self.reads_to_filter = set()  # Set of read names to filter out
        
        # Statistics
        self.stats = FilterStats()
        
        logging.info(f"Initialized Sieve")
        logging.info(f"Input FASTA: {self.input_fasta}")
        logging.info(f"Read list: {self.read_list_file}")
        logging.info(f"Output FASTA: {self.output_fasta}")
        logging.info(f"Filtering classes: {self.ctcr_classes}")
        logging.info(f"Using seqkit: {self.use_seqkit}")
    
    def _setup_logging(self):
        """Configure logging settings"""
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
    
    def _check_seqkit_available(self) -> bool:
        """Check if seqkit is available in system PATH"""
        try:
            result = subprocess.run(['seqkit', 'version'], 
                                  capture_output=True, 
                                  text=True,
                                  timeout=5)
            return result.returncode == 0
        except (subprocess.SubprocessError, FileNotFoundError):
            logging.warning("seqkit not found, using pure Python implementation")
            return False
    
    def load_read_list(self):
        """
        Load read names from CSV file that match the CtcR classes to filter
        
        Expected CSV format:
        readName,readClass
        m84074_250307_021057_s1/100272738/ccs,CtcR-inversion
        ...
        """
        logging.info(f"Loading read list from {self.read_list_file}...")
        
        # Try to detect delimiter if it's tab-separated
        if self.delimiter == 't':
            self.delimiter = '\t'
        
        # Track statistics for debugging
        class_counts = {}
        duplicate_reads = set()
        all_reads_seen = set()
        
        with open(self.read_list_file, 'r') as f:
            # Try to detect the actual delimiter
            first_line = f.readline().strip()
            f.seek(0)  # Reset to beginning
            
            # Auto-detect delimiter if needed
            if ',' in first_line and '\t' not in first_line:
                actual_delimiter = ','
            elif '\t' in first_line and ',' not in first_line:
                actual_delimiter = '\t'
            else:
                actual_delimiter = self.delimiter
            
            reader = csv.DictReader(f, delimiter=actual_delimiter)
            
            # Check if required columns exist
            if 'readName' not in reader.fieldnames or 'readClass' not in reader.fieldnames:
                # Try case-insensitive match
                fieldnames_lower = [fn.lower() for fn in reader.fieldnames]
                if 'readname' not in fieldnames_lower or 'readclass' not in fieldnames_lower:
                    raise ValueError(f"CSV file must contain 'readName' and 'readClass' columns. "
                                   f"Found columns: {reader.fieldnames}")
                
                # Map to correct column names
                readname_col = reader.fieldnames[fieldnames_lower.index('readname')]
                readclass_col = reader.fieldnames[fieldnames_lower.index('readclass')]
            else:
                readname_col = 'readName'
                readclass_col = 'readClass'
            
            # Load reads to filter
            for row in reader:
                self.stats.csv_total_reads += 1
                read_name = row[readname_col].strip()
                read_class = row[readclass_col].strip()
                
                # Track all occurrences
                if read_name in all_reads_seen:
                    duplicate_reads.add(read_name)
                all_reads_seen.add(read_name)
                
                # Count by class
                class_counts[read_class] = class_counts.get(read_class, 0) + 1
                
                # Check if this read should be filtered
                if read_class in self.ctcr_classes:
                    self.reads_to_filter.add(read_name)
                    self.stats.csv_ctcr_reads += 1
        
        # Print detailed statistics
        logging.info(f"Loaded {self.stats.csv_total_reads} total rows from CSV")
        logging.info(f"Class distribution:")
        for cls, count in sorted(class_counts.items()):
            if cls in self.ctcr_classes:
                logging.info(f"  {cls}: {count} (will filter)")
            else:
                logging.info(f"  {cls}: {count}")
        
        logging.info(f"Total CtcR classifications: {self.stats.csv_ctcr_reads}")
        logging.info(f"Unique reads to filter: {len(self.reads_to_filter)}")
        
        if duplicate_reads:
            logging.warning(f"Found {len(duplicate_reads)} duplicate read names in CSV")
            logging.debug(f"First 5 duplicates: {list(duplicate_reads)[:5]}")
    
    def filter_fasta(self):
        """
        Filter FASTA file using pure Python implementation
        Removes sequences whose IDs are in the filter set
        """
        logging.info(f"Filtering {self.input_fasta}...")
        
        # Track what we actually filter for debugging
        filtered_ids = []
        not_in_csv_but_filtered = []
        fasta_read_ids = set()
        
        with open(self.input_fasta, 'r') as infile, \
             open(self.output_fasta, 'w') as outfile:
            
            current_header = None
            current_sequence = []
            write_current = True
            current_read_id = None
            
            for line in infile:
                line = line.rstrip('\n')
                
                if line.startswith('>'):
                    # Process previous sequence if exists
                    if current_header is not None:
                        self.stats.total_reads += 1
                        fasta_read_ids.add(current_read_id)
                        
                        if write_current:
                            # Write retained sequence
                            outfile.write(f"{current_header}\n")
                            outfile.write('\n'.join(current_sequence) + '\n')
                            self.stats.retained_reads += 1
                        else:
                            # Track filtered sequence
                            self.stats.filtered_reads += 1
                            filtered_ids.append(current_read_id)
                    
                    # Start new sequence
                    current_header = line
                    current_sequence = []
                    
                    # Extract read name from header
                    # Header format: >readName optional_description
                    current_read_id = line[1:].split()[0]  # Remove '>' and get first part
                    
                    # Check if this read should be filtered
                    write_current = current_read_id not in self.reads_to_filter
                    
                elif line:  # Non-empty sequence line
                    current_sequence.append(line)
            
            # Process last sequence
            if current_header is not None:
                self.stats.total_reads += 1
                fasta_read_ids.add(current_read_id)
                
                if write_current:
                    outfile.write(f"{current_header}\n")
                    outfile.write('\n'.join(current_sequence) + '\n')
                    self.stats.retained_reads += 1
                else:
                    self.stats.filtered_reads += 1
                    filtered_ids.append(current_read_id)
        
        # Diagnostic information
        logging.info(f"FASTA processing complete:")
        logging.info(f"  Total unique read IDs in FASTA: {len(fasta_read_ids)}")
        logging.info(f"  Sequences filtered: {len(filtered_ids)}")
        
        # Check for mismatches
        reads_in_csv_not_in_fasta = self.reads_to_filter - fasta_read_ids
        if reads_in_csv_not_in_fasta:
            logging.warning(f"Found {len(reads_in_csv_not_in_fasta)} reads in CSV filter list but NOT in FASTA")
            logging.debug(f"First 5 examples: {list(reads_in_csv_not_in_fasta)[:5]}")
        
        # Check if we have duplicates in FASTA
        if len(filtered_ids) != len(set(filtered_ids)):
            duplicates_in_filtered = len(filtered_ids) - len(set(filtered_ids))
            logging.warning(f"Found {duplicates_in_filtered} duplicate sequences in FASTA that were filtered")
    
    def filter_with_seqkit(self):
        """
        Alternative: Use seqkit for filtering (faster for very large files)
        """
        if not self.use_seqkit:
            return self.filter_fasta()
        
        logging.info("Using seqkit for filtering...")
        
        try:
            # Create temporary file with read names to filter
            temp_file = Path("temp_reads_to_filter.txt")
            with open(temp_file, 'w') as f:
                for read_name in self.reads_to_filter:
                    f.write(f"{read_name}\n")
            
            # Use seqkit grep with -v flag to exclude matches
            cmd = f"seqkit grep -v -f {temp_file} {self.input_fasta} -o {self.output_fasta}"
            
            result = subprocess.run(cmd, shell=True, 
                                  capture_output=True, 
                                  text=True,
                                  timeout=300)
            
            if result.returncode != 0:
                logging.warning(f"seqkit filtering failed: {result.stderr}")
                logging.info("Falling back to Python implementation...")
                temp_file.unlink(missing_ok=True)
                return self.filter_fasta()
            
            # Clean up
            temp_file.unlink(missing_ok=True)
            
            # Get statistics
            self._calculate_stats_from_files()
            
        except Exception as e:
            logging.warning(f"seqkit filtering failed: {e}")
            logging.info("Falling back to Python implementation...")
            return self.filter_fasta()
    
    def _calculate_stats_from_files(self):
        """Calculate statistics by counting sequences in files"""
        # Count input sequences
        self.stats.total_reads = self._count_sequences(self.input_fasta)
        
        # Count output sequences
        self.stats.retained_reads = self._count_sequences(self.output_fasta)
        
        # Calculate filtered
        self.stats.filtered_reads = self.stats.total_reads - self.stats.retained_reads
    
    def _count_sequences(self, fasta_file: Path) -> int:
        """Count sequences in a FASTA file"""
        if self.use_seqkit:
            count = self.get_seqkit_count(fasta_file)
            if count is not None:
                return count
        
        # Fallback to Python counting
        count = 0
        with open(fasta_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    count += 1
        return count
    
    def get_seqkit_count(self, file_path: Path) -> Optional[int]:
        """
        Get sequence count using seqkit
        
        Args:
            file_path: Path to FASTA file
            
        Returns:
            Number of sequences or None if seqkit fails
        """
        try:
            cmd = f"seqkit stats {file_path} -T"
            result = subprocess.run(cmd, shell=True, 
                                  capture_output=True, 
                                  text=True,
                                  timeout=30)
            
            if result.returncode == 0:
                # Parse seqkit output (tab-delimited)
                lines = result.stdout.strip().split('\n')
                if len(lines) > 1:
                    # Second line contains stats, 4th column is num_seqs
                    fields = lines[1].split('\t')
                    if len(fields) >= 4:
                        return int(fields[3].replace(',', ''))
        except (subprocess.SubprocessError, ValueError) as e:
            logging.debug(f"Failed to get seqkit count: {e}")
        
        return None
    
    def print_summary(self):
        """Print filtering summary"""
        print("\n" + "="*60)
        print("FILTERING SUMMARY")
        print("="*60)
        print(f"Input FASTA:       {self.input_fasta}")
        print(f"Read list CSV:     {self.read_list_file}")
        print(f"Output FASTA:      {self.output_fasta}")
        print("-"*60)
        print(f"CSV Statistics:")
        print(f"  Total reads:     {self.stats.csv_total_reads:,}")
        print(f"  CtcR reads:      {self.stats.csv_ctcr_reads:,}")
        print("-"*60)
        print(f"FASTA Statistics:")
        print(f"  Total sequences: {self.stats.total_reads:,}")
        print(f"  Filtered out:    {self.stats.filtered_reads:,} "
              f"({self.stats.filtered_percentage:.2f}%)")
        print(f"  Retained:        {self.stats.retained_reads:,} "
              f"({self.stats.retained_percentage:.2f}%)")
        print("="*60)
    
    def run(self):
        """Execute the filtering pipeline"""
        logging.info("Starting FASTA filtering pipeline...")
        
        # Load reads to filter from CSV
        self.load_read_list()
        
        # Perform filtering
        if self.use_seqkit:
            self.filter_with_seqkit()
        else:
            self.filter_fasta()
        
        # Display summary
        self.print_summary()
        
        logging.info(f"Filtering complete. Output saved to {self.output_fasta}")
        
        return self.stats


def main():
    """Command-line interface"""
    parser = argparse.ArgumentParser(
        description='Filter FASTA sequences based on CSV classification file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Expected CSV format:
  readName,readClass
  read_id_1,CtcR-perfect
  read_id_2,CtcR-inversion
  read_id_3,other_class
  ...

The script will filter out (remove) sequences from the FASTA file
that are classified as CtcR-perfect, CtcR-inversion, or CtcR-hybrid
in the CSV file.

Examples:
  # Basic usage (filter out CtcR sequences)
  %(prog)s -i input.fasta -l readlist.csv -o filtered.fasta
  
  # Use tab-delimited CSV
  %(prog)s -i input.fasta -l readlist.tsv -o filtered.fasta -d t
  
  # Filter specific classes only
  %(prog)s -i input.fasta -l readlist.csv -o filtered.fasta \\
           --classes CtcR-perfect CtcR-hybrid
  
  # Keep only CtcR sequences (inverse filter)
  %(prog)s -i input.fasta -l readlist.csv -o filtered.fasta --inverse
        """
    )
    
    parser.add_argument('-i', '--input', 
                       required=True,
                       help='Input FASTA file')
    
    parser.add_argument('-l', '--list',
                       required=True,
                       help='CSV file with readName and readClass columns')
    
    parser.add_argument('-o', '--output', 
                       required=True,
                       help='Output FASTA file (filtered sequences)')
    
    parser.add_argument('-d', '--delimiter',
                       default=',',
                       choices=[',', 't', 'c'],
                       help='Delimiter in CSV file: "," for comma (default), "t" for tab')
    
    parser.add_argument('--classes',
                       nargs='+',
                       help='Specific classes to filter (default: CtcR-perfect, CtcR-inversion, CtcR-hybrid)')
    
    parser.add_argument('--inverse',
                       action='store_true',
                       help='Inverse filter: keep only sequences with specified classes')
    
    parser.add_argument('--no-seqkit',
                       action='store_true',
                       help='Disable seqkit usage, use pure Python implementation')
    
    parser.add_argument('-v', '--verbose',
                       action='store_true',
                       help='Enable verbose logging')
    
    parser.add_argument('--diagnose',
                       action='store_true',
                       help='Run diagnostic mode to check ID matching')
    
    args = parser.parse_args()
    
    # Set logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Convert delimiter
    delimiter = args.delimiter
    if delimiter == 't':
        delimiter = '\t'
    elif delimiter == 'c':
        delimiter = ','
    
    try:
        # Prepare classes to filter
        if args.classes:
            ctcr_classes = set(args.classes)
        else:
            ctcr_classes = None  # Use defaults
        
        # Create filter instance
        sieve = Sieve(
            input_fasta=args.input,
            read_list_file=args.list,
            output_fasta=args.output,
            delimiter=delimiter,
            ctcr_classes=ctcr_classes,
            use_seqkit=not args.no_seqkit
        )
        
        # Diagnostic mode
        if args.diagnose:
            logging.info("Running in diagnostic mode...")
            sieve.load_read_list()
            
            # Check first few FASTA headers
            logging.info("\nChecking FASTA headers vs CSV read names...")
            with open(sieve.input_fasta, 'r') as f:
                fasta_ids = []
                for line in f:
                    if line.startswith('>'):
                        fasta_id = line[1:].split()[0]
                        fasta_ids.append(fasta_id)
                        if len(fasta_ids) <= 5:
                            logging.info(f"FASTA ID example {len(fasta_ids)}: '{fasta_id}'")
            
            # Show examples from CSV
            csv_examples = list(sieve.reads_to_filter)[:5]
            for i, csv_id in enumerate(csv_examples, 1):
                logging.info(f"CSV ID example {i}: '{csv_id}'")
            
            # Check for exact matches
            matches = 0
            for fasta_id in fasta_ids[:100]:  # Check first 100
                if fasta_id in sieve.reads_to_filter:
                    matches += 1
            
            logging.info(f"\nOut of first 100 FASTA sequences, {matches} have matching IDs in CSV filter list")
            
            if matches == 0:
                logging.warning("No matches found! ID format might be different between CSV and FASTA")
                logging.info("Please check if the readName format in CSV matches the sequence IDs in FASTA")
            
            sys.exit(0)
        
        # Handle inverse filtering
        if args.inverse:
            # For inverse mode, we want to keep only the CtcR sequences
            # So we invert the filter set
            sieve.load_read_list()
            
            # Get all read IDs from FASTA
            all_reads = set()
            with open(sieve.input_fasta, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        read_id = line[1:].split()[0]
                        all_reads.add(read_id)
            
            # Invert: filter out everything NOT in the CtcR list
            sieve.reads_to_filter = all_reads - sieve.reads_to_filter
            
            logging.info("Inverse mode: keeping only sequences WITH specified classes")
            
            # Now run filtering with inverted set
            if sieve.use_seqkit:
                sieve.filter_with_seqkit()
            else:
                sieve.filter_fasta()
            
            sieve.print_summary()
        else:
            # Normal mode
            stats = sieve.run()
        
        # Exit with appropriate code
        sys.exit(0)
        
    except Exception as e:
        logging.error(f"Error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
