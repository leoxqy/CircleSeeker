"""
Sieve: FASTA Read Filtering and Processing Module.

This module provides functionality for filtering and processing FASTA sequence
reads based on a list of read identifiers. It efficiently processes large FASTA
files while maintaining detailed statistics about the filtering process.

Key features:
- FASTA file filtering based on read lists
- Support for both tab and comma delimiters
- Read count tracking and statistics
- Memory-efficient processing
- Comprehensive logging system

Typical usage:
    sieve = Sieve(input_fasta, read_list_file, output_fasta)
    sieve.run()

Copyright (c) 2024 CircleSeeker Team
"""

import argparse
import subprocess
import os
import logging

class Sieve():
    def __init__(self, input_fasta, read_list_file, output_fasta, delimiter = 't'):
        # logging settings
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        
        self.input_fasta = input_fasta
        self.read_list_file = read_list_file
        self.output_fasta = output_fasta
        self.delimiter ='\t' if delimiter == 't' else ','
        self.total_reads = 0
        self.unmatched_reads = 0
    
    @staticmethod
    def run_command(command):
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            raise Exception(f"Command failed: {command}\nError: {stderr.decode()}")
        return stdout.decode()

    def extract_unmatched_reads(self):
        # Extract matched read names from the file and remove duplicates
        delimiter = self.delimiter
        command = f"cut -d'{delimiter}' -f1 {self.read_list_file} | sort | uniq > matched_reads.txt"
        self.run_command(command)
        
        # Use seqkit to extract unmatched reads
        command = f"seqkit grep -v -f matched_reads.txt {self.input_fasta} > {self.output_fasta}"
        self.run_command(command)
        
        # Clean up temporary file
        os.remove("matched_reads.txt")

    def count_reads(self, file_path):
        command = f"seqkit stats {file_path} | awk 'NR>1 {{print $4}}'"
        output = self.run_command(command).strip()
        # Remove commas from the output before converting to int
        count = int(output.replace(',', ''))
        return count

    def run(self):
        # Count total reads in input file
        logging.info("Counting total reads in input file...")
        self.total_reads = self.count_reads(self.input_fasta)
        
        # Extract unmatched reads
        logging.info("Extracting unmatched reads...")
        self.extract_unmatched_reads()
        
        # Count unmatched reads
        logging.info("Counting unmatched reads...")
        self.unmatched_reads = self.count_reads(self.output_fasta)
        
        logging.info("Processing complete. Results:")
        logging.info(f"1. Unmatched reads file: {self.output_fasta}")
        logging.info(f"2. Total number of reads in the input fasta file: {self.total_reads}")
        logging.info(f"3. Number of unmatched reads: {self.unmatched_reads}")
        logging.info(f"4. Percentage of unmatched reads: {(self.unmatched_reads / self.total_reads) * 100:.2f}%")

def main():
    parser = argparse.ArgumentParser(description='Extract unmatched reads from a FASTA file')
    parser.add_argument('-i', '--input', required=True, help='Input FASTA file')
    parser.add_argument('-l', '--list', required=True, help='Read list file')
    parser.add_argument('-o', '--output', required=True, help='Output FASTA file')
    parser.add_argument('-d', '--delimiter', default='t', choices=['t', 'c'],
                      help='Delimiter in read list file: "t" for tab (default), "c" for comma')

    args = parser.parse_args()

    sieve = Sieve(args.input, args.list, args.output, args.delimiter)
    sieve.run()

if __name__ == '__main__':
    main()