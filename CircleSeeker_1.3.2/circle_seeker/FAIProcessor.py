#!/usr/bin/env python3

"""
FAIProcessor: FASTA Index Processing for Circular DNA Analysis.

This module manages FASTA index (FAI) file processing for different types of
circular DNA sequences (UECC, MECC, CECC). It provides utilities for handling
and processing various types of circular DNA FASTA files and their indices.

Key features:
- Support for multiple circular DNA types (UECC, MECC, CECC)
- Automated FAI file generation and processing
- Efficient file handling and validation
- Comprehensive error handling and logging
- Support for inferred circular DNA sequences

Typical usage:
    processor = FAIProcessor(ueccdna="ucc.fa", meccdna="mcc.fa")
    processor.process()

Copyright (c) 2024 CircleSeeker Team
"""

import argparse
import logging
import os
import subprocess
import sys

class FAIProcessor:
    def __init__(self, ueccdna=None, meccdna=None, ceccdna=None, mceccdna=None, inferred_ueccdna=None):
        """
        Initialize FAIProcessor with file paths
        
        Args:
            ueccdna (str): Path to UeccDNA FASTA file
            meccdna (str): Path to MeccDNA FASTA file
            ceccdna (str): Path to CeccDNA FASTA file
            mceccdna (str): Path to MCeccDNA FASTA file
            inferred_ueccdna (str): Path to Inferred UeccDNA FASTA file
        """
        self.ueccdna = ueccdna
        self.meccdna = meccdna
        self.ceccdna = ceccdna
        self.mceccdna = mceccdna
        self.inferred_ueccdna = inferred_ueccdna
        self.setup_logging()

    def setup_logging(self):
        """Configure logging settings"""
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[logging.StreamHandler()]
        )

    def is_file_empty(self, filepath):
        """Check if file is empty"""
        return os.path.getsize(filepath) == 0

    def get_csv_file(self, fasta_file):
        """Get corresponding CSV file path"""
        return fasta_file.replace('.fasta', '.csv')  # 修改这里：.fa -> .fasta

    def remove_related_files(self, fasta_file, check_csv=False):
        """Remove FASTA file and its related files"""
        try:
            # Remove FASTA file
            if os.path.exists(fasta_file):
                os.remove(fasta_file)
                logging.info(f"Deleted empty FASTA file: {fasta_file}")
            
            # Remove related CSV file
            if check_csv:
                csv_file = self.get_csv_file(fasta_file)
                if os.path.exists(csv_file):
                    os.remove(csv_file)
                    logging.info(f"Deleted related CSV file: {csv_file}")
                
            # Remove index file (if exists)
            fai_file = fasta_file + '.fai'
            if os.path.exists(fai_file):
                os.remove(fai_file)
                logging.info(f"Deleted index file: {fai_file}")
                
        except Exception as e:
            logging.error(f"Error while deleting files: {str(e)}")

    def check_csv_exists(self, fasta_file):
        """Check if corresponding CSV file exists"""
        csv_file = self.get_csv_file(fasta_file)
        if not os.path.exists(csv_file):
            logging.error(f"Corresponding CSV file not found: {csv_file}")
            return False
        return True

    def process_fasta_file(self, fasta_file, require_csv=False):
        """Process single FASTA file"""
        if not os.path.exists(fasta_file):
            logging.warning(f"File does not exist: {fasta_file}")
            return False
        
        logging.info(f"Starting to process file: {fasta_file}")
        
        # Check if CSV file exists if required
        if require_csv and not self.check_csv_exists(fasta_file):
            return False
        
        # Check if file is empty
        if self.is_file_empty(fasta_file):
            logging.warning(f"File is empty: {fasta_file}")
            self.remove_related_files(fasta_file, check_csv=require_csv)
            return False
        
        # Execute samtools faidx command
        try:
            result = subprocess.run(
                ['samtools', 'faidx', fasta_file], 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE,
                text=True
            )
            
            if result.returncode == 0:
                logging.info(f"Successfully created index for file: {fasta_file}")
                return True
            else:
                logging.error(f"Error creating index for file: {fasta_file}")
                logging.error(f"Error message: {result.stderr}")
                return False
                
        except Exception as e:
            logging.error(f"Error executing samtools command: {str(e)}")
            return False

    def process_all_files(self):
        """Process all FASTA files"""
        logging.info("Starting FASTA file processing script")

        # Process all files
        # UeccDNA and MeccDNA don't require CSV check
        if self.ueccdna:
            self.process_fasta_file(self.ueccdna)
        if self.meccdna:
            self.process_fasta_file(self.meccdna)
        if self.inferred_ueccdna:
            self.process_fasta_file(self.inferred_ueccdna)
        
        # CeccDNA and MCeccDNA require CSV check
        if self.ceccdna:
            self.process_fasta_file(self.ceccdna, require_csv=True)
        if self.mceccdna:
            self.process_fasta_file(self.mceccdna, require_csv=True)

        logging.info("Script execution completed")

def main():
    # Create argument parser
    parser = argparse.ArgumentParser(description='Process FASTA files and create indices')
    parser.add_argument('-uc', '--ueccdna', type=str, required=True, help='UeccDNA FASTA file path')
    parser.add_argument('-ui', '--inferred_ueccdna', type=str, required=True, help='Inferred UeccDNA FASTA file path')
    parser.add_argument('-m', '--meccdna', type=str, required=True, help='MeccDNA FASTA file path')
    parser.add_argument('-c', '--ceccdna', type=str, required=True, help='CeccDNA FASTA file path')
    parser.add_argument('-mc', '--mceccdna', type=str, required=True, help='MCeccDNA FASTA file path')
    args = parser.parse_args()

    # Initialize and run processor
    processor = FAIProcessor(
        ueccdna=args.ueccdna,
        meccdna=args.meccdna,
        ceccdna=args.ceccdna,
        mceccdna=args.mceccdna,
        inferred_ueccdna=args.inferred_ueccdna
    )
    processor.process_all_files()

if __name__ == "__main__":
    main()