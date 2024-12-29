#!/usr/bin/env python3
# coding: utf-8

"""
Tamer: Sequence Processing and Standardization Module.

This module handles the processing and standardization of circular DNA sequences,
managing both UECC and MECC data types. It provides functionality for sequence
extraction, naming convention application, and data format integration.

Key features:
- Sequence extraction and processing
- Standardized naming convention application
- FASTA and CSV file integration
- Efficient sequence dictionary management
- Comprehensive logging system

Typical usage:
    tamer = Tamer(input_csv, input_fasta, output_csv)
    tamer.run()

Copyright (c) 2024 CircleSeeker Team
"""

import argparse
import pandas as pd
from Bio import SeqIO
import logging

class Tamer:
    def __init__(self, input_tecc_csv, input_fasta, output_csv):
        self.input_csv = input_tecc_csv
        self.input_fasta = input_fasta
        self.output_csv = output_csv
        self.df = None
        self.fasta_dict = None
        
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    def read_fasta(self):
        self.fasta_dict = SeqIO.to_dict(SeqIO.parse(self.input_fasta, "fasta"))
        logging.info(f"Read FASTA file: {self.input_fasta}")

    def process_sequences(self):
        # Add new columns for eName and eSeq
        self.df['eName'] = ''
        self.df['eSeq'] = ''
        
        # Process all sequences
        for idx, row in self.df.iterrows():
            if row['eClass'] == 'Uecc':
                self.df.loc[idx, 'eName'] = f"{row['eChr']}-{row['eStart']}-{row['eEnd']}"
            elif row['eClass'] == 'Mecc':
                self.df.loc[idx, 'eName'] = (
                    f"{row['qname']}|Second|{row['eLength']}|{row['eRepeatNum']}|circular"
                )
            
            # Extract sequence
            if row['qname'] in self.fasta_dict:
                seq = self.fasta_dict[row['qname']][int(row['qstart']):int(row['qend'])]
                self.df.loc[idx, 'eSeq'] = str(seq.seq)
            else:
                logging.warning(f"Sequence {row['qname']} not found in FASTA file. Skipping.")
        
        logging.info(f"Processed {len(self.df)} sequences")

    def run(self):
        try:
            self.df = pd.read_csv(self.input_csv)
            logging.info(f"Read CSV file: {self.input_csv}")
            
            self.read_fasta()
            self.process_sequences()
            
            # Save the enhanced CSV file
            self.df.to_csv(self.output_csv, index=False)
            logging.info(f"Wrote enhanced CSV file to {self.output_csv}")

        except Exception as e:
            logging.error(f"An error occurred: {e}")

def main():
    parser = argparse.ArgumentParser(
        description='Process TECC CSV and FASTA files to generate enhanced CSV with sequence information'
    )
    
    parser.add_argument('-i', '--input-csv', required=True,
                      help='Input TECC CSV file')
    
    parser.add_argument('-f', '--fasta', required=True,
                      help='Input FASTA file')
    
    parser.add_argument('-o', '--output-csv', required=True,
                      help='Output enhanced CSV file')

    args = parser.parse_args()

    # Initialize and run Tamer
    tamer = Tamer(
        input_tecc_csv=args.input_csv,
        input_fasta=args.fasta,
        output_csv=args.output_csv
    )
    tamer.run()

if __name__ == '__main__':
    main()