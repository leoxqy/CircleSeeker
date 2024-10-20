#!/usr/bin/env python3
# coding: utf-8

import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import logging
import os

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class TreatOther:
    def __init__(self, xecc_fai, xecc_fa, cecc_csv, 
                 cecc_fa, uecc_csv, xecc_csv_output, 
                 cecc_csv_output, 
                 uecc_inferred_csv_output, xecc_fa_output, 
                 cecc_fa_output):
        self.xecc_fai = xecc_fai
        self.xecc_fa = xecc_fa
        self.cecc_csv = cecc_csv
        self.cecc_fa = cecc_fa
        self.uecc_csv = uecc_csv
        self.xecc_csv_output = xecc_csv_output
        self.xecc_fa_output = xecc_fa_output
        self.cecc_csv_output = cecc_csv_output
        self.cecc_fa_output = cecc_fa_output
        self.uecc_inferred_csv_output = uecc_inferred_csv_output
        
    def process_xecc(self):
        logging.info("Processing Xecc data...")
        
        # Read fai file
        df = pd.read_csv(self.xecc_fai, sep='\t', header=None, names=['query_id', 'eLength', 'offset', 'linebases', 'linewidth'])
        
        # Add required columns
        df['eClass'] = 'Xecc'
        df['eState'] = 'Confirmed-eccDNA'
        df['eReads'] = df['query_id'].str.split('|').str[0]
        
        # Generate eName
        df['eName'] = 'Xecc|Confirmed|' + (df.index + 1).astype(str)
        
        # Select and order columns
        df_output = df[['eName', 'eLength', 'eClass', 'eState', 'eReads']]
        df_output.to_csv(self.xecc_csv_output, index=False)
        
        # Update fasta file
        name_map = dict(zip(df['query_id'], df['eName']))
        with open(self.xecc_fa_output, 'w') as out_fa:
            for record in SeqIO.parse(self.xecc_fa, "fasta"):
                if record.id in name_map:
                    record.id = name_map[record.id]
                    record.description = ""
                SeqIO.write(record, out_fa, "fasta")
        
        logging.info(f"Xecc processing complete. Output files: {self.xecc_csv_output}, {self.xecc_fa_output}")

    def process_cecc(self):
        logging.info("Processing Cecc data...")
        
        # Read CSV file
        df = pd.read_csv(self.cecc_csv)
        
        # Select columns
        columns = ['query_id', 'consLen', 'copyNum', 'LocNum', 'identity', 'alignment_length', 'subject_id', 's_start', 's_end', 'strand', 'readName']
        
        # Create a new DataFrame with selected columns
        df_processed = df[columns].copy()
        
        # Rename columns
        df_processed.columns = ['query_id', 'eLength', 'eRepeatNum', 'LocNum', 'Identity', 'Alignment_length', 'Chr', 'Start', 'End', 'Strand', 'eReads']
        
        # Add required columns
        df_processed['eClass'] = 'Cecc'
        df_processed['eState'] = 'Confirmed-eccDNA'
        
        # Generate eName
        df_processed['eName'] = df_processed.groupby('query_id').ngroup().apply(lambda x: f"Cecc|Confirmed|{x+1}")
        
        # Select and order columns
        columns_order = ['eName', 'eLength', 'eRepeatNum', 'eClass', 'eState', 'eReads', 'LocNum', 'Identity', 'Alignment_length', 'Chr', 'Start', 'End', 'Strand']
        df_output = df_processed[columns_order]
        df_output.to_csv(self.cecc_csv_output, index=False)
        
        # Update fasta file
        name_map = dict(zip(df_processed['query_id'], df_processed['eName']))
        with open(self.cecc_fa_output, 'w') as out_fa:
            for record in SeqIO.parse(self.cecc_fa, "fasta"):
                if record.id in name_map:
                    record.id = name_map[record.id]
                    record.description = ""
                SeqIO.write(record, out_fa, "fasta")
        
        logging.info(f"Cecc processing complete. Output files: {self.cecc_csv_output}, {self.cecc_fa_output}")

    def process_uecc(self):
        logging.info("Processing Uecc data...")
        
        # Read CSV file
        df = pd.read_csv(self.uecc_csv)
        
        # Select columns and create a new DataFrame
        columns = ['Chromosome', 'Start', 'End', 'Coverage', 'Num_Reads', 'Avg_Read_Length', 'Num_Split_Reads', 'Avg_Split_Read_Length', 'Score']
        df_processed = df[columns].copy()
        
        # Rename columns
        df_processed.columns = ['eChr', 'eStart', 'eEnd', 'Coverage', 'Num_Reads', 'Avg_Read_Length', 'Num_Split_Reads', 'Avg_Split_Read_Length', 'Score']
        
        # Add required columns
        df_processed['eLength'] = df_processed['eEnd'] - df_processed['eStart'] + 1
        df_processed['eName'] = df_processed.apply(lambda row: f"Uecc|{row['eChr']}-{row['eStart']}-{row['eEnd']}", axis=1)
        df_processed['eClass'] = 'Uecc'
        df_processed['eState'] = 'Inferred-eccDNA'
        
        # Add filter condition: eLength >= 100
        df_processed = df_processed[df_processed['eLength'] >= 100]
        
        # Select and order columns
        columns_order = ['eName', 'eChr', 'eStart', 'eEnd', 'eLength', 'Num_Reads', 'Avg_Read_Length', 'Num_Split_Reads', 'Avg_Split_Read_Length', 'Score']
        df_output = df_processed[columns_order]
        df_output.to_csv(self.uecc_inferred_csv_output, index=False)
        
        logging.info(f"Uecc processing complete. Output file: {self.uecc_inferred_csv_output}")

    def process_eccDNA(self):
        self.process_xecc()
        self.process_cecc()
        self.process_uecc()

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Process Xecc, Cecc, and Uecc DNA data.")
    parser.add_argument("--xecc_fai", required=True, help="Path to xecc.fa.fai file")
    parser.add_argument("--xecc_fa", required=True, help="Path to xecc.fa file")
    parser.add_argument("--cecc_csv", required=True, help="Path to CeccDNA.csv file")
    parser.add_argument("--cecc_fa", required=True, help="Path to cecc.fa file")
    parser.add_argument("--uecc_csv", required=True, help="Path to circular_dna_results_scored.csv file")
    parser.add_argument("--prefix", required=True, help="Prefix for output files")
    args = parser.parse_args()

    treat_other = TreatOther(
        args.xecc_fai, args.xecc_fa, args.cecc_csv, args.cecc_fa, args.uecc_csv,
        args.prefix
    )
    xecc_csv, cecc_csv, xecc_fa, cecc_fa, uecc_inferred_csv = treat_other.process_eccDNA()
    print(f"Processing complete. Results saved:")
    print(f"XeccDNA CSV: {xecc_csv}")
    print(f"XeccDNA FASTA: {xecc_fa}")
    print(f"CeccDNA CSV: {cecc_csv}")
    print(f"CeccDNA FASTA: {cecc_fa}")
    print(f"UeccDNA inferred CSV: {uecc_inferred_csv}")
