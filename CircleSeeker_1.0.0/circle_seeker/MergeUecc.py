#!/usr/bin/env python3
# coding: utf-8

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import logging
from collections import OrderedDict
import os

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class MergeUecc:
    def __init__(self, ueccdna_part1, tecc_analysis, uecc_part1, uecc_part2, output_csv, output_fasta):
        self.ueccdna_part1 = ueccdna_part1
        self.tecc_analysis = tecc_analysis
        self.uecc_part1 = uecc_part1
        self.uecc_part2 = uecc_part2
        self.output_csv = output_csv
        self.output_fasta = output_fasta

    def read_csv(self, file_path):
        return pd.read_csv(file_path)

    def process_dataframes(self, df1, df2):
        df1_processed = df1[['eChr', 'eStart', 'eEnd', 'eStrand', 'eLength', 'eRepeatNum', 'MatDegree', 'eClass', 'eReads', 'query_id']]
        df2_processed = df2[df2['eClass'] == 'Uecc'][['eChr', 'eStart', 'eEnd', 'eStrand', 'eLength', 'eRepeatNum', 'MatDegree', 'eClass', 'eReads']]
        df2_processed['query_id'] = df2_processed.apply(lambda row: f"{row['eChr']}-{row['eStart']}-{row['eEnd']}", axis=1)
        combined_df = pd.concat([df1_processed, df2_processed], ignore_index=True)
        combined_df['eName'] = combined_df.apply(lambda row: f"{row['eClass']}|{row['eChr']}-{row['eStart']}-{row['eEnd']}", axis=1)
        combined_df['eState'] = 'Confirmed-eccDNA'
        return combined_df

    def merge_fasta_files(self, file1, file2, output_file):
        records = list(SeqIO.parse(file1, "fasta")) + list(SeqIO.parse(file2, "fasta"))
        SeqIO.write(records, output_file, "fasta")

    def update_fasta_names(self, fasta_file, name_map, output_file):
        records = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id in name_map:
                record.id = name_map[record.id]
                record.description = ""
            records.append(record)
        SeqIO.write(records, output_file, "fasta")

    def remove_duplicates(self, input_file, output_file):
        seen_sequences = OrderedDict()
        duplicates = 0
        for record in SeqIO.parse(input_file, "fasta"):
            sequence = str(record.seq)
            if sequence not in seen_sequences:
                seen_sequences[sequence] = record
            else:
                duplicates += 1

        with open(output_file, "w") as output_handle:
            SeqIO.write(seen_sequences.values(), output_handle, "fasta")
        
        logging.info(f"Removed {duplicates} duplicate sequences")
        return len(seen_sequences)

    def process_final_dataframe(self, df):
        grouped = df.groupby('eName').agg({
            'eChr': 'first',
            'eStart': 'first',
            'eEnd': 'first',
            'eStrand': lambda x: max(set(x), key=list(x).count) if len(set(x)) > 1 else x.iloc[0],
            'eLength': 'first',
            'eRepeatNum': 'sum',
            'MatDegree': 'mean',
            'eClass': 'first',
            'eState': 'first',
            'eReads': lambda x: ','.join(set(x))
        }).reset_index()

        grouped['eReadNum'] = grouped['eReads'].apply(lambda x: len(set(x.split(','))))
        grouped['MatDegree'] = grouped['MatDegree'].round(2)

        columns_order = ['eName', 'eChr', 'eStart', 'eEnd', 'eStrand', 'eLength', 'eRepeatNum', 
                         'MatDegree', 'eClass', 'eState', 'eReadNum', 'eReads']
        return grouped[columns_order]

    def run_merge_uecc(self):
        logging.info("Reading CSV files...")
        df1 = self.read_csv(self.ueccdna_part1)
        df2 = self.read_csv(self.tecc_analysis)
        
        logging.info("Processing DataFrames...")
        combined_df = self.process_dataframes(df1, df2)
        
        logging.info("Processing final DataFrame...")
        final_df = self.process_final_dataframe(combined_df)
        
        logging.info("Saving processed CSV...")
        final_df.to_csv(self.output_csv, index=False)
        
        logging.info("Merging FASTA files...")
        self.merge_fasta_files(self.uecc_part1, self.uecc_part2, "temp_merged.fa")
        
        logging.info("Updating FASTA sequence names...")
        name_map = dict(zip(combined_df['query_id'], combined_df['eName']))
        self.update_fasta_names("temp_merged.fa", name_map, "temp_renamed.fa")
        
        logging.info("Removing duplicate sequences...")
        unique_count = self.remove_duplicates("temp_renamed.fa", self.output_fasta)
        logging.info(f"Final FASTA file contains {unique_count} unique sequences")
        
        # Clean up temporary files
        os.remove("temp_merged.fa")
        os.remove("temp_renamed.fa")
        
        logging.info("Process completed successfully.")

def main():
    parser = argparse.ArgumentParser(description="Process eccDNA data, generate output files, and remove duplicates.")
    parser.add_argument("--ueccdna_part1", required=True, help="Input UeccDNA_part1.csv file")
    parser.add_argument("--tecc_analysis", required=True, help="Input tecc_analysis_results.csv file")
    parser.add_argument("--uecc_part1", required=True, help="Input uecc_part1.fa file")
    parser.add_argument("--uecc_part2", required=True, help="Input uecc_part2.fa file")
    parser.add_argument("--output_csv", required=True, help="Output CSV file")
    parser.add_argument("--output_fasta", required=True, help="Output FASTA file (deduplicated)")
    args = parser.parse_args()

    merger = MergeUecc(args.ueccdna_part1, args.tecc_analysis, args.uecc_part1, args.uecc_part2, args.output_csv, args.output_fasta)
    merger.run_merge_uecc()

if __name__ == "__main__":
    main()