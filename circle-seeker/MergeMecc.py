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

class MergeMecc:
    def __init__(self, meccdna_part1, tecc_analysis, mecc_part1, mecc_part2, output_csv, output_fasta):
        self.meccdna_part1 = meccdna_part1
        self.tecc_analysis = tecc_analysis
        self.mecc_part1 = mecc_part1
        self.mecc_part2 = mecc_part2
        self.output_csv = output_csv
        self.output_fasta = output_fasta

    def read_csv(self, file_path):
        return pd.read_csv(file_path)

    def process_meccdna_part1(self, df):
        rows = []
        for _, row in df.iterrows():
            align_regions = row['AlignRegion'].split(';')
            if row['consLen'] < 100:
                continue
            for region in align_regions:
                mat_degree, location = region.split('|')
                chr, start_end = location.split('-', 1)
                start, end = map(int, start_end.split('-'))
                rows.append({
                    'eChr': chr,
                    'eStart': start,
                    'eEnd': end,
                    'eLength': row['consLen'],
                    'eRepeatNum': row['copyNum'],
                    'MatDegree': float(mat_degree),
                    'eClass': 'Mecc',
                    'eReads': row['readName'],
                    'query_id': row['query_id']
                })
        return pd.DataFrame(rows)

    def process_tecc_analysis(self, df):
        mecc_df = df[df['eClass'] == 'Mecc'].copy()
        mecc_df = mecc_df[mecc_df['eLength'] >= 100]
        
        def get_query_id(group):
            longest_row = group.loc[group['eLength'].idxmax()]
            return f"{group.name}|Second|{longest_row['eLength']}|{longest_row['eRepeatNum']}|circular"
        
        query_ids = mecc_df.groupby('eReads').apply(get_query_id)
        mecc_df['query_id'] = mecc_df['eReads'].map(query_ids)
        
        return mecc_df[['eChr', 'eStart', 'eEnd', 'eLength', 'eRepeatNum', 'MatDegree', 'eClass', 'eReads', 'query_id']]

    def process_dataframes(self, df1, df2):
        combined_df = pd.concat([df1, df2], ignore_index=True)
        combined_df['eName'] = combined_df.groupby('eReads').ngroup().apply(lambda x: f"Mecc|Confirmed|{x+1}")
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
            'eLength': 'first',
            'eRepeatNum': 'sum',
            'MatDegree': 'mean',
            'eClass': 'first',
            'eState': 'first',
            'eReads': lambda x: ','.join(set(x))
        }).reset_index()

        grouped['eReadNum'] = grouped['eReads'].apply(lambda x: len(set(x.split(','))))
        grouped['MatDegree'] = grouped['MatDegree'].round(2)

        columns_order = ['eName', 'eChr', 'eStart', 'eEnd', 'eLength', 'eRepeatNum', 
                         'MatDegree', 'eClass', 'eState', 'eReadNum', 'eReads']
        return grouped[columns_order]

    def run_merge_mecc(self):
        logging.info("Reading CSV files...")
        df1 = self.read_csv(self.meccdna_part1)
        df2 = self.read_csv(self.tecc_analysis)
        
        logging.info("Processing MeccDNA_part1.csv...")
        df1_processed = self.process_meccdna_part1(df1)
        
        logging.info("Processing tecc_analysis_results.csv...")
        df2_processed = self.process_tecc_analysis(df2)
        
        logging.info("Merging DataFrames...")
        combined_df = self.process_dataframes(df1_processed, df2_processed)
        
        logging.info("Processing final DataFrame...")
        final_df = self.process_final_dataframe(combined_df)
        
        logging.info("Saving processed CSV...")
        final_df.to_csv(self.output_csv, index=False)
        
        logging.info("Merging FASTA files...")
        self.merge_fasta_files(self.mecc_part1, self.mecc_part2, "temp_merged.fa")
        
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
    parser = argparse.ArgumentParser(description="Process Mecc DNA data, generate output files, and remove duplicates.")
    parser.add_argument("--meccdna_part1", required=True, help="Input MeccDNA_part1.csv file")
    parser.add_argument("--tecc_analysis", required=True, help="Input tecc_analysis_results.csv file")
    parser.add_argument("--mecc_part1", required=True, help="Input mecc_part1.fa file")
    parser.add_argument("--mecc_part2", required=True, help="Input mecc_part2.fa file")
    parser.add_argument("--output_csv", required=True, help="Output CSV file")
    parser.add_argument("--output_fasta", required=True, help="Output FASTA file (deduplicated)")
    args = parser.parse_args()

    merger = MergeMecc(args.meccdna_part1, args.tecc_analysis, args.mecc_part1, args.mecc_part2, args.output_csv, args.output_fasta)
    merger.run_merge_mecc()

if __name__ == "__main__":
    main()