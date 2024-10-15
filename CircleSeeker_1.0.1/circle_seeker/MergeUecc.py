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
        return pd.read_csv(file_path, dtype={'eStart': int, 'eEnd': int})

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
        seen_ids = OrderedDict()
        duplicates = 0
        for record in SeqIO.parse(input_file, "fasta"):
            if record.id not in seen_ids:
                seen_ids[record.id] = record
            else:
                duplicates += 1

        with open(output_file, "w") as output_handle:
            SeqIO.write(seen_ids.values(), output_handle, "fasta")

        logging.info(f"Removed {duplicates} duplicate sequences by ID")
        return len(seen_ids)

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

    def sort_dataframe(self, df, sort_by='start'):
        if sort_by == 'start':
            return df.sort_values(by=['eStart', 'eChr', 'eEnd']).reset_index(drop=True)
        else:  # sort_by == 'end'
            return df.sort_values(by=['eEnd', 'eChr', 'eStart']).reset_index(drop=True)

    def assign_group_ids(self, df, start_diff=50, end_diff=50):
        df['prev_eStart'] = df.groupby('eChr')['eStart'].shift(1)
        df['prev_eEnd'] = df.groupby('eChr')['eEnd'].shift(1)
        
        df['eStart_diff'] = (df['eStart'] - df['prev_eStart']).abs()
        df['eEnd_diff'] = (df['eEnd'] - df['prev_eEnd']).abs()
        
        df['new_group'] = (df['eStart_diff'] > start_diff) | (df['eEnd_diff'] > end_diff) | (df['eStart_diff'].isnull())
        
        df['group_id'] = df.groupby('eChr')['new_group'].cumsum()
        
        df = df.drop(columns=['prev_eStart', 'prev_eEnd', 'eStart_diff', 'eEnd_diff', 'new_group'])
        
        return df

    def merge_groups(self, df):
        df_sorted = df.sort_values(['eChr', 'group_id', 'MatDegree', 'eLength'], ascending=[True, True, False, False])
        
        aggregation_functions = {
            'eName': 'first',
            'eChr': 'first',
            'eStart': 'first',
            'eEnd': 'first',
            'eStrand': 'first',
            'eLength': 'first',
            'eRepeatNum': 'sum',
            'MatDegree': 'max',
            'eClass': 'first',
            'eState': 'first',
            'eReadNum': 'sum',
            'eReads': lambda x: ','.join(x.dropna().astype(str).tolist())
        }
        
        merged_df = df_sorted.groupby(['eChr', 'group_id']).agg(aggregation_functions).reset_index(drop=True)
        
        return merged_df

    def process_csv(self, df, sort_by):
        df_sorted = self.sort_dataframe(df, sort_by=sort_by)
        df_grouped = self.assign_group_ids(df_sorted, start_diff=50, end_diff=50)
        merged_df = self.merge_groups(df_grouped)
        return merged_df.drop(columns=['group_id'], errors='ignore')

    def filter_fasta_by_csv(self, csv_file, input_fasta, output_fasta):
        # Read the CSV file to get the list of valid IDs
        df = pd.read_csv(csv_file)
        valid_ids = set(df['eName'])

        # Filter the FASTA file
        filtered_records = []
        for record in SeqIO.parse(input_fasta, "fasta"):
            if record.id in valid_ids:
                filtered_records.append(record)

        # Write the filtered sequences to the output file
        SeqIO.write(filtered_records, output_fasta, "fasta")
        
        logging.info(f"Filtered FASTA file contains {len(filtered_records)} sequences")
        return len(filtered_records)

    def run_merge_uecc(self):
        logging.info("Reading CSV files...")
        df1 = self.read_csv(self.ueccdna_part1)
        df2 = self.read_csv(self.tecc_analysis)
        
        logging.info("Processing DataFrames...")
        combined_df = self.process_dataframes(df1, df2)
        
        logging.info("Processing final DataFrame...")
        final_df = self.process_final_dataframe(combined_df)
        
        logging.info("Processing CSV with start sorting...")
        final_df = self.process_csv(final_df, 'start')
        
        logging.info("Processing CSV with end sorting...")
        final_df = self.process_csv(final_df, 'end')
        
        logging.info("Saving processed CSV...")
        final_df.to_csv(self.output_csv, index=False)
        
        logging.info("Merging FASTA files...")
        self.merge_fasta_files(self.uecc_part1, self.uecc_part2, "temp_merged.fa")
        
        logging.info("Updating FASTA sequence names...")
        name_map = dict(zip(combined_df['query_id'], combined_df['eName']))
        self.update_fasta_names("temp_merged.fa", name_map, "temp_renamed.fa")
        
        logging.info("Removing duplicate sequences...")
        unique_count = self.remove_duplicates("temp_renamed.fa", "temp_deduplicated.fa")
        logging.info(f"Deduplicated FASTA file contains {unique_count} unique sequences")

        logging.info("Filtering FASTA sequences based on CSV entries...")
        final_count = self.filter_fasta_by_csv(self.output_csv, "temp_deduplicated.fa", self.output_fasta)
        logging.info(f"Final filtered FASTA file contains {final_count} sequences")
        
        # Clean up temporary files
        os.remove("temp_merged.fa")
        os.remove("temp_renamed.fa")
        os.remove("temp_deduplicated.fa")
        
        logging.info("Process completed successfully.")

def main():
    parser = argparse.ArgumentParser(description="Process eccDNA data, generate output files, and remove duplicates.")
    parser.add_argument("--ueccdna_part1", required=True, help="Input UeccDNA_part1.csv file")
    parser.add_argument("--tecc_analysis", required=True, help="Input tecc_analysis_results.csv file")
    parser.add_argument("--uecc_part1", required=True, help="Input uecc_part1.fa file")
    parser.add_argument("--uecc_part2", required=True, help="Input uecc_part2.fa file")
    parser.add_argument("--output_csv", required=True, help="Output CSV file")
    parser.add_argument("--output_fasta", required=True, help="Output FASTA file (filtered)")
    args = parser.parse_args()

    merger = MergeUecc(args.ueccdna_part1, args.tecc_analysis, args.uecc_part1, args.uecc_part2, args.output_csv, args.output_fasta)
    merger.run_merge_uecc()

if __name__ == "__main__":
    main()
