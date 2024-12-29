#!/usr/bin/env python3
# coding: utf-8

"""
MergeUecc: Unique eccDNA (UECC) Merger.

This module manages the merging of Unique eccDNA (UECC)
data from multiple sources. It provides functionality for combining and processing
UECC data while ensuring data integrity and proper sequence handling.

Key features:
- UECC data merging from multiple sources
- Sequence validation and deduplication
- FASTA file generation for merged sequences
- CSV output for detailed UECC information
- Ordered dictionary implementation for sequence tracking

Typical usage:
    merger = MergeUecc(uecc_part1, uecc_part2, output_csv, output_fasta)
    merger.run_merge_uecc()

Version: 1.0.1
License: GNU General Public License v2
Copyright (c) 2024 CircleSeeker Team
"""

import pandas as pd
import argparse
import logging
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os
from collections import OrderedDict

class MergeUecc:
    def __init__(self, uecc_part1, uecc_part2, output_csv, output_fasta):
        # Set logging configuration
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        self.logger = logging.getLogger(__name__)
        
        self.uecc_part1 = uecc_part1
        self.uecc_part2 = uecc_part2
        self.output_csv = output_csv
        self.output_fasta = output_fasta

    def read_csv(self, file_path):
        return pd.read_csv(file_path, dtype={'eStart': int, 'eEnd': int})

    def process_dataframes(self, df1, df2):
        # Process first dataframe (uecc_part1)
        df1_processed = df1[['eChr', 'eStart', 'eEnd', 'eStrand', 'eLength', 
                             'eRepeatNum', 'MatDegree', 'eClass', 'eReads', 'eSeq']]
        
        # Process second dataframe (uecc_part2)
        # Filter for Uecc class and select relevant columns
        df2_processed = df2[df2['eClass'] == 'Uecc'][['eChr', 'eStart', 'eEnd', 'eStrand', 
                                                       'eLength', 'eRepeatNum', 'MatDegree', 
                                                       'eClass', 'eReads', 'eSeq']]
        
        # Combine dataframes
        combined_df = pd.concat([df1_processed, df2_processed], ignore_index=True)
        
        # Generate eName
        combined_df['eName'] = combined_df.apply(
            lambda row: f"{row['eClass']}|{row['eChr']}-{row['eStart']}-{row['eEnd']}", 
            axis=1
        )
        
        combined_df['eState'] = 'Confirmed-eccDNA'
        return combined_df

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
            'eReads': lambda x: ','.join(set(x)),
            'eSeq': 'first'  # Keep the first sequence for each group
        }).reset_index()

        grouped['eReadNum'] = grouped['eReads'].apply(lambda x: len(set(x.split(','))))
        grouped['MatDegree'] = grouped['MatDegree'].round(2)

        columns_order = ['eName', 'eChr', 'eStart', 'eEnd', 'eStrand', 'eLength', 
                         'eRepeatNum', 'MatDegree', 'eClass', 'eState', 'eReadNum', 
                         'eReads', 'eSeq']
        return grouped[columns_order]

    # The following are new methods for final merge processing
    def sort_dataframe(self, df, sort_by='start'):
        if sort_by == 'start':
            return df.sort_values(by=['eChr', 'eStart', 'eEnd']).reset_index(drop=True)
        else:  # sort_by == 'end'
            return df.sort_values(by=['eChr', 'eEnd', 'eStart']).reset_index(drop=True)

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
            'eReads': lambda x: ','.join(set(','.join(x.dropna().astype(str).tolist()).split(','))),
            'eSeq': 'first'
        }
        
        merged_df = df_sorted.groupby(['eChr', 'group_id']).agg(aggregation_functions).reset_index(drop=True)
        
        # After merging, eReads may contain duplicate reads, need to deduplicate again
        merged_df['eReads'] = merged_df['eReads'].apply(lambda r: ','.join(sorted(set(r.split(','))))) 
        merged_df['eReadNum'] = merged_df['eReads'].apply(lambda x: len(set(x.split(','))))

        return merged_df

    def process_csv(self, df, sort_by):
        df_sorted = self.sort_dataframe(df, sort_by=sort_by)
        df_grouped = self.assign_group_ids(df_sorted, start_diff=50, end_diff=50)
        merged_df = self.merge_groups(df_grouped)
        return merged_df.drop(columns=['group_id'], errors='ignore')

    def create_fasta_from_df(self, df, output_file):
        # Create SeqRecord objects from DataFrame
        records = []
        for _, row in df.iterrows():
            seq_record = SeqRecord(
                Seq(row['eSeq']),
                id=row['eName'],
                description=""
            )
            records.append(seq_record)
        
        # Write to FASTA file
        SeqIO.write(records, output_file, "fasta")
        self.logger.info(f"Written {len(records)} sequences to FASTA file")

    def run_merge_uecc(self):
        self.logger.info("Reading CSV files...")
        df1 = self.read_csv(self.uecc_part1)
        df2 = self.read_csv(self.uecc_part2)
        
        self.logger.info("Processing DataFrames...")
        combined_df = self.process_dataframes(df1, df2)
        
        self.logger.info("Processing final DataFrame...")
        final_df = self.process_final_dataframe(combined_df)
        
        # Add merge rule processing
        self.logger.info("Processing CSV with start sorting and merging...")
        final_df = self.process_csv(final_df, 'start')
        
        self.logger.info("Processing CSV with end sorting and merging...")
        final_df = self.process_csv(final_df, 'end')

        self.logger.info("Saving processed CSV...")
        final_df.to_csv(self.output_csv, index=False)
        
        self.logger.info("Creating FASTA file from sequences...")
        self.create_fasta_from_df(final_df, self.output_fasta)
        
        self.logger.info("Process completed successfully.")

def main():
    parser = argparse.ArgumentParser(description="Process eccDNA data and generate output files.")
    parser.add_argument("--uecc_part1", required=True, help="Input Uecc_part1.csv file")
    parser.add_argument("--uecc_part2", required=True, help="Input uecc_part2.csv file")
    parser.add_argument("--output_csv", required=True, help="Output CSV file")
    parser.add_argument("--output_fasta", required=True, help="Output FASTA file")
    args = parser.parse_args()

    merger = MergeUecc(
        args.uecc_part1,
        args.uecc_part2,
        args.output_csv,
        args.output_fasta
    )
    merger.run_merge_uecc()

if __name__ == "__main__":
    main()