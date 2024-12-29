#!/usr/bin/env python3
# coding: utf-8

"""
MergeMecc: Multiple alignment eccDNA (MECC) Merger.

This module handles the merging of MECC data from multiple sources. It processes and combines MECC
data parts while maintaining sequence integrity and generating comprehensive outputs.

Key features:
- MECC data merging from multiple sources
- Sequence validation and processing
- FASTA file generation for merged sequences
- CSV output for detailed MECC information
- Robust error handling and logging

Typical usage:
    merger = MergeMecc(meccdna_part1, meccdna_part2, output_csv, output_fasta)
    merger.run_merge_mecc()

Version: 1.3.3
License: GNU General Public License v2
Copyright (c) 2024 CircleSeeker Team
"""

import argparse
import logging
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

class MergeMecc:
    def __init__(self, meccdna_part1, meccdna_part2, output_csv, output_fasta):
        # Set up logging
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        self.logger = logging.getLogger(__name__)
        
        self.meccdna_part1 = meccdna_part1
        self.meccdna_part2 = meccdna_part2
        self.output_csv = output_csv
        self.output_fasta = output_fasta

    def run_merge_mecc(self):
        """Execute the complete merge process"""
        # Read input files
        self.logger.info("Reading input files...")
        df1 = self.read_csv(self.meccdna_part1)
        df2 = self.read_csv(self.meccdna_part2)
        
        # Process data
        self.logger.info("Processing MeccDNA part1 data...")
        df1_processed = self.process_meccdna_part1(df1)
        
        self.logger.info("Processing MeccDNA part2 data...")
        df2_processed = self.process_meccdna_part2(df2)
        
        # Merge dataframes
        self.logger.info("Combining processed data...")
        combined_df = pd.concat([df1_processed, df2_processed], ignore_index=True)
        
        # Handle empty data cases
        if combined_df.empty:
            self.logger.warning("No valid data after combining. Creating empty output files.")
            pd.DataFrame(columns=['eName', 'MapNum', 'eRepeatNum', 'eState', 'eChr', 'eStart', 
                                'eEnd', 'eLength', 'MatDegree', 'eClass', 'eReadNum', 'eReads', 
                                'eSeq']).to_csv(self.output_csv, index=False)
            with open(self.output_fasta, 'w') as f:
                pass
            return
        
        if 'query_id' not in combined_df.columns:
            self.logger.error("Required column 'query_id' not found in processed data.")
            return
        
        # Process query groups
        self.logger.info("Processing query groups...")
        query_groups = self.process_query_groups(combined_df)
        
        self.logger.info("Merging similar feature groups...")
        merged_query_groups = self.merge_same_feature_groups(query_groups)
        
        self.logger.info("Integrating data...")
        integrated_df = self.integrate_data(merged_query_groups, combined_df)
        
        self.logger.info("Adding names and states...")
        final_df = self.add_name_and_state(integrated_df)
        
        # Prepare and save final results
        self.logger.info("Preparing final results...")
        final_results = self.prepare_final_dataframe(final_df)
        final_results.to_csv(self.output_csv, index=False)
        
        # Create FASTA file from final results
        self.logger.info("Creating FASTA file from sequences...")
        self.create_fasta_from_df(final_results, self.output_fasta)
        
        self.logger.info("Processing completed successfully.")

    def read_csv(self, file_path):
        """Read CSV file and handle empty file cases"""
        try:
            return pd.read_csv(file_path)
        except (FileNotFoundError, pd.errors.EmptyDataError):
            self.logger.warning(f"File {file_path} is empty or doesn't exist. Creating empty DataFrame.")
            return pd.DataFrame()

    def process_meccdna_part1(self, df):
        """Process MeccDNA part1 data"""
        if df.empty:
            self.logger.warning("MeccDNA part1 DataFrame is empty. Returning empty DataFrame.")
            return pd.DataFrame(columns=['eChr', 'eStart', 'eEnd', 'eLength', 'eRepeatNum', 
                                      'MatDegree', 'eClass', 'eReads', 'query_id', 'eSeq'])
        
        rows = []
        for _, row in df.iterrows():
            # Skip records with consLen <= 100
            if row['consLen'] <= 100:
                continue
            
            align_regions = row['AlignRegion'].split(';')
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
                    'query_id': row['query_id'],
                    'eSeq': row['eSeq']
                })
        return pd.DataFrame(rows)

    def process_meccdna_part2(self, df):
        """Process MeccDNA part2 data"""
        empty_df = pd.DataFrame(columns=['eChr', 'eStart', 'eEnd', 'eLength', 'eRepeatNum', 
                                      'MatDegree', 'eClass', 'eReads', 'query_id', 'eSeq'])
        
        if df.empty:
            self.logger.warning("MeccDNA part2 DataFrame is empty. Returning empty DataFrame.")
            return empty_df

        # Filter for Mecc class and length > 100
        mecc_df = df[(df['eClass'] == 'Mecc') & (df['eLength'] > 100)].copy()
        if mecc_df.empty:
            self.logger.warning("No valid Mecc records found in MeccDNA part2. Returning empty DataFrame.")
            return empty_df
        
        mecc_df['query_id'] = mecc_df['qname']
        
        return mecc_df[['eChr', 'eStart', 'eEnd', 'eLength', 'eRepeatNum', 'MatDegree', 
                       'eClass', 'eReads', 'query_id', 'eSeq']]

    def process_query_groups(self, df):
        """Process query groups and feature values"""
        df_processed = df.copy()
        df_processed['location_feature'] = df_processed.apply(
            lambda row: f"{row['eChr']}-{row['eStart']}-{row['eEnd']}", 
            axis=1
        )
        
        query_groups = df_processed.groupby('query_id').agg({
            'location_feature': list,
            'eRepeatNum': 'first',
            'eReads': 'first',
        }).reset_index()
        
        query_groups['feature_count'] = query_groups['location_feature'].apply(len)
        valid_groups = query_groups[query_groups['feature_count'] >= 1].copy()
        valid_groups['location_feature'] = valid_groups['location_feature'].apply(sorted)
        
        return valid_groups

    def merge_same_feature_groups(self, query_groups):
        """Merge query groups with the same feature values"""
        query_groups['feature_str'] = query_groups['location_feature'].apply(lambda x: ';'.join(sorted(x)))
        feature_groups = query_groups.groupby('feature_str').agg(list).reset_index()
        
        merged_groups = []
        for _, group in feature_groups.iterrows():
            merged_groups.append({
                'query_id': group['query_id'][0],
                'location_feature': query_groups[query_groups['query_id'] == group['query_id'][0]]['location_feature'].iloc[0],
                'eRepeatNum': sum(group['eRepeatNum']),
                'eReads': ';'.join(group['eReads'])
            })
        
        return pd.DataFrame(merged_groups)

    def integrate_data(self, merged_df, raw_combined_df):
        """Integrate merged_df and raw_combined_df"""
        # Create a mapping from query_id to merged info
        merged_info_map = merged_df.set_index('query_id')[['eRepeatNum', 'eReads']].to_dict('index')
        
        # Filter valid raw data
        valid_query_ids = set(merged_df['query_id'])
        filtered_raw = raw_combined_df[raw_combined_df['query_id'].isin(valid_query_ids)].copy()
        
        # Use vectorized operations to update data
        filtered_raw['eRepeatNum'] = filtered_raw['query_id'].map(
            lambda x: merged_info_map[x]['eRepeatNum'])
        filtered_raw['eReads'] = filtered_raw['query_id'].map(
            lambda x: merged_info_map[x]['eReads'])
        
        # Select needed columns
        result = filtered_raw[['eChr', 'eStart', 'eEnd', 'eLength', 'MatDegree', 
                             'eClass', 'eRepeatNum', 'eReads', 'query_id', 'eSeq']]
        
        return result

    def add_name_and_state(self, df):
        """Add eName and eState"""
        unique_query_ids = sorted(df['query_id'].unique())
        query_id_to_num = {qid: idx+1 for idx, qid in enumerate(unique_query_ids)}
        
        df['eName'] = df['query_id'].map(lambda x: f"Mecc|Confirmed|{query_id_to_num[x]}")
        df['eState'] = 'Confirmed-eccDNA'
        
        df['eName_num'] = df['eName'].str.extract(r'(\d+)').astype(int)
        df = df.sort_values('eName_num')
        df = df.drop('eName_num', axis=1)
        
        return df

    def create_fasta_from_df(self, df, output_file):
        """Create FASTA file from DataFrame and ensure sequences are not duplicated"""
        # Group by eName and take the first sequence record
        unique_df = df.drop_duplicates(subset=['eName'], keep='first')
        
        records = []
        for _, row in unique_df.iterrows():
            if pd.notna(row['eSeq']):  # Ensure sequence is not NaN
                seq_record = SeqRecord(
                    Seq(row['eSeq']),
                    id=row['eName'],
                    description=""
                )
                records.append(seq_record)
        
        # Write to FASTA file
        SeqIO.write(records, output_file, "fasta")
        self.logger.info(f"Written {len(records)} unique sequences to FASTA file")
        if len(records) < len(df):
            self.logger.info(f"Removed {len(df) - len(records)} duplicate sequences")

    def prepare_final_dataframe(self, df):
        """Prepare final dataframe, add MapNum column"""
        final_df = df.copy()
        
        # Calculate eReadNum
        final_df['eReadNum'] = final_df['eReads'].apply(lambda x: len(str(x).split(';')))
        
        # Calculate MapNum (number of rows per query_id)
        map_num_series = df.groupby('query_id').size()
        final_df['MapNum'] = final_df['query_id'].map(map_num_series)
        
        # Set column order
        columns_order = [
            'eName', 'MapNum', 'eRepeatNum', 'eState', 'eChr', 'eStart', 'eEnd',
            'eLength', 'MatDegree', 'eClass', 'eReadNum', 'eReads', 'eSeq'
        ]
        
        # Sort
        final_df['sort_key'] = final_df['eName'].str.extract(r'(\d+)').astype(int)
        final_df = final_df.sort_values('sort_key')
        final_df = final_df.drop('sort_key', axis=1)
        
        # Select final columns
        return final_df[columns_order]

def main():
    parser = argparse.ArgumentParser(description="Process Mecc DNA data and generate output files.")
    parser.add_argument("--meccdna_part1", required=True, help="Input MeccDNA_part1.csv file")
    parser.add_argument("--meccdna_part2", required=True, help="Input MeccDNA_part2.csv file")
    parser.add_argument("--output_csv", required=True, help="Output CSV file")
    parser.add_argument("--output_fasta", required=True, help="Output FASTA file")
    args = parser.parse_args()

    merger = MergeMecc(
        args.meccdna_part1,
        args.meccdna_part2,
        args.output_csv,
        args.output_fasta
    )
    merger.run_merge_mecc()

if __name__ == "__main__":
    main()