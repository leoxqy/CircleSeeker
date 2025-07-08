#!/usr/bin/env python3
# coding: utf-8

"""
MergeUeccInf: Unique eccDNA (UECC) Inference Analysis and Merging Module.

This module provides functionality for analyzing and merging inferred Unique
eccDNA (UECC) data. It processes both confirmed and
inferred UECC instances, generating comprehensive analysis outputs.

Key features:
- Analysis of inferred UECC instances
- Merging of confirmed and inferred UECCs
- Shared count analysis and tracking
- FASTA output generation
- Temporary file management
- Comprehensive result reporting

Typical usage:
    analyzer = UeccAnalyzer(inf_process_csv, confirmed_csv, inferred_output,
                            shared_count_output, fasta_output, merged_output)
    analyzer.run()

Version: 1.0.1
License: GNU General Public License v2
Copyright (c) 2024 CircleSeeker Team
"""

import sys
import csv

# ========== Increase CSV field size limit at script start ========== #
try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    max_int = sys.maxsize
    while True:
        try:
            csv.field_size_limit(max_int)
            break
        except OverflowError:
            max_int = int(max_int / 10)

# Import other modules after setting field size limit
import pandas as pd
import logging
import argparse
import os

class UeccAnalyzer:
    def __init__(
        self,
        inf_process_csv,
        confirmed_csv,
        inferred_output,
        shared_count_output,
        fasta_output,
        merged_output,
        keep_tmp=False
    ):
        self.inf_process_csv = inf_process_csv
        self.confirmed_csv = confirmed_csv
        self.inferred_output = inferred_output
        self.shared_count_output = shared_count_output
        self.fasta_output = fasta_output
        self.merged_output = merged_output
        self.keep_tmp = keep_tmp
        
        # Set log level based on keep_tmp parameter
        log_level = logging.DEBUG if self.keep_tmp else logging.INFO
        logging.basicConfig(
            level=log_level,
            format='%(asctime)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        self.logger = logging.getLogger(__name__)
        
        self.logger.info(
            f"Initialized UeccAnalyzer with input files: "
            f"{inf_process_csv}, {confirmed_csv}"
        )

    def process_inferred(self):
        """
        1. Read and process inferred eccDNA data (inf_process_csv, column names already modified).
        2. Filter entries with length >= 100.
        3. Normalize and rename to eChr, eStart, eEnd, eLength, etc., generate eName, eClass, eState.
        4. Output the filtered CSV.
        """
        self.logger.info("Processing Inferred Uecc data...")
        df = pd.read_csv(self.inf_process_csv)

        # Column names in inf_process_csv
        #   chrom, start, end, length, num_reads, high_cov, weighted_cov,
        #   split_read_count, split_read_ratio, Mean_MAPQ, High_MAPQ_ratio,
        #   score, Avg_Read_Length, Avg_Split_Read_Length, reads_list, eSeq
        columns_to_keep = [
            'chrom',
            'start',
            'end',
            'length',
            'num_reads',
            'Avg_Read_Length',
            'split_read_count',
            'Avg_Split_Read_Length',
            'score',
            'eSeq'
        ]
        
        # Keep only required columns
        df_processed = df[columns_to_keep].copy()

        # Rename to standardized fields
        df_processed.rename(
            columns={
                'chrom': 'eChr',
                'start': 'eStart',
                'end': 'eEnd',
                'length': 'eLength',
                'num_reads': 'Num_Reads',
                'split_read_count': 'Num_Split_Reads',
                'Avg_Split_Read_Length': 'Avg_Split_Read_Length',
                'score': 'Score'
            },
            inplace=True
        )

        # Generate eName
        df_processed['eName'] = df_processed.apply(
            lambda row: f"Uecc|{row['eChr']}-{row['eStart']}-{row['eEnd']}", axis=1
        )
        # Fixed columns
        df_processed['eClass'] = 'Uecc'
        df_processed['eState'] = 'Inferred-eccDNA'

        # Filter: length >= 100
        df_filtered = df_processed[df_processed['eLength'] >= 100].copy()

        # Required column order
        final_columns = [
            'eName',
            'eChr',
            'eStart',
            'eEnd',
            'eLength',
            'Num_Reads',
            'Avg_Read_Length',
            'Num_Split_Reads',
            'Avg_Split_Read_Length',
            'Score',
            'eSeq',
            'eClass',
            'eState'
        ]
        df_filtered = df_filtered[final_columns]

        df_filtered.to_csv(self.inferred_output, index=False)
        self.logger.info(
            f"Inferred Uecc processing complete. Output: {self.inferred_output}"
        )
        return self.inferred_output

    def compare_and_count(self):
        """
        1. Compare process_inferred output (inferred_output) with original column names from confirmed results (confirmed_csv).
        2. Count inferred records that match within ±30bp range.
        3. Output to shared_count_output.
        """
        self.logger.info("Comparing Inferred and Confirmed Uecc data...")

        # Read inferred results (already renamed)
        df_inferred = pd.read_csv(self.inferred_output)
        
        # Read confirmed results (not yet renamed), assuming minimum columns: Chromosome,Start,End,Length
        df_confirmed = pd.read_csv(self.confirmed_csv)

        # Only rename necessary columns
        df_confirmed.rename(
            columns={
                'Chromosome': 'eChr',
                'Start': 'eStart',
                'End': 'eEnd',
                'Length': 'eLength'
            },
            inplace=True
        )
        # If the confirmed file contains other columns (Score, eSeq, etc.), add them similarly:
        # df_confirmed.rename(columns={'Score': 'Score', 'eSeq': 'eSeq'}, inplace=True)
        # ...

        # For merge_and_output convenience, also generate eName in confirmed table
        if 'eName' not in df_confirmed.columns:
            df_confirmed['eName'] = df_confirmed.apply(
                lambda row: f"Uecc|{row['eChr']}-{row['eStart']}-{row['eEnd']}", axis=1
            )

        # Group by chromosome to speed up comparison
        inferred_by_chr = {}
        for chr_name, subdf in df_inferred.groupby('eChr'):
            inferred_by_chr[chr_name] = subdf

        matched_inferred_set = set()

        for _, row in df_confirmed.iterrows():
            c_chr = row['eChr']
            c_start = row['eStart']
            c_end = row['eEnd']

            if c_chr in inferred_by_chr:
                sub_inferred = inferred_by_chr[c_chr]
                matched = sub_inferred[
                    (sub_inferred['eStart'].sub(c_start).abs() <= 30) &
                    (sub_inferred['eEnd'].sub(c_end).abs() <= 30)
                ]
                if not matched.empty:
                    # Add matched inferred eName to set to ensure each inferred is counted only once
                    matched_inferred_set.update(matched['eName'])

        count = len(matched_inferred_set)
        pd.DataFrame({'SharedCount': [count]}).to_csv(self.shared_count_output, index=False)
        self.logger.info(f"Shared count: {count}. Output: {self.shared_count_output}")
        return count

    def generate_fasta(self):
        """
        Generate FASTA from inferred results (inferred_output).
        """
        self.logger.info("Generating FASTA from Inferred Uecc data...")
        df_inferred = pd.read_csv(self.inferred_output)

        with open(self.fasta_output, 'w') as f:
            for _, row in df_inferred.iterrows():
                eName = row['eName']
                eSeq = row['eSeq']
                f.write(f">{eName}\n{eSeq}\n")

        self.logger.info(f"FASTA generated: {self.fasta_output}")
        return self.fasta_output

    def merge_and_output(self):
        """
        1. Merge (already renamed) inferred_output with (renamed) confirmed_csv.
        2. Use eName as index, if confirmed has non-empty values, overwrite inferred empty values.
        3. Output to merged_output.
        """
        self.logger.info("Merging Inferred and Confirmed Uecc data for final output...")

        df_inferred = pd.read_csv(self.inferred_output)
        df_confirmed = pd.read_csv(self.confirmed_csv)
        
        # Rename confirmed table again to ensure field alignment
        df_confirmed.rename(
            columns={
                'Chromosome': 'eChr',
                'Start': 'eStart',
                'End': 'eEnd',
                'Length': 'eLength'
            },
            inplace=True
        )
        # If there are other columns, add them here
        # ...

        # Ensure confirmed table has eName
        if 'eName' not in df_confirmed.columns:
            df_confirmed['eName'] = df_confirmed.apply(
                lambda row: f"Uecc|{row['eChr']}-{row['eStart']}-{row['eEnd']}", axis=1
            )

        df_inferred.set_index('eName', inplace=True)
        df_confirmed.set_index('eName', inplace=True)

        # Merge indices
        all_names = df_inferred.index.union(df_confirmed.index)
        merged_df = pd.DataFrame(index=all_names)

        # Columns to merge (adjust according to needs)
        columns_to_merge = [
            'eChr', 'eStart', 'eEnd', 'eLength', 
            'Num_Reads', 'Avg_Read_Length', 
            'Num_Split_Reads', 'Avg_Split_Read_Length',
            'Score', 'eSeq', 'eState'
        ]

        for col in columns_to_merge:
            inferred_vals = df_inferred[col] if col in df_inferred.columns else None
            confirmed_vals = df_confirmed[col] if col in df_confirmed.columns else None

            if inferred_vals is not None and confirmed_vals is not None:
                # First put inferred, then use confirmed non-empty values to overwrite
                merged_df[col] = inferred_vals.reindex(all_names)
                confirmed_nonnull = confirmed_vals.dropna()
                merged_df.loc[confirmed_nonnull.index, col] = confirmed_nonnull
            elif confirmed_vals is not None:
                merged_df[col] = confirmed_vals.reindex(all_names)
            else:
                merged_df[col] = inferred_vals.reindex(all_names)

        # Restore eName
        merged_df.reset_index(inplace=True)
        merged_df.rename(columns={'index': 'eName'}, inplace=True)

        # Data type conversion & fill empty values
        for col in ['eStart', 'eEnd', 'eLength', 'Num_Reads', 'Num_Split_Reads']:
            if col in merged_df.columns:
                merged_df[col] = merged_df[col].fillna(0).astype(int)

        if 'eState' in merged_df.columns:
            merged_df['eState'] = merged_df['eState'].fillna('Inferred-eccDNA')
        if 'eSeq' in merged_df.columns:
            merged_df['eSeq'] = merged_df['eSeq'].fillna('')

        # Output required columns
        final_cols = ['eName'] + columns_to_merge
        final_cols = [c for c in final_cols if c in merged_df.columns]
        merged_df = merged_df[final_cols]

        merged_df.to_csv(self.merged_output, index=False)
        self.logger.info(f"Merged output generated: {self.merged_output}")

    def run(self):
        # Execute four main functions in order
        self.process_inferred()
        self.compare_and_count()
        self.generate_fasta()
        self.merge_and_output()
        self.logger.info("All tasks completed.")


def main():
    parser = argparse.ArgumentParser(
        description="Analyze Uecc DNA data with merges, logging included in class."
    )
    parser.add_argument(
        "--inf_process_csv",
        required=True,
        help="Path to the input CSV (e.g. demo_Final.Inf.process.Uecc.csv) "
             "with columns: chrom,start,end,length,num_reads,high_cov,weighted_cov,"
             "split_read_count,split_read_ratio,Mean_MAPQ,High_MAPQ_ratio,score,"
             "Avg_Read_Length,Avg_Split_Read_Length,reads_list,eSeq"
    )
    parser.add_argument(
        "--confirmed_csv",
        required=True,
        help="Path to the confirmed Uecc CSV (e.g. demo.Final.Confirmed.Uecc.csv), "
             "with original columns such as Chromosome,Start,End,Length, etc."
    )
    parser.add_argument(
        "--inferred_output",
        required=True,
        help="Path to output the processed (filtered & renamed) inferred result CSV."
    )
    parser.add_argument(
        "--shared_count_output",
        required=True,
        help="Path to output the CSV which contains a single row 'SharedCount' for the overlap count."
    )
    parser.add_argument(
        "--fasta_output",
        required=True,
        help="Path to output the FASTA file generated from the inferred result."
    )
    parser.add_argument(
        "--merged_output",
        required=True,
        help="Path to output the final merged CSV."
    )
    parser.add_argument(
        "--keep-tmp",
        action="store_true",
        help="Keep temporary files and show debug logs"
    )

    args = parser.parse_args()

    analyzer = UeccAnalyzer(
        inf_process_csv=args.inf_process_csv,
        confirmed_csv=args.confirmed_csv,
        inferred_output=args.inferred_output,
        shared_count_output=args.shared_count_output,
        fasta_output=args.fasta_output,
        merged_output=args.merged_output,
        keep_tmp=args.keep_tmp
    )
    analyzer.run()

if __name__ == "__main__":
    main()