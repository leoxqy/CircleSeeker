#!/usr/bin/env python
# coding: utf-8

"""
Ringmaster: Circular DNA Classification and Analysis Orchestrator.

This module serves as the central coordinator for circular DNA classification and
analysis. It processes BLAST results to identify and categorize different types
of circular DNA structures (UECC, MECC, CECC) with parallel processing capabilities.

Key features:
- Multi-threaded BLAST result processing
- Classification of multiple circular DNA types
- Progress tracking with tqdm
- Comprehensive logging system
- Flexible output directory management

Typical usage:
    ringmaster = Ringmaster(blast_results_file, circular_seq_fasta,
                          Uecc_output_csv, Mecc_output_csv, Cecc_output_csv)
    ringmaster.run()

Version: 1.3.3
License: GNU General Public License v2
Copyright (c) 2024 CircleSeeker Team
"""

import pandas as pd
import numpy as np
import argparse
from Bio import SeqIO
import os
import time
from multiprocessing import Pool
from tqdm import tqdm
import logging

class Ringmaster:
    def __init__(self, blast_results_file, circular_seq_fasta, Uecc_output_csv, Mecc_output_csv, Cecc_output_csv, process_xecc=False, num_threads=1, prefix=None, output_dir=None):
        # Configure logging
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        self.logger = logging.getLogger(__name__)

        self.blast_results_file = blast_results_file
        self.circular_seq_fasta = circular_seq_fasta
        self.process_xecc = process_xecc
        self.num_threads = num_threads

        # Set prefix and output directory
        self.prefix = prefix if prefix else ""
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            self.output_dir = output_dir
        else:
            self.output_dir = os.getcwd()

        # Update output file paths to the specified directory (if provided)
        self.Uecc_output_csv = os.path.join(self.output_dir, os.path.basename(Uecc_output_csv))
        self.Mecc_output_csv = os.path.join(self.output_dir, os.path.basename(Mecc_output_csv))
        self.Cecc_output_csv = os.path.join(self.output_dir, os.path.basename(Cecc_output_csv))

        # Xecc output file
        xecc_filename = f"{self.prefix}_XeccDNA.fa" if self.prefix else "XeccDNA.fa"
        self.xecc_output_fa = os.path.join(self.output_dir, xecc_filename)

    def read_blast_results(self):
        columns = ['query_id', 'subject_id', 'identity', 'alignment_length', 'mismatches',
                'gap_opens', 'q_start', 'q_end', 's_start', 's_end', 'evalue', 'bit_score']
        df = pd.read_csv(self.blast_results_file, sep='\t', header=None, names=columns)

        df['strand'] = np.where(df['s_start'] < df['s_end'], '+', '-')
        df.loc[df['strand'] == '-', ['s_start', 's_end']] = df.loc[df['strand'] == '-', ['s_end', 's_start']].values

        return df

    def process_query_id(self, df):
        split_cols = df['query_id'].str.split('|', expand=True)
        df['readName'] = split_cols[0]
        df['repN'] = split_cols[1]
        df['consLen'] = split_cols[2].astype(int)
        df['copyNum'] = split_cols[3].astype(float)
        return df

    def add_calculated_columns(self, df):
        df['Rlength'] = df['s_end'] - df['s_start'] + 1
        df['gap_Length'] = df['consLen'] - df['Rlength']
        df['Gap_Percentage'] = ((df['gap_Length'].abs() / df['consLen']) * 100).round(2)
        return df

    def filter_low_gap_percentage(self,df):
        return df[(df['Gap_Percentage'] <= 20) | (df['gap_Length'].abs() <= 100)].copy()

    def process_group(self, group):
        if len(group) == 1:
            return 'Uecc', group.index[0]

        start_diffs = np.abs(group['s_start'].values[:, None] - group['s_start'].values)
        end_diffs = np.abs(group['s_end'].values[:, None] - group['s_end'].values)

        close_pairs = (start_diffs <= 5) | (end_diffs <= 5)
        np.fill_diagonal(close_pairs, False)

        if close_pairs.any():
            min_gap_index = group['Gap_Percentage'].idxmin()
            return 'Uecc', min_gap_index
        else:
            return 'Mecc', None

    def process_Uecc(self, df):
        Uecc_df = df[df['eClass'] == 'Uecc'].copy()
        columns = ['query_id', 'subject_id', 's_start', 's_end', 'strand', 'consLen', 'copyNum', 'Gap_Percentage', 'eClass', 'gap_Length', 'readName', 'q_start', 'q_end']
        Uecc_df = Uecc_df[columns].copy()

        column_mapping = {
            'subject_id': 'eChr',
            'consLen': 'eLength',
            'copyNum': 'eRepeatNum',
            'strand': 'eStrand',
            'readName': 'eReads'
        }
        Uecc_df.rename(columns=column_mapping, inplace=True)

        Uecc_df['eStart'] = Uecc_df['s_start'] - Uecc_df['gap_Length']
        Uecc_df['eEnd'] = Uecc_df['s_end']

        mask = Uecc_df['eStart'] < 1
        Uecc_df.loc[mask, 'eStart'] = Uecc_df.loc[mask, 's_start']
        Uecc_df.loc[mask, 'eEnd'] = Uecc_df.loc[mask, 's_end'] + Uecc_df.loc[mask, 'gap_Length']

        Uecc_df['MatDegree'] = 100 - Uecc_df['Gap_Percentage']
        Uecc_df['eStart'] = Uecc_df['eStart'].astype(int)
        Uecc_df['eEnd'] = Uecc_df['eEnd'].astype(int)
        Uecc_df['eName'] = Uecc_df['eChr'] + "-" + Uecc_df['eStart'].astype(str) + "-" + Uecc_df['eEnd'].astype(str)

        result_columns = ['eName', 'eChr', 'eStart', 'eEnd', 'eStrand', 'eLength', 'eRepeatNum', 'MatDegree', 'eClass', 'eReads', 'query_id', 'q_start', 'q_end']
        return Uecc_df[result_columns].copy()

    def generate_align_region(self, row):
        if row['s_start'] - row['gap_Length'] < 1:
            return f"{row['identity']}|{row['subject_id']}-{row['s_start']}-{row['s_end'] + row['gap_Length']}"
        else:
            return f"{row['identity']}|{row['subject_id']}-{row['s_start'] - row['gap_Length']}-{row['s_end']}"

    def process_mecc_group(self, group):
        q_start_mode = group['q_start'].mode().iloc[0]
        q_start_row = group[group['q_start'] == q_start_mode].iloc[0]

        q_end_mode = group['q_end'].mode().iloc[0]
        q_end_row = group[group['q_end'] == q_end_mode].iloc[0]

        align_regions = group.apply(self.generate_align_region, axis=1).tolist()

        return pd.Series({
            'query_id': group['query_id'].iloc[0],
            'q_start': q_start_row['q_start'],
            'q_end': q_end_row['q_end'],
            'copyNum': group['copyNum'].iloc[0],
            'consLen': group['consLen'].iloc[0],
            'readName': group['readName'].iloc[0],
            'eClass': group['eClass'].iloc[0],
            'AlignRegion': ';'.join(align_regions),
            'MatchNum': len(group)
        })

    def process_mecc(self, df):
        mecc_df = df[df['eClass'] == 'Mecc'].copy()
        result_df = mecc_df.groupby('query_id').apply(self.process_mecc_group).reset_index(drop=True)

        result_df['gap_Length_2'] = result_df['consLen'] - (result_df['q_end'] - result_df['q_start'] + 1)
        result_df['tStart'] = result_df['q_start'] - result_df['gap_Length_2']
        result_df['tEnd'] = result_df['q_end']

        mask = result_df['tStart'] < 1
        result_df.loc[mask, 'tStart'] = result_df.loc[mask, 'q_start']
        result_df.loc[mask, 'tEnd'] = result_df.loc[mask, 'q_end'] + result_df.loc[mask, 'gap_Length_2']

        result_df['tStart'] = result_df['tStart'].astype(int)
        result_df['tEnd'] = result_df['tEnd'].astype(int)

        return result_df

    def process_other_blast_rows(self, df):
        other_df = df[~df['eClass'].isin(['Uecc', 'Mecc'])].copy()
        group_sizes = other_df.groupby('query_id').size()
        query_ids_to_remove = group_sizes[group_sizes > 100].index
        other_df = other_df[~other_df['query_id'].isin(query_ids_to_remove)]
        return other_df

    def calculate_score(self, row1, row2):
        q_diff = min(abs(row1['q_start'] - row2['q_start']), abs(row1['q_end'] - row2['q_end']))
        subject_score = 100 if row1['subject_id'] == row2['subject_id'] else 0
        identity_score = min(row1['identity'], row2['identity'])
        s_diff = min(abs(row1['s_start'] - row2['s_start']), abs(row1['s_end'] - row2['s_end']))
        s_diff_score = 1 / (s_diff + 1)
        return q_diff * 0.4 + subject_score * 0.3 + identity_score * 0.2 + s_diff_score * 0.1

    def analyze_group(self, group):
        group = group.sort_values('q_start')
        Prelists = []
        Prelist_summaries = []
        processed = set()

        query_id = group['query_id'].iloc[0]
        cons_len = group['consLen'].iloc[0]

        for i, row in group.iterrows():
            if i in processed:
                continue

            current_chain = [row]
            processed.add(i)

            while True:
                potential_next = group[(group['q_start'] >= current_chain[-1]['q_end'] - 10) &
                                    (group['q_start'] <= current_chain[-1]['q_end'] + 10) &
                                    (~group.index.isin(processed))]

                if potential_next.empty:
                    break

                if len(potential_next) > 1:
                    scores = potential_next.apply(lambda x: self.calculate_score(current_chain[-1], x), axis=1)
                    next_row = potential_next.loc[scores.idxmax()]
                else:
                    next_row = potential_next.iloc[0]

                current_chain.append(next_row)
                processed.add(next_row.name)

            Prelist = pd.DataFrame(current_chain)
            start_q_start = Prelist['q_start'].min()
            end_q_end = Prelist['q_end'].max()
            pre_length = end_q_end - start_q_start + 1

            if pre_length > cons_len and len(Prelist) > 1:
                Prelist = pd.concat([Prelist.iloc[1:], Prelist.iloc[:1]]).reset_index(drop=True)
                Prelist['Prelist_ID'] = len(Prelists) + 1
                Prelists.append(Prelist)

                summary = {
                    'query_id': query_id,
                    'Prelist_ID': len(Prelists),
                    'start_q_start': start_q_start,
                    'end_q_end': end_q_end,
                    'PreLength': pre_length,
                    'consLen': cons_len,
                    'num_rows': len(Prelist)
                }
                Prelist_summaries.append(summary)

        if Prelists:
            return pd.concat(Prelists, ignore_index=True), pd.DataFrame(Prelist_summaries)
        else:
            return pd.DataFrame(), pd.DataFrame()

    def process_detailed_results(self, df):
        grouped = df.groupby(['query_id', 'Prelist_ID'])
        results = []

        for (query_id, prelist_id), group in grouped:
            group['Repeat1'] = group.duplicated(subset=['s_start'], keep=False)
            group['Repeat2'] = group.duplicated(subset=['s_end'], keep=False)
            repeated_rows = group[(group['Repeat1'] | group['Repeat2'])].copy()

            if not repeated_rows.empty:
                repeated_rows.loc[:, 'LocNum'] = range(1, len(repeated_rows) + 1)
                repeated_rows.loc[:, 's_length'] = repeated_rows['s_end'] - repeated_rows['s_start'] + 1

                selected_columns = ['Prelist_ID', 'query_id', 'q_start', 'q_end', 'consLen', 'LocNum',
                                  'identity', 'alignment_length', 'subject_id', 's_start', 's_end',
                                  's_length', 'strand', 'readName', 'repN', 'copyNum']

                results.append(repeated_rows[selected_columns])

        if results:
            return pd.concat(results, ignore_index=True)
        else:
            return pd.DataFrame()

    def process_detailed_results_further(self, df):
        grouped = df.groupby(['query_id', 'Prelist_ID'])
        results = []

        for (query_id, prelist_id), group in grouped:
            if len(group) <= 2:
                continue

            group = group.iloc[:-2]
            group = group.drop_duplicates(subset=['s_start', 's_end'])
            total_s_length = group['s_length'].sum()
            cons_len = group['consLen'].iloc[0]
            ratio = total_s_length / cons_len

            if 0.8 <= ratio <= 1.2:
                group['ratio'] = ratio
                results.append(group)

        if results:
            final_df = pd.concat(results, ignore_index=True)
            prelist_counts = final_df.groupby('query_id')['Prelist_ID'].nunique()
            final_df['eClass'] = final_df.apply(lambda row: 'Cecc' if prelist_counts[row['query_id']] == 1 else 'Ceccm', axis=1)
            return final_df
        else:
            return pd.DataFrame()

    def process_sequences(self, circular_seq_fasta, Uecc_df, mecc_df, Cecc_df):
        Uecc_df['eSeq'] = ''
        mecc_df['eSeq'] = ''
        Cecc_df['eSeq'] = ''

        Uecc_info = dict(zip(Uecc_df['query_id'], zip(Uecc_df.index, Uecc_df['q_start'], Uecc_df['q_end'])))
        mecc_info = dict(zip(mecc_df['query_id'], zip(mecc_df.index, mecc_df['q_start'], mecc_df['q_end'])))
        Cecc_info = dict(zip(Cecc_df['query_id'], zip(Cecc_df.index, Cecc_df['q_start'], Cecc_df['q_end'])))

        for record in SeqIO.parse(circular_seq_fasta, "fasta"):
            if record.id in Uecc_info:
                idx, start, end = Uecc_info[record.id]
                sequence = str(record.seq[start-1:end])
                Uecc_df.at[idx, 'eSeq'] = sequence
            elif record.id in mecc_info:
                idx, start, end = mecc_info[record.id]
                sequence = str(record.seq[start-1:end])
                mecc_df.at[idx, 'eSeq'] = sequence
            elif record.id in Cecc_info:
                idx, start, end = Cecc_info[record.id]
                sequence = str(record.seq[start-1:end])
                Cecc_df.at[idx, 'eSeq'] = sequence

        return Uecc_df, mecc_df, Cecc_df

    def process_xecc_sequences(self, processed_ids=None):
        if processed_ids is None:
            processed_ids = set()

        xecc_sequences = []

        for record in SeqIO.parse(self.circular_seq_fasta, "fasta"):
            if record.id not in processed_ids:
                half_length = len(record.seq) // 2
                sub_record = record[:half_length]
                sub_record.id = record.id
                sub_record.description = record.description
                xecc_sequences.append(sub_record)

        if xecc_sequences:
            SeqIO.write(xecc_sequences, self.xecc_output_fa, "fasta")
            self.logger.info(f"Wrote {len(xecc_sequences)} XeccDNA sequences to {self.xecc_output_fa}")

        return len(xecc_sequences)

    def split_into_chunks(self, df, chunk_size):
        groups = list(df.groupby('query_id'))
        chunks = []
        current_chunk = []
        current_size = 0

        for name, group in groups:
            group_size = len(group)
            if current_size + group_size > chunk_size and current_chunk:
                chunks.append(pd.concat(current_chunk))
                current_chunk = []
                current_size = 0
            current_chunk.append(group)
            current_size += group_size

        if current_chunk:
            chunks.append(pd.concat(current_chunk))

        return chunks

    def process_chunk(self, chunk_data):
        chunk_df = chunk_data

        df = chunk_df.copy()
        df = self.add_calculated_columns(df)
        low_gap_df = self.filter_low_gap_percentage(df)

        low_gap_df['occurrence_count'] = low_gap_df.groupby('query_id')['query_id'].transform('count')
        low_gap_df['eClass'] = ''
        low_gap_df.loc[low_gap_df['occurrence_count'] == 1, 'eClass'] = 'Uecc'

        duplicate_mask = low_gap_df['occurrence_count'] > 1
        duplicate_groups = low_gap_df[duplicate_mask].groupby('query_id')

        results = []
        for name, group in duplicate_groups:
            eclass, keep_index = self.process_group(group)
            if eclass == 'Uecc':
                results.append((name, eclass, [keep_index]))
            else:
                results.append((name, eclass, group.index.tolist()))

        for name, eclass, indices in results:
            low_gap_df.loc[low_gap_df['query_id'] == name, 'eClass'] = eclass
            if eclass == 'Uecc':
                low_gap_df = low_gap_df[~((low_gap_df['query_id'] == name) & (~low_gap_df.index.isin(indices)))]

        low_gap_df = low_gap_df.drop('occurrence_count', axis=1)

        eClass_mapping = low_gap_df[['query_id', 'eClass']].drop_duplicates()
        df = df.merge(eClass_mapping, on='query_id', how='left')

        return df, low_gap_df

    def run(self):
        self.logger.info("Starting Ringmaster analysis")
        start_time = time.time()

        df = self.read_blast_results()
        df = self.process_query_id(df)

        chunk_size = 20000
        chunks = self.split_into_chunks(df, chunk_size)
        self.logger.info(f"Processing {len(chunks)} data chunks")

        with Pool(processes=self.num_threads) as pool:
            results = pool.map(self.process_chunk, chunks)

        dfs, low_gap_dfs = zip(*results)
        df = pd.concat(dfs, ignore_index=True)
        low_gap_df = pd.concat(low_gap_dfs, ignore_index=True)

        Uecc_results = self.process_Uecc(low_gap_df)
        mecc_results = self.process_mecc(low_gap_df)

        other_results = self.process_other_blast_rows(df)
        if other_results.empty:
            self.logger.warning("No other BLAST results found")
            Cecc_results = pd.DataFrame()
        else:
            self.logger.info(f"Processing {len(other_results)} other BLAST results")
            detailed_results_list = []
            summary_results_list = []

            grouped = other_results.groupby('query_id')

            for name, group in tqdm(grouped, total=len(grouped), desc="Analyzing groups"):
                detailed, summary = self.analyze_group(group)
                if not detailed.empty:
                    detailed_results_list.append(detailed)
                    summary_results_list.append(summary)

            if detailed_results_list:
                detailed_results = pd.concat(detailed_results_list, ignore_index=True)
                processed_detailed_results = self.process_detailed_results(detailed_results)
                if not processed_detailed_results.empty:
                    Cecc_results = self.process_detailed_results_further(processed_detailed_results)
                else:
                    Cecc_results = pd.DataFrame()
            else:
                Cecc_results = pd.DataFrame()

        Uecc_results, mecc_results, Cecc_results = self.process_sequences(
            self.circular_seq_fasta, 
            Uecc_results, 
            mecc_results, 
            Cecc_results
        )

        Uecc_results.to_csv(self.Uecc_output_csv, index=False)
        mecc_results.to_csv(self.Mecc_output_csv, index=False)
        if not Cecc_results.empty:
            Cecc_results.to_csv(self.Cecc_output_csv, index=False)

        processed_ids = set()
        processed_ids.update(Uecc_results['query_id'].unique())
        processed_ids.update(mecc_results['query_id'].unique())
        if not Cecc_results.empty:
            processed_ids.update(Cecc_results['query_id'].unique())

        if self.process_xecc:
            num_xecc = self.process_xecc_sequences(processed_ids)
            self.logger.info(f"Processed {num_xecc} XeccDNA sequences")

        execution_time = time.time() - start_time
        self.logger.info(f"Analysis completed in {execution_time:.2f} seconds")

        print(f"Processing complete. Output files:")
        print(f" - Uecc: {self.Uecc_output_csv}")
        print(f" - Mecc: {self.Mecc_output_csv}")
        print(f" - Cecc: {self.Cecc_output_csv if not Cecc_results.empty else 'No Cecc results'}")
        if self.process_xecc:
            print(f"Additional output: {self.xecc_output_fa}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process BLAST results and generate output files')
    parser.add_argument('-i', '--input', required=True, help='Input BLAST results file path')
    parser.add_argument('-fa', '--input_fasta', required=True, help='Input FASTA file path')
    parser.add_argument('-u', '--Uecc_output_csv', default='UeccDNA.csv', help='Uecc output CSV file path')
    parser.add_argument('-m', '--Mecc_output_csv', default='MeccDNA.csv', help='Mecc output CSV file path')
    parser.add_argument('-c', '--Cecc_output_csv', default='CeccDNA.csv', help='Cecc output CSV file path')
    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of CPU cores to use')
    parser.add_argument('-X', '--process_xecc', action='store_true', help='Generate XeccDNA.fa file for unclassified sequences')
    parser.add_argument('--prefix', default=None, help='Prefix (sample name) for output files, especially Xecc')
    parser.add_argument('--output_dir', default=None, help='Directory to save output files')

    args = parser.parse_args()

    ringmaster = Ringmaster(
        args.input,
        args.input_fasta,
        args.Uecc_output_csv,
        args.Mecc_output_csv,
        args.Cecc_output_csv,
        args.process_xecc,
        args.threads,
        prefix=args.prefix,
        output_dir=args.output_dir
    )
    ringmaster.run()
