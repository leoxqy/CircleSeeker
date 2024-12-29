"""
MergeCecc: Chimeric eccDNA (CECC) Processing Module.

This module handles the processing and merging of Chimeric eccDNA (CECC) data. It provides functionality for processing both standard CECC and
Multiple alignment Chimeric CECC (MCECC) sequences, including sequence analysis and FASTA file generation.

Key features:
- CECC and MCECC sequence processing
- FASTA file generation for both CECC and MCECC
- CSV data processing and validation
- Sequence record management using BioPython
- Comprehensive logging system

Typical usage:
    processor = CeccProcessor(args)
    processor.run()

Version: 1.0.1
License: GNU General Public License v2
Copyright (c) 2024 CircleSeeker Team
"""

import pandas as pd
import logging
import argparse
import sys
import os
from typing import List, Dict, Tuple
from collections import defaultdict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

class CeccProcessor:
    def __init__(self, args):
        self.input_csv = args.input_csv
        self.output_cecc_fasta = args.output_cecc_fasta
        self.output_mcecc_fasta = args.output_mcecc_fasta
        self.final_cecc_csv = args.final_cecc_csv
        self.final_mcecc_csv = args.final_mcecc_csv
        
        # Set up logger
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        self.logger.addHandler(console_handler)
        
        self.cecc_id_mapping = {}
        self.mcecc_id_mapping = {}

    def _fill_sequences(self, df: pd.DataFrame) -> pd.DataFrame:
        """Fill missing sequences based on query_id groups."""
        self.logger.info("Filling missing sequences...")
        
        # Create temporary column for sorting
        df['has_seq'] = df['eSeq'].notna()
        
        # Sort dataframe by query_id and sequence presence
        df_sorted = df.sort_values(['query_id', 'has_seq'], ascending=[True, False])
        
        # Group by query_id and forward fill the eSeq
        df_filled = df_sorted.groupby('query_id', as_index=False).apply(
            lambda x: x.fillna(method='ffill')
        )
        
        # Remove temporary column
        df_filled = df_filled.drop('has_seq', axis=1)
        
        # Verify all sequences are filled
        missing_seq = df_filled['eSeq'].isna().sum()
        if missing_seq > 0:
            self.logger.warning(f"Found {missing_seq} rows with missing sequences after filling")
        
        return df_filled

    def _split_by_eclass(self, df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Split data into Cecc and Ceccm based on 'eClass'."""
        cecc_df = df[df['eClass'] == 'Cecc']
        ceccm_df = df[df['eClass'] == 'Ceccm']
        return cecc_df, ceccm_df

    def _create_location_signature(self, group_df: pd.DataFrame) -> tuple:
        """Create a location signature for each group."""
        group_df = group_df.sort_values(['subject_id', 's_start', 's_end'])
        positions = []
        current_subject = None
        current_positions = []

        for _, row in group_df.iterrows():
            if current_subject != row['subject_id']:
                if current_positions:
                    positions.append((current_subject, tuple(sorted(current_positions))))
                current_subject = row['subject_id']
                current_positions = []
            current_positions.append((row['s_start'], row['s_end']))

        if current_positions:
            positions.append((current_subject, tuple(sorted(current_positions))))

        positions.sort(key=lambda x: x[0])
        return tuple(positions)

    def _merge_groups(self, groups: List[pd.DataFrame]) -> pd.DataFrame:
        """Merge groups with the same location signature."""
        if not groups:
            return pd.DataFrame()

        base_query_id = min(g['query_id'].iloc[0] for g in groups)
        base_group = next(g for g in groups if g['query_id'].iloc[0] == base_query_id)

        merged_query_ids = ';'.join(sorted(set(g['query_id'].iloc[0] for g in groups)))
        all_reads = set()
        for g in groups:
            reads = g['readName'].iloc[0].split('|')
            all_reads.update(reads)
        merged_readnames = '|'.join(sorted(all_reads))
        total_copynum = sum(float(g['copyNum'].iloc[0]) for g in groups)

        merged_data = base_group.copy()
        merged_data['Pquery_id'] = base_query_id
        merged_data['Mquery_id'] = merged_query_ids
        merged_data['readName'] = merged_readnames
        merged_data['copyNum'] = total_copynum

        return merged_data

    def _generate_map_id(self, df: pd.DataFrame) -> pd.DataFrame:
        """Generate continuous MapID for the DataFrame."""
        if df.empty:
            return df
            
        if 'Prelist_ID' not in df.columns:
            self.logger.error("Missing Prelist_ID column")
            raise ValueError("Missing Prelist_ID column")

        result_df = df.copy()

        try:
            prelist_ids = []
            seen_ids = set()

            for pid in result_df['Prelist_ID']:
                if pid not in seen_ids:
                    prelist_ids.append(pid)
                    seen_ids.add(pid)

            id_mapping = {pid: idx + 1 for idx, pid in enumerate(prelist_ids)}
            result_df['MapID'] = result_df['Prelist_ID'].map(id_mapping)
            result_df = result_df.drop('Prelist_ID', axis=1)

            return result_df

        except Exception as e:
            self.logger.error(f"Error generating MapID: {str(e)}")
            raise

    def process_cecc(self, cecc_df: pd.DataFrame) -> pd.DataFrame:
        """Process Cecc data using location-based strategy."""
        if cecc_df.empty:
            return pd.DataFrame()
            
        location_dict = defaultdict(list)

        for query_id, group in cecc_df.groupby('query_id'):
            loc_signature = self._create_location_signature(group)
            location_dict[loc_signature].append(group)

        merged_results = []
        for signature, groups in location_dict.items():
            merged_group = self._merge_groups(groups)
            merged_results.append(merged_group)

        if merged_results:
            final_df = pd.concat(merged_results, ignore_index=True)
            final_df['ReadsNum'] = final_df['readName'].apply(lambda x: len(x.split('|')))

            chimera_counts = {}
            for pquery_id in final_df['Pquery_id'].unique():
                rows = cecc_df[cecc_df['query_id'].str.startswith(pquery_id)]
                chimera_counts[pquery_id] = len(rows)

            final_df['ChimeraNum'] = final_df['Pquery_id'].map(chimera_counts)

            columns_order = [
                'Pquery_id', 'Mquery_id', 'consLen', 'LocNum', 'ChimeraNum', 'ReadsNum',
                'identity', 'alignment_length', 'subject_id', 's_start', 's_end',
                's_length', 'strand', 'readName', 'copyNum', 'eClass', 'eSeq'
            ]
            final_df = final_df[columns_order]

            return final_df
        else:
            return pd.DataFrame()

    def process_ceccm(self, ceccm_df: pd.DataFrame) -> pd.DataFrame:
        """Process Ceccm data and generate features."""
        if ceccm_df.empty:
            return pd.DataFrame()

        columns_to_drop = ['q_start', 'q_end', 'ratio']
        ceccm_df = ceccm_df.drop(columns=[col for col in columns_to_drop if col in ceccm_df.columns])

        features_dict = self.process_ceccm_features(ceccm_df)
        merged_df = self.find_and_merge_queries(features_dict, ceccm_df)

        if not merged_df.empty:
            chimera_counts = {}
            for query_id in merged_df['query_id'].unique():
                unique_prelist_ids = len(ceccm_df[ceccm_df['query_id'] == query_id]['Prelist_ID'].unique())
                chimera_counts[query_id] = unique_prelist_ids

            merged_df['ChimeraNum'] = merged_df['query_id'].map(chimera_counts)
            merged_df['ReadsNum'] = merged_df['readName'].apply(lambda x: len(x.split('|')))
            merged_df = self._generate_map_id(merged_df)

        return merged_df

    def process_ceccm_features(self, ceccm_df: pd.DataFrame) -> Dict[str, List[Dict]]:
        """Process Ceccm data and group by query_id and Prelist_ID to generate features."""
        query_features = {}

        for query_id, qgroup in ceccm_df.groupby('query_id'):
            features_list = []

            for prelist_id, ngroup in qgroup.groupby('Prelist_ID'):
                ngroup_sorted = ngroup.sort_values(['subject_id', 's_start', 's_end'])

                feature = {
                    'chr': '|'.join(ngroup_sorted['subject_id'].astype(str)),
                    'start': '|'.join(ngroup_sorted['s_start'].astype(str)),
                    'end': '|'.join(ngroup_sorted['s_end'].astype(str))
                }

                features_list.append(feature)

            query_features[query_id] = features_list

        return query_features

    def find_and_merge_queries(self, features_dict: Dict[str, List[Dict]], df: pd.DataFrame) -> pd.DataFrame:
        """Find and merge query_id groups with the same features."""
        if df.empty:
            return pd.DataFrame()

        feature_to_queries = defaultdict(set)

        for query_id, features in features_dict.items():
            for feature in features:
                feature_tuple = (
                    feature['chr'],
                    feature['start'],
                    feature['end']
                )
                feature_to_queries[feature_tuple].add(query_id)

        connected_queries = defaultdict(set)

        for queries in feature_to_queries.values():
            if len(queries) > 1:
                queries = list(queries)
                for i in range(len(queries)):
                    connected_queries[queries[i]].update(queries)

        merge_groups = []
        while connected_queries:
            start_query = next(iter(connected_queries))
            group = set()
            to_check = {start_query}

            while to_check:
                query = to_check.pop()
                if query not in group:
                    group.add(query)
                    connected = connected_queries.get(query, set())
                    to_check.update(connected - group)

            if len(group) > 1:
                merge_groups.append(list(group))

            for query in group:
                connected_queries.pop(query, None)

        merged_results = []
        processed_queries = set()

        for group in merge_groups:
            max_features_query = max(group, key=lambda x: len(features_dict[x]))
            base_data = df[df['query_id'] == max_features_query].copy()

            all_reads = set()
            for qid in group:
                reads = df[df['query_id'] == qid]['readName'].iloc[0].split('|')
                all_reads.update(reads)
            merged_readnames = '|'.join(sorted(all_reads))

            merged_query_ids = ';'.join(sorted(group))
            total_copynum = sum(
                float(df[df['query_id'] == qid]['copyNum'].iloc[0])
                for qid in group
            )

            base_data['Mquery_id'] = merged_query_ids
            base_data['readName'] = merged_readnames
            base_data['copyNum'] = total_copynum

            merged_results.append(base_data)
            processed_queries.update(group)

        unmerged_queries = set(features_dict.keys()) - processed_queries
        for query_id in unmerged_queries:
            single_data = df[df['query_id'] == query_id].copy()
            single_data['Mquery_id'] = query_id
            merged_results.append(single_data)

        if merged_results:
            final_df = pd.concat(merged_results, ignore_index=True)
            return final_df
        else:
            return pd.DataFrame()

    def _generate_fasta_from_df(self, df: pd.DataFrame, output_file: str, id_column: str):
        """Generate FASTA file from DataFrame."""
        if df.empty:
            with open(output_file, 'w') as f:
                pass
            return

        records = []
        for _, row in df.iterrows():
            if pd.notna(row['eSeq']):
                record = SeqRecord(
                    Seq(row['eSeq']),
                    id=str(row[id_column]),
                    description=""
                )
                records.append(record)

        # Write unique sequences only
        seen_ids = set()
        unique_records = []
        for record in records:
            if record.id not in seen_ids:
                seen_ids.add(record.id)
                unique_records.append(record)

        SeqIO.write(unique_records, output_file, "fasta")

    def rename_ids(self, cecc_df: pd.DataFrame, mcecc_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Rename IDs in dataframes."""
        try:
            cecc_df_renamed = self.process_cecc_ids(cecc_df)
            mcecc_df_renamed = self.process_mcecc_ids(mcecc_df)
            return cecc_df_renamed, mcecc_df_renamed
        except Exception as e:
            self.logger.error(f"An error occurred during ID renaming: {str(e)}")
            raise

    def process_cecc_ids(self, cecc_df: pd.DataFrame) -> pd.DataFrame:
        """Process Cecc dataframe and create ID mappings."""
        if cecc_df.empty:
            return cecc_df
            
        cecc_pquery_ids = sorted(cecc_df['Pquery_id'].unique())
        self.cecc_id_mapping = {
            old_id: self._create_new_id("Cecc", idx+1)
            for idx, old_id in enumerate(cecc_pquery_ids)
        }
        cecc_df['Pquery_id'] = cecc_df['Pquery_id'].map(self.cecc_id_mapping)
        return cecc_df

    def process_mcecc_ids(self, mcecc_df: pd.DataFrame) -> pd.DataFrame:
        """Process MCecc dataframe and create ID mappings."""
        if mcecc_df.empty:
            return mcecc_df
            
        mcecc_query_ids = sorted(mcecc_df['query_id'].unique())
        self.mcecc_id_mapping = {
            old_id: self._create_new_id("MCecc", idx+1)
            for idx, old_id in enumerate(mcecc_query_ids)
        }
        mcecc_df['query_id'] = mcecc_df['query_id'].map(self.mcecc_id_mapping)
        return mcecc_df

    def _create_new_id(self, prefix: str, index: int) -> str:
        """Create a new ID."""
        return f"{prefix}|Confirmed|{index}"

    def rename_cecc_columns(self, cecc_df: pd.DataFrame) -> pd.DataFrame:
        """Process Cecc dataframe to rename columns."""
        if cecc_df.empty:
            return pd.DataFrame(columns=[
                'eName', 'ChimeraNum', 'eLength', 'eRepeatNum', 'eClass', 'eState',
                'LocNum', 'Identity', 'Alignment_length', 'Chr', 'Start', 'End',
                'Strand', 'ReadsNum', 'eReads', 'eSeq'
            ])

        new_df = pd.DataFrame()
        new_df['eName'] = cecc_df['Pquery_id']
        new_df['ChimeraNum'] = cecc_df['ChimeraNum']
        new_df['eLength'] = cecc_df['consLen']
        new_df['eRepeatNum'] = cecc_df['copyNum']
        new_df['eClass'] = 'Cecc'
        new_df['eState'] = 'Confirmed-eccDNA'
        new_df['LocNum'] = cecc_df['LocNum']
        new_df['Identity'] = cecc_df['identity']
        new_df['Alignment_length'] = cecc_df['alignment_length']
        new_df['Chr'] = cecc_df['subject_id']
        new_df['Start'] = cecc_df['s_start']
        new_df['End'] = cecc_df['s_end']
        new_df['Strand'] = cecc_df['strand']
        new_df['ReadsNum'] = cecc_df['ReadsNum']
        new_df['eReads'] = cecc_df['readName']
        new_df['eSeq'] = cecc_df['eSeq']

        cecc_columns = [
            'eName', 'ChimeraNum', 'eLength', 'eRepeatNum', 'eClass', 'eState',
            'LocNum', 'Identity', 'Alignment_length', 'Chr', 'Start', 'End',
            'Strand', 'ReadsNum', 'eReads', 'eSeq'
        ]

        return new_df[cecc_columns]

    def rename_mcecc_columns(self, mcecc_df: pd.DataFrame) -> pd.DataFrame:
        """Process MCecc dataframe to rename columns."""
        if mcecc_df.empty:
            return pd.DataFrame(columns=[
                'eName', 'MapID', 'ChimeraNum', 'eLength', 'eRepeatNum', 'eClass',
                'eState', 'LocNum', 'Identity', 'Alignment_length', 'Chr', 'Start',
                'End', 'Strand', 'ReadsNum', 'eReads', 'eSeq'
            ])

        new_df = pd.DataFrame()
        new_df['eName'] = mcecc_df['query_id']
        new_df['MapID'] = mcecc_df['MapID']
        new_df['ChimeraNum'] = mcecc_df['ChimeraNum']
        new_df['eLength'] = mcecc_df['consLen']
        new_df['eRepeatNum'] = mcecc_df['copyNum']
        new_df['eClass'] = 'MCecc'
        new_df['eState'] = 'Confirmed-eccDNA'
        new_df['LocNum'] = mcecc_df['LocNum']
        new_df['Identity'] = mcecc_df['identity']
        new_df['Alignment_length'] = mcecc_df['alignment_length']
        new_df['Chr'] = mcecc_df['subject_id']
        new_df['Start'] = mcecc_df['s_start']
        new_df['End'] = mcecc_df['s_end']
        new_df['Strand'] = mcecc_df['strand']
        new_df['ReadsNum'] = mcecc_df['ReadsNum']
        new_df['eReads'] = mcecc_df['readName']
        new_df['eSeq'] = mcecc_df['eSeq']

        mcecc_columns = [
            'eName', 'MapID', 'ChimeraNum', 'eLength', 'eRepeatNum', 'eClass',
            'eState', 'LocNum', 'Identity', 'Alignment_length', 'Chr', 'Start',
            'End', 'Strand', 'ReadsNum', 'eReads', 'eSeq'
        ]

        return new_df[mcecc_columns]

    def run(self):
        """Run the complete processing workflow."""
        try:
            self.logger.info(f"Reading input CSV: {self.input_csv}")
            
            # Read and fill sequences
            df = pd.read_csv(self.input_csv)
            df = self._fill_sequences(df)
            
            # Split and process data
            cecc_df, ceccm_df = self._split_by_eclass(df)
            
            self.logger.info(f"Found {len(cecc_df)} Cecc records and {len(ceccm_df)} MCecc records")
            
            # Process Cecc data
            if not cecc_df.empty:
                cecc_df_processed = self.process_cecc(cecc_df)
                self.logger.info(f"Processed {len(cecc_df_processed)} Cecc records")
            else:
                cecc_df_processed = pd.DataFrame()
                self.logger.warning("No Cecc data found in input file")

            # Process MCecc data
            if not ceccm_df.empty:
                mcecc_df_processed = self.process_ceccm(ceccm_df)
                self.logger.info(f"Processed {len(mcecc_df_processed)} MCecc records")
            else:
                mcecc_df_processed = pd.DataFrame()
                self.logger.warning("No MCecc data found in input file")

            # Generate FASTA files from processed data
            if not cecc_df_processed.empty:
                cecc_df_renamed, mcecc_df_renamed = self.rename_ids(cecc_df_processed, mcecc_df_processed)
                self._generate_fasta_from_df(cecc_df_renamed, self.output_cecc_fasta, 'Pquery_id')
                self.logger.info(f"Generated Cecc FASTA file: {self.output_cecc_fasta}")
            
            if not mcecc_df_processed.empty:
                self._generate_fasta_from_df(mcecc_df_renamed, self.output_mcecc_fasta, 'query_id')
                self.logger.info(f"Generated MCecc FASTA file: {self.output_mcecc_fasta}")

            # Generate final CSV files
            final_cecc_df = self.rename_cecc_columns(cecc_df_renamed)
            final_cecc_df.to_csv(self.final_cecc_csv, index=False)
            self.logger.info(f"Written Cecc results to {self.final_cecc_csv}")

            final_mcecc_df = self.rename_mcecc_columns(mcecc_df_renamed)
            final_mcecc_df.to_csv(self.final_mcecc_csv, index=False)
            self.logger.info(f"Written MCecc results to {self.final_mcecc_csv}")

            self.logger.info("Processing completed successfully")

        except Exception as e:
            self.logger.error(f"Error during processing: {str(e)}")
            raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process Cecc and MCecc data with sequences.')
    parser.add_argument('-i', '--input_csv', required=True, help='Input CSV file with sequences')
    parser.add_argument('-fa_cecc', '--output_cecc_fasta', required=True, help='Output Cecc FASTA file')
    parser.add_argument('-fa_mcecc', '--output_mcecc_fasta', required=True, help='Output MCecc FASTA file')
    parser.add_argument('-final_cecc_csv', '--final_cecc_csv', required=True, help='Final Cecc CSV file')
    parser.add_argument('-final_mcecc_csv', '--final_mcecc_csv', required=True, help='Final MCecc CSV file')

    args = parser.parse_args()

    processor = CeccProcessor(args)
    processor.run()