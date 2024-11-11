# cecc_processor.py

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import logging
import argparse
import sys
import os
import time
from typing import List, Dict, Tuple
from collections import defaultdict

class CeccProcessor:
    def __init__(self, args):
        self.input_csv = args.input_csv
        self.fasta_input = args.fasta_input
        self.output_cecc_fasta = args.output_cecc_fasta
        self.output_mcecc_fasta = args.output_mcecc_fasta
        self.final_cecc_csv = args.final_cecc_csv
        self.final_mcecc_csv = args.final_mcecc_csv

        # print(f"Input CSV: {self.input_csv}")
        # print(f"Fasta Input: {self.fasta_input}")
        # print(f"Output Cecc Fasta: {self.output_cecc_fasta}")
        # print(f"Output Mcecc Fasta: {self.output_mcecc_fasta}")
        # print(f"Final Cecc CSV: {self.final_cecc_csv}")
        # print(f"Final Mcecc CSV: {self.final_mcecc_csv}")
        
        # Set up logger
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        self.logger.addHandler(console_handler)
        
        self.cecc_id_mapping = {}
        self.mcecc_id_mapping = {}

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
                's_length', 'strand', 'readName', 'copyNum', 'eClass'
            ]
            final_df = final_df[columns_order]

            return final_df
        else:
            return pd.DataFrame()

    def process_ceccm(self, ceccm_df: pd.DataFrame) -> pd.DataFrame:
        """Process Ceccm data and generate features."""
        if ceccm_df.empty:
            return pd.DataFrame()

        columns_to_drop = ['q_start', 'q_end', 'ratio', 'eClass']
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

    def extract_fasta_sequences(self, cecc_df_processed, mcecc_df_processed):
        """Extract sequences from FASTA file based on processed dataframes."""
        try:
            cecc_ids = self._get_query_ids(cecc_df_processed, is_cecc=True)
            mcecc_ids = self._get_query_ids(mcecc_df_processed, is_cecc=False)
            
            # Create empty FASTA files if no data
            if not cecc_ids:
                with open(self.output_cecc_fasta, 'w') as f:
                    pass
            else:
                self._extract_and_write_sequences(cecc_ids, self.output_cecc_fasta)
                
            if not mcecc_ids:
                with open(self.output_mcecc_fasta, 'w') as f:
                    pass
            else:
                self._extract_and_write_sequences(mcecc_ids, self.output_mcecc_fasta)
                
        except Exception as e:
            self.logger.error(f"An error occurred during FASTA extraction: {str(e)}")
            raise

    def _get_query_ids(self, df: pd.DataFrame, is_cecc: bool = True) -> set:
        """Get query_ids from the DataFrame."""
        query_ids = set()
        if df.empty:
            return query_ids
            
        if is_cecc:
            if 'Pquery_id' not in df.columns:
                raise ValueError(f"Column Pquery_id not found in DataFrame")
            query_ids.update(df['Pquery_id'])
        else:
            if 'query_id' not in df.columns:
                raise ValueError(f"Column query_id not found in DataFrame")
            query_ids.update(df['query_id'])
        return query_ids

    def _extract_and_write_sequences(self, query_ids: set, output_file: str):
        """Extract and write sequences for specified query_ids."""
        if not query_ids:
            # Create empty file if no sequences to extract
            with open(output_file, 'w') as f:
                pass
            return
            
        with open(output_file, 'w') as output_handle:
            for record in SeqIO.parse(self.fasta_input, "fasta"):
                if record.id in query_ids:
                    SeqIO.write(record, output_handle, "fasta")

    def rename_ids(self, cecc_df, mcecc_df):
        """Rename IDs in dataframes and FASTA files."""
        try:
            cecc_df_renamed = self.process_cecc_ids(cecc_df)
            mcecc_df_renamed = self.process_mcecc_ids(mcecc_df)
            self.process_fasta_files()
            return cecc_df_renamed, mcecc_df_renamed
        except Exception as e:
            self.logger.error(f"An error occurred during ID renaming: {str(e)}")
            raise

    def process_cecc_ids(self, cecc_df):
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

    def process_mcecc_ids(self, mcecc_df):
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

    def process_fasta_files(self):
        """Process FASTA files."""
        if self.cecc_id_mapping:
            self._process_single_fasta(
                self.output_cecc_fasta,
                self.output_cecc_fasta,
                self.cecc_id_mapping
            )

        if self.mcecc_id_mapping:
            self._process_single_fasta(
                self.output_mcecc_fasta,
                self.output_mcecc_fasta,
                self.mcecc_id_mapping
            )

    def _process_single_fasta(self, input_file: str, output_file: str, id_mapping: dict):
        """Process a single FASTA file."""
        if not os.path.exists(input_file):
            self.logger.warning(f"Input file does not exist: {input_file}")
            return

        new_records = []
        for record in SeqIO.parse(input_file, "fasta"):
            if record.id in id_mapping:
                new_record = SeqRecord(
                    seq=record.seq,
                    id=id_mapping[record.id],
                    description=""
                )
                new_records.append(new_record)

        SeqIO.write(new_records, output_file, "fasta")

    def rename_cecc_columns(self, cecc_df):
        """Process Cecc dataframe to rename columns."""
        if cecc_df.empty:
            return pd.DataFrame(columns=[
                'eName', 'ChimeraNum', 'eLength', 'eRepeatNum', 'eClass', 'eState',
                'LocNum', 'Identity', 'Alignment_length', 'Chr', 'Start', 'End',
                'Strand', 'ReadsNum', 'eReads'
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

        cecc_columns = [
            'eName', 'ChimeraNum', 'eLength', 'eRepeatNum', 'eClass', 'eState',
            'LocNum', 'Identity', 'Alignment_length', 'Chr', 'Start', 'End',
            'Strand', 'ReadsNum', 'eReads'
        ]

        return new_df[cecc_columns]

    def rename_mcecc_columns(self, mcecc_df):
        """Process MCecc dataframe to rename columns."""
        if mcecc_df.empty:
            return pd.DataFrame(columns=[
                'eName', 'MapID', 'ChimeraNum', 'eLength', 'eRepeatNum', 'eClass',
                'eState', 'LocNum', 'Identity', 'Alignment_length', 'Chr', 'Start',
                'End', 'Strand', 'ReadsNum', 'eReads'
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

        mcecc_columns = [
            'eName', 'MapID', 'ChimeraNum', 'eLength', 'eRepeatNum', 'eClass',
            'eState', 'LocNum', 'Identity', 'Alignment_length', 'Chr', 'Start',
            'End', 'Strand', 'ReadsNum', 'eReads'
        ]

        return new_df[mcecc_columns]

    def run(self):
        """Run the complete processing workflow."""
        try:
            self.logger.info(f"Input: CSV={self.input_csv}, FASTA={self.fasta_input}")
            
            # Step 1: Process Cecc and MCecc data
            df = pd.read_csv(self.input_csv)
            cecc_df, ceccm_df = self._split_by_eclass(df)
            
            # 添加数据存在性检查的日志
            self.logger.info(f"Found {len(cecc_df)} Cecc records and {len(ceccm_df)} MCecc records")
            
            # 处理Cecc数据
            if not cecc_df.empty:
                cecc_df_processed = self.process_cecc(cecc_df)
                self.logger.info(f"Processed {len(cecc_df_processed)} Cecc records")
            else:
                cecc_df_processed = pd.DataFrame()
                self.logger.warning("No Cecc data found in input file")

            # 处理MCecc数据
            if not ceccm_df.empty:
                mcecc_df_processed = self.process_ceccm(ceccm_df)
                self.logger.info(f"Processed {len(mcecc_df_processed)} MCecc records")
            else:
                mcecc_df_processed = pd.DataFrame()
                self.logger.warning("No MCecc data found in input file")

            # Step 2: Extract sequences from FASTA file
            if not cecc_df_processed.empty or not mcecc_df_processed.empty:
                self.extract_fasta_sequences(cecc_df_processed, mcecc_df_processed)

            # Step 3: Rename IDs in dataframes and FASTA files
            if not cecc_df_processed.empty or not mcecc_df_processed.empty:
                cecc_df_renamed, mcecc_df_renamed = self.rename_ids(cecc_df_processed, mcecc_df_processed)
            else:
                cecc_df_renamed = pd.DataFrame()
                mcecc_df_renamed = pd.DataFrame()

            # Step 4: Rename columns in dataframes and write files
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
    parser = argparse.ArgumentParser(description='Process Cecc and MCecc data in one step.')
    parser.add_argument('-i', '--input_csv', required=True, help='Input CSV file')
    parser.add_argument('-fa', '--fasta_input', required=True, help='Input FASTA file')
    parser.add_argument('-fa_cecc', '--output_cecc_fasta', required=True, help='Output Cecc FASTA file')
    parser.add_argument('-fa_mcecc', '--output_mcecc_fasta', required=True, help='Output MCecc FASTA file')
    parser.add_argument('-final_cecc_csv', '--final_cecc_csv', required=True, help='Final Cecc CSV file')
    parser.add_argument('-final_mcecc_csv', '--final_mcecc_csv', required=True, help='Final MCecc CSV file')

    args = parser.parse_args()

    processor = CeccProcessor(args)
    processor.run()