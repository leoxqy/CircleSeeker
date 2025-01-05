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

Version: 1.0.1
License: GNU General Public License v2
Copyright (c) 2024 CircleSeeker Team
"""

import argparse
import logging
import pandas as pd
from collections import defaultdict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

class MergeMecc:
    def __init__(self, meccdna_part1, meccdna_part2, output_csv, output_fasta):
        # logging
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
        """Main workflow"""
        # 1) Read input
        self.logger.info("Reading input files...")
        df1 = self.read_csv(self.meccdna_part1)
        df2 = self.read_csv(self.meccdna_part2)

        # 2) Basic processing
        self.logger.info("Processing MeccDNA part1 data...")
        df1_processed = self.process_meccdna_part1(df1)

        self.logger.info("Processing MeccDNA part2 data...")
        df2_processed = self.process_meccdna_part2(df2)

        # 3) Combine
        self.logger.info("Combining data...")
        combined_df = pd.concat([df1_processed, df2_processed], ignore_index=True)
        if combined_df.empty:
            self.logger.warning("No valid data after combining. Creating empty output files...")
            self.write_empty_output()
            return

        if 'query_id' not in combined_df.columns:
            self.logger.error("Required column 'query_id' missing.")
            return

        self.logger.info("Force recalc eLength,eEnd from eSeq...")
        combined_df['eLength'] = combined_df['eSeq'].apply(
            lambda s: len(str(s)) if pd.notnull(s) else 0
        )
        combined_df['eEnd'] = combined_df['eStart'] + combined_df['eLength'] - 1

        self.logger.info("Generating approximate location_feature by dividing start/end by 50 & using integer division.")
        combined_df['location_feature'] = combined_df.apply(
            lambda row: self.generate_approx_location_feature(
                row['eChr'], row['eStart'], row['eEnd']
            ),
            axis=1
        )

        # 6) partial intersection merging with union-find
        self.logger.info("Grouping query_ids by union-find, if they share any location_feature.")
        query_groups = self.build_query_groups(combined_df)
        merged_groups = self.union_find_merge(query_groups)

        # 7) integrate => keep the rep query_id
        self.logger.info("Integrating merged info back to combined_df...")
        integrated_df = self.integrate_data(merged_groups, combined_df)

        # 8) add name/state, finalize & output
        self.logger.info("Adding name/state to final records...")
        final_df = self.add_name_and_state(integrated_df)

        self.logger.info("Preparing final DataFrame & saving CSV...")
        final_results = self.prepare_final_dataframe(final_df)
        final_results.to_csv(self.output_csv, index=False)

        self.logger.info("Creating FASTA from final results...")
        self.create_fasta_from_df(final_results, self.output_fasta)

        self.logger.info("Processing completed.")

    def generate_approx_location_feature(self, chr_, start, end):
        """
        使用 (start//50, end//50) 作为 approximate coordinate bin（取整而非四舍五入）
        """
        start_approx = start // 50
        end_approx   = end // 50
        
        return f"{chr_}-{start_approx}-{end_approx}"

    def build_query_groups(self, df):
        """
        group by query_id => collect all location_feature into a list
        sum up eRepeatNum, combine eReads
        also store feature_count
        """
        grouped = df.groupby('query_id').agg({
            'location_feature': list,
            'eRepeatNum': 'sum',   # or first, up to you
            'eReads': lambda x: ';'.join(set(';'.join(x).split(';'))),
        }).reset_index()

        grouped['feature_count'] = grouped['location_feature'].apply(len)
        return grouped

    def union_find_merge(self, query_groups):
        """
        union-find => if 2 qid share any location_feature => same cluster
        pick the one with largest feature_count in that cluster as rep
        sum eRepeatNum, merge eReads
        """
        if query_groups.empty:
            return query_groups

        from collections import defaultdict
        location_map = defaultdict(list)
        for _, row in query_groups.iterrows():
            qid = row['query_id']
            for loc in row['location_feature']:
                location_map[loc].append(qid)

        all_qids = query_groups['query_id'].unique().tolist()
        parent = {q: q for q in all_qids}

        def find(x):
            if parent[x] != x:
                parent[x] = find(parent[x])
            return parent[x]

        def union(a, b):
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[rb] = ra

        # link
        for loc, qlist in location_map.items():
            if len(qlist) > 1:
                base = qlist[0]
                for other in qlist[1:]:
                    union(base, other)

        # cluster
        clusters = defaultdict(list)
        rowdicts = query_groups.to_dict('records')
        for row in rowdicts:
            root = find(row['query_id'])
            clusters[root].append(row)

        merged_rows = []
        for root_id, items in clusters.items():
            all_loc = set()
            sum_repeat = 0
            all_reads = set()
            max_feat_count = 0
            rep_qid = None

            for it in items:
                all_loc.update(it['location_feature'])
                sum_repeat += it['eRepeatNum']
                all_reads.update(it['eReads'].split(';'))
                if it['feature_count'] > max_feat_count:
                    max_feat_count = it['feature_count']
                    rep_qid = it['query_id']

            merged_rows.append({
                'query_id': rep_qid,
                'location_feature': sorted(all_loc),
                'eRepeatNum': sum_repeat,
                'eReads': ';'.join(sorted(all_reads)),
                'feature_count': len(all_loc)
            })
        return pd.DataFrame(merged_rows)

    def integrate_data(self, merged_df, raw_df):
        """Keep the rep query_id, update eRepeatNum,eReads"""
        if merged_df.empty:
            return raw_df.iloc[0:0]

        merged_map = merged_df.set_index('query_id')[['eRepeatNum','eReads']].to_dict('index')
        valid_qids = set(merged_df['query_id'])
        filtered = raw_df[raw_df['query_id'].isin(valid_qids)].copy()

        filtered['eRepeatNum'] = filtered['query_id'].map(lambda x: merged_map[x]['eRepeatNum'])
        filtered['eReads']     = filtered['query_id'].map(lambda x: merged_map[x]['eReads'])

        return filtered[[
            'eChr','eStart','eEnd','eLength','MatDegree',
            'eClass','eRepeatNum','eReads','query_id','eSeq'
        ]]

    def read_csv(self, file_path):
        try:
            return pd.read_csv(file_path)
        except (FileNotFoundError, pd.errors.EmptyDataError):
            self.logger.warning(f"File {file_path} not found or empty.")
            return pd.DataFrame()

    def process_meccdna_part1(self, df):
        if df.empty:
            return pd.DataFrame(columns=[
                'eChr','eStart','eEnd','eLength','eRepeatNum',
                'MatDegree','eClass','eReads','query_id','eSeq'
            ])

        rows=[]
        for _,row in df.iterrows():
            if row['consLen'] <= 100:
                continue
            if 'AlignRegion' not in row or pd.isnull(row['AlignRegion']):
                continue
            for region in row['AlignRegion'].split(';'):
                mat_degree, loc = region.split('|')
                chr_, se = loc.split('-',1)
                st, ed = map(int, se.split('-'))
                rows.append({
                    'eChr': chr_,
                    'eStart': st,
                    'eEnd': ed,
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
        empty_df = pd.DataFrame(columns=[
            'eChr','eStart','eEnd','eLength','eRepeatNum',
            'MatDegree','eClass','eReads','query_id','eSeq'
        ])
        if df.empty:
            return empty_df
        mecc_df = df[(df['eClass'] == 'Mecc') & (df['eLength'] > 100)].copy()
        if mecc_df.empty:
            return empty_df
        if 'qname' in mecc_df.columns:
            mecc_df['query_id'] = mecc_df['qname']

        return mecc_df[[
            'eChr','eStart','eEnd','eLength','eRepeatNum','MatDegree',
            'eClass','eReads','query_id','eSeq'
        ]]

    def add_name_and_state(self, df):
        unique_qids = sorted(df['query_id'].unique())
        qid_to_idx = {q: i+1 for i,q in enumerate(unique_qids)}
        df['eName'] = df['query_id'].map(lambda x: f"Mecc|Confirmed|{qid_to_idx[x]}")
        df['eState'] = 'Confirmed-eccDNA'
        df['tmp_sort'] = df['eName'].str.extract(r'(\d+)').astype(int)
        df = df.sort_values('tmp_sort').drop('tmp_sort', axis=1)
        return df

    def prepare_final_dataframe(self, df):
        df2 = df.copy()
        df2['eReadNum'] = df2['eReads'].apply(lambda s: len(str(s).split(';')))
        mapnum = df2.groupby('query_id').size()
        df2['MapNum'] = df2['query_id'].map(mapnum)
        cols = [
            'eName','MapNum','eRepeatNum','eState','eChr','eStart',
            'eEnd','eLength','MatDegree','eClass','eReadNum','eReads','eSeq'
        ]
        df2['sort_key'] = df2['eName'].str.extract(r'(\d+)').astype(int)
        df2 = df2.sort_values('sort_key').drop('sort_key', axis=1)
        return df2[cols]

    def create_fasta_from_df(self, df, output_fasta):
        unique_df = df.drop_duplicates(subset=['eName'], keep='first')
        records = []
        for _, row in unique_df.iterrows():
            if pd.notnull(row['eSeq']):
                rcd = SeqRecord(Seq(row['eSeq']), id=row['eName'], description="")
                records.append(rcd)
        SeqIO.write(records, output_fasta, "fasta")
        self.logger.info(f"Written {len(records)} sequences to FASTA.")
        if len(records) < len(df):
            self.logger.info(f"Removed {len(df) - len(records)} duplicates.")

    def write_empty_output(self):
        pd.DataFrame(columns=[
            'eName','MapNum','eRepeatNum','eState','eChr','eStart',
            'eEnd','eLength','MatDegree','eClass','eReadNum','eReads','eSeq'
        ]).to_csv(self.output_csv, index=False)
        with open(self.output_fasta, 'w'):
            pass

def main():
    parser = argparse.ArgumentParser(description="Force eSeq-based eLength,eEnd, approximate location by start//50.")
    parser.add_argument("--meccdna_part1", required=True)
    parser.add_argument("--meccdna_part2", required=True)
    parser.add_argument("--output_csv", required=True)
    parser.add_argument("--output_fasta", required=True)
    args = parser.parse_args()

    merger = MergeMecc(
        args.meccdna_part1,
        args.meccdna_part2,
        args.output_csv,
        args.output_fasta
    )
    merger.run_merge_mecc()

if __name__=="__main__":
    main()