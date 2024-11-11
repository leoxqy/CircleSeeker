#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import argparse
import logging

class Juggler:
    def __init__(self, input_blast, input_fai, output_num_linr, output_tecc, output_read_class):
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        self.logger = logging.getLogger(__name__)
        
        self.input_blast = input_blast
        self.input_fai = input_fai
        self.output_num_linr = output_num_linr
        self.output_tecc = output_tecc
        self.output_read_class = output_read_class

    def read_and_merge_data(self):
        self.logger.info("Reading and merging input files")
        # Read BLAST results file
        df1 = pd.read_csv(self.input_blast, sep=r'\s+', header=None,
                          names=['qname', 'tname', 'identity', 'alnlen', 'mismatch', 'gapopen',
                                 'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits'])

        # Read FASTA index file
        df2 = pd.read_csv(self.input_fai, sep=r'\s+', header=None,
                          names=['qname', 'length', 'offset', 'linebases', 'linewidth'])

        # Merge two DataFrames using 'qname' as the merge key
        merged_df = pd.merge(df1, df2, on='qname', how='left')

        # Add Strand column
        merged_df['Strand'] = np.where(merged_df['tstart'] < merged_df['tend'], '+', '-')

        # Ensure tstart < tend
        merged_df.loc[merged_df['tstart'] > merged_df['tend'], ['tstart', 'tend']] = \
            merged_df.loc[merged_df['tstart'] > merged_df['tend'], ['tend', 'tstart']].values

        # Calculate eff_Ratio and round to two decimal places
        merged_df['eff_Ratio'] = (merged_df['alnlen'] / merged_df['length'] * 100).round(2)

        # Count qname frequency
        qname_counts = merged_df['qname'].value_counts().reset_index()
        qname_counts.columns = ['qname', 'frequency']

        # Merge frequency information back to the original DataFrame
        merged_df = pd.merge(merged_df, qname_counts, on='qname', how='left')

        self.logger.info(f"Processed {len(merged_df)} records")
        return merged_df

    def process_single_occurrence(self, merged_df):
        # Filter rows where qname occurs only once
        single_occurrence = merged_df[merged_df['frequency'] == 1]

        # In rows that occur only once, retain rows where eff_Ratio > 98
        filtered_single_df = single_occurrence[single_occurrence['eff_Ratio'] > 98]

        # Calculate the number of retained rows
        kept_rows = len(filtered_single_df)

        # Output Num_LinR to CSV file
        num_linr_df = pd.DataFrame({'Num_LinR': [kept_rows]})
        num_linr_df.to_csv(self.output_num_linr, index=False)

        return filtered_single_df

    def calculate_score(self, row1, row2):
        """Calculate the score between two rows"""
        q_diff = min(abs(row1['qstart'] - row2['qstart']), abs(row1['qend'] - row2['qend']))
        subject_score = 100 if row1['tname'] == row2['tname'] else 0
        identity_score = min(row1['identity'], row2['identity'])
        s_diff = min(abs(row1['tstart'] - row2['tstart']), abs(row1['tend'] - row2['tend']))
        s_diff_score = 1 / (s_diff + 1)  # Avoid division by zero

        # Comprehensive score, weights can be adjusted as needed
        return q_diff * 0.4 + subject_score * 0.3 + identity_score * 0.2 + s_diff_score * 0.1

    def analyze_group(self, group):
        """Analyze each qname group, build alignment chains"""
        group = group.sort_values('qstart')
        PreCtcRlists = []
        processed = set()
        list_number = 1

        for i, row in group.iterrows():
            if i in processed:
                continue

            current_chain = [row]
            processed.add(i)

            while True:
                potential_next = group[
                    (group['qstart'] >= current_chain[-1]['qend'] - 10) &
                    (group['qstart'] <= current_chain[-1]['qend'] + 10) &
                    (~group.index.isin(processed))
                ]

                if potential_next.empty:
                    break

                if len(potential_next) > 1:
                    scores = potential_next.apply(lambda x: self.calculate_score(current_chain[-1], x), axis=1)
                    next_row = potential_next.loc[scores.idxmax()]
                else:
                    next_row = potential_next.iloc[0]

                current_chain.append(next_row)
                processed.add(next_row.name)

            # Add PreCtcRlist number to each row in the current chain
            chain_df = pd.DataFrame(current_chain)
            chain_df['PreCtcRlist_number'] = list_number
            PreCtcRlists.append(chain_df)
            list_number += 1

        return pd.concat(PreCtcRlists, ignore_index=True)

    def process_multi_occurrence(self, merged_df):
        # Filter rows where qname occurs more than 2 times and less than 100 times
        multi_occurrence_df = merged_df[(merged_df['frequency'] > 2) & (merged_df['frequency'] < 100)]

        # Group by qname and apply analyze_group function
        prectcrlist_df = multi_occurrence_df.groupby('qname', group_keys=False).apply(self.analyze_group).reset_index(drop=True)

        return prectcrlist_df

    def merge_overlapping_lists(self, prectcrlist_df):
        def merge_lists(group):
            prectcrlists = group.groupby('PreCtcRlist_number')
            merged = []
            merged_numbers = set()

            for num1, list1 in prectcrlists:
                if num1 in merged_numbers or len(list1) < 3:
                    continue

                merged_list = list1.copy()
                merged_numbers.add(num1)

                for num2, list2 in prectcrlists:
                    if num2 in merged_numbers or len(list2) < 3:
                        continue

                    if set(list1['tstart']).intersection(set(list2['tstart'])) or \
                       set(list1['tend']).intersection(set(list2['tend'])):
                        merged_list = pd.concat([merged_list, list2])
                        merged_numbers.add(num2)

                if len(merged_list) >= 3:
                    merged_list['PreCtcRlist_number'] = f"{num1}m" if len(merged_numbers) > 1 else num1
                    merged.append(merged_list)

            return pd.concat(merged) if merged else pd.DataFrame()

        # First, filter out all PreCtcRlist shorter than 3 rows
        prectcrlist_df = prectcrlist_df.groupby(['qname', 'PreCtcRlist_number']).filter(lambda x: len(x) >= 3)

        # Filter out qname with more than one PreCtcRlist_number
        multi_list_qnames = prectcrlist_df.groupby('qname')['PreCtcRlist_number'].nunique()
        multi_list_qnames = multi_list_qnames[multi_list_qnames > 1].index

        # Process these qname
        processed_df = prectcrlist_df[prectcrlist_df['qname'].isin(multi_list_qnames)].groupby('qname').apply(merge_lists)
        processed_df = processed_df.reset_index(drop=True)

        # Add rows with only one PreCtcRlist_number back to the result
        single_list_df = prectcrlist_df[~prectcrlist_df['qname'].isin(multi_list_qnames)]
        final_prectcrlist_df = pd.concat([processed_df, single_list_df]).reset_index(drop=True)

        return final_prectcrlist_df

    def process_group(self, group):
        # Ensure tname is consistent
        if len(group['tname'].unique()) > 1:
            return None

        # Get original sequence length
        original_length = group['length'].iloc[0]

        # Check for duplicates in tstart and tend columns
        def check_repeats(series, tolerance=30):
            sorted_vals = sorted(series)
            repeats = []
            for i in range(len(sorted_vals) - 1):
                if abs(sorted_vals[i] - sorted_vals[i + 1]) <= tolerance:
                    repeats.extend([sorted_vals[i], sorted_vals[i + 1]])
            return list(set(repeats))

        tstart_repeats = check_repeats(group['tstart'])
        tend_repeats = check_repeats(group['tend'])

        # Mark rows with duplicates
        group['Repeat1'] = group['tstart'].apply(lambda x: x in tstart_repeats)
        group['Repeat2'] = group['tend'].apply(lambda x: x in tend_repeats)

        # Check if there are rows marked as both Repeat1 and Repeat2
        repeat_rows = group[group['Repeat1'] & group['Repeat2']]
        if not repeat_rows.empty:
            # Calculate eLength
            eLength = max(repeat_rows['tend']) - min(repeat_rows['tstart']) + 1

            # Check if eLength is greater than original length
            if eLength > original_length:
                return None  # Return None if eLength is greater than original length

            # Calculate the number of repeats for Repeat1 and Repeat2
            repeat1_count = sum(group['Repeat1'])
            repeat2_count = sum(group['Repeat2'])

            # Take the minimum of Repeat1 and Repeat2 as eRepeatNum
            eRepeatNum = min(repeat1_count, repeat2_count)

            # Calculate the most common Strand value
            most_common_strand = group['Strand'].mode().iloc[0]

            # Check if Strand is consistent
            inversion_flag = 1 if group['Strand'].nunique() > 1 else 0

            # Get qstart and qend from the first row marked as both Repeat1 and Repeat2
            first_repeat_row = repeat_rows.iloc[0]
            qstart = first_repeat_row['qstart']
            qend = first_repeat_row['qend']

            # Return a DataFrame
            return pd.DataFrame({
                'qname': [group['qname'].iloc[0]],
                'PreCtcRlist_number': [group['PreCtcRlist_number'].iloc[0]],
                'eChr': [group['tname'].iloc[0]],
                'eStart': [min(repeat_rows['tstart'])],
                'eEnd': [max(repeat_rows['tend'])],
                'eStrand': [most_common_strand],
                'eLength': [eLength],
                'eRepeatNum': [eRepeatNum],
                'MatDegree': [round(repeat_rows['identity'].mean(), 2)],
                'eClass': ['Tecc'],
                'inversion': [inversion_flag],
                'eReads': [group['qname'].iloc[0]],
                'qstart': [qstart],  # New column
                'qend': [qend]  # New column
            })
        return None

    def process_final_prectcrlist(self, final_prectcrlist_df):
        # Group by qname and PreCtcRlist_number and apply processing function
        result_dfs = []
        for (qname, prectc), group in final_prectcrlist_df.groupby(['qname', 'PreCtcRlist_number']):
            result = self.process_group(group)
            if result is not None:
                result_dfs.append(result)

        # Merge all results
        if result_dfs:
            output_df = pd.concat(result_dfs, ignore_index=True)
        else:
            output_df = pd.DataFrame()

        return output_df

    def generate_read_classification(self, output_df):
        read_classification = []
        grouped = output_df.groupby('qname')

        for qname, group in grouped:
            prectc_numbers = group['PreCtcRlist_number'].unique()
            eLengths = group['eLength'].unique()

            if len(prectc_numbers) > 1:
                # Check if eLength is within 50bp
                max_length = max(eLengths)
                min_length = min(eLengths)
                if max_length - min_length <= 50:
                    read_classification.append({'readName': qname, 'readClass': 'CtcR-perfect'})
                    output_df.loc[output_df['qname'] == qname, 'eClass'] = 'Mecc'
                else:
                    read_classification.append({'readName': qname, 'readClass': 'CtcR-multiple'})
                    output_df.loc[output_df['qname'] == qname, 'eClass'] = 'Uecc'
            elif any('m' in str(num) for num in prectc_numbers):
                read_classification.append({'readName': qname, 'readClass': 'CtcR-hybrid'})
                output_df.loc[output_df['qname'] == qname, 'eClass'] = 'Uecc'
            else:  # Only one PreCtcRlist_number
                inversion = group['inversion'].iloc[0]
                if inversion == 0:
                    read_classification.append({'readName': qname, 'readClass': 'CtcR-perfect'})
                else:
                    read_classification.append({'readName': qname, 'readClass': 'CtcR-inversion'})
                output_df.loc[output_df['qname'] == qname, 'eClass'] = 'Uecc'

        # Convert read_classification to DataFrame
        read_df = pd.DataFrame(read_classification)

        # Save Read_df to CSV file, using comma as separator
        read_df.to_csv(self.output_read_class, index=False)

        return output_df

    def filter_mecc_groups(self, output_df):
        # Filter rows where eClass is Mecc
        mecc_df = output_df[output_df['eClass'] == 'Mecc']
        
        # Group by qname
        grouped = mecc_df.groupby('qname')
        
        # Store groups to keep
        keep_groups = []
        
        for qname, group in grouped:
            eRepeatNums = group['eRepeatNum'].unique()
            if max(eRepeatNums) - min(eRepeatNums) < 2:
                keep_groups.append(qname)
        
        # Keep Mecc groups that meet the condition and all non-Mecc rows
        filtered_df = pd.concat([
            output_df[output_df['qname'].isin(keep_groups)],
            output_df[output_df['eClass'] != 'Mecc']
        ])
        
        return filtered_df

    def run(self):
        self.logger.info("Starting Juggler pipeline")
        
        merged_df = self.read_and_merge_data()
        filtered_single_df = self.process_single_occurrence(merged_df)
        self.logger.info(f"Found {len(filtered_single_df)} single occurrence sequences")

        prectcrlist_df = self.process_multi_occurrence(merged_df)
        final_prectcrlist_df = self.merge_overlapping_lists(prectcrlist_df)
        output_df = self.process_final_prectcrlist(final_prectcrlist_df)
        output_df = self.generate_read_classification(output_df)
        output_df = self.filter_mecc_groups(output_df)
        
        output_df.to_csv(self.output_tecc, index=False)
        self.logger.info("Pipeline completed successfully")
