import argparse
import time
import pandas as pd
from itertools import combinations
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import time
from tqdm import tqdm
import logging

class Carousel:
    def __init__(self, input_file, output_file, read_list_file, circular_fasta_file):
        self.input_file = input_file
        self.output_file = output_file
        self.read_list_file = read_list_file
        self.circular_fasta_file = circular_fasta_file
        # Setup logging
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        self.logger = logging.getLogger(__name__)

    def run(self):
        self.process()

    def read_and_preprocess_data(self):
        self.logger.info("Reading and preprocessing data")
        # Define column names
        columns = [
            'readName', 'repN', 'copyNum', 'readLen', 'start', 'end', 
            'consLen', 'aveMatch', 'fullLen', 'subPos', 'consSeq'
        ]
        
        # Read CSV file and preprocess data
        df = pd.read_csv(self.input_file, sep='\t', header=None, names=columns)
                
        # Convert 'subPos' to string and replace commas with pipes
        df['subPos'] = df['subPos'].astype(str).replace(',', '|', regex=True)
        
        # Calculate effective length percentage
        df['Effective_Length'] = ((df['end'] - df['start'] + 1) / df['readLen'] * 100).round(2)
        
        # Filter rows with average match > 99
        df = df[df['aveMatch'] > 99]
        
        # Remove fullLen column
        df = df.drop(columns=['fullLen'])
        
        self.logger.info(f"Filtering reads with more than 10 repeated regions")
        # Group by readName and filter out groups with more than 5 items
        group_sizes = df.groupby('readName').size()
        total_groups_before = len(group_sizes)
        valid_readnames = group_sizes[group_sizes <= 10].index
        df = df[df['readName'].isin(valid_readnames)]
        
        total_groups_after = len(df['readName'].unique())
        
        # Prepare and print the filtered report
        filtered_report = f"""
        Filtering Report:
        - Initial reads: {total_groups_before}
        - Reads removed: {total_groups_before - total_groups_after}
        - Final reads: {total_groups_after}
        """
        self.logger.info(filtered_report)
        
        return df

    # OLD VERSION v1.0.2
    # def read_and_preprocess_data(self):
    #     # Define column names
    #     columns = [
    #         'readName', 'repN', 'copyNum', 'readLen', 'start', 'end', 
    #         'consLen', 'aveMatch', 'fullLen', 'subPos', 'consSeq'
    #     ]
        
    #     # Read CSV file and preprocess data
    #     df = pd.read_csv(self.input_file, sep='\t', header=None, names=columns)
    #     df['subPos'] = df['subPos'].str.replace(',', '|')  # Replace commas with pipes in subPos column
    #     df['Effective_Length'] = ((df['end'] - df['start'] + 1) / df['readLen'] * 100).round(2)  # Calculate effective length percentage
    #     df = df[df['aveMatch'] > 99]  # Filter rows with average match > 99
    #     df = df.drop(columns=['fullLen'])  # Remove fullLen column
        
    #     return df
    
    @staticmethod
    def combine_rows_data(rows):
        # Combine multiple rows into one
        new_row = rows.iloc[0].copy()
        new_row['repN'] = 'repM0'
        new_row['copyNum'] = rows['copyNum'].sum()
        new_row['start'] = rows['start'].min()
        new_row['end'] = rows['end'].max()
        min_start_row = rows.loc[rows['start'].idxmin()]
        new_row['consLen'] = min_start_row['consLen']
        new_row['aveMatch'] = min_start_row['aveMatch']
        new_row['subPos'] = '|'.join(rows['subPos'])
        new_row['consSeq'] = min_start_row['consSeq']
        new_row['Effective_Length'] = min(100, rows['Effective_Length'].sum())
        return new_row

    @staticmethod
    def combine_rows(group):
        # Combine rows within a group
        group = group.reset_index(drop=True)
        
        if len(group) == 2:
            # Handle the case with only two rows
            cons_len_diff = abs(group['consLen'].iloc[0] - group['consLen'].iloc[1])
            if cons_len_diff <= 15:
                return pd.DataFrame([Carousel.combine_rows_data(group)])
            else:
                total_effective_length = group['Effective_Length'].sum()
                max_single_effective_length = group['Effective_Length'].max()
                
                if abs(total_effective_length - 100) <= abs(max_single_effective_length - 100):
                    return group
                else:
                    return group[group['Effective_Length'] == max_single_effective_length]
        
        if group['Effective_Length'].sum() <= 102:
            return pd.DataFrame([Carousel.combine_rows_data(group)])
        
        # Find the best combination
        best_combination = None
        best_total = 0
        for r in range(1, len(group) + 1):
            for combo in combinations(range(len(group)), r):
                subset = group.iloc[list(combo)]
                subset_total = subset['Effective_Length'].sum()
                if 98 <= subset_total <= 102:
                    cons_len_range = subset['consLen'].max() - subset['consLen'].min()
                    if cons_len_range <= 15:
                        return pd.DataFrame([Carousel.combine_rows_data(subset)])
                elif subset_total < 98 and subset_total > best_total:
                    best_combination = subset
                    best_total = subset_total
        
        if best_combination is not None:
            return pd.DataFrame([Carousel.combine_rows_data(best_combination)])
        return group


    def process_repeated_readnames(self, df):
         # Find repeated readNames
        repeated_readnames = df['readName'].value_counts()[df['readName'].value_counts() > 1].index
        repeated_df = df[df['readName'].isin(repeated_readnames)]
        self.logger.info(f"Processing {len(repeated_readnames)} repeated readNames")
        
        # Process repeated readNames
        processed_repeated = []
        
        # Create a tqdm progress bar
        with tqdm(total=len(repeated_readnames), desc="Processing repeated readNames") as pbar:
            for name, group in repeated_df.groupby('readName'):
                processed_group = self.combine_rows(group)
                processed_repeated.append(processed_group)
                pbar.update(1)
        
        processed_repeated = pd.concat(processed_repeated, ignore_index=True)

        # Combine processed repeated and single occurrence readNames
        single_occurrence = df[~df['readName'].isin(repeated_readnames)]
        result_df = pd.concat([single_occurrence, processed_repeated], ignore_index=True)

        # print("process_repeated_readnames completed")
        return result_df

    @staticmethod
    def classify_readnames(df):
        # Classify readNames
        read_list_df = pd.DataFrame(columns=['readName', 'readClass'])
        
        repeated_readnames = df['readName'].value_counts()[df['readName'].value_counts() > 1].index
        read_list_df = pd.concat([read_list_df, 
                                  pd.DataFrame({'readName': repeated_readnames, 'readClass': 'CtcR-multiple'})], 
                                 ignore_index=True)
        
        single_occurrence = df[~df['readName'].isin(repeated_readnames)]
        single_repM = single_occurrence[single_occurrence['repN'].str.startswith('repM', na=False)]
        read_list_df = pd.concat([read_list_df, 
                                  pd.DataFrame({'readName': single_repM['readName'], 'readClass': 'CtcR-inversion'})], 
                                 ignore_index=True)
        
        remaining = df[~df['readName'].isin(read_list_df['readName'])]
        
        ctcr_p = remaining[remaining['Effective_Length'] >= 99]
        read_list_df = pd.concat([read_list_df, 
                                  pd.DataFrame({'readName': ctcr_p['readName'], 'readClass': 'CtcR-perfect'})], 
                                 ignore_index=True)
        
        ctcr_i = remaining[(remaining['Effective_Length'] >= 70) & (remaining['Effective_Length'] < 99)]
        read_list_df = pd.concat([read_list_df, 
                                  pd.DataFrame({'readName': ctcr_i['readName'], 'readClass': 'CtcR-hybrid'})], 
                                 ignore_index=True)
        
        return read_list_df

    @staticmethod
    def circularize_sequences(sequences):
        # Circularize sequences
        circular_sequences = []
        for seq in sequences:
            circular_seq = seq.seq + seq.seq
            circular_record = SeqRecord(circular_seq, 
                                        id=seq.id + "|circular", 
                                        description="")
            circular_sequences.append(circular_record)
        # print(f"Circularized {len(circular_sequences)} sequences")
        return circular_sequences

    @staticmethod
    def write_fasta(sequences, output_file):
        # Write FASTA file
        SeqIO.write(sequences, output_file, "fasta")
        # print(f"Wrote {len(sequences)} sequences to {output_file}")

    def process(self):
        start_time = time.time()
        self.logger.info("Starting Carousel processing")
        
        # Read and preprocess data
        start_time = time.time()
        df = self.read_and_preprocess_data()
        # print(f"Read and preprocess data completed in {time.time() - start_time:.2f} seconds")
        
        # Process repeated readNames
        start_time = time.time()
        result_df = self.process_repeated_readnames(df)
        # print(f"Process repeated readNames completed in {time.time() - start_time:.2f} seconds")
        
        # Discard decimal part of copyNum
        start_time = time.time()
        result_df['copyNum'] = result_df['copyNum'].astype(int)
        # print(f"Discard decimal part of copyNum completed in {time.time() - start_time:.2f} seconds")
        
        # Add new column and adjust it to be the first column, including copyNum
        start_time = time.time()
        result_df['readName|repN|consLen|copyNum'] = (
            result_df['readName'] + '|' + 
            result_df['repN'] + '|' + 
            result_df['consLen'].astype(str) + '|' + 
            result_df['copyNum'].astype(str)
        )
        columns = ['readName|repN|consLen|copyNum'] + [col for col in result_df.columns if col != 'readName|repN|consLen|copyNum']
        result_df = result_df[columns]
        # print(f"Add new column and adjust columns completed in {time.time() - start_time:.2f} seconds")
        
        # Save processed results
        start_time = time.time()
        result_df.to_csv(self.output_file, index=False)
        # print(f"Save processed results completed in {time.time() - start_time:.2f} seconds")
        
        # Classify readNames
        start_time = time.time()
        read_list_df = self.classify_readnames(result_df)
        # print(f"Classify readNames completed in {time.time() - start_time:.2f} seconds")
        
        # Save read_list results
        start_time = time.time()
        read_list_df.to_csv(self.read_list_file, index=False)
        # print(f"Save read_list results completed in {time.time() - start_time:.2f} seconds")
        
        # Generate circularized FASTA file
        start_time = time.time()
        sequences = [SeqRecord(Seq(row['consSeq']), id=row['readName|repN|consLen|copyNum'], description="") for _, row in result_df.iterrows()]
        circular_sequences = self.circularize_sequences(sequences)
        # print(f"Generate circularized sequences completed in {time.time() - start_time:.2f} seconds")
        
        # Write circularized sequences
        start_time = time.time()
        self.write_fasta(circular_sequences, self.circular_fasta_file)
        # print(f"Write circularized sequences completed in {time.time() - start_time:.2f} seconds")
        
        # Print processing result information
        print(f"Carousel Processing complete")
        # print(f"Total reads: {len(result_df)}")
        # print(f"Number of classified reads: {len(read_list_df)}")
        
        print(f"Total processing time: {time.time() - start_time:.2f} seconds")
        
        self.logger.info(f"Processing completed in {time.time() - start_time:.2f} seconds")
        
def main():
    parser = argparse.ArgumentParser(description="Process eccDNA data using Carousel")
    parser.add_argument("-i", "--input", required=True, help="Input CSV file")
    parser.add_argument("-o", "--output", required=True, help="Output CSV file")
    parser.add_argument("-rl", "--read_list", required=True, help="Output read list CSV file")
    parser.add_argument("-cf", "--circular_fasta", required=True, help="Output circular FASTA file")

    args = parser.parse_args()

    carousel = Carousel(args.input, args.output, args.read_list, args.circular_fasta)
    carousel.process()

if __name__ == "__main__":
    main()
