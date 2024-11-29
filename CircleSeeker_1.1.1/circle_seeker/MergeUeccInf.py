#!/usr/bin/env python3
# coding: utf-8

import pandas as pd
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class UeccProcessor:
    def __init__(self, uecc_csv, uecc_output):
        self.uecc_csv = uecc_csv
        self.uecc_output = uecc_output
        
    def process_uecc(self):
        logging.info("Processing Uecc data...")
        
        # Read CSV file
        df = pd.read_csv(self.uecc_csv)
        
        # Select columns and create a new DataFrame
        columns = ['Chromosome', 'Start', 'End', 'Coverage', 'Num_Reads', 
                  'Avg_Read_Length', 'Num_Split_Reads', 'Avg_Split_Read_Length', 'Score']
        df_processed = df[columns].copy()
        
        # Rename columns
        df_processed.columns = ['eChr', 'eStart', 'eEnd', 'Coverage', 'Num_Reads',
                              'Avg_Read_Length', 'Num_Split_Reads', 'Avg_Split_Read_Length', 'Score']
        
        # Add required columns
        df_processed['eLength'] = df_processed['eEnd'] - df_processed['eStart'] + 1
        df_processed['eName'] = df_processed.apply(
            lambda row: f"Uecc|{row['eChr']}-{row['eStart']}-{row['eEnd']}", axis=1
        )
        df_processed['eClass'] = 'Uecc'
        df_processed['eState'] = 'Inferred-eccDNA'
        
        # Add filter condition: eLength >= 100
        df_processed = df_processed[df_processed['eLength'] >= 100]
        
        # Select and order columns
        columns_order = ['eName', 'eChr', 'eStart', 'eEnd', 'eLength', 'Num_Reads',
                        'Avg_Read_Length', 'Num_Split_Reads', 'Avg_Split_Read_Length', 'Score']
        df_output = df_processed[columns_order]
        df_output.to_csv(self.uecc_output, index=False)
        
        logging.info(f"Uecc processing complete. Output file: {self.uecc_output}")
        return self.uecc_output

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Process Uecc DNA data.")
    parser.add_argument("--uecc_csv", required=True, 
                        help="Path to circular_dna_results_scored.csv file")
    parser.add_argument("--output", required=True, 
                        help="Path to output file")
    args = parser.parse_args()

    processor = UeccProcessor(args.uecc_csv, args.output)
    output_file = processor.process_uecc()
    print(f"Processing complete. Results saved to: {output_file}")
