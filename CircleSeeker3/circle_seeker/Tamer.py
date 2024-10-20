#!/usr/bin/env python3
# coding: utf-8

import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import logging

class Tamer:
    def __init__(self, input_tecc_csv, input_fasta, output_uecc, output_mecc):
        self.input_csv = input_tecc_csv
        self.input_fasta = input_fasta
        self.output_uecc = output_uecc
        self.output_mecc = output_mecc
        self.df = None
        self.fasta_dict = None
        
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    def read_fasta(self):
        self.fasta_dict = SeqIO.to_dict(SeqIO.parse(self.input_fasta, "fasta"))
        logging.info(f"Read FASTA file: {self.input_fasta}")

    def write_fasta(self, sequences, file_path):
        with open(file_path, 'w') as f:
            SeqIO.write(sequences, f, "fasta")
        logging.info(f"Wrote {len(sequences)} sequences to {file_path}")

    def process_uecc(self):
        uecc_sequences = []
        for _, row in self.df[self.df['eClass'] == 'Uecc'].iterrows():
            if row['qname'] in self.fasta_dict:
                seq = self.fasta_dict[row['qname']][int(row['qstart']):int(row['qend'])]
                seq_id = f"{row['eChr']}-{row['eStart']}-{row['eEnd']}"
                uecc_sequences.append(SeqRecord(seq.seq, id=seq_id, description=""))
            else:
                logging.warning(f"Sequence {row['qname']} not found in FASTA file. Skipping.")
        return uecc_sequences

    def process_mecc(self):
        mecc_sequences = []
        mecc_df = self.df[self.df['eClass'] == 'Mecc']
        
        for qname, group in mecc_df.groupby('qname'):
            if qname in self.fasta_dict:
                longest_row = group.loc[group['eLength'].idxmax()]
                seq = self.fasta_dict[qname][int(longest_row['qstart']):int(longest_row['qend'])]
                seq_id = f"{qname}|Second|{longest_row['eLength']}|{longest_row['eRepeatNum']}|circular"
                mecc_sequences.append(SeqRecord(seq.seq, id=seq_id, description=""))
            else:
                logging.warning(f"Sequence {qname} not found in FASTA file. Skipping.")
        return mecc_sequences

    def run(self):
        try:
            self.df = pd.read_csv(self.input_csv)
            logging.info(f"Read CSV file: {self.input_csv}")
            
            self.read_fasta()
            
            uecc_sequences = self.process_uecc()
            self.write_fasta(uecc_sequences, self.output_uecc)
            
            mecc_sequences = self.process_mecc()
            self.write_fasta(mecc_sequences, self.output_mecc)

        except Exception as e:
            logging.error(f"An error occurred: {e}")

