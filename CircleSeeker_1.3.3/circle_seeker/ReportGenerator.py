#!/usr/bin/env python3
# coding: utf-8

"""
ReportGenerator: Comprehensive Report Generation for CircleSeeker Results.

This module handles the generation of detailed reports and summaries for circular DNA
analysis results. It processes various input files to create human-readable reports
in multiple formats (text, CSV, and HTML).

Features:
- Aggregates results from multiple analysis stages
- Generates statistical summaries of circular DNA findings
- Creates detailed HTML reports with interactive visualizations
- Supports multiple output formats for downstream analysis

Typical usage:
    generator = ReportGenerator(fai, ctcr1, ctcr2, output="report.txt")
    generator.generate_report()

Version: 1.3.3
License: GNU General Public License v2
Copyright (c) 2024 CircleSeeker Team
"""

import pandas as pd
import logging
import argparse
import subprocess
import tempfile
import os
import sys
import shutil
from typing import List, Dict
from collections import Counter

class ReportGenerator:
    def __init__(self, fai, ctcr1, ctcr2, linr, uecc_fai, mecc_fai, cecc_fai=None, mcecc_fai=None,
                 uecc_inferred_fai=None, shared_count_csv=None,
                 output="report.txt", summary_output="summary.csv", html_output="report.html",
                 sample_name="Unknown", keep_tmp=False):
        self.fai = fai
        self.ctcr1 = ctcr1
        self.ctcr2 = ctcr2
        self.linr = linr
        self.uecc_fai = uecc_fai
        self.mecc_fai = mecc_fai
        self.cecc_fai = cecc_fai
        self.mcecc_fai = mcecc_fai
        self.uecc_inferred_fai = uecc_inferred_fai
        self.shared_count_csv = shared_count_csv
        self.output = output
        self.summary_output = summary_output
        self.html_output = html_output
        self.sample_name = sample_name
        self.keep_tmp = keep_tmp

        # Set log level based on keep_tmp parameter
        log_level = logging.DEBUG if self.keep_tmp else logging.INFO
        logging.basicConfig(
            level=log_level,
            format='%(asctime)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        self.logger = logging.getLogger(__name__)
        
        self.logger.info("Initialized ReportGenerator")

    def count_reads(self):
        try:
            with open(self.fai, 'r') as f:
                return sum(1 for line in f)
        except Exception as e:
            self.logger.error(f"Error reading FAI file: {str(e)}")
            return 0

    def count_ctcr(self):
        try:
            df1 = pd.read_csv(self.ctcr1)
            df2 = pd.read_csv(self.ctcr2)

            if 'readName' not in df1.columns and 'qName' in df1.columns:
                df1 = df1.rename(columns={'qName': 'readName'})
            if 'readName' not in df2.columns and 'qName' in df2.columns:
                df2 = df2.rename(columns={'qName': 'readName'})

            combined = pd.concat([df1[['readName']], df2[['readName']]], ignore_index=True)
            return len(combined['readName'].unique())
        except Exception as e:
            self.logger.error(f"Error counting CtcR: {str(e)}")
            return 0

    def classify_ctcr(self):
        try:
            df1 = pd.read_csv(self.ctcr1)
            df2 = pd.read_csv(self.ctcr2)
            
            if 'readClass' not in df1.columns and 'class' in df1.columns:
                df1 = df1.rename(columns={'class': 'readClass'})
            if 'readClass' not in df2.columns and 'class' in df2.columns:
                df2 = df2.rename(columns={'class': 'readClass'})
            
            combined = pd.concat([df1[['readClass']], df2[['readClass']]], ignore_index=True)
            class_counts = Counter(combined['readClass'])

            expected_classes = ['CtcR-multiple', 'CtcR-inversion', 'CtcR-perfect', 'CtcR-hybrid']
            for class_name in expected_classes:
                if class_name not in class_counts:
                    class_counts[class_name] = 0
            return class_counts
        except Exception as e:
            self.logger.error(f"Error classifying CtcR: {str(e)}")
            return Counter()

    def count_linr(self):
        try:
            df = pd.read_csv(self.linr)
            if 'Num_LinR' in df.columns:
                return df['Num_LinR'].iloc[0]
            elif 'num_LinR' in df.columns:
                return df['num_LinR'].iloc[0]
            else:
                self.logger.warning("Warning: Column 'Num_LinR' not found in LinR file")
                return 0
        except Exception as e:
            self.logger.error(f"Error counting LinR: {str(e)}")
            return 0

    def read_fai_info(self, fai_file):
        if fai_file is None or not os.path.exists(fai_file):
            return pd.DataFrame()

        names = []
        lengths = []
        with open(fai_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    names.append(parts[0])
                    lengths.append(int(parts[1]))

        return pd.DataFrame({
            'Name': names,
            'Length': lengths
        })

    def read_fai_as_lengths(self, fai_file):
        if fai_file is None or not os.path.exists(fai_file):
            return pd.Series(dtype=int)
            
        lengths = []
        with open(fai_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    lengths.append(int(parts[1]))
        return pd.Series(lengths, dtype=int)

    def read_csv_if_exists(self, file_path):
        if file_path and os.path.exists(file_path):
            return pd.read_csv(file_path)
        return pd.DataFrame()

    def create_summary_csv(self):
        dataframes = []

        # Confirmed Uecc
        uecc_df = self.read_fai_info(self.uecc_fai)
        if not uecc_df.empty:
            uecc_df['Class'] = 'Uecc'
            uecc_df['Status'] = 'Confirmed'
            dataframes.append(uecc_df)

        # Confirmed Mecc
        mecc_df = self.read_fai_info(self.mecc_fai)
        if not mecc_df.empty:
            mecc_df['Class'] = 'Mecc'
            mecc_df['Status'] = 'Confirmed'
            dataframes.append(mecc_df)

        # Confirmed Cecc
        if self.cecc_fai:
            cecc_df = self.read_fai_info(self.cecc_fai)
            if not cecc_df.empty:
                cecc_df['Class'] = 'Cecc'
                cecc_df['Status'] = 'Confirmed'
                dataframes.append(cecc_df)

        # Confirmed MCecc
        if self.mcecc_fai:
            mcecc_df = self.read_fai_info(self.mcecc_fai)
            if not mcecc_df.empty:
                mcecc_df['Class'] = 'MCecc'
                mcecc_df['Status'] = 'Confirmed'
                dataframes.append(mcecc_df)

        # Inferred Uecc from FAI
        if self.uecc_inferred_fai:
            inferred_df = self.read_fai_info(self.uecc_inferred_fai)
            if not inferred_df.empty:
                inferred_df['Class'] = 'Uecc'
                inferred_df['Status'] = 'Inferred'
                dataframes.append(inferred_df)

        if dataframes:
            all_eccdna = pd.concat(dataframes, ignore_index=True)
            all_eccdna = all_eccdna[['Class', 'Name', 'Length', 'Status']]
            all_eccdna.to_csv(self.summary_output, index=False)
            self.logger.info(f"Summary CSV has been generated and saved to {self.summary_output}")
            return all_eccdna
        else:
            self.logger.warning("No data available for summary CSV")
            return pd.DataFrame(columns=['Class', 'Name', 'Length', 'Status'])

    def get_empty_stats(self):
        return {
            'count': 0,
            'min_length': 0,
            'max_length': 0,
            'mean_length': 0,
            'median_length': 0,
            'mode_length': 0
        }

    def analyze_eccdna_fai(self, fai_file):
        if fai_file is None:
            return self.get_empty_stats()

        lengths = self.read_fai_as_lengths(fai_file)
        if lengths.empty:
            return self.get_empty_stats()

        return {
            'count': len(lengths),
            'min_length': lengths.min(),
            'max_length': lengths.max(),
            'mean_length': lengths.mean(),
            'median_length': lengths.median(),
            'mode_length': lengths.mode().iloc[0] if not lengths.empty else 0
        }

    def generate_tabular_report(self):
        try:
            total_reads = self.count_reads()
            linr_count = self.count_linr()
            ctcr_count = self.count_ctcr()
            ctcr_classes = self.classify_ctcr()

            if total_reads == 0:
                self.logger.warning("Warning: No reads found in FAI file")

            other_reads = total_reads - linr_count - ctcr_count

            read_stats = pd.DataFrame({
                'Category': ['Total Reads', 'Linear Reads (LinR)', 'Concatemeric tandem copies Reads (CtcR)', 'Other Reads'],
                'Count': [total_reads, linr_count, ctcr_count, other_reads],
                'Percentage': [
                    '100%',
                    f'{linr_count/total_reads:.2%}' if total_reads>0 else '0%',
                    f'{ctcr_count/total_reads:.2%}' if total_reads>0 else '0%',
                    f'{other_reads/total_reads:.2%}' if total_reads>0 else '0%'
                ]
            })

            uecc = self.analyze_eccdna_fai(self.uecc_fai)
            mecc = self.analyze_eccdna_fai(self.mecc_fai)
            cecc = self.analyze_eccdna_fai(self.cecc_fai)
            mcecc = self.analyze_eccdna_fai(self.mcecc_fai)
            inferred = self.analyze_eccdna_fai(self.uecc_inferred_fai)

            confirmed_eccdna = uecc['count'] + mecc['count'] + cecc['count'] + mcecc['count']

            N = 0
            if self.shared_count_csv and os.path.exists(self.shared_count_csv):
                df_shared = pd.read_csv(self.shared_count_csv)
                if 'SharedCount' in df_shared.columns and not df_shared.empty:
                    N = df_shared['SharedCount'].iloc[0]

            total_inferred = inferred['count']
            total_eccdna = confirmed_eccdna + total_inferred - N

            X = 0.00
            if total_inferred > 0:
                X = (N / total_inferred)*100

            # Overlap information
            overlap_info = f"Between Confirmed-UeccDNA and Inferred-UeccDNA, there are {N} overlapping eccDNAs, accounting for {X:.2f}% of Inferred-UeccDNA."

            ctcr_stats = pd.DataFrame({
                'Category': ['CtcR-multiple', 'CtcR-inversion', 'CtcR-perfect', 'CtcR-hybrid'],
                'Count': [
                    ctcr_classes['CtcR-multiple'],
                    ctcr_classes['CtcR-inversion'],
                    ctcr_classes['CtcR-perfect'],
                    ctcr_classes['CtcR-hybrid']
                ],
                'Percentage': [
                    f'{ctcr_classes["CtcR-multiple"]/ctcr_count:.2%}' if ctcr_count>0 else '0%',
                    f'{ctcr_classes["CtcR-inversion"]/ctcr_count:.2%}' if ctcr_count>0 else '0%',
                    f'{ctcr_classes["CtcR-perfect"]/ctcr_count:.2%}' if ctcr_count>0 else '0%',
                    f'{ctcr_classes["CtcR-hybrid"]/ctcr_count:.2%}' if ctcr_count>0 else '0%'
                ]
            })

            eccdna_types = []
            eccdna_counts = []
            min_lengths = []
            max_lengths = []
            mean_lengths = []
            median_lengths = []
            mode_lengths = []

            def append_stats(name, stats):
                if stats['count'] > 0:
                    eccdna_types.append(name)
                    eccdna_counts.append(stats['count'])
                    min_lengths.append(stats['min_length'])
                    max_lengths.append(stats['max_length'])
                    mean_lengths.append(f"{stats['mean_length']:.2f}")
                    median_lengths.append(f"{stats['median_length']:.2f}")
                    mode_lengths.append(stats['mode_length'])

            append_stats('UeccDNA', uecc)
            append_stats('MeccDNA', mecc)
            append_stats('CeccDNA', cecc)
            append_stats('MCeccDNA', mcecc)
            append_stats('Inferred UeccDNA', inferred)

            eccdna_stats = pd.DataFrame({
                'Type': eccdna_types,
                'Count': eccdna_counts,
                'Min Length (bp)': min_lengths,
                'Max Length (bp)': max_lengths,
                'Mean Length (bp)': mean_lengths,
                'Median Length (bp)': median_lengths,
                'Mode Length (bp)': mode_lengths
            })

            all_confirmed_lengths = pd.concat([
                self.read_fai_as_lengths(self.uecc_fai),
                self.read_fai_as_lengths(self.mecc_fai),
                self.read_fai_as_lengths(self.cecc_fai) if self.cecc_fai else pd.Series(dtype=int),
                self.read_fai_as_lengths(self.mcecc_fai) if self.mcecc_fai else pd.Series(dtype=int)
            ])

            inferred_lengths = self.read_fai_as_lengths(self.uecc_inferred_fai) if self.uecc_inferred_fai else pd.Series(dtype=int)
            all_eccdna_lengths = pd.concat([all_confirmed_lengths, inferred_lengths])

            overall_stats = pd.DataFrame({
                'Category': ['All Confirmed eccDNA', 'All eccDNA (Confirmed + Inferred)'],
                'Count': [confirmed_eccdna, total_eccdna],
                'Min Length (bp)': [
                    all_confirmed_lengths.min() if not all_confirmed_lengths.empty else 0,
                    all_eccdna_lengths.min() if not all_eccdna_lengths.empty else 0
                ],
                'Max Length (bp)': [
                    all_confirmed_lengths.max() if not all_confirmed_lengths.empty else 0,
                    all_eccdna_lengths.max() if not all_eccdna_lengths.empty else 0
                ],
                'Mean Length (bp)': [
                    f"{all_confirmed_lengths.mean():.2f}" if not all_confirmed_lengths.empty else "0.00",
                    f"{all_eccdna_lengths.mean():.2f}" if not all_eccdna_lengths.empty else "0.00"
                ],
                'Median Length (bp)': [
                    f"{all_confirmed_lengths.median():.2f}" if not all_confirmed_lengths.empty else "0.00",
                    f"{all_eccdna_lengths.median():.2f}" if not all_eccdna_lengths.empty else "0.00"
                ],
                'Mode Length (bp)': [
                    all_confirmed_lengths.mode().iloc[0] if not all_confirmed_lengths.empty else 0,
                    all_eccdna_lengths.mode().iloc[0] if not all_eccdna_lengths.empty else 0
                ]
            })

            return {
                'read_stats': read_stats,
                'ctcr_stats': ctcr_stats,
                'eccdna_stats': eccdna_stats,
                'overall_stats': overall_stats,
                'overlap_info': overlap_info
            }

        except Exception as e:
            self.logger.error(f"Error generating report: {str(e)}")
            raise

    def save_text_report(self, stats):
        with open(self.output, 'w') as f:
            f.write(f"eccDNA Analysis Report of {self.sample_name} Sample\n")
            f.write("=====================\n\n")

            f.write("1. Read Statistics\n")
            f.write(stats['read_stats'].to_string(index=False))
            f.write("\n\n")

            f.write("2. CtcR Classification\n")
            f.write(stats['ctcr_stats'].to_string(index=False))
            f.write("\n\n")

            if not stats['eccdna_stats'].empty:
                f.write("3. eccDNA Statistics\n")
                f.write(stats['eccdna_stats'].to_string(index=False))
                f.write("\n\n")

            f.write("4. Overall Statistics\n")
            f.write(stats['overall_stats'].to_string(index=False))
            f.write("\n\n")

            # Overlap information
            f.write("5. Overlap Information\n")
            f.write(stats['overlap_info'] + "\n")

        self.logger.info(f"Text report has been generated and saved to {self.output}")

    def save_html_report(self, stats):
        html_template = """<!DOCTYPE html>
        <html>
        <head>
            <title>eccDNA Analysis Report of {sample_name} Sample</title>
            <style>
                body {{
                    font-family: Arial, sans-serif;
                    margin: 20px;
                    color: #333;
                    max-width: 1200px;
                    margin: auto;
                    padding: 20px;
                }}
                h1 {{
                    color: #2c3e50;
                    border-bottom: 2px solid #2c3e50;
                    padding-bottom: 10px;
                    text-align: center;
                }}
                h2 {{
                    color: #34495e;
                    margin-top: 30px;
                    background: #f8f9fa;
                    padding: 10px;
                    border-radius: 5px;
                }}
                table {{
                    border-collapse: collapse;
                    width: 100%;
                    margin: 20px 0;
                    box-shadow: 0 1px 3px rgba(0,0,0,0.2);
                }}
                th, td {{
                    border: 1px solid #ddd;
                    padding: 12px;
                    text-align: left;
                }}
                th {{
                    background-color: #f5f6fa;
                    font-weight: bold;
                }}
                tr:nth-child(even) {{
                    background-color: #f9f9f9;
                }}
                tr:hover {{
                    background-color: #f5f5f5;
                }}
                .container {{
                    background: white;
                    padding: 20px;
                    border-radius: 8px;
                    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                    margin-bottom: 20px;
                }}
                .footer {{
                    text-align: center;
                    margin-top: 30px;
                    padding: 20px;
                    color: #666;
                    font-size: 0.9em;
                }}
            </style>
        </head>
        <body>
            <div class="container">
                <h1>eccDNA Analysis Report of {sample_name} Sample</h1>
                
                <h2>1. Read Statistics</h2>
                {read_stats}
                
                <h2>2. CtcR Classification</h2>
                {ctcr_stats}
                
                <h2>3. eccDNA Statistics</h2>
                {eccdna_stats}
                
                <h2>4. Overall Statistics</h2>
                {overall_stats}

                <h2>5. Overlap Information</h2>
                <p>{overlap_info}</p>

            </div>
            <div class="footer">
                Generated by eccDNA Analysis Tool CircleSeeker
            </div>
        </body>
        </html>"""

        eccdna_stats_html = stats['eccdna_stats'].to_html(index=False, classes='dataframe') if not stats['eccdna_stats'].empty else "<p>No eccDNA Stats Available</p>"

        html_content = html_template.format(
            sample_name=self.sample_name,
            read_stats=stats['read_stats'].to_html(index=False, classes='dataframe'),
            ctcr_stats=stats['ctcr_stats'].to_html(index=False, classes='dataframe'),
            eccdna_stats=eccdna_stats_html,
            overall_stats=stats['overall_stats'].to_html(index=False, classes='dataframe'),
            overlap_info=stats['overlap_info']
        )
        
        with open(self.html_output, 'w') as f:
            f.write(html_content)
            
        self.logger.info(f"HTML report has been generated and saved to {self.html_output}")

    def run(self):
        try:
            self.logger.info("Starting report generation...")
            stats = self.generate_tabular_report()
            self.logger.info("Generating text report...")
            self.save_text_report(stats)
            self.logger.info("Generating HTML report...")
            self.save_html_report(stats)
            self.logger.info("Generating summary CSV...")
            self.create_summary_csv()
            self.logger.info("All reports generated successfully!")
        except Exception as e:
            self.logger.error(f"Error during report generation: {str(e)}")
            raise

def main():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.StreamHandler()]
    )
    logger = logging.getLogger(__name__)

    parser = argparse.ArgumentParser(description="Generate eccDNA analysis report")
    # Required parameters
    parser.add_argument("--fai", required=True, help="Path to All_HIFI_reads.fasta.fai")
    parser.add_argument("--ctcr1", required=True, help="Path to CtcR.List_patr1.csv")
    parser.add_argument("--ctcr2", required=True, help="Path to CtcR.List_patr2.csv")
    parser.add_argument("--linr", required=True, help="Path to Num_LinR.csv")
    parser.add_argument("--uecc", required=True, help="Path to Merge.Uecc.fa.fai")
    parser.add_argument("--mecc", required=True, help="Path to Merge.Mecc.fa.fai")
    parser.add_argument("--sample_name", required=True, help="Sample name for the report title")

    # Optional parameters
    parser.add_argument("--cecc", help="Path to Cecc.Confirmed.fa.fai")
    parser.add_argument("--mcecc", help="Path to Mecc.Confirmed.fa.fai")
    parser.add_argument("--uecc_inferred_fai", help="Path to Uecc.Inferred.fa.fai (Inferred UeccDNA)")
    parser.add_argument("--shared_count_csv", help="Path to demo.Final.SharedCount.csv")
    parser.add_argument("--output", default="report.txt", help="Path to output report file")
    parser.add_argument("--summary_output", default="summary.csv", help="Path to output summary CSV file")
    parser.add_argument("--html_output", default="report.html", help="Path to output HTML report file")
    parser.add_argument("--keep-tmp", action="store_true", help="Keep temporary files and show debug logs")

    try:
        args = parser.parse_args()

        report_generator = ReportGenerator(
            fai=args.fai,
            ctcr1=args.ctcr1,
            ctcr2=args.ctcr2,
            linr=args.linr,
            uecc_fai=args.uecc,
            mecc_fai=args.mecc,
            cecc_fai=args.cecc,
            mcecc_fai=args.mcecc,
            uecc_inferred_fai=args.uecc_inferred_fai,
            shared_count_csv=args.shared_count_csv,
            output=args.output,
            summary_output=args.summary_output,
            html_output=args.html_output,
            sample_name=args.sample_name,
            keep_tmp=args.keep_tmp
        )

        report_generator.run()

    except Exception as e:
        logger.error(f"Error: {str(e)}")
        raise

if __name__ == "__main__":
    main()