import pandas as pd
import numpy as np
from collections import Counter
import argparse
import os

class ReportGenerator:
    def __init__(self, fai, ctcr1, ctcr2, linr, uecc, mecc, xecc, cecc, uecc_inferred, output):
        self.fai = fai
        self.ctcr1 = ctcr1
        self.ctcr2 = ctcr2
        self.linr = linr
        self.uecc = uecc
        self.mecc = mecc
        self.xecc = xecc
        self.cecc = cecc
        self.uecc_inferred = uecc_inferred
        self.output = output

    def count_reads(self):
        with open(self.fai, 'r') as f:
            return sum(1 for line in f)

    def count_ctcr(self):
        df1 = pd.read_csv(self.ctcr1)
        df2 = pd.read_csv(self.ctcr2)
        combined = pd.concat([df1, df2])
        return len(combined['readName'].unique())

    def classify_ctcr(self):
        df1 = pd.read_csv(self.ctcr1)
        df2 = pd.read_csv(self.ctcr2)
        combined = pd.concat([df1, df2])
        return Counter(combined['readClass'])

    def count_linr(self):
        df = pd.read_csv(self.linr)
        return df['Num_LinR'].iloc[0]

    def analyze_eccdna(self, file):
        df = pd.read_csv(file)
        count = len(df['eName'].unique())
        lengths = df['eLength']
        return {
            'count': count,
            'min_length': lengths.min(),
            'max_length': lengths.max(),
            'mean_length': lengths.mean(),
            'median_length': lengths.median(),
            'mode_length': lengths.mode().iloc[0]
        }

    def generate_report(self):
        total_reads = self.count_reads()
        linr_count = self.count_linr()
        ctcr_count = self.count_ctcr()
        ctcr_classes = self.classify_ctcr()

        uecc = self.analyze_eccdna(self.uecc)
        mecc = self.analyze_eccdna(self.mecc)
        xecc = self.analyze_eccdna(self.xecc)
        cecc = self.analyze_eccdna(self.cecc)
        uecc_inferred = self.analyze_eccdna(self.uecc_inferred)

        confirmed_eccdna = uecc['count'] + mecc['count'] + xecc['count'] + cecc['count']
        total_eccdna = confirmed_eccdna + uecc_inferred['count']

        all_confirmed_lengths = pd.concat([
            pd.read_csv(self.uecc)['eLength'],
            pd.read_csv(self.mecc)['eLength'],
            pd.read_csv(self.xecc)['eLength'],
            pd.read_csv(self.cecc)['eLength']
        ])

        all_eccdna_lengths = pd.concat([all_confirmed_lengths, pd.read_csv(self.uecc_inferred)['eLength']])

        report = f"""
eccDNA Analysis Report
======================

1. Read Statistics:
   Total Reads: {total_reads}
   Linear Reads (LinR): {linr_count} ({linr_count/total_reads:.2%})
   Circular Reads (CtcR): {ctcr_count} ({ctcr_count/total_reads:.2%})

2. CtcR Classification:
   CtcR-multiple: {ctcr_classes['CtcR-multiple']} ({ctcr_classes['CtcR-multiple']/ctcr_count:.2%})
   CtcR-inversion: {ctcr_classes['CtcR-inversion']} ({ctcr_classes['CtcR-inversion']/ctcr_count:.2%})
   CtcR-perfect: {ctcr_classes['CtcR-perfect']} ({ctcr_classes['CtcR-perfect']/ctcr_count:.2%})
   CtcR-hybrid: {ctcr_classes.get('CtcR-hybrid', 0)} ({ctcr_classes.get('CtcR-hybrid', 0)/ctcr_count:.2%})

3. Confirmed eccDNA Analysis:
   Total Confirmed eccDNA: {confirmed_eccdna}

   UeccDNA: {uecc['count']}
     Shortest: {uecc['min_length']} bp
     Longest: {uecc['max_length']} bp
     Average: {uecc['mean_length']:.2f} bp
     Median: {uecc['median_length']:.2f} bp
     Mode: {uecc['mode_length']} bp

   MeccDNA: {mecc['count']}
     Shortest: {mecc['min_length']} bp
     Longest: {mecc['max_length']} bp
     Average: {mecc['mean_length']:.2f} bp
     Median: {mecc['median_length']:.2f} bp
     Mode: {mecc['mode_length']} bp

   XeccDNA: {xecc['count']}
     Shortest: {xecc['min_length']} bp
     Longest: {xecc['max_length']} bp
     Average: {xecc['mean_length']:.2f} bp
     Median: {xecc['median_length']:.2f} bp
     Mode: {xecc['mode_length']} bp

   CeccDNA: {cecc['count']}
     Shortest: {cecc['min_length']} bp
     Longest: {cecc['max_length']} bp
     Average: {cecc['mean_length']:.2f} bp
     Median: {cecc['median_length']:.2f} bp
     Mode: {cecc['mode_length']} bp

4. Inferred eccDNA Analysis:
   Total Inferred UeccDNA: {uecc_inferred['count']}
     Shortest: {uecc_inferred['min_length']} bp
     Longest: {uecc_inferred['max_length']} bp
     Average: {uecc_inferred['mean_length']:.2f} bp
     Median: {uecc_inferred['median_length']:.2f} bp
     Mode: {uecc_inferred['mode_length']} bp

5. Overall eccDNA Statistics:
   Total eccDNA (Confirmed + Inferred): {total_eccdna}

   All Confirmed eccDNA:
     Shortest: {all_confirmed_lengths.min()} bp
     Longest: {all_confirmed_lengths.max()} bp
     Average: {all_confirmed_lengths.mean():.2f} bp
     Median: {all_confirmed_lengths.median():.2f} bp
     Mode: {all_confirmed_lengths.mode().iloc[0]} bp

   All eccDNA (Confirmed + Inferred):
     Shortest: {all_eccdna_lengths.min()} bp
     Longest: {all_eccdna_lengths.max()} bp
     Average: {all_eccdna_lengths.mean():.2f} bp
     Median: {all_eccdna_lengths.median():.2f} bp
     Mode: {all_eccdna_lengths.mode().iloc[0]} bp
"""

        return report

    def save_report(self, report):
        with open(self.output, "w") as f:
            f.write(report)
        print(f"Report has been generated and saved to {self.output}")

    def run(self):
        report = self.generate_report()
        self.save_report(report)

def main():
    parser = argparse.ArgumentParser(description="Generate eccDNA analysis report")
    parser.add_argument("--fai", required=True, help="Path to All_HIFI_reads.fasta.fai")
    parser.add_argument("--ctcr1", required=True, help="Path to CtcR.List_patr1.csv")
    parser.add_argument("--ctcr2", required=True, help="Path to CtcR.List_patr2.csv")
    parser.add_argument("--linr", required=True, help="Path to Num_LinR.csv")
    parser.add_argument("--uecc", required=True, help="Path to Merge.Uecc.csv")
    parser.add_argument("--mecc", required=True, help="Path to Merge.Mecc.csv")
    parser.add_argument("--xecc", required=True, help="Path to Xecc.Confirmed.csv")
    parser.add_argument("--cecc", required=True, help="Path to Cecc.Confirmed.csv")
    parser.add_argument("--uecc_inferred", required=True, help="Path to Uecc.Inferred.csv")
    parser.add_argument("--output", required=True, help="Path to output report file")

    args = parser.parse_args()

    report_generator = ReportGenerator(
        args.fai, args.ctcr1, args.ctcr2, args.linr, args.uecc, args.mecc,
        args.xecc, args.cecc, args.uecc_inferred, args.output
    )
    report_generator.run()

if __name__ == "__main__":
    main()