#!/usr/bin/env python3
# coding: utf-8

"""
EccDNA Analysis Report Generator
Analyzes eccDNA FASTA files and generates comprehensive HTML reports
"""

import pandas as pd
import numpy as np
from Bio import SeqIO
import argparse
import logging
import os
from collections import Counter
from typing import Dict, List, Set, Tuple
import sys

class EccDNAAnalyzer:
    def __init__(self, 
                 uecc_fasta: str,
                 mecc_fasta: str, 
                 cecc_fasta: str,
                 inferred_fasta: str,
                 reads_fasta: str,
                 processed_csv: str,
                 output_html: str = "eccdna_report.html",
                 sample_name: str = "Unknown"):
        
        self.uecc_fasta = uecc_fasta
        self.mecc_fasta = mecc_fasta
        self.cecc_fasta = cecc_fasta
        self.inferred_fasta = inferred_fasta
        self.reads_fasta = reads_fasta
        self.processed_csv = processed_csv
        self.output_html = output_html
        self.sample_name = sample_name
        
        # Setup logging
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        self.logger = logging.getLogger(__name__)
        
    def read_fasta_sequences(self, fasta_file: str) -> Dict[str, int]:
        """Read FASTA file and return dictionary of sequence IDs and lengths"""
        sequences = {}
        try:
            for record in SeqIO.parse(fasta_file, "fasta"):
                sequences[str(record.id)] = len(record.seq)
            self.logger.info(f"Read {len(sequences)} sequences from {fasta_file}")
        except Exception as e:
            self.logger.error(f"Error reading {fasta_file}: {str(e)}")
        return sequences
    
    def calculate_statistics(self, lengths: List[int]) -> Dict:
        """Calculate statistical measures for sequence lengths"""
        if not lengths:
            return {
                'count': 0,
                'min_length': 0,
                'max_length': 0,
                'mean_length': 0,
                'median_length': 0,
                'mode_length': 0,
                'std_dev': 0
            }
        
        lengths_array = np.array(lengths)
        mode_result = Counter(lengths).most_common(1)
        
        return {
            'count': len(lengths),
            'min_length': int(np.min(lengths_array)),
            'max_length': int(np.max(lengths_array)),
            'mean_length': float(np.mean(lengths_array)),
            'median_length': float(np.median(lengths_array)),
            'mode_length': mode_result[0][0] if mode_result else 0,
            'std_dev': float(np.std(lengths_array))
        }
    
    def analyze_fasta_file(self, fasta_file: str, label: str) -> Tuple[Dict, Dict]:
        """Analyze a FASTA file and return statistics and sequences"""
        sequences = self.read_fasta_sequences(fasta_file)
        lengths = list(sequences.values())
        stats = self.calculate_statistics(lengths)
        stats['label'] = label
        
        self.logger.info(f"{label}: {stats['count']} sequences, "
                        f"length range: {stats['min_length']}-{stats['max_length']} bp")
        
        return stats, sequences
    
    def find_overlaps(self, inferred_fasta: str, confirmed_fasta: str, 
                     label: str) -> Set[str]:
        """
        Find overlapping sequences based on sequence content
        Note: This compares actual sequence content, not just IDs
        """
        # Read inferred sequences
        inferred_seq_dict = {}
        for record in SeqIO.parse(inferred_fasta, "fasta"):
            seq_str = str(record.seq).upper()  # Normalize to uppercase
            if seq_str not in inferred_seq_dict:
                inferred_seq_dict[seq_str] = []
            inferred_seq_dict[seq_str].append(record.id)
        
        # Find overlaps with confirmed sequences
        overlapping_ids = set()
        for record in SeqIO.parse(confirmed_fasta, "fasta"):
            seq_str = str(record.seq).upper()
            if seq_str in inferred_seq_dict:
                # Add all inferred IDs that match this sequence
                overlapping_ids.update(inferred_seq_dict[seq_str])
        
        self.logger.info(f"Found {len(overlapping_ids)} sequence overlaps between inferred and {label}")
        return overlapping_ids
    
    def analyze_reads_classification(self) -> Dict:
        """Analyze read classification from CSV file"""
        try:
            # Read the CSV file
            df = pd.read_csv(self.processed_csv)
            
            # Get all reads from FASTA
            all_reads = self.read_fasta_sequences(self.reads_fasta)
            total_reads = len(all_reads)
            
            # Count CtcR reads by class
            ctcr_counts = df['readClass'].value_counts().to_dict()
            total_ctcr = len(df)
            
            # Calculate other reads
            ctcr_read_names = set(df['readName'].values)
            all_read_names = set(all_reads.keys())
            other_reads = all_read_names - ctcr_read_names
            other_count = len(other_reads)
            
            # Ensure all CtcR classes are present
            for ctcr_class in ['CtcR-perfect', 'CtcR-hybrid', 'CtcR-inversion']:
                if ctcr_class not in ctcr_counts:
                    ctcr_counts[ctcr_class] = 0
            
            return {
                'total_reads': total_reads,
                'total_ctcr': total_ctcr,
                'other_reads': other_count,
                'ctcr_perfect': ctcr_counts.get('CtcR-perfect', 0),
                'ctcr_hybrid': ctcr_counts.get('CtcR-hybrid', 0),
                'ctcr_inversion': ctcr_counts.get('CtcR-inversion', 0)
            }
            
        except Exception as e:
            self.logger.error(f"Error analyzing reads classification: {str(e)}")
            return {
                'total_reads': 0,
                'total_ctcr': 0,
                'other_reads': 0,
                'ctcr_perfect': 0,
                'ctcr_hybrid': 0,
                'ctcr_inversion': 0
            }
    
    def generate_html_report(self, results: Dict):
        """Generate HTML report with all analysis results"""
        
        html_template = """<!DOCTYPE html>
<html>
<head>
    <title>eccDNA Analysis Report - {sample_name}</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            margin: 0;
            padding: 0;
            background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
            min-height: 100vh;
        }}
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            padding: 20px;
        }}
        h1 {{
            color: #2c3e50;
            text-align: center;
            padding: 30px;
            background: white;
            border-radius: 10px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
            margin-bottom: 30px;
        }}
        h2 {{
            color: #34495e;
            margin-top: 40px;
            padding: 15px;
            background: linear-gradient(90deg, #667eea 0%, #764ba2 100%);
            color: white;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }}
        .stat-card {{
            background: white;
            padding: 20px;
            border-radius: 10px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
            transition: transform 0.3s ease;
        }}
        .stat-card:hover {{
            transform: translateY(-5px);
            box-shadow: 0 4px 12px rgba(0,0,0,0.15);
        }}
        .stat-card h3 {{
            color: #667eea;
            margin-top: 0;
            font-size: 1.2em;
            border-bottom: 2px solid #667eea;
            padding-bottom: 10px;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            background: white;
            border-radius: 10px;
            overflow: hidden;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
            margin: 20px 0;
        }}
        th {{
            background: linear-gradient(90deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 15px;
            text-align: left;
            font-weight: 600;
        }}
        td {{
            padding: 12px 15px;
            border-bottom: 1px solid #ecf0f1;
        }}
        tr:hover {{
            background-color: #f8f9fa;
        }}
        .overlap-section {{
            background: white;
            padding: 25px;
            border-radius: 10px;
            margin: 20px 0;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }}
        .overlap-grid {{
            display: grid;
            grid-template-columns: repeat(3, 1fr);
            gap: 20px;
            margin-top: 20px;
        }}
        .overlap-card {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            text-align: center;
        }}
        .overlap-card .number {{
            font-size: 2.5em;
            font-weight: bold;
            margin: 10px 0;
        }}
        .overlap-card .label {{
            font-size: 1.1em;
            opacity: 0.95;
        }}
        .percentage {{
            color: #27ae60;
            font-weight: bold;
        }}
        .footer {{
            text-align: center;
            margin-top: 50px;
            padding: 20px;
            color: #7f8c8d;
            background: white;
            border-radius: 10px;
        }}
        .summary-box {{
            background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%);
            color: white;
            padding: 20px;
            border-radius: 10px;
            margin: 20px 0;
            font-size: 1.1em;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 eccDNA Analysis Report - {sample_name}</h1>
        
        <h2>📊 1. Read Classification Statistics</h2>
        <div class="stats-grid">
            <div class="stat-card">
                <h3>Total Reads</h3>
                <p style="font-size: 2em; color: #3498db; margin: 10px 0;">{total_reads:,}</p>
            </div>
            <div class="stat-card">
                <h3>CtcR Reads</h3>
                <p style="font-size: 2em; color: #e74c3c; margin: 10px 0;">{total_ctcr:,}</p>
                <p style="color: #7f8c8d;">({ctcr_percentage:.2f}% of total)</p>
            </div>
            <div class="stat-card">
                <h3>Other Reads</h3>
                <p style="font-size: 2em; color: #95a5a6; margin: 10px 0;">{other_reads:,}</p>
                <p style="color: #7f8c8d;">({other_percentage:.2f}% of total)</p>
            </div>
        </div>
        
        <h2>🔍 2. CtcR Subtype Distribution</h2>
        <table>
            <thead>
                <tr>
                    <th>CtcR Subtype</th>
                    <th>Count</th>
                    <th>Percentage of CtcR</th>
                    <th>Percentage of Total Reads</th>
                </tr>
            </thead>
            <tbody>
                <tr>
                    <td><strong>CtcR-perfect</strong></td>
                    <td>{ctcr_perfect:,}</td>
                    <td>{ctcr_perfect_pct:.2f}%</td>
                    <td>{ctcr_perfect_total_pct:.2f}%</td>
                </tr>
                <tr>
                    <td><strong>CtcR-hybrid</strong></td>
                    <td>{ctcr_hybrid:,}</td>
                    <td>{ctcr_hybrid_pct:.2f}%</td>
                    <td>{ctcr_hybrid_total_pct:.2f}%</td>
                </tr>
                <tr>
                    <td><strong>CtcR-inversion</strong></td>
                    <td>{ctcr_inversion:,}</td>
                    <td>{ctcr_inversion_pct:.2f}%</td>
                    <td>{ctcr_inversion_total_pct:.2f}%</td>
                </tr>
            </tbody>
        </table>
        
        <h2>📈 3. eccDNA Sequence Statistics</h2>
        <table>
            <thead>
                <tr>
                    <th>eccDNA Type</th>
                    <th>Count</th>
                    <th>Min Length (bp)</th>
                    <th>Max Length (bp)</th>
                    <th>Mean Length (bp)</th>
                    <th>Median Length (bp)</th>
                    <th>Mode Length (bp)</th>
                    <th>Std Dev</th>
                </tr>
            </thead>
            <tbody>
                {eccdna_stats_rows}
            </tbody>
        </table>
        
        <div class="summary-box" style="margin-top: 20px;">
            <strong>Note on All eccDNA Statistics:</strong> The "All eccDNA (Combined)" row represents the aggregate statistics 
            of all eccDNA sequences from UeccDNA, MeccDNA, CeccDNA, and Inferred eccDNA combined. 
            This provides an overall view of the eccDNA size distribution in the sample.
        </div>
        
        <h2>🔄 4. Overlap Analysis Between Inferred and Confirmed eccDNA</h2>
        <div class="overlap-section">
            <div class="summary-box">
                <strong>Total Inferred eccDNA sequences:</strong> {total_inferred:,}
            </div>
            <div class="overlap-grid">
                <div class="overlap-card">
                    <div class="label">Overlaps with UeccDNA</div>
                    <div class="number">{uecc_overlaps:,}</div>
                    <div class="label">{uecc_overlap_pct:.2f}% of inferred</div>
                </div>
                <div class="overlap-card">
                    <div class="label">Overlaps with MeccDNA</div>
                    <div class="number">{mecc_overlaps:,}</div>
                    <div class="label">{mecc_overlap_pct:.2f}% of inferred</div>
                </div>
                <div class="overlap-card">
                    <div class="label">Overlaps with CeccDNA</div>
                    <div class="number">{cecc_overlaps:,}</div>
                    <div class="label">{cecc_overlap_pct:.2f}% of inferred</div>
                </div>
            </div>
            <div class="summary-box" style="margin-top: 20px; background: linear-gradient(135deg, #4facfe 0%, #00f2fe 100%);">
                <strong>Summary:</strong> Among {total_inferred:,} inferred eccDNA sequences, 
                {uecc_overlaps:,} share identical sequences with UeccDNA ({uecc_overlap_pct:.2f}%), 
                {mecc_overlaps:,} share identical sequences with MeccDNA ({mecc_overlap_pct:.2f}%), 
                and {cecc_overlaps:,} share identical sequences with CeccDNA ({cecc_overlap_pct:.2f}%).
                <br><br>
                <em>Note: Overlaps are determined by comparing actual sequence content, not sequence IDs. 
                Sequences are considered overlapping if they have identical nucleotide content.</em>
            </div>
        </div>
        
        <div class="footer">
            <p>Generated by eccDNA Analysis Tool <strong>CircleSeeker2</strong></p>
            <p>Report generated for sample: <strong>{sample_name}</strong></p>
        </div>
    </div>
</body>
</html>"""
        
        # Generate eccDNA stats rows
        eccdna_stats_rows = ""
        for stats in results['eccdna_stats']:
            # Add special styling for the combined stats row
            if stats['label'] == "All eccDNA (Combined)":
                eccdna_stats_rows += f"""
                <tr style="border-top: 3px solid #667eea; font-weight: bold; background-color: #f0f4ff;">
                    <td><strong>{stats['label']}</strong></td>
                    <td>{stats['count']:,}</td>
                    <td>{stats['min_length']:,}</td>
                    <td>{stats['max_length']:,}</td>
                    <td>{stats['mean_length']:.2f}</td>
                    <td>{stats['median_length']:.2f}</td>
                    <td>{stats['mode_length']:,}</td>
                    <td>{stats['std_dev']:.2f}</td>
                </tr>"""
            else:
                eccdna_stats_rows += f"""
                <tr>
                    <td><strong>{stats['label']}</strong></td>
                    <td>{stats['count']:,}</td>
                    <td>{stats['min_length']:,}</td>
                    <td>{stats['max_length']:,}</td>
                    <td>{stats['mean_length']:.2f}</td>
                    <td>{stats['median_length']:.2f}</td>
                    <td>{stats['mode_length']:,}</td>
                    <td>{stats['std_dev']:.2f}</td>
                </tr>"""
        
        # Calculate percentages for reads
        reads = results['reads_classification']
        total = reads['total_reads'] if reads['total_reads'] > 0 else 1
        ctcr_total = reads['total_ctcr'] if reads['total_ctcr'] > 0 else 1
        
        # Calculate percentages for overlaps
        inferred_total = results['inferred_stats']['count'] if results['inferred_stats']['count'] > 0 else 1
        
        # Fill in the template
        html_content = html_template.format(
            sample_name=self.sample_name,
            total_reads=reads['total_reads'],
            total_ctcr=reads['total_ctcr'],
            other_reads=reads['other_reads'],
            ctcr_percentage=(reads['total_ctcr'] / total * 100),
            other_percentage=(reads['other_reads'] / total * 100),
            ctcr_perfect=reads['ctcr_perfect'],
            ctcr_hybrid=reads['ctcr_hybrid'],
            ctcr_inversion=reads['ctcr_inversion'],
            ctcr_perfect_pct=(reads['ctcr_perfect'] / ctcr_total * 100),
            ctcr_hybrid_pct=(reads['ctcr_hybrid'] / ctcr_total * 100),
            ctcr_inversion_pct=(reads['ctcr_inversion'] / ctcr_total * 100),
            ctcr_perfect_total_pct=(reads['ctcr_perfect'] / total * 100),
            ctcr_hybrid_total_pct=(reads['ctcr_hybrid'] / total * 100),
            ctcr_inversion_total_pct=(reads['ctcr_inversion'] / total * 100),
            eccdna_stats_rows=eccdna_stats_rows,
            total_inferred=results['inferred_stats']['count'],
            uecc_overlaps=len(results['overlaps']['uecc']),
            mecc_overlaps=len(results['overlaps']['mecc']),
            cecc_overlaps=len(results['overlaps']['cecc']),
            uecc_overlap_pct=(len(results['overlaps']['uecc']) / inferred_total * 100),
            mecc_overlap_pct=(len(results['overlaps']['mecc']) / inferred_total * 100),
            cecc_overlap_pct=(len(results['overlaps']['cecc']) / inferred_total * 100)
        )
        
        # Write HTML file
        with open(self.output_html, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        self.logger.info(f"HTML report generated: {self.output_html}")
    
    def run(self):
        """Main analysis workflow"""
        self.logger.info("Starting eccDNA analysis...")
        
        # Analyze confirmed eccDNA files
        uecc_stats, uecc_seqs = self.analyze_fasta_file(self.uecc_fasta, "UeccDNA")
        mecc_stats, mecc_seqs = self.analyze_fasta_file(self.mecc_fasta, "MeccDNA")
        cecc_stats, cecc_seqs = self.analyze_fasta_file(self.cecc_fasta, "CeccDNA")
        
        # Analyze inferred eccDNA
        inferred_stats, inferred_seqs = self.analyze_fasta_file(self.inferred_fasta, "Inferred eccDNA")
        
        # Find overlaps (based on sequence content)
        uecc_overlaps = self.find_overlaps(self.inferred_fasta, self.uecc_fasta, "UeccDNA")
        mecc_overlaps = self.find_overlaps(self.inferred_fasta, self.mecc_fasta, "MeccDNA")
        cecc_overlaps = self.find_overlaps(self.inferred_fasta, self.cecc_fasta, "CeccDNA")
        
        # Calculate combined statistics for all eccDNA
        all_eccdna_lengths = []
        all_eccdna_lengths.extend(list(uecc_seqs.values()))
        all_eccdna_lengths.extend(list(mecc_seqs.values()))
        all_eccdna_lengths.extend(list(cecc_seqs.values()))
        all_eccdna_lengths.extend(list(inferred_seqs.values()))
        
        all_stats = self.calculate_statistics(all_eccdna_lengths)
        all_stats['label'] = "All eccDNA (Combined)"
        
        # Analyze reads classification
        reads_classification = self.analyze_reads_classification()
        
        # Compile results
        results = {
            'eccdna_stats': [uecc_stats, mecc_stats, cecc_stats, inferred_stats, all_stats],
            'inferred_stats': inferred_stats,
            'overlaps': {
                'uecc': uecc_overlaps,
                'mecc': mecc_overlaps,
                'cecc': cecc_overlaps
            },
            'reads_classification': reads_classification
        }
        
        # Generate HTML report
        self.generate_html_report(results)
        
        # Print summary to console
        self.print_summary(results)
        
        self.logger.info("Analysis complete!")
        
    def print_summary(self, results):
        """Print analysis summary to console"""
        print("\n" + "="*60)
        print(f"eccDNA Analysis Summary - {self.sample_name}")
        print("="*60)
        
        print("\n📊 Read Classification:")
        reads = results['reads_classification']
        print(f"  Total Reads: {reads['total_reads']:,}")
        print(f"  CtcR Reads: {reads['total_ctcr']:,} ({reads['total_ctcr']/reads['total_reads']*100:.2f}%)")
        print(f"    - CtcR-perfect: {reads['ctcr_perfect']:,}")
        print(f"    - CtcR-hybrid: {reads['ctcr_hybrid']:,}")
        print(f"    - CtcR-inversion: {reads['ctcr_inversion']:,}")
        print(f"  Other Reads: {reads['other_reads']:,} ({reads['other_reads']/reads['total_reads']*100:.2f}%)")
        
        print("\n🧬 eccDNA Statistics:")
        for stats in results['eccdna_stats']:
            print(f"\n  {stats['label']}:")
            print(f"    Count: {stats['count']:,}")
            print(f"    Length range: {stats['min_length']:,} - {stats['max_length']:,} bp")
            print(f"    Mean length: {stats['mean_length']:.2f} bp")
            
        print("\n🔄 Overlap Analysis:")
        inferred_total = results['inferred_stats']['count']
        print(f"  Total Inferred eccDNA: {inferred_total:,}")
        print(f"  Overlaps with UeccDNA: {len(results['overlaps']['uecc']):,} ({len(results['overlaps']['uecc'])/inferred_total*100:.2f}%)")
        print(f"  Overlaps with MeccDNA: {len(results['overlaps']['mecc']):,} ({len(results['overlaps']['mecc'])/inferred_total*100:.2f}%)")
        print(f"  Overlaps with CeccDNA: {len(results['overlaps']['cecc']):,} ({len(results['overlaps']['cecc'])/inferred_total*100:.2f}%)")
        
        print("\n" + "="*60)

def main():
    parser = argparse.ArgumentParser(
        description="Analyze eccDNA FASTA files and generate comprehensive reports",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python eccdna_analyzer.py \\
    --uecc step6_UeccDNA_pre.fasta \\
    --mecc step6_MeccDNA_pre.fasta \\
    --cecc step6_CeccDNA_pre.fasta \\
    --inferred step13.renamed.fasta \\
    --reads step1_Test.HeLa_10k.fasta \\
    --csv step2_processed.csv \\
    --sample "HeLa_10k" \\
    --output report.html
        """
    )
    
    # Required arguments
    parser.add_argument("--uecc", required=True, 
                       help="Path to UeccDNA FASTA file (e.g., step6_UeccDNA_pre.fasta)")
    parser.add_argument("--mecc", required=True,
                       help="Path to MeccDNA FASTA file (e.g., step6_MeccDNA_pre.fasta)")
    parser.add_argument("--cecc", required=True,
                       help="Path to CeccDNA FASTA file (e.g., step6_CeccDNA_pre.fasta)")
    parser.add_argument("--inferred", required=True,
                       help="Path to inferred eccDNA FASTA file (e.g., step13.renamed.fasta)")
    parser.add_argument("--reads", required=True,
                       help="Path to reads FASTA file (e.g., step1_Test.HeLa_10k.fasta)")
    parser.add_argument("--csv", required=True,
                       help="Path to processed CSV file with read classifications (e.g., step2_processed.csv)")
    
    # Optional arguments
    parser.add_argument("--sample", default="Unknown",
                       help="Sample name for the report (default: Unknown)")
    parser.add_argument("--output", default="eccdna_report.html",
                       help="Output HTML report filename (default: eccdna_report.html)")
    
    args = parser.parse_args()
    
    # Check if all input files exist
    input_files = [args.uecc, args.mecc, args.cecc, args.inferred, args.reads, args.csv]
    for file_path in input_files:
        if not os.path.exists(file_path):
            print(f"Error: Input file '{file_path}' not found!")
            sys.exit(1)
    
    # Create analyzer instance and run
    analyzer = EccDNAAnalyzer(
        uecc_fasta=args.uecc,
        mecc_fasta=args.mecc,
        cecc_fasta=args.cecc,
        inferred_fasta=args.inferred,
        reads_fasta=args.reads,
        processed_csv=args.csv,
        output_html=args.output,
        sample_name=args.sample
    )
    
    try:
        analyzer.run()
        print(f"\n✅ Success! HTML report saved to: {args.output}")
    except Exception as e:
        logging.error(f"Analysis failed: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
