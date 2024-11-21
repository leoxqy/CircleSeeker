import pandas as pd
import numpy as np
from collections import Counter
import argparse
import os
import logging

class ReportGenerator:
    def __init__(self, fai, ctcr1, ctcr2, linr, uecc_fai, mecc_fai, cecc_fai=None, mcecc_fai=None, 
                 uecc_inferred=None, output="report.txt", summary_output="summary.csv", html_output="report.html"):
        self.fai = fai
        self.ctcr1 = ctcr1
        self.ctcr2 = ctcr2
        self.linr = linr
        self.uecc_fai = uecc_fai
        self.mecc_fai = mecc_fai
        self.cecc_fai = cecc_fai
        self.mcecc_fai = mcecc_fai
        self.uecc_inferred = uecc_inferred
        self.output = output
        self.summary_output = summary_output
        self.html_output = html_output
        
        # 简化logger初始化，只获取logger实例
        self.logger = logging.getLogger(__name__)

    def count_reads(self):
        """Count total number of reads from FAI file"""
        try:
            with open(self.fai, 'r') as f:
                return sum(1 for line in f)
        except Exception as e:
            self.logger.error(f"Error reading FAI file: {str(e)}")
            return 0

    def count_ctcr(self):
        """Count unique circular reads from two CSV files"""
        try:
            # Read both CSV files
            df1 = pd.read_csv(self.ctcr1)
            df2 = pd.read_csv(self.ctcr2)
            
            # For first file, ensure 'readName' column exists
            if 'readName' not in df1.columns and 'qName' in df1.columns:
                df1 = df1.rename(columns={'qName': 'readName'})
                
            # For second file, ensure 'readName' column exists
            if 'readName' not in df2.columns and 'qName' in df2.columns:
                df2 = df2.rename(columns={'qName': 'readName'})
            
            # Combine dataframes and count unique reads
            combined = pd.concat([df1[['readName']], df2[['readName']]])
            return len(combined['readName'].unique())
        except Exception as e:
            self.logger.error(f"Error counting CtcR: {str(e)}")
            return 0

    def classify_ctcr(self):
        """Classify circular reads from two CSV files"""
        try:
            # Read both CSV files
            df1 = pd.read_csv(self.ctcr1)
            df2 = pd.read_csv(self.ctcr2)
            
            # For first file, ensure 'readClass' column exists
            if 'readClass' not in df1.columns and 'class' in df1.columns:
                df1 = df1.rename(columns={'class': 'readClass'})
                
            # For second file, ensure 'readClass' column exists
            if 'readClass' not in df2.columns and 'class' in df2.columns:
                df2 = df2.rename(columns={'class': 'readClass'})
            
            # Combine dataframes and count classes
            combined = pd.concat([df1[['readClass']], df2[['readClass']]])
            class_counts = Counter(combined['readClass'])
            
            # Ensure all expected classes are represented
            expected_classes = ['CtcR-multiple', 'CtcR-inversion', 'CtcR-perfect', 'CtcR-hybrid']
            for class_name in expected_classes:
                if class_name not in class_counts:
                    class_counts[class_name] = 0
                    
            return class_counts
        except Exception as e:
            self.logger.error(f"Error classifying CtcR: {str(e)}")
            return Counter()

    def count_linr(self):
        """Count linear reads from CSV file"""
        try:
            df = pd.read_csv(self.linr)
            # Check if column exists with either name
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
        """Read FAI file and return both names and lengths"""
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
        """Read FAI file and return sequence lengths as a pandas Series"""
        if fai_file is None or not os.path.exists(fai_file):
            return pd.Series()
            
        lengths = []
        with open(fai_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    lengths.append(int(parts[1]))
        return pd.Series(lengths)

    def read_csv_if_exists(self, file_path):
        """Read CSV file if it exists, otherwise return empty DataFrame"""
        if file_path and os.path.exists(file_path):
            return pd.read_csv(file_path)
        return pd.DataFrame()

    def create_summary_csv(self):
        """Create a summary CSV with all eccDNA information"""
        dataframes = []

        # Read confirmed eccDNA from FAI files
        uecc_df = self.read_fai_info(self.uecc_fai)
        if not uecc_df.empty:
            uecc_df['Class'] = 'Uecc'
            uecc_df['Status'] = 'Confirmed'
            dataframes.append(uecc_df)

        mecc_df = self.read_fai_info(self.mecc_fai)
        if not mecc_df.empty:
            mecc_df['Class'] = 'Mecc'
            mecc_df['Status'] = 'Confirmed'
            dataframes.append(mecc_df)

        if self.cecc_fai:
            cecc_df = self.read_fai_info(self.cecc_fai)
            if not cecc_df.empty:
                cecc_df['Class'] = 'Cecc'
                cecc_df['Status'] = 'Confirmed'
                dataframes.append(cecc_df)

        if self.mcecc_fai:
            mcecc_df = self.read_fai_info(self.mcecc_fai)
            if not mcecc_df.empty:
                mcecc_df['Class'] = 'MCecc'
                mcecc_df['Status'] = 'Confirmed'
                dataframes.append(mcecc_df)

        # Read inferred eccDNA from CSV if it exists
        if self.uecc_inferred:
            inferred_df = self.read_csv_if_exists(self.uecc_inferred)
            if not inferred_df.empty:
                inferred_df = pd.DataFrame({
                    'Name': inferred_df['eName'],
                    'Length': inferred_df['eLength'],
                    'Class': 'Uecc',
                    'Status': 'Inferred'
                })
                dataframes.append(inferred_df)

        # Combine all DataFrames if any exist
        if dataframes:
            all_eccdna = pd.concat(dataframes, ignore_index=True)
            # Reorder columns to match requested format
            all_eccdna = all_eccdna[['Class', 'Name', 'Length', 'Status']]
            # Save to CSV
            all_eccdna.to_csv(self.summary_output, index=False)
            self.logger.info(f"Summary CSV has been generated and saved to {self.summary_output}")
            return all_eccdna
        else:
            self.logger.warning("No data available for summary CSV")
            return pd.DataFrame(columns=['Class', 'Name', 'Length', 'Status'])

    def get_empty_stats(self):
        """Return empty statistics dictionary"""
        return {
            'count': 0,
            'min_length': 0,
            'max_length': 0,
            'mean_length': 0,
            'median_length': 0,
            'mode_length': 0
        }

    def analyze_eccdna_fai(self, fai_file):
        """Analyze eccDNA using FAI file"""
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

    def analyze_inferred_eccdna(self, file):
        """Analyze inferred eccDNA from CSV"""
        if file is None:
            return self.get_empty_stats()
            
        df = self.read_csv_if_exists(file)
        if df.empty:
            return self.get_empty_stats()
            
        return {
            'count': len(df['eName'].unique()) if 'eName' in df else 0,
            'min_length': df['eLength'].min() if 'eLength' in df else 0,
            'max_length': df['eLength'].max() if 'eLength' in df else 0,
            'mean_length': df['eLength'].mean() if 'eLength' in df else 0,
            'median_length': df['eLength'].median() if 'eLength' in df else 0,
            'mode_length': df['eLength'].mode().iloc[0] if 'eLength' in df and not df.empty else 0
        }

    def generate_tabular_report(self):
        """Generate report in tabular format with error handling"""
        try:
            total_reads = self.count_reads()
            linr_count = self.count_linr()
            ctcr_count = self.count_ctcr()
            ctcr_classes = self.classify_ctcr()

            if total_reads == 0:
                self.logger.warning("Warning: No reads found in FAI file")

            # Calculate other reads
            other_reads = total_reads - linr_count - ctcr_count

            # Create DataFrames for each section with Other Reads included
            read_stats = pd.DataFrame({
                'Category': ['Total Reads', 'Linear Reads (LinR)', 'Concatemeric tandem copies Reads (CtcR)', 'Other Reads'],
                'Count': [total_reads, linr_count, ctcr_count, other_reads],
                'Percentage': [
                    '100%',
                    f'{linr_count/total_reads:.2%}',
                    f'{ctcr_count/total_reads:.2%}',
                    f'{other_reads/total_reads:.2%}'
                ]
            })

            # Use FAI files for confirmed eccDNA
            uecc = self.analyze_eccdna_fai(self.uecc_fai)
            mecc = self.analyze_eccdna_fai(self.mecc_fai)
            cecc = self.analyze_eccdna_fai(self.cecc_fai)
            mcecc = self.analyze_eccdna_fai(self.mcecc_fai)
            uecc_inferred = self.analyze_inferred_eccdna(self.uecc_inferred)

            confirmed_eccdna = uecc['count'] + mecc['count'] + cecc['count'] + mcecc['count']
            total_eccdna = confirmed_eccdna + (uecc_inferred['count'] if self.uecc_inferred else 0)

            ctcr_stats = pd.DataFrame({
                'Category': ['CtcR-multiple', 'CtcR-inversion', 'CtcR-perfect', 'CtcR-hybrid'],
                'Count': [
                    ctcr_classes['CtcR-multiple'],
                    ctcr_classes['CtcR-inversion'],
                    ctcr_classes['CtcR-perfect'],
                    ctcr_classes.get('CtcR-hybrid', 0)
                ],
                'Percentage': [
                    f'{ctcr_classes["CtcR-multiple"]/ctcr_count:.2%}',
                    f'{ctcr_classes["CtcR-inversion"]/ctcr_count:.2%}',
                    f'{ctcr_classes["CtcR-perfect"]/ctcr_count:.2%}',
                    f'{ctcr_classes.get("CtcR-hybrid", 0)/ctcr_count:.2%}'
                ]
            })

            # Create list of eccDNA types and their stats
            eccdna_types = []
            eccdna_counts = []
            min_lengths = []
            max_lengths = []
            mean_lengths = []
            median_lengths = []
            mode_lengths = []

            # Add confirmed types if they have data
            if uecc['count'] > 0:
                eccdna_types.append('UeccDNA')
                eccdna_counts.append(uecc['count'])
                min_lengths.append(uecc['min_length'])
                max_lengths.append(uecc['max_length'])
                mean_lengths.append(f"{uecc['mean_length']:.2f}")
                median_lengths.append(f"{uecc['median_length']:.2f}")
                mode_lengths.append(uecc['mode_length'])

            if mecc['count'] > 0:
                eccdna_types.append('MeccDNA')
                eccdna_counts.append(mecc['count'])
                min_lengths.append(mecc['min_length'])
                max_lengths.append(mecc['max_length'])
                mean_lengths.append(f"{mecc['mean_length']:.2f}")
                median_lengths.append(f"{mecc['median_length']:.2f}")
                mode_lengths.append(mecc['mode_length'])

            if cecc['count'] > 0:
                eccdna_types.append('CeccDNA')
                eccdna_counts.append(cecc['count'])
                min_lengths.append(cecc['min_length'])
                max_lengths.append(cecc['max_length'])
                mean_lengths.append(f"{cecc['mean_length']:.2f}")
                median_lengths.append(f"{cecc['median_length']:.2f}")
                mode_lengths.append(cecc['mode_length'])

            if mcecc['count'] > 0:
                eccdna_types.append('MCeccDNA')
                eccdna_counts.append(mcecc['count'])
                min_lengths.append(mcecc['min_length'])
                max_lengths.append(mcecc['max_length'])
                mean_lengths.append(f"{mcecc['mean_length']:.2f}")
                median_lengths.append(f"{mcecc['median_length']:.2f}")
                mode_lengths.append(mcecc['mode_length'])

            # Add inferred type if it exists and has data
            if self.uecc_inferred and uecc_inferred['count'] > 0:
                eccdna_types.append('Inferred UeccDNA')
                eccdna_counts.append(uecc_inferred['count'])
                min_lengths.append(uecc_inferred['min_length'])
                max_lengths.append(uecc_inferred['max_length'])
                mean_lengths.append(f"{uecc_inferred['mean_length']:.2f}")
                median_lengths.append(f"{uecc_inferred['median_length']:.2f}")
                mode_lengths.append(uecc_inferred['mode_length'])

            eccdna_stats = pd.DataFrame({
                'Type': eccdna_types,
                'Count': eccdna_counts,
                'Min Length (bp)': min_lengths,
                'Max Length (bp)': max_lengths,
                'Mean Length (bp)': mean_lengths,
                'Median Length (bp)': median_lengths,
                'Mode Length (bp)': mode_lengths
            })

            # Combine lengths from FAI files
            all_confirmed_lengths = pd.concat([
                self.read_fai_as_lengths(self.uecc_fai),
                self.read_fai_as_lengths(self.mecc_fai),
                self.read_fai_as_lengths(self.cecc_fai) if self.cecc_fai else pd.Series(),
                self.read_fai_as_lengths(self.mcecc_fai) if self.mcecc_fai else pd.Series()
            ])

            # Add inferred lengths if they exist
            if self.uecc_inferred:
                inferred_df = self.read_csv_if_exists(self.uecc_inferred)
                if not inferred_df.empty:
                    inferred_lengths = inferred_df['eLength']
                    all_eccdna_lengths = pd.concat([all_confirmed_lengths, inferred_lengths])
                else:
                    all_eccdna_lengths = all_confirmed_lengths
            else:
                all_eccdna_lengths = all_confirmed_lengths

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
                'overall_stats': overall_stats
            }

        except Exception as e:
            self.logger.error(f"Error generating report: {str(e)}")
            raise

    def save_text_report(self, stats):
        """Save report in text format with tables"""
        with open(self.output, 'w') as f:
            f.write("eccDNA Analysis Report\n")
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
            
        self.logger.info(f"Text report has been generated and saved to {self.output}")

    def save_html_report(self, stats):
        """Save report in HTML format"""
        html_template = """<!DOCTYPE html>
        <html>
        <head>
            <title>eccDNA Analysis Report</title>
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
                <h1>eccDNA Analysis Report</h1>
                
                <h2>1. Read Statistics</h2>
                {read_stats}
                
                <h2>2. CtcR Classification</h2>
                {ctcr_stats}
                
                <h2>3. eccDNA Statistics</h2>
                {eccdna_stats}
                
                <h2>4. Overall Statistics</h2>
                {overall_stats}
            </div>
            <div class="footer">
                Generated by eccDNA Analysis Tool CircleSeeker
            </div>
        </body>
        </html>"""
                
        # Format the template with the tables
        html_content = html_template.format(
            read_stats=stats['read_stats'].to_html(index=False, classes='dataframe'),
            ctcr_stats=stats['ctcr_stats'].to_html(index=False, classes='dataframe'),
            eccdna_stats=stats['eccdna_stats'].to_html(index=False, classes='dataframe'),
            overall_stats=stats['overall_stats'].to_html(index=False, classes='dataframe')
        )
        
        # Save the HTML file
        with open(self.html_output, 'w') as f:
            f.write(html_content)
            
        self.logger.info(f"HTML report has been generated and saved to {self.html_output}")


    def run(self):
        """Execute all report generation steps"""
        try:
            self.logger.info("Starting report generation...")
            
            # Generate all statistics
            stats = self.generate_tabular_report()
            
            # Save text report
            self.logger.info("Generating text report...")
            self.save_text_report(stats)
            
            # Save HTML report
            self.logger.info("Generating HTML report...")
            self.save_html_report(stats)
            
            # Create summary CSV
            self.logger.info("Generating summary CSV...")
            self.create_summary_csv()
            
            self.logger.info("All reports generated successfully!")
            
        except Exception as e:
            self.logger.error(f"Error during report generation: {str(e)}")
            raise

def main():
    # 只在main函数中配置一次logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.StreamHandler()]  # 只添加一个控制台处理器
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
    
    # Optional parameters
    parser.add_argument("--cecc", help="Path to Cecc.Confirmed.fa.fai")
    parser.add_argument("--mcecc", help="Path to Mecc.Confirmed.fa.fai")
    parser.add_argument("--uecc_inferred", help="Path to Uecc.Inferred.csv")
    parser.add_argument("--output", default="report.txt", help="Path to output report file")
    parser.add_argument("--summary_output", default="summary.csv", help="Path to output summary CSV file")
    parser.add_argument("--html_output", default="report.html", help="Path to output HTML report file")

    try:
        args = parser.parse_args()

        # Create report generator instance
        report_generator = ReportGenerator(
            args.fai, args.ctcr1, args.ctcr2, args.linr,
            args.uecc, args.mecc, args.cecc, args.mcecc,
            args.uecc_inferred, args.output, args.summary_output,
            args.html_output
        )

        # Run the report generation
        report_generator.run()

    except Exception as e:
        logger.error(f"Error: {str(e)}")
        raise

if __name__ == "__main__":
    main()
