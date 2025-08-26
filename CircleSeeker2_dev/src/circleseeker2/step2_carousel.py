#!/usr/bin/env python3
"""
step2_carousel.py - Enhanced Circular DNA Analysis Pipeline

This script processes tandem repeat sequences to identify and classify circular DNA.
It implements intelligent overlap handling, consistency analysis, and multi-pattern recognition.

Usage:
    python step2_carousel.py -i input.txt -o output.csv -f circular.fasta

Author: CircleSeeker Team
Date: 2024
Version: 2.1.0 (Performance optimized and bug fixes)
"""

import argparse
import sys
import logging
import pandas as pd
import numpy as np
import networkx as nx
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
try:
    from intervaltree import IntervalTree
    INTERVALTREE_AVAILABLE = True
except ImportError:
    INTERVALTREE_AVAILABLE = False
    logging.warning("intervaltree not installed, falling back to O(n^2) overlap detection")


class Carousel:
    """Enhanced Carousel for circular DNA analysis"""
    
    def __init__(self, input_file, output_file, circular_fasta, logger=None):
        self.input_file = input_file
        self.output_file = output_file
        self.circular_fasta = circular_fasta
        
        # Use provided logger or create a new one
        if logger is None:
            self.logger = logging.getLogger(__name__)
            if not self.logger.handlers:
                handler = logging.StreamHandler()
                formatter = logging.Formatter('[%(asctime)s] %(levelname)s - %(message)s')
                handler.setFormatter(formatter)
                self.logger.addHandler(handler)
                self.logger.setLevel(logging.INFO)
        else:
            self.logger = logger
    
    def read_and_preprocess_data(self):
        """Read and preprocess tandem repeat data"""
        self.logger.info("Reading and preprocessing data...")
        
        # Define column names
        columns = ['readName', 'repN', 'copyNum', 'readLen', 'start', 'end', 
                   'consLen', 'aveMatch', 'fullLen', 'subPos', 'consSeq']
        
        # Read TSV file
        df = pd.read_csv(self.input_file, sep='\t', header=None, names=columns)
        
        # Convert subPos to string and replace commas with pipes
        df['subPos'] = df['subPos'].astype(str).str.replace(',', '|')
        
        # Calculate Effective_Length
        df['Effective_Length'] = ((df['end'] - df['start'] + 1) / df['readLen']) * 100
        
        # Drop fullLen column
        df = df.drop(columns=['fullLen'])
        
        # Filter low quality (aveMatch < 99)
        df = df[df['aveMatch'] >= 99]
        
        # Reset index
        df = df.reset_index(drop=True)
        
        self.logger.info(f"Loaded {len(df)} regions from {df['readName'].nunique()} reads")
        
        return df
    
    def separate_simple_complex(self, df):
        """Separate reads into simple and complex categories"""
        self.logger.info("Separating simple and complex reads...")
        
        # Count occurrences of each readName
        read_counts = df['readName'].value_counts()
        
        # Simple reads: appear only once
        simple_reads = read_counts[read_counts == 1].index
        # Complex reads: appear multiple times
        complex_reads = read_counts[read_counts > 1].index
        
        # Separate data
        df_simple = df[df['readName'].isin(simple_reads)].copy()
        df_complex = df[df['readName'].isin(complex_reads)].copy()
        
        self.logger.info(f"Simple reads: {len(df_simple)}, Complex reads: {df_complex['readName'].nunique()}")
        
        return df_simple, df_complex
    
    def classify_simple_reads(self, df_simple):
        """Classify simple reads based on effective length"""
        self.logger.info("Classifying simple reads...")
        
        # Create classification column
        df_simple['classification'] = 'Other'  # Default to Other
        
        # Classify based on effective length
        # CtcR-perfect: >=99%
        df_simple.loc[df_simple['Effective_Length'] >= 99, 'classification'] = 'CtcR-perfect'
        
        # CtcR-hybrid: 70-99%
        df_simple.loc[(df_simple['Effective_Length'] >= 70) & 
                      (df_simple['Effective_Length'] < 99), 'classification'] = 'CtcR-hybrid'
        
        return df_simple
    
    def process_complex_group_with_graph_optimized(self, group_df):
        """Process single complex read group using optimized overlap detection"""
        # Reset index for operation
        group_df = group_df.reset_index(drop=True)
        
        # Build overlap graph
        G = nx.Graph()
        
        # Add nodes
        for idx in group_df.index:
            G.add_node(idx)
        
        if INTERVALTREE_AVAILABLE and len(group_df) > 10:
            # Use interval tree for better performance on larger groups
            tree = IntervalTree()
            for idx, row in group_df.iterrows():
                tree[row['start']:row['end']+1] = idx
            
            # Find overlaps efficiently
            processed_pairs = set()
            for idx, row in group_df.iterrows():
                overlaps = tree[row['start']:row['end']+1]
                for interval in overlaps:
                    other_idx = interval.data
                    if other_idx != idx:
                        pair = tuple(sorted([idx, other_idx]))
                        if pair not in processed_pairs:
                            other_row = group_df.loc[other_idx]
                            overlap = min(row['end'], other_row['end']) - max(row['start'], other_row['start']) + 1
                            
                            if overlap > 50:
                                G.add_edge(idx, other_idx, overlap=overlap)
                                processed_pairs.add(pair)
        else:
            # Fallback to O(n^2) method for small groups or when intervaltree not available
            for i in group_df.index:
                for j in group_df.index:
                    if i < j:
                        # Check for overlap (fixed calculation with +1)
                        overlap = min(group_df.loc[i, 'end'], group_df.loc[j, 'end']) - \
                                 max(group_df.loc[i, 'start'], group_df.loc[j, 'start']) + 1
                        
                        # Significant overlap if > 50bp
                        if overlap > 50:
                            G.add_edge(i, j, overlap=overlap)
        
        # Find all connected components
        selected_indices = []
        
        for component in nx.connected_components(G):
            if len(component) == 1:
                # No overlap, keep directly
                selected_indices.extend(component)
            else:
                # Has overlap, keep the longest coverage
                component_df = group_df.loc[list(component)]
                component_df = component_df.copy()
                component_df['coverage_length'] = component_df['end'] - component_df['start'] + 1
                best_idx = component_df['coverage_length'].idxmax()
                selected_indices.append(best_idx)
        
        # Return processed data
        processed_df = group_df.loc[selected_indices].copy()
        
        return processed_df
    
    def process_all_complex_reads(self, df_complex):
        """Process all complex reads"""
        self.logger.info("Processing complex reads with optimized algorithm...")
        
        processed_groups = []
        
        # Process by readName groups
        for readName, group in df_complex.groupby('readName'):
            # Process this group
            processed_group = self.process_complex_group_with_graph_optimized(group)
            
            # Calculate total effective length for the group
            group_total_effective_length = processed_group['Effective_Length'].sum()
            
            # Add group total effective length to each row
            processed_group['group_total_effective_length'] = group_total_effective_length
            
            # Save processed group
            processed_groups.append(processed_group)
        
        # Merge all processed groups
        df_complex_processed = pd.concat(processed_groups, ignore_index=True)
        
        self.logger.info(f"Processed: {len(df_complex)} -> {len(df_complex_processed)} regions")
        
        return df_complex_processed
    
    def analyze_consLen_consistency(self, df_complex_processed):
        """Analyze consLen consistency within each group"""
        self.logger.info("Analyzing consLen consistency...")
        
        consistency_results = []
        
        # Analyze by readName groups
        for readName, group in df_complex_processed.groupby('readName'):
            consLen_values = group['consLen'].values
            
            # Calculate statistics
            mean_consLen = np.mean(consLen_values)
            std_consLen = np.std(consLen_values)
            min_consLen = np.min(consLen_values)
            max_consLen = np.max(consLen_values)
            
            # Coefficient of variation (CV)
            cv = (std_consLen / mean_consLen * 100) if mean_consLen > 0 else 0
            
            # Max/min ratio (with zero protection)
            ratio_max_min = max_consLen / min_consLen if min_consLen > 0 else np.inf
            
            # Check for integer multiple relationship (with zero protection)
            if min_consLen > 0:
                base_unit = min_consLen
                is_integer_multiple = True
                multiples = []
                for consLen in consLen_values:
                    multiple = consLen / base_unit
                    multiples.append(multiple)
                    if abs(multiple - round(multiple)) > 0.05:  # 5% tolerance
                        is_integer_multiple = False
            else:
                is_integer_multiple = False
                multiples = None
            
            consistency_results.append({
                'readName': readName,
                'num_regions': len(group),
                'consLen_list': list(consLen_values),
                'mean_consLen': mean_consLen,
                'std_consLen': std_consLen,
                'cv_percent': cv,
                'min_consLen': min_consLen,
                'max_consLen': max_consLen,
                'ratio_max_min': ratio_max_min,
                'is_integer_multiple': is_integer_multiple,
                'multiples': multiples if is_integer_multiple else None,
                'group_total_effective_length': group['group_total_effective_length'].iloc[0]
            })
        
        # Convert to DataFrame
        consistency_df = pd.DataFrame(consistency_results)
        
        # Add consistency classification
        def classify_consistency(row):
            if row['cv_percent'] < 1:  # CV < 1%
                return 'highly_consistent'
            elif row['cv_percent'] < 5:  # CV < 5%
                return 'consistent'
            elif row['is_integer_multiple']:
                return 'integer_multiple'
            else:
                return 'variable'
        
        consistency_df['consistency_type'] = consistency_df.apply(classify_consistency, axis=1)
        
        # Statistics
        highly_consistent = len(consistency_df[consistency_df['consistency_type'] == 'highly_consistent'])
        self.logger.info(f"Highly consistent groups: {highly_consistent}")
        
        return consistency_df
    
    def process_and_classify_complex_reads(self, df_complex_processed, consistency_analysis):
        """Process and classify complex reads"""
        self.logger.info("Classifying complex reads...")
        
        # Get highly consistent groups with multiple records
        highly_consistent_multi = consistency_analysis[
            (consistency_analysis['consistency_type'] == 'highly_consistent') & 
            (consistency_analysis['num_regions'] > 1)
        ]['readName'].tolist()
        
        processed_results = []
        
        # Process by readName groups
        for readName, group in df_complex_processed.groupby('readName'):
            
            if readName in highly_consistent_multi:
                # Highly consistent multi-record groups: merge to single record
                first_row = group.iloc[0].copy()
                
                # Sum Effective_Length and copyNum (with proper rounding)
                first_row['Effective_Length'] = group['Effective_Length'].sum()
                first_row['copyNum'] = np.round(group['copyNum'].sum())
                
                # Classify
                if first_row['Effective_Length'] >= 70:
                    first_row['classification'] = 'CtcR-inversion'
                else:
                    first_row['classification'] = 'Other'
                
                processed_results.append(first_row)
                
            else:
                # Other groups
                if len(group) == 1:
                    # Single record, use simple read rules
                    row = group.iloc[0].copy()
                    eff_length = row['Effective_Length']
                    
                    if eff_length >= 99:
                        row['classification'] = 'CtcR-perfect'
                    elif eff_length >= 70:
                        row['classification'] = 'CtcR-hybrid'
                    else:
                        row['classification'] = 'Other'
                    
                    processed_results.append(row)
                    
                else:
                    # Multiple records but not highly consistent, mark all as hybrid
                    for _, row in group.iterrows():
                        row_copy = row.copy()
                        row_copy['classification'] = 'CtcR-hybrid'
                        processed_results.append(row_copy)
        
        # Merge results
        df_complex_final = pd.DataFrame(processed_results)
        
        return df_complex_final
    
    def create_main_results_table(self, df_simple_classified, df_complex_final):
        """Create main results table"""
        # Merge simple and complex read results
        df_main = pd.concat([df_simple_classified, df_complex_final], ignore_index=True)
        
        # Convert copyNum to integer with proper rounding
        df_main['copyNum'] = np.round(df_main['copyNum']).astype(int)
        
        # Create unique ID column (with proper type conversion)
        df_main['unique_id'] = (
            df_main['readName'].astype(str) + '|' + 
            df_main['repN'].astype(str) + '|' + 
            df_main['consLen'].astype(str) + '|' + 
            df_main['copyNum'].astype(str)
        )
        
        # Move unique_id to first column
        cols = df_main.columns.tolist()
        cols = ['unique_id'] + [col for col in cols if col != 'unique_id']
        df_main = df_main[cols]
        
        return df_main
    
    def circularize_sequences(self, df_main):
        """Generate circularized sequences"""
        circular_sequences = []
        
        for _, row in df_main.iterrows():
            # Get sequence
            seq = row['consSeq']
            
            # Circularize: concatenate sequence with itself
            circular_seq = seq + seq
            
            # Create SeqRecord object
            record = SeqRecord(
                Seq(circular_seq),
                id=row['unique_id'] + '|circular',
                description=""
            )
            
            circular_sequences.append(record)
        
        return circular_sequences
    
    def write_fasta(self, sequences, output_file):
        """Write sequences to FASTA file"""
        SeqIO.write(sequences, output_file, "fasta")
    
    def create_readname_classification(self, df_main):
        """
        Create readName-based classification table
        Each readName appears only once with its classification
        """
        # Get readName and classification
        df_readname = df_main[['readName', 'classification']].copy()
        
        # For reads with multiple records, determine final classification
        # Priority: CtcR-perfect > CtcR-inversion > CtcR-hybrid > Other
        priority_map = {
            'CtcR-perfect': 1,
            'CtcR-inversion': 2,
            'CtcR-hybrid': 3,
            'Other': 4
        }
        
        # Add priority column
        df_readname['priority'] = df_readname['classification'].map(priority_map)
        
        # Sort by readName and priority
        df_readname = df_readname.sort_values(['readName', 'priority'])
        
        # Keep the highest priority classification for each readName
        df_readname = df_readname.drop_duplicates(subset=['readName'], keep='first')
        
        # Drop priority column
        df_readname = df_readname.drop(columns=['priority'])
        
        # Rename columns to match expected format
        df_readname.columns = ['readName', 'readClass']
        
        return df_readname
    
    def process(self):
        """Main processing pipeline"""
        self.logger.info("=" * 60)
        self.logger.info("Carousel - Enhanced Circular DNA Analysis")
        self.logger.info("=" * 60)
        
        # 1. Read and preprocess
        df = self.read_and_preprocess_data()
        
        # 2. Separate simple and complex reads
        df_simple, df_complex = self.separate_simple_complex(df)
        
        # 3. Process simple reads
        df_simple = self.classify_simple_reads(df_simple)
        
        # 4. Process complex reads
        if not df_complex.empty:
            # Process overlaps with optimized algorithm
            df_complex_processed = self.process_all_complex_reads(df_complex)
            
            # Analyze consistency
            consistency_analysis = self.analyze_consLen_consistency(df_complex_processed)
            
            # Classify complex reads
            df_complex_final = self.process_and_classify_complex_reads(
                df_complex_processed, consistency_analysis
            )
        else:
            df_complex_final = pd.DataFrame()
        
        # 5. Create main results table
        df_main = self.create_main_results_table(df_simple, df_complex_final)
        
        # 6. Create readName-based classification (one entry per readName)
        df_classification = self.create_readname_classification(df_main)
        
        # 7. Save classification results
        df_classification.to_csv(self.output_file, index=False)
        self.logger.info(f"Classification results saved to '{self.output_file}'")
        
        # 8. Generate circular FASTA
        circular_sequences = self.circularize_sequences(df_main)
        self.write_fasta(circular_sequences, self.circular_fasta)
        self.logger.info(f"Circular sequences saved to '{self.circular_fasta}'")
        
        # 9. Display statistics
        self.logger.info("=" * 60)
        self.logger.info("Final Classification Statistics:")
        self.logger.info("-" * 60)
        
        # Statistics from df_classification (unique readNames)
        for class_name, count in df_classification['readClass'].value_counts().items():
            self.logger.info(f"{class_name}: {count}")
        
        self.logger.info("-" * 60)
        self.logger.info(f"Total unique reads: {len(df_classification)}")
        self.logger.info(f"Total records in main table: {len(df_main)}")
        self.logger.info("=" * 60)
        
        # Additional detailed statistics
        if len(df_main) > len(df_classification):
            self.logger.info("Note: Some reads have multiple records in the main table")
            self.logger.info("The classification file contains one entry per unique readName")
        
        return df_main, df_classification


def setup_logging(level=logging.INFO):
    """Setup logging configuration"""
    logging.basicConfig(
        level=level,
        format='[%(asctime)s] %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    return logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description="Carousel - Enhanced Circular DNA Analysis Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Output Files:
  1. Classification CSV: Contains readName and readClass columns
     - Each readName appears only once
     - For reads with multiple regions, the highest priority classification is used
     - Priority: CtcR-perfect > CtcR-inversion > CtcR-hybrid > Other
  
  2. Circular FASTA: Contains circularized sequences
     - May contain multiple sequences for complex reads
     - Sequence IDs include region information
        
Examples:
  %(prog)s -i tandem_repeats.txt -o classification.csv -f circular.fasta
  %(prog)s -i input.txt -o output.csv -f seqs.fasta -v  # verbose mode
        """
    )
    
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input TSV file with tandem repeat data"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output CSV file with read classifications (readName,readClass format)"
    )
    parser.add_argument(
        "-f", "--fasta",
        required=True,
        help="Output FASTA file with circular sequences"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose output"
    )
    
    args = parser.parse_args()
    
    # Setup logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = setup_logging(log_level)
    
    # Check for optional dependencies
    if not INTERVALTREE_AVAILABLE:
        logger.info("Note: Install 'intervaltree' package for better performance with large datasets")
        logger.info("  pip install intervaltree")
    
    # Create and run Carousel
    carousel = Carousel(
        input_file=args.input,
        output_file=args.output,
        circular_fasta=args.fasta,
        logger=logger
    )
    
    try:
        df_main, df_classification = carousel.process()
        logger.info("Analysis completed successfully!")
        
        # Verify output format
        if 'readName' not in df_classification.columns or 'readClass' not in df_classification.columns:
            logger.error("Output format error: Missing required columns")
            sys.exit(1)
            
        # Check for duplicates
        duplicate_count = df_classification['readName'].duplicated().sum()
        if duplicate_count > 0:
            logger.warning(f"Warning: Found {duplicate_count} duplicate readNames in output")
        
    except Exception as e:
        logger.error(f"Error during processing: {e}")
        import traceback
        logger.debug(traceback.format_exc())
        sys.exit(1)


if __name__ == "__main__":
    main()
