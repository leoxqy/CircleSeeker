"""
Gatekeeper - Simple eccDNA Classification Module

This module performs initial classification of eccDNA from BLAST results:
- Uecc: Unique/simple circular DNA 
- Mecc: Multiple-copy repeat circular DNA
- Unclassified: All alignments from queries not classified as Uecc/Mecc

Migrated from step4_gatekeeper.py to the new CircleSeeker2 architecture.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Tuple, List, Set, Optional
import pandas as pd
import numpy as np


class GatekeeperClassifier:
    """Gatekeeper classifier for eccDNA"""
    
    # BLAST column names
    BLAST_COLUMNS = [
        'query_id', 'subject_id', 'identity', 'alignment_length',
        'mismatches', 'gap_opens', 'q_start', 'q_end', 's_start',
        's_end', 'evalue', 'bit_score', 'sstrand'
    ]
    
    def __init__(self, gap_threshold: float = 10.0, min_full_length_coverage: float = 95.0,
                 logger: Optional[logging.Logger] = None):
        """
        Initialize Gatekeeper classifier
        
        Args:
            gap_threshold: Maximum gap percentage for quality filtering (default: 10%)
            min_full_length_coverage: Minimum coverage for full-length classification (default: 95%)
            logger: Optional logger instance
        """
        self.gap_threshold = gap_threshold
        self.min_full_length_coverage = min_full_length_coverage
        self.stats = {}
        
        # Setup logger
        self.logger = logger or logging.getLogger(self.__class__.__name__)
    
    def read_blast_results(self, blast_file: Path) -> pd.DataFrame:
        """
        Read and preprocess BLAST results
        
        Args:
            blast_file: Path to BLAST results file
            
        Returns:
            Preprocessed DataFrame
        """
        self.logger.info(f"Reading BLAST results from {blast_file}")
        
        try:
            df = pd.read_csv(blast_file, sep='\t', header=None, names=self.BLAST_COLUMNS)
        except pd.errors.EmptyDataError:
            self.logger.error(f"BLAST file is empty: {blast_file}")
            raise
        except Exception as e:
            self.logger.error(f"Failed to read BLAST file: {e}")
            raise
            
        # Process strand information
        df['strand'] = df['sstrand'].apply(lambda x: '+' if x == 'plus' else '-')
        neg_strand_mask = df['strand'] == '-'
        df.loc[neg_strand_mask, ['s_start', 's_end']] = df.loc[neg_strand_mask, ['s_end', 's_start']].values
        
        # Parse query IDs with proper type conversion
        try:
            split_cols = df['query_id'].str.split('|', expand=True)
            df['readName'] = split_cols[0].astype(str)
            df['consLen'] = pd.to_numeric(split_cols[2], errors='coerce').fillna(0).astype(int)
            df['copyNum'] = pd.to_numeric(split_cols[3], errors='coerce').fillna(0).astype(float)
            
            # Validate consLen to avoid division by zero later
            invalid_conslen = df['consLen'] <= 0
            if invalid_conslen.any():
                self.logger.warning(f"Found {invalid_conslen.sum()} entries with invalid consLen (<=0)")
                # Filter out invalid entries
                df = df[~invalid_conslen].copy()
                
        except Exception as e:
            self.logger.error(f"Failed to parse query IDs: {e}")
            raise
            
        # Calculate derived columns with safe division
        df['Rlength'] = df['s_end'] - df['s_start'] + 1
        df['gap_Length'] = df['consLen'] - df['Rlength']
        
        # Safe division for Gap_Percentage
        df['Gap_Percentage'] = 0.0
        valid_mask = df['consLen'] > 0
        df.loc[valid_mask, 'Gap_Percentage'] = (
            (df.loc[valid_mask, 'gap_Length'].abs() / df.loc[valid_mask, 'consLen']) * 100
        ).round(2)
        
        self.logger.info(f"Loaded {len(df):,} alignments from {df['query_id'].nunique():,} queries")
        
        # Store statistics
        self.stats['total_alignments'] = len(df)
        self.stats['total_queries'] = df['query_id'].nunique()
        
        return df
    
    def classify_uecc_mecc(self, df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, Set[str]]:
        """
        Classify Uecc and Mecc based on high-quality alignments
        
        Args:
            df: Full DataFrame
            
        Returns:
            Tuple of (uecc_df, mecc_df, classified_query_ids)
        """
        self.logger.info("=" * 60)
        self.logger.info("Step 1: Classify using high-quality alignments")
        
        # Filter by quality
        high_quality = df[df['Gap_Percentage'] <= self.gap_threshold].copy()
        self.logger.info(f"High-quality alignments: {len(high_quality):,} / {len(df):,}")
        
        if high_quality.empty:
            self.logger.warning("No high-quality alignments found")
            return pd.DataFrame(), pd.DataFrame(), set()
        
        # Process each query based on high-quality alignments
        uecc_list = []
        mecc_list = []
        classified_queries = set()
        
        for query_id, group in high_quality.groupby('query_id'):
            if len(group) == 1:
                # Single alignment -> Uecc
                group_copy = group.copy()
                group_copy['eClass'] = 'Uecc'
                group_copy['classification_reason'] = 'Single alignment'
                uecc_list.append(group_copy)
                classified_queries.add(query_id)
            else:
                # Multiple alignments - check for overlaps
                processed_group = self._process_overlaps_for_query(group)
                
                if len(processed_group) == 1:
                    # After overlap removal -> Uecc
                    processed_copy = processed_group.copy()
                    processed_copy['eClass'] = 'Uecc'
                    processed_copy['classification_reason'] = 'Single after overlap removal'
                    uecc_list.append(processed_copy)
                    classified_queries.add(query_id)
                else:
                    # Check if it's Mecc (full-length repeats)
                    if self._is_full_length_repeat(processed_group):
                        processed_copy = processed_group.copy()
                        processed_copy['eClass'] = 'Mecc'
                        mecc_list.append(processed_copy)
                        classified_queries.add(query_id)
                    # If not Mecc, it remains unclassified
        
        # Combine results
        uecc_df = pd.concat(uecc_list, ignore_index=False) if uecc_list else pd.DataFrame()
        mecc_df = pd.concat(mecc_list, ignore_index=False) if mecc_list else pd.DataFrame()
        
        self.logger.info(f"Uecc: {len(classified_queries & set(uecc_df['query_id'].unique() if not uecc_df.empty else []))} queries")
        self.logger.info(f"Mecc: {len(classified_queries & set(mecc_df['query_id'].unique() if not mecc_df.empty else []))} queries")
        
        return uecc_df, mecc_df, classified_queries
    
    def _process_overlaps_for_query(self, group: pd.DataFrame) -> pd.DataFrame:
        """
        Process overlaps for a single query's alignments
        
        Args:
            group: DataFrame for single query
            
        Returns:
            Processed DataFrame with overlaps resolved
        """
        kept_alignments = []
        
        for chr_id, chr_group in group.groupby('subject_id'):
            overlapping_indices = self._find_overlaps_sweepline(chr_group)
            
            if overlapping_indices:
                # Select best from overlapping group
                overlap_df = chr_group.loc[list(overlapping_indices)]
                best_idx = overlap_df['Gap_Percentage'].idxmin()
                kept_alignments.append(chr_group.loc[[best_idx]])
                
                # Keep non-overlapping
                non_overlap = chr_group[~chr_group.index.isin(overlapping_indices)]
                if len(non_overlap) > 0:
                    kept_alignments.append(non_overlap)
            else:
                kept_alignments.append(chr_group)
        
        if kept_alignments:
            return pd.concat(kept_alignments, ignore_index=False)
        return pd.DataFrame()
    
    def _find_overlaps_sweepline(self, group: pd.DataFrame) -> Set[int]:
        """
        Find overlapping intervals using sweep-line algorithm
        
        Args:
            group: DataFrame with alignments on same chromosome
            
        Returns:
            Set of overlapping indices
        """
        if len(group) <= 1:
            return set()
        
        events = []
        for idx, row in group.iterrows():
            events.append((row['s_start'], 0, idx))  # 0 for start event
            events.append((row['s_end'], 1, idx))    # 1 for end event
        
        # Sort by position, then by event type (starts before ends)
        events.sort(key=lambda x: (x[0], x[1]))
        
        active_intervals = set()
        overlapping_indices = set()
        
        for pos, event_type, idx in events:
            if event_type == 0:  # start event
                if active_intervals:
                    overlapping_indices.add(idx)
                    overlapping_indices.update(active_intervals)
                active_intervals.add(idx)
            else:  # end event
                active_intervals.discard(idx)
        
        return overlapping_indices
    
    def _is_full_length_repeat(self, group: pd.DataFrame) -> bool:
        """
        Check if alignments represent full-length repeats (Mecc)
        
        Args:
            group: DataFrame with multiple alignments for one query
            
        Returns:
            True if at least 2 full-length copies
        """
        # Filter out entries with zero or invalid consLen
        valid_group = group[group['consLen'] > 0]
        
        if valid_group.empty or len(valid_group) < 2:
            return False
        
        # Calculate coverage for each alignment
        coverages = (valid_group['Rlength'] / valid_group['consLen']) * 100
        full_length_count = (coverages >= self.min_full_length_coverage).sum()
        
        return full_length_count >= 2
    
    def extract_unclassified(self, df_original: pd.DataFrame, classified_queries: Set[str]) -> pd.DataFrame:
        """
        Extract all alignments from unclassified queries
        
        Args:
            df_original: Original complete DataFrame
            classified_queries: Set of query IDs that were classified
            
        Returns:
            DataFrame with all alignments from unclassified queries
        """
        self.logger.info("=" * 60)
        self.logger.info("Step 2: Extract unclassified queries")
        
        all_queries = set(df_original['query_id'].unique())
        unclassified_queries = all_queries - classified_queries
        
        self.logger.info(f"Total queries: {len(all_queries):,}")
        self.logger.info(f"Classified queries: {len(classified_queries):,}")
        self.logger.info(f"Unclassified queries: {len(unclassified_queries):,}")
        
        # Extract ALL alignments for unclassified queries
        unclassified_df = df_original[df_original['query_id'].isin(unclassified_queries)].copy()
        
        if not unclassified_df.empty:
            # Add reason column for downstream analysis
            unclassified_df['unclass_reason'] = 'Not_classified'
            
            # Add quality indicator
            unclassified_df['quality_category'] = unclassified_df['Gap_Percentage'].apply(
                lambda x: 'High_quality' if x <= self.gap_threshold else 'Low_quality'
            )
            
            # Statistics
            hq_count = (unclassified_df['quality_category'] == 'High_quality').sum()
            lq_count = (unclassified_df['quality_category'] == 'Low_quality').sum()
            
            self.logger.info(f"Unclassified alignments: {len(unclassified_df):,}")
            self.logger.info(f"  - High quality (Gap <= {self.gap_threshold}%): {hq_count:,}")
            self.logger.info(f"  - Low quality (Gap > {self.gap_threshold}%): {lq_count:,}")
        
        return unclassified_df
    
    def format_outputs(self, uecc_df: pd.DataFrame, mecc_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Format output DataFrames with appropriate columns
        
        Args:
            uecc_df: Uecc DataFrame
            mecc_df: Mecc DataFrame
            
        Returns:
            Tuple of formatted (uecc_df, mecc_df)
        """
        # Format Uecc
        if not uecc_df.empty:
            uecc_formatted = uecc_df.copy()
            uecc_formatted['MatDegree'] = 100 - uecc_formatted['Gap_Percentage']
            uecc_formatted['eName'] = (
                uecc_formatted['subject_id'].astype(str) + "-" + 
                uecc_formatted['s_start'].astype(str) + "-" + 
                uecc_formatted['s_end'].astype(str)
            )
            
            # Rename columns
            uecc_formatted.rename(columns={
                'subject_id': 'eChr',
                's_start': 'eStart',
                's_end': 'eEnd',
                'strand': 'eStrand',
                'consLen': 'eLength',
                'copyNum': 'eRepeatNum',
                'readName': 'eReads'
            }, inplace=True)
        else:
            uecc_formatted = uecc_df
        
        # Format Mecc
        if not mecc_df.empty:
            mecc_formatted = mecc_df.copy()
            mecc_formatted.rename(columns={
                'subject_id': 'eChr',
                's_start': 'eStart', 
                's_end': 'eEnd',
                'strand': 'eStrand',
                'consLen': 'eLength'
            }, inplace=True)
        else:
            mecc_formatted = mecc_df
        
        return uecc_formatted, mecc_formatted
    
    def run(self, blast_file: Path) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        Run complete Gatekeeper classification
        
        Args:
            blast_file: Path to BLAST results
            
        Returns:
            Tuple of (uecc_df, mecc_df, unclassified_df)
        """
        self.logger.info("=" * 60)
        self.logger.info("Starting Gatekeeper classification")
        self.logger.info("=" * 60)
        
        # Read and preprocess BLAST results
        df_original = self.read_blast_results(blast_file)
        
        # Check if we have valid data after preprocessing
        if df_original.empty:
            self.logger.warning("No valid data after preprocessing")
            return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
        
        # Step 1: Classify Uecc and Mecc based on high-quality alignments
        uecc_df, mecc_df, classified_queries = self.classify_uecc_mecc(df_original)
        
        # Step 2: Extract ALL alignments from unclassified queries
        unclassified_df = self.extract_unclassified(df_original, classified_queries)
        
        # Format outputs
        uecc_df, mecc_df = self.format_outputs(uecc_df, mecc_df)
        
        # Final statistics
        uecc_count = uecc_df['query_id'].nunique() if not uecc_df.empty else 0
        mecc_count = mecc_df['query_id'].nunique() if not mecc_df.empty else 0
        unclassified_count = unclassified_df['query_id'].nunique() if not unclassified_df.empty else 0
        
        # Quality breakdown for unclassified
        unclass_hq_queries = 0
        unclass_lq_only_queries = 0
        if not unclassified_df.empty:
            for query_id, group in unclassified_df.groupby('query_id'):
                if (group['quality_category'] == 'High_quality').any():
                    unclass_hq_queries += 1
                else:
                    unclass_lq_only_queries += 1
        
        self.logger.info("=" * 60)
        self.logger.info("Classification Summary:")
        self.logger.info(f"  Uecc: {uecc_count:,} queries")
        self.logger.info(f"  Mecc: {mecc_count:,} queries")
        self.logger.info(f"  Unclassified: {unclassified_count:,} queries")
        if unclassified_count > 0:
            self.logger.info(f"    - With high-quality alignments: {unclass_hq_queries:,}")
            self.logger.info(f"    - Only low-quality alignments: {unclass_lq_only_queries:,}")
        self.logger.info("=" * 60)
        
        return uecc_df, mecc_df, unclassified_df