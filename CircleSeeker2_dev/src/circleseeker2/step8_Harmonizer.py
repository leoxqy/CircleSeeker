#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
eccDNA_cluster_dedup_no_priority.py

Independent deduplication for each eccDNA type using their respective CD-HIT clusters.
No cross-type priority judgments - each type is processed completely independently.

Key features:
- Each type (Uecc/Mecc/Cecc) uses its own .clstr file for deduplication
- Uecc: keeps 1 representative row per cluster, aggregates metadata
- Mecc/Cecc: keeps ALL rows of representative ID, aggregates metadata from cluster members
- Outputs 11 standardized files maintaining compatibility with existing pipeline

Author: Your Team
Date: 2025
"""

import argparse
import logging
import sys
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set
import pandas as pd
import numpy as np

logger = logging.getLogger("eccDNA_dedup")

# ================== Configuration ==================
FIXED_FRONT_COLUMNS = ['eccDNA_id', 'eChr', 'eStart0', 'eEnd0', 'eStrand', 'eLength', 'MatDegree', 'copyNum']

# ================== Logging Setup ==================
def setup_logging(level: str) -> None:
    """Configure logging with specified level."""
    lv = getattr(logging, level.upper(), None)
    if not isinstance(lv, int):
        raise ValueError(f"Invalid log level: {level}")
    logging.basicConfig(
        level=lv,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%H:%M:%S"
    )

# ================== CD-HIT Cluster Parser ==================
class CDHitClusters:
    """Parse and manage CD-HIT cluster information."""
    
    def __init__(self) -> None:
        self.members: Dict[str, List[str]] = {}  # cluster_id -> member_ids
        self.rep_of: Dict[str, str] = {}  # cluster_id -> representative_id
    
    def parse(self, clstr_file: Path) -> None:
        """Parse CD-HIT .clstr file."""
        logger.info(f"Parsing CD-HIT cluster file: {clstr_file}")
        
        if not clstr_file.exists():
            raise FileNotFoundError(f"Cluster file not found: {clstr_file}")
        
        current_cluster_id: Optional[str] = None
        current_members: List[str] = []
        current_rep: Optional[str] = None
        
        with open(clstr_file, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                # New cluster header
                if line.startswith(">Cluster"):
                    # Save previous cluster if exists
                    if current_cluster_id is not None and current_members:
                        self.members[current_cluster_id] = current_members.copy()
                        if current_rep:
                            self.rep_of[current_cluster_id] = current_rep
                    
                    # Parse cluster ID
                    tokens = line.split()
                    try:
                        current_cluster_id = str(int(tokens[1]))
                    except (IndexError, ValueError):
                        current_cluster_id = line.split(">Cluster", 1)[1].strip()
                    
                    current_members = []
                    current_rep = None
                    continue
                
                # Member line - extract ID
                # Pattern: "0       12634nt, >U23805... *"
                match = re.search(r'>\s*([^\.>\s][^\.]*?)\.\.\.', line)
                if not match:
                    match = re.search(r'>\s*([^\s,]+)', line)
                
                if match:
                    member_id = match.group(1)
                    current_members.append(member_id)
                    
                    # Check if this is the representative (marked with *)
                    if line.endswith("*"):
                        current_rep = member_id
        
        # Save last cluster
        if current_cluster_id is not None and current_members:
            self.members[current_cluster_id] = current_members.copy()
            if current_rep:
                self.rep_of[current_cluster_id] = current_rep
        
        logger.info(f"Parsed {len(self.members)} clusters, {len(self.rep_of)} with representatives")
    
    def get_id_to_cluster_map(self) -> Dict[str, str]:
        """Create mapping from member ID to cluster ID."""
        mapping = {}
        for cluster_id, member_list in self.members.items():
            for member_id in member_list:
                mapping[member_id] = cluster_id
        return mapping

# ================== Helper Functions ==================
def to_numeric_safe(series, default=0):
    """Safely convert series to numeric, filling NaN with default."""
    try:
        return pd.to_numeric(series, errors='coerce').fillna(default)
    except Exception:
        return pd.Series([default] * len(series))

def format_two_decimals(value) -> str:
    """Format value to 2 decimal places."""
    try:
        if pd.isna(value):
            return ""
        return f"{float(value):.2f}"
    except Exception:
        return ""

def concat_unique_semicolon(series: pd.Series) -> str:
    """Concatenate unique values with semicolon."""
    unique_vals = series.dropna().astype(str).unique()
    return ';'.join(unique_vals) if len(unique_vals) else ""

def merge_read_lists(series: pd.Series) -> str:
    """Merge semicolon-separated read lists, removing duplicates."""
    all_reads = set()
    for reads_str in series.dropna().astype(str):
        if reads_str:
            tokens = [t.strip() for t in reads_str.split(';') if t.strip()]
            all_reads.update(tokens)
    return ';'.join(sorted(all_reads)) if all_reads else ""

def majority_vote(series: pd.Series) -> Optional[str]:
    """Return the most common value in series."""
    clean_series = series.dropna().astype(str)
    if clean_series.empty:
        return None
    return clean_series.value_counts().idxmax()

def count_reads_from_string(s: str) -> int:
    """Count unique reads from semicolon-separated string."""
    if pd.isna(s) or not isinstance(s, str) or not s.strip():
        return 0
    return len({x.strip() for x in s.split(';') if x.strip()})

def ensure_int64_column(df: pd.DataFrame, col: str) -> None:
    """Ensure column is Int64 type."""
    if col not in df.columns:
        df[col] = pd.Series([pd.NA] * len(df), dtype="Int64")
    else:
        try:
            df[col] = pd.to_numeric(df[col], errors='coerce').astype("Int64")
        except Exception:
            df[col] = pd.Series([pd.NA] * len(df), dtype="Int64")

# ================== Data Processing Functions ==================
def normalize_coordinates(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize coordinate columns across different input formats."""
    df = df.copy()
    
    # Handle chromosome column
    if 'eChr' not in df.columns:
        if 'subject_id' in df.columns:
            df['eChr'] = df['subject_id']
        elif 'chr' in df.columns:
            df['eChr'] = df['chr']
    
    # Handle start coordinate (0-based)
    if 'eStart0' not in df.columns:
        if 'eStart' in df.columns:
            df['eStart0'] = to_numeric_safe(df['eStart'], 0) - 1
        elif 's_start' in df.columns:
            df['eStart0'] = to_numeric_safe(df['s_start'], 0) - 1
        elif 'start' in df.columns:
            df['eStart0'] = to_numeric_safe(df['start'], 0) - 1
    
    # Handle end coordinate
    if 'eEnd0' not in df.columns:
        if 'eEnd' in df.columns:
            df['eEnd0'] = to_numeric_safe(df['eEnd'], 0)
        elif 's_end' in df.columns:
            df['eEnd0'] = to_numeric_safe(df['s_end'], 0)
        elif 'end' in df.columns:
            df['eEnd0'] = to_numeric_safe(df['end'], 0)
    
    ensure_int64_column(df, 'eStart0')
    ensure_int64_column(df, 'eEnd0')
    
    # Handle strand
    if 'eStrand' not in df.columns:
        if 'strand' in df.columns:
            df['eStrand'] = df['strand']
        elif 'sstrand' in df.columns:
            df['eStrand'] = df['sstrand']
    
    # Handle length
    if 'eLength' not in df.columns:
        if 'consLen' in df.columns:
            df['eLength'] = to_numeric_safe(df['consLen'], 0)
        else:
            df['eLength'] = (pd.to_numeric(df['eEnd0'], errors='coerce') - 
                           pd.to_numeric(df['eStart0'], errors='coerce'))
    ensure_int64_column(df, 'eLength')
    
    # Handle copyNum
    if 'copyNum' not in df.columns:
        df['copyNum'] = pd.NA
    
    # Handle MatDegree
    if 'MatDegree' not in df.columns:
        if 'identity' in df.columns:
            df['MatDegree'] = to_numeric_safe(df['identity'], 0)
        elif 'Gap_Percentage' in df.columns:
            df['MatDegree'] = 100 - to_numeric_safe(df['Gap_Percentage'], 0)
        else:
            df['MatDegree'] = pd.NA
    
    return df

def reorder_columns_for_output(df: pd.DataFrame) -> pd.DataFrame:
    """Reorder columns: fixed front, others, then reads_count and readName at end."""
    df = df.copy()
    
    # Ensure required columns exist
    if 'readName' not in df.columns:
        df['readName'] = ""
    if 'reads_count' not in df.columns:
        df['reads_count'] = pd.NA
    
    # Separate columns
    front_cols = [c for c in FIXED_FRONT_COLUMNS if c in df.columns]
    all_cols = set(df.columns)
    middle_cols = [c for c in df.columns if c not in front_cols and c not in ('reads_count', 'readName')]
    
    # Reorder
    ordered_cols = front_cols + middle_cols + ['reads_count', 'readName']
    return df[ordered_cols]

# ================== Type-Specific Processing ==================
def process_uecc(df: pd.DataFrame, clusters: CDHitClusters, dtype: str) -> pd.DataFrame:
    """Process Uecc data: one representative row per cluster."""
    if df.empty:
        return df
    
    df = df.copy()
    id_to_cluster = clusters.get_id_to_cluster_map()
    df['cluster_id'] = df['eccDNA_id'].map(id_to_cluster)
    
    # Handle singletons
    singleton_mask = df['cluster_id'].isna()
    df.loc[singleton_mask, 'cluster_id'] = 'singleton:' + df.loc[singleton_mask, 'eccDNA_id'].astype(str)
    
    # Select representative for each cluster
    rep_map = clusters.rep_of
    
    # Calculate metrics for selection
    if 'eLength' in df.columns:
        metric = df.groupby(['cluster_id', 'eccDNA_id'])['eLength'].max()
    else:
        metric = df.groupby(['cluster_id', 'eccDNA_id']).size()
    
    metric_df = metric.reset_index().rename(columns={metric.name: 'metric'})
    
    # Add representative flag
    rep_df = pd.DataFrame({'cluster_id': list(rep_map.keys()), 'rep_id': list(rep_map.values())})
    metric_df = metric_df.merge(rep_df, on='cluster_id', how='left')
    metric_df['is_rep'] = (metric_df['eccDNA_id'] == metric_df['rep_id']).astype(int)
    
    # Sort to prioritize: representative > metric
    metric_df = metric_df.sort_values(['cluster_id', 'is_rep', 'metric'], 
                                     ascending=[True, False, False])
    
    # Get chosen representative for each cluster
    chosen = metric_df.groupby('cluster_id').head(1)[['cluster_id', 'eccDNA_id']]
    chosen = chosen.rename(columns={'eccDNA_id': 'chosen_id'})
    
    # Aggregate metadata from all cluster members
    agg_dict = {}
    
    if 'eReads' in df.columns:
        agg_dict['eReads'] = df.groupby('cluster_id')['eReads'].apply(merge_read_lists)
    
    if 'eRepeatNum' in df.columns:
        # Sum eRepeatNum for unique IDs in cluster
        unique_by_id = df[['cluster_id', 'eccDNA_id', 'eRepeatNum']].dropna().drop_duplicates()
        agg_dict['eRepeatNum'] = unique_by_id.groupby('cluster_id')['eRepeatNum'].apply(
            lambda s: to_numeric_safe(s, 0).sum()
        )
    
    if 'eStrand' in df.columns:
        agg_dict['eStrand'] = df.groupby('cluster_id')['eStrand'].apply(majority_vote)
    
    if agg_dict:
        agg_df = pd.DataFrame(agg_dict)
    else:
        agg_df = pd.DataFrame(index=df['cluster_id'].unique())
    
    # Join with chosen representatives
    result = df.merge(chosen, on='cluster_id', how='left')
    result = result[result['eccDNA_id'] == result['chosen_id']].drop_duplicates(subset=['cluster_id', 'eccDNA_id'])
    
    # Merge aggregated data
    result = result.merge(agg_df, on='cluster_id', how='left', suffixes=('', '_agg'))
    
    # Use aggregated values where available
    for col in ['eReads', 'eRepeatNum', 'eStrand']:
        col_agg = col + '_agg'
        if col_agg in result.columns:
            result[col] = result[col_agg].where(result[col_agg].notna(), result.get(col))
            result = result.drop(columns=[col_agg])
    
    # Add merge tracking
    cluster_counts = df.groupby('cluster_id')['eccDNA_id'].nunique()
    result['num_merged'] = result['cluster_id'].map(cluster_counts).fillna(1).astype(int)
    
    merged_ids = df.groupby('cluster_id')['eccDNA_id'].apply(
        lambda s: ';'.join(sorted(s.unique()))
    )
    result['merged_from_ids'] = result['cluster_id'].map(merged_ids).fillna(result['eccDNA_id'])
    
    result['orig_eccDNA_id'] = result['eccDNA_id']
    result['data_type'] = dtype
    
    return result

def process_mecc_cecc(df: pd.DataFrame, clusters: CDHitClusters, dtype: str) -> pd.DataFrame:
    """Process Mecc/Cecc data: keep all rows of representative ID."""
    if df.empty:
        return df
    
    full_df = df.copy()
    id_to_cluster = clusters.get_id_to_cluster_map()
    full_df['cluster_id'] = full_df['eccDNA_id'].map(id_to_cluster)
    
    # Handle singletons
    singleton_mask = full_df['cluster_id'].isna()
    full_df.loc[singleton_mask, 'cluster_id'] = 'singleton:' + full_df.loc[singleton_mask, 'eccDNA_id'].astype(str)
    
    # Determine metric column for selection
    if dtype == 'Mecc':
        metric_col = 'copyNum' if 'copyNum' in full_df.columns else None
        tiebreak_col = 'eLength' if 'eLength' in full_df.columns else None
    else:  # Cecc
        metric_col = 'consLen' if 'consLen' in full_df.columns else 'eLength'
        tiebreak_col = None
    
    # Calculate metrics
    if metric_col and metric_col in full_df.columns:
        metric = full_df.groupby(['cluster_id', 'eccDNA_id'])[metric_col].max()
    else:
        metric = full_df.groupby(['cluster_id', 'eccDNA_id']).size()
    
    metric_df = metric.reset_index().rename(columns={metric.name: 'metric'})
    
    # Add tiebreaker if needed
    if tiebreak_col and tiebreak_col in full_df.columns:
        tie = full_df.groupby(['cluster_id', 'eccDNA_id'])[tiebreak_col].max()
        tie_df = tie.reset_index().rename(columns={tie.name: 'tie'})
        metric_df = metric_df.merge(tie_df, on=['cluster_id', 'eccDNA_id'], how='left')
    else:
        metric_df['tie'] = 0
    
    # Add representative flag
    rep_map = clusters.rep_of
    rep_df = pd.DataFrame({'cluster_id': list(rep_map.keys()), 'rep_id': list(rep_map.values())})
    metric_df = metric_df.merge(rep_df, on='cluster_id', how='left')
    metric_df['is_rep'] = (metric_df['eccDNA_id'] == metric_df['rep_id']).astype(int)
    
    # Sort to select best representative
    metric_df = metric_df.sort_values(['cluster_id', 'is_rep', 'metric', 'tie'], 
                                     ascending=[True, False, False, False])
    
    # Get chosen representative for each cluster
    chosen = metric_df.groupby('cluster_id').head(1)[['cluster_id', 'eccDNA_id']]
    chosen = chosen.rename(columns={'eccDNA_id': 'chosen_id'})
    
    # Filter to keep only rows of chosen representatives
    result = full_df.merge(chosen, on='cluster_id', how='left')
    result = result[result['eccDNA_id'] == result['chosen_id']].copy()
    
    # Aggregate metadata from all cluster members
    agg_dict = {}
    
    if 'readName' in full_df.columns:
        agg_dict['readName'] = full_df.groupby('cluster_id')['readName'].apply(concat_unique_semicolon)
    
    if 'eReads' in full_df.columns:
        agg_dict['eReads'] = full_df.groupby('cluster_id')['eReads'].apply(merge_read_lists)
    
    if 'copyNum' in full_df.columns:
        # Sum copyNum for unique IDs
        unique_by_id = full_df[['cluster_id', 'eccDNA_id', 'copyNum']].dropna().drop_duplicates()
        agg_dict['copyNum'] = unique_by_id.groupby('cluster_id')['copyNum'].apply(
            lambda s: to_numeric_safe(s, 0).sum()
        )
    
    # Handle strand columns
    for strand_col in ['strand', 'eStrand', 'sstrand']:
        if strand_col in full_df.columns:
            agg_dict[strand_col] = full_df.groupby('cluster_id')[strand_col].apply(majority_vote)
    
    if agg_dict:
        agg_df = pd.DataFrame(agg_dict)
        result = result.merge(agg_df, on='cluster_id', how='left', suffixes=('', '_agg'))
        
        # Use aggregated values
        for col in agg_dict.keys():
            col_agg = col + '_agg'
            if col_agg in result.columns:
                result[col] = result[col_agg].where(pd.notna(result[col_agg]), result.get(col))
                result = result.drop(columns=[col_agg])
    
    # Add merge tracking
    cluster_counts = full_df.groupby('cluster_id')['eccDNA_id'].nunique()
    result['num_merged'] = result['cluster_id'].map(cluster_counts).fillna(1).astype(int)
    
    merged_ids = full_df.groupby('cluster_id')['eccDNA_id'].apply(
        lambda s: ';'.join(sorted(s.unique()))
    )
    result['merged_from_ids'] = result['cluster_id'].map(merged_ids).fillna(result['eccDNA_id'])
    
    result['orig_eccDNA_id'] = result['eccDNA_id']
    result['data_type'] = dtype
    
    return result

def renumber_eccdna_ids(df: pd.DataFrame, dtype: str) -> pd.DataFrame:
    """Renumber eccDNA IDs sequentially by type."""
    if df.empty:
        return df
    
    df = df.copy()
    unique_clusters = sorted(df['cluster_id'].astype(str).unique())
    
    # Create mapping: cluster_id -> new eccDNA_id
    id_mapping = {}
    for i, cluster_id in enumerate(unique_clusters):
        new_id = f"{dtype}DNA{i+1}"
        id_mapping[cluster_id] = new_id
    
    df['eccDNA_id'] = df['cluster_id'].astype(str).map(id_mapping)
    
    return df

# ================== Output Generation ==================
def write_uecc_outputs(df: pd.DataFrame, output_dir: Path, prefix: Optional[str], drop_seq: bool = False) -> None:
    """Generate Uecc output files."""
    if df.empty:
        logger.info("Uecc data empty, skipping outputs")
        return
    
    df = normalize_coordinates(df)
    
    # Sort by eccDNA_id using natural sort for consistent output
    # Extract numeric part for sorting
    df['type_prefix'] = df['eccDNA_id'].str.extract(r'([A-Za-z]+)', expand=False)
    df['id_number'] = df['eccDNA_id'].str.extract(r'(\d+)', expand=False).astype(int)
    df = df.sort_values(['type_prefix', 'id_number'])
    df = df.drop(columns=['type_prefix', 'id_number'])
    
    # Derive readName and reads_count
    if 'eReads' in df.columns:
        readName_series = df['eReads'].fillna("")
    else:
        readName_series = pd.Series([""] * len(df), index=df.index)
    
    reads_count_series = readName_series.apply(count_reads_from_string)
    
    # Handle copyNum: use existing or derive from eRepeatNum
    if 'copyNum' in df.columns and df['copyNum'].notna().any():
        copy_num = df['copyNum']
    else:
        copy_num = df.get('eRepeatNum', pd.Series([pd.NA] * len(df), index=df.index))
    
    # Build core dataframe
    core = pd.DataFrame({
        'eccDNA_id': df['eccDNA_id'],
        'eChr': df.get('eChr', pd.NA),
        'eStart0': df['eStart0'],
        'eEnd0': df['eEnd0'],
        'eStrand': df.get('eStrand', pd.NA),
        'eLength': df['eLength'],
        'MatDegree': df.get('MatDegree', pd.NA).apply(format_two_decimals),
        'copyNum': copy_num,
        'eClass': 'Uecc',
        'num_merged': df.get('num_merged', 1),
        'merged_from_ids': df.get('merged_from_ids', df['eccDNA_id']),
        'reads_count': reads_count_series,
        'readName': readName_series
    })
    
    core = reorder_columns_for_output(core)
    
    # Write core CSV (drop eSeq if requested)
    core_df = core.copy()
    if drop_seq and 'eSeq' in core_df.columns:
        core_df = core_df.drop(columns=['eSeq'])
    
    output_file = output_dir / f"{prefix}_UeccDNA.core.csv" if prefix else output_dir / "UeccDNA.core.csv"
    core_df.to_csv(output_file, index=False)
    logger.info(f"Wrote {output_file}")
    
    # Write BED file
    bed_repeat = df['eRepeatNum'] if 'eRepeatNum' in df.columns else copy_num
    bed = pd.DataFrame({
        'chrom': core['eChr'],
        'chromStart': core['eStart0'],
        'chromEnd': core['eEnd0'],
        'name': core['eccDNA_id'],
        'score': core['reads_count'],
        'strand': core['eStrand'],
        'eLength': core['eLength'],
        'eClass': 'Uecc',
        'eRepeatNum': bed_repeat,
        'MatDegree': core['MatDegree']
    })
    
    bed_file = output_dir / f"{prefix}_UeccDNA.bed" if prefix else output_dir / "UeccDNA.bed"
    bed.to_csv(bed_file, sep='\t', header=False, index=False)
    logger.info(f"Wrote {bed_file}")
    
    # Write FASTA if sequence available
    if 'eSeq' in df.columns:
        fa_file = output_dir / f"{prefix}_UeccDNA.fa" if prefix else output_dir / "UeccDNA.fa"
        seq_map = df.set_index('eccDNA_id')['eSeq'].to_dict()
        copy_map = pd.Series(copy_num.values, index=df['eccDNA_id']).to_dict()
        
        with open(fa_file, 'w') as f:
            for _, row in core.iterrows():
                seq = seq_map.get(row['eccDNA_id'])
                repeats = copy_map.get(row['eccDNA_id'])
                repeats_str = str(int(repeats)) if pd.notna(repeats) else 'NA'
                reads_str = str(int(row['reads_count'])) if pd.notna(row['reads_count']) else 'NA'
                
                header = f">{row['eccDNA_id']}|{row['eChr']}:{row['eStart0']}-{row['eEnd0']}({row['eStrand']})|length={row['eLength']}|repeats={repeats_str}|reads={reads_str}"
                f.write(header + "\n")
                
                if isinstance(seq, str) and seq:
                    for i in range(0, len(seq), 60):
                        f.write(seq[i:i+60] + "\n")
                else:
                    f.write("N\n")
        
        logger.info(f"Wrote {fa_file}")
    else:
        logger.warning("No eSeq column in Uecc data, skipping FASTA output")

def write_mecc_outputs(df: pd.DataFrame, output_dir: Path, prefix: Optional[str], drop_seq: bool = False) -> None:
    """Generate Mecc output files."""
    if df.empty:
        logger.info("Mecc data empty, skipping outputs")
        return
    
    df = normalize_coordinates(df).copy()
    
    # Sort by eccDNA_id first (natural sort), then by position for multi-site entries
    # Extract numeric part for sorting
    df['type_prefix'] = df['eccDNA_id'].str.extract(r'([A-Za-z]+)', expand=False)
    df['id_number'] = df['eccDNA_id'].str.extract(r'(\d+)', expand=False).astype(int)
    df = df.sort_values(['type_prefix', 'id_number', 'eStart0'])
    df = df.drop(columns=['type_prefix', 'id_number'])
    
    # Add hit indexing
    df['hit_index'] = df.groupby('eccDNA_id').cumcount() + 1
    hit_counts = df.groupby('eccDNA_id').size().rename('hit_count')
    df = df.merge(hit_counts, on='eccDNA_id', how='left')
    
    # Calculate reads_count per eccDNA
    if 'readName' in df.columns:
        reads_count = df.groupby('eccDNA_id')['readName'].apply(
            lambda s: count_reads_from_string(';'.join(s.dropna().astype(str)))
        ).rename('reads_count')
    else:
        reads_count = df.groupby('eccDNA_id').size().rename('reads_count')
    
    # Build core dataframe (excluding certain columns)
    exclude_cols = ['identity', 'alignment_length', 'evalue', 'bit_score']
    
    core = pd.DataFrame({
        'eccDNA_id': df['eccDNA_id'],
        'eChr': df.get('eChr', pd.NA),
        'eStart0': df['eStart0'],
        'eEnd0': df['eEnd0'],
        'eStrand': df.get('eStrand', pd.NA),
        'eLength': df['eLength'],
        'MatDegree': df.get('MatDegree', pd.NA).apply(format_two_decimals),
        'copyNum': df.get('copyNum', pd.NA),
        'eClass': 'Mecc',
        'hit_index': df['hit_index'],
        'hit_count': df['hit_count'],
        'num_merged': df.get('num_merged', 1),
        'merged_from_ids': df.get('merged_from_ids', df['eccDNA_id'])
    })
    
    core = core.merge(reads_count, on='eccDNA_id', how='left')
    core['readName'] = df.get('readName', pd.NA)
    
    core = reorder_columns_for_output(core)
    
    # Write core CSV (drop eSeq if requested)
    core_df = core.copy()
    if drop_seq and 'eSeq' in df.columns:
        # Note: eSeq might be in original df but not in core
        pass
    
    output_file = output_dir / f"{prefix}_MeccSites.core.csv" if prefix else output_dir / "MeccSites.core.csv"
    core_df.to_csv(output_file, index=False)
    logger.info(f"Wrote {output_file}")
    
    # Write Sites BED (score=1)
    bed = pd.DataFrame({
        'chrom': core['eChr'],
        'chromStart': core['eStart0'],
        'chromEnd': core['eEnd0'],
        'name': core['eccDNA_id'] + "|site" + core['hit_index'].astype(str) + "/" + core['hit_count'].astype(str),
        'score': 1,
        'strand': core['eStrand'],
        'eLength': core['eLength'],
        'eClass': 'Mecc',
        'copyNum': core['copyNum']
    })
    
    bed_file = output_dir / f"{prefix}_MeccSites.bed" if prefix else output_dir / "MeccSites.bed"
    bed.to_csv(bed_file, sep='\t', header=False, index=False)
    logger.info(f"Wrote {bed_file}")
    
    # Write BestSite BED
    if 'bit_score' in df.columns:
        best_idx = df.groupby('eccDNA_id')['bit_score'].idxmax()
        best_df = df.loc[best_idx].copy()
    else:
        best_df = df.groupby('eccDNA_id').head(1).copy()
    
    best_df = best_df.merge(reads_count, on='eccDNA_id', how='left')
    
    best_bed = pd.DataFrame({
        'chrom': best_df['eChr'],
        'chromStart': best_df['eStart0'],
        'chromEnd': best_df['eEnd0'],
        'name': best_df['eccDNA_id'],
        'score': best_df['reads_count'].fillna(1).astype(int),
        'strand': best_df.get('eStrand', pd.NA),
        'eLength': best_df.get('eLength', pd.NA),
        'eClass': 'Mecc',
        'copyNum': best_df.get('copyNum', pd.NA)
    })
    
    best_bed_file = output_dir / f"{prefix}_MeccBestSite.bed" if prefix else output_dir / "MeccBestSite.bed"
    best_bed.to_csv(best_bed_file, sep='\t', header=False, index=False)
    logger.info(f"Wrote {best_bed_file}")
    
    # Write FASTA
    fa_file = output_dir / f"{prefix}_Mecc.fa" if prefix else output_dir / "Mecc.fa"
    
    with open(fa_file, 'w') as f:
        meta = pd.DataFrame({'eccDNA_id': df['eccDNA_id'].unique()})
        meta = meta.merge(df.groupby('eccDNA_id')['eLength'].max().rename('length'), on='eccDNA_id', how='left')
        meta = meta.merge(df.groupby('eccDNA_id')['copyNum'].max().rename('copies'), on='eccDNA_id', how='left')
        meta = meta.merge(reads_count, on='eccDNA_id', how='left')
        meta = meta.merge(df.groupby('eccDNA_id')['hit_count'].max().rename('sites'), on='eccDNA_id', how='left')
        
        if 'eSeq' in df.columns:
            seq_pick = df.dropna(subset=['eSeq']).groupby('eccDNA_id').head(1)[['eccDNA_id', 'eSeq']]
            meta = meta.merge(seq_pick, on='eccDNA_id', how='left')
        
        for _, row in meta.iterrows():
            reads_n = int(row['reads_count']) if pd.notna(row['reads_count']) else 1
            sites_n = int(row['sites']) if pd.notna(row['sites']) else 1
            length_n = int(row['length']) if pd.notna(row['length']) else 'NA'
            copies_n = row['copies'] if pd.notna(row['copies']) else 'NA'
            
            header = f">{row['eccDNA_id']}|multi_loci:{sites_n}_sites|length={length_n}|copies={copies_n}|reads={reads_n}"
            f.write(header + "\n")
            
            seq = row.get('eSeq')
            if isinstance(seq, str) and seq:
                for i in range(0, len(seq), 60):
                    f.write(seq[i:i+60] + "\n")
            else:
                f.write("N\n")
    
    logger.info(f"Wrote {fa_file}")

def write_cecc_outputs(df: pd.DataFrame, output_dir: Path, prefix: Optional[str], drop_seq: bool = False) -> None:
    """Generate Cecc output files."""
    if df.empty:
        logger.info("Cecc data empty, skipping outputs")
        return
    
    df = normalize_coordinates(df).copy()
    
    # Sort by eccDNA_id first using natural sort
    # Extract numeric part for sorting
    df['type_prefix'] = df['eccDNA_id'].str.extract(r'([A-Za-z]+)', expand=False)
    df['id_number'] = df['eccDNA_id'].str.extract(r'(\d+)', expand=False).astype(int)
    df = df.sort_values(['type_prefix', 'id_number'])
    df = df.drop(columns=['type_prefix', 'id_number'])
    
    # Handle segment indexing
    if 'segment_order' in df.columns:
        try:
            df['seg_index'] = pd.to_numeric(df['segment_order'], errors='coerce').astype('Int64')
        except Exception:
            df['seg_index'] = df.groupby('eccDNA_id')['eStart0'].rank(method='first').astype(int)
    else:
        df['seg_index'] = df.groupby('eccDNA_id')['eStart0'].rank(method='first').astype(int)
    
    seg_total = df.groupby('eccDNA_id').size().rename('seg_total')
    df = df.merge(seg_total, on='eccDNA_id', how='left')
    
    # Assign junction roles
    if 'junction_role' not in df.columns:
        def assign_role(row):
            if row['seg_index'] == 1:
                return 'head'
            elif row['seg_index'] == row['seg_total']:
                return 'tail'
            else:
                return 'body'
        df['junction_role'] = df.apply(assign_role, axis=1)
    
    # Calculate reads_count
    if 'readName' in df.columns:
        reads_per_id = df.groupby('eccDNA_id')['readName'].apply(
            lambda s: count_reads_from_string(';'.join(s.dropna().astype(str)))
        ).rename('reads_count')
    else:
        reads_per_id = df.groupby('eccDNA_id').size().rename('reads_count')
    
    df = df.merge(reads_per_id, on='eccDNA_id', how='left')
    
    # Build core dataframe
    core = pd.DataFrame({
        'eccDNA_id': df['eccDNA_id'],
        'eChr': df.get('eChr', pd.NA),
        'eStart0': df['eStart0'],
        'eEnd0': df['eEnd0'],
        'eStrand': df.get('eStrand', pd.NA),
        'seg_index': df['seg_index'],
        'seg_total': df['seg_total'],
        'junction_role': df.get('junction_role', pd.NA),
        'eLength': df['eLength'],
        'MatDegree': df.get('MatDegree', pd.NA).apply(format_two_decimals),
        'copyNum': df.get('copyNum', pd.NA),
        'num_merged': df.get('num_merged', 1),
        'merged_from_ids': df.get('merged_from_ids', df['eccDNA_id']),
        'reads_count': df['reads_count'],
        'readName': df.get('readName', pd.NA)
    })
    
    core = reorder_columns_for_output(core)
    
    # Write core CSV
    core_df = core.copy()
    
    output_file = output_dir / f"{prefix}_CeccSegments.core.csv" if prefix else output_dir / "CeccSegments.core.csv"
    core_df.to_csv(output_file, index=False)
    logger.info(f"Wrote {output_file}")
    
    # Write Segments BED
    name_col = (core['eccDNA_id'] + "|seg" + core['seg_index'].astype(str) + 
                "/" + core['seg_total'].astype(str) + "|" + core['junction_role'].astype(str))
    
    bed = pd.DataFrame({
        'chrom': core['eChr'],
        'chromStart': core['eStart0'],
        'chromEnd': core['eEnd0'],
        'name': name_col,
        'score': core['reads_count'],
        'strand': core['eStrand'],
        'eLength': core['eLength'],
        'eClass': 'Cecc',
        'copyNum': core['copyNum'],
        'MatDegree': core['MatDegree']
    })
    
    bed_file = output_dir / f"{prefix}_CeccSegments.bed" if prefix else output_dir / "CeccSegments.bed"
    bed.to_csv(bed_file, sep='\t', header=False, index=False)
    logger.info(f"Wrote {bed_file}")
    
    # Write Junctions BEDPE
    junction_rows = []
    for eid, sub in core.sort_values(['eccDNA_id', 'seg_index']).groupby('eccDNA_id'):
        rows = list(sub.to_dict('records'))
        for i in range(len(rows) - 1):
            a, b = rows[i], rows[i + 1]
            junction_rows.append({
                'chrom1': a['eChr'],
                'start1': int(a['eEnd0']) - 1 if pd.notna(a['eEnd0']) else pd.NA,
                'end1': a['eEnd0'],
                'chrom2': b['eChr'],
                'start2': b['eStart0'],
                'end2': int(b['eStart0']) + 1 if pd.notna(b['eStart0']) else pd.NA,
                'name': f"{eid}|seg{a['seg_index']}->seg{b['seg_index']}",
                'score': a['reads_count'],
                'strand1': a['eStrand'],
                'strand2': b['eStrand']
            })
    
    if junction_rows:
        bedpe = pd.DataFrame(junction_rows)
        bedpe_file = output_dir / f"{prefix}_CeccJunctions.bedpe" if prefix else output_dir / "CeccJunctions.bedpe"
        bedpe.to_csv(bedpe_file, sep='\t', header=False, index=False)
        logger.info(f"Wrote {bedpe_file}")
    else:
        logger.info("No junctions to output for Cecc")
    
    # Write FASTA
    fa_file = output_dir / f"{prefix}_Cecc.fa" if prefix else output_dir / "Cecc.fa"
    
    with open(fa_file, 'w') as f:
        seg_n = core.groupby('eccDNA_id')['seg_total'].max().rename('N')
        copies = core.groupby('eccDNA_id')['copyNum'].max().rename('copies')
        reads = core.groupby('eccDNA_id')['reads_count'].max().rename('reads')
        length = core.groupby('eccDNA_id')['eLength'].max().rename('length')
        
        junctions = core.sort_values(['eccDNA_id', 'seg_index']).groupby('eccDNA_id')['eChr'].apply(
            lambda s: '-'.join(s.astype(str).tolist())
        ).rename('junctions')
        
        meta = pd.concat([seg_n, copies, reads, length, junctions], axis=1).reset_index()
        
        if 'eSeq' in df.columns:
            seq_pick = df.dropna(subset=['eSeq']).groupby('eccDNA_id').head(1)[['eccDNA_id', 'eSeq']]
            meta = meta.merge(seq_pick, on='eccDNA_id', how='left')
        
        for _, row in meta.iterrows():
            segments_n = int(row['N']) if pd.notna(row['N']) else 1
            length_n = int(row['length']) if pd.notna(row['length']) else 'NA'
            copies_n = row['copies'] if pd.notna(row['copies']) else 'NA'
            reads_n = int(row['reads']) if pd.notna(row['reads']) else 'NA'
            
            header = f">{row['eccDNA_id']}|segments:{segments_n}|junctions:{row['junctions']}|length={length_n}|copies={copies_n}|reads={reads_n}"
            f.write(header + "\n")
            
            seq = row.get('eSeq')
            if isinstance(seq, str) and seq:
                for i in range(0, len(seq), 60):
                    f.write(seq[i:i+60] + "\n")
            else:
                f.write("N\n")
    
    logger.info(f"Wrote {fa_file}")

# ================== Main Processing ==================
def process_all_types(
    uecc_csv: Optional[Path],
    uecc_clstr: Optional[Path],
    mecc_csv: Optional[Path],
    mecc_clstr: Optional[Path],
    cecc_csv: Optional[Path],
    cecc_clstr: Optional[Path],
    output_dir: Path,
    prefix: Optional[str],
    drop_seq: bool
) -> int:
    """Process all eccDNA types independently and generate outputs."""
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    results = {}
    
    # Process Uecc
    if uecc_csv and uecc_clstr:
        logger.info("Processing Uecc data...")
        if not uecc_csv.exists():
            logger.warning(f"Uecc CSV not found: {uecc_csv}")
        elif not uecc_clstr.exists():
            logger.warning(f"Uecc cluster file not found: {uecc_clstr}")
        else:
            uecc_clusters = CDHitClusters()
            uecc_clusters.parse(uecc_clstr)
            
            uecc_df = pd.read_csv(uecc_csv)
            if 'eccDNA_id' not in uecc_df.columns:
                raise ValueError(f"Uecc CSV missing 'eccDNA_id' column: {uecc_csv}")
            
            processed_uecc = process_uecc(uecc_df, uecc_clusters, 'Uecc')
            processed_uecc = renumber_eccdna_ids(processed_uecc, 'Uecc')
            
            # Store original eSeq if present before potentially dropping
            uecc_has_seq = 'eSeq' in processed_uecc.columns
            
            results['Uecc'] = processed_uecc
            write_uecc_outputs(processed_uecc, output_dir, prefix, drop_seq)
    
    # Process Mecc
    if mecc_csv and mecc_clstr:
        logger.info("Processing Mecc data...")
        if not mecc_csv.exists():
            logger.warning(f"Mecc CSV not found: {mecc_csv}")
        elif not mecc_clstr.exists():
            logger.warning(f"Mecc cluster file not found: {mecc_clstr}")
        else:
            mecc_clusters = CDHitClusters()
            mecc_clusters.parse(mecc_clstr)
            
            mecc_df = pd.read_csv(mecc_csv)
            if 'eccDNA_id' not in mecc_df.columns:
                raise ValueError(f"Mecc CSV missing 'eccDNA_id' column: {mecc_csv}")
            
            processed_mecc = process_mecc_cecc(mecc_df, mecc_clusters, 'Mecc')
            processed_mecc = renumber_eccdna_ids(processed_mecc, 'Mecc')
            
            # Store original eSeq if present before potentially dropping
            mecc_has_seq = 'eSeq' in processed_mecc.columns
            
            results['Mecc'] = processed_mecc
            write_mecc_outputs(processed_mecc, output_dir, prefix, drop_seq)
    
    # Process Cecc
    if cecc_csv and cecc_clstr:
        logger.info("Processing Cecc data...")
        if not cecc_csv.exists():
            logger.warning(f"Cecc CSV not found: {cecc_csv}")
        elif not cecc_clstr.exists():
            logger.warning(f"Cecc cluster file not found: {cecc_clstr}")
        else:
            cecc_clusters = CDHitClusters()
            cecc_clusters.parse(cecc_clstr)
            
            cecc_df = pd.read_csv(cecc_csv)
            if 'eccDNA_id' not in cecc_df.columns:
                raise ValueError(f"Cecc CSV missing 'eccDNA_id' column: {cecc_csv}")
            
            processed_cecc = process_mecc_cecc(cecc_df, cecc_clusters, 'Cecc')
            processed_cecc = renumber_eccdna_ids(processed_cecc, 'Cecc')
            
            # Store original eSeq if present before potentially dropping
            cecc_has_seq = 'eSeq' in processed_cecc.columns
            
            results['Cecc'] = processed_cecc
            write_cecc_outputs(processed_cecc, output_dir, prefix, drop_seq)
    
    # Print summary
    print("\n" + "=" * 50)
    print("PROCESSING COMPLETE")
    print("=" * 50)
    
    for dtype in ['Uecc', 'Mecc', 'Cecc']:
        if dtype in results:
            df = results[dtype]
            unique_ids = df['eccDNA_id'].nunique()
            total_rows = len(df)
            merged_clusters = df[df['num_merged'] > 1]['cluster_id'].nunique()
            print(f"{dtype}: {unique_ids} unique eccDNAs, {total_rows} total rows, {merged_clusters} merged clusters")
        else:
            print(f"{dtype}: Not processed")
    
    print("=" * 50)
    
    return 0

# ================== CLI ==================
def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Deduplicate eccDNA by type using respective CD-HIT clusters (no cross-type priority)"
    )
    
    # Input files
    parser.add_argument("-iU", "--input-uecc", type=Path, help="Uecc CSV file")
    parser.add_argument("-cU", "--cluster-uecc", type=Path, help="Uecc cluster file (.clstr)")
    parser.add_argument("-iM", "--input-mecc", type=Path, help="Mecc CSV file")
    parser.add_argument("-cM", "--cluster-mecc", type=Path, help="Mecc cluster file (.clstr)")
    parser.add_argument("-iC", "--input-cecc", type=Path, help="Cecc CSV file")
    parser.add_argument("-cC", "--cluster-cecc", type=Path, help="Cecc cluster file (.clstr)")
    
    # Output options
    parser.add_argument("-o", "--output-dir", type=Path, default=Path("merged_output"),
                       help="Output directory (default: merged_output)")
    parser.add_argument("-p", "--prefix", type=str, default=None,
                       help="Output filename prefix")
    parser.add_argument("--drop-seq", action="store_true",
                       help="Drop eSeq column (FASTA will use N if missing)")
    
    # Logging
    parser.add_argument("-l", "--log-level", default="INFO",
                       choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                       help="Logging level (default: INFO)")
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(args.log_level)
    
    # Validate inputs
    if not any([args.input_uecc, args.input_mecc, args.input_cecc]):
        logger.error("At least one input CSV file must be provided")
        sys.exit(1)
    
    # Check paired inputs
    if args.input_uecc and not args.cluster_uecc:
        logger.error("Uecc CSV provided but no cluster file (-cU)")
        sys.exit(1)
    if args.input_mecc and not args.cluster_mecc:
        logger.error("Mecc CSV provided but no cluster file (-cM)")
        sys.exit(1)
    if args.input_cecc and not args.cluster_cecc:
        logger.error("Cecc CSV provided but no cluster file (-cC)")
        sys.exit(1)
    
    try:
        sys.exit(process_all_types(
            uecc_csv=args.input_uecc,
            uecc_clstr=args.cluster_uecc,
            mecc_csv=args.input_mecc,
            mecc_clstr=args.cluster_mecc,
            cecc_csv=args.input_cecc,
            cecc_clstr=args.cluster_cecc,
            output_dir=args.output_dir,
            prefix=args.prefix,
            drop_seq=args.drop_seq
        ))
    except Exception as e:
        logger.error(f"Processing failed: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
