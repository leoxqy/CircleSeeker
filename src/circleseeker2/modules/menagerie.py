"""
Menagerie - eccDNA Sequence Processing and Classification

This module processes and classifies eccDNA sequences into U/M/C types:
- UeccDNA: Unique/simple circular DNA sequences
- MeccDNA: Multiple-copy repeat circular DNA sequences  
- CeccDNA: Complex circular DNA sequences
- XeccDNA: Unclassified sequences (optional)

Key functions:
- Loads FASTA sequences and extracts circular segments
- Assigns simplified IDs (U1, M1, C1, etc.)
- Exports separate FASTA files for each eccDNA type
- Calculates maturation degrees and validates sequences

Migrated from step6_menagerie.py to the new CircleSeeker2 architecture.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Set
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from circleseeker2.exceptions import PipelineError


class MenagerieConfig:
    """Configuration for Menagerie processor."""
    
    def __init__(self):
        self.process_xecc = False
        self.validate_sequences = True
        self.add_numbering = True
        self.export_fasta = True
        self.aggregate_uecc = True
        self.calculate_mat_degree = True
        
        # XeccDNA ID formatting
        self.xecc_prefix_queryid = True
        self.xecc_id_separator = "__"


class Menagerie:
    """eccDNA sequence processor and classifier."""
    
    def __init__(self, config: Optional[MenagerieConfig] = None, 
                 logger: Optional[logging.Logger] = None):
        """Initialize Menagerie processor."""
        self.config = config or MenagerieConfig()
        self.logger = logger or logging.getLogger(self.__class__.__name__)
        self.fasta_sequences: Dict[str, str] = {}
        
        # FASTA records for each type
        self.fasta_records_u: List[SeqRecord] = []
        self.fasta_records_m: List[SeqRecord] = []
        self.fasta_records_c: List[SeqRecord] = []
        
        # Counters
        self.uecc_counter = 0
        self.mecc_counter = 0
        self.cecc_counter = 0
    
    @staticmethod
    def _sanitize_fasta_id(id_str: str) -> str:
        """Sanitize FASTA ID to be valid."""
        id_str = re.sub(r"\\s+", "_", str(id_str))
        id_str = re.sub(r"[^A-Za-z0-9._-]", "_", id_str)
        return id_str
    
    def load_fasta(self, fasta_file: Path) -> None:
        """Load FASTA sequences into memory."""
        self.logger.info(f"Loading FASTA file: {fasta_file}")
        try:
            self.fasta_sequences = {
                record.id: str(record.seq)
                for record in SeqIO.parse(str(fasta_file), "fasta")
            }
            self.logger.info(f"Loaded {len(self.fasta_sequences):,} sequences")
        except Exception as e:
            raise PipelineError(f"Failed to load FASTA file: {e}")
    
    def _compute_eSeq_column(self, df: pd.DataFrame) -> pd.Series:
        """Compute eSeq column efficiently using vectorized operations."""
        if 'query_id' not in df.columns:
            return pd.Series([''] * len(df), index=df.index)
        
        # Find length and start columns
        len_col = next((c for c in ['eLength', 'consLen'] if c in df.columns), None)
        q_start_col = 'q_start' if 'q_start' in df.columns else None
        
        if not len_col or not q_start_col:
            return pd.Series([''] * len(df), index=df.index)
        
        # Get unique query data
        first_rows = df[['query_id', q_start_col, len_col]].dropna().drop_duplicates(
            subset=['query_id'], keep='first'
        )
        
        # Build sequence cache
        cache: Dict[str, str] = {}
        
        qids = first_rows['query_id'].to_numpy()
        q_starts = pd.to_numeric(first_rows[q_start_col], errors='coerce').fillna(0).astype(np.int64).to_numpy()
        cons_lens = pd.to_numeric(first_rows[len_col], errors='coerce').fillna(0).astype(np.int64).to_numpy()
        
        # Extract sequences
        for qid, q_start, cons_len in zip(qids, q_starts, cons_lens):
            # Try direct lookup first
            seq_str = self.fasta_sequences.get(qid)
            
            # If not found, try alternative lookups
            if not seq_str:
                # Try looking up by base read name (before |rep...)
                base_id = qid.split('|')[0] if '|' in qid else qid
                seq_str = self.fasta_sequences.get(base_id)
                
                # If still not found, try fuzzy match by base ID
                if not seq_str:
                    # Look for any key that starts with the base ID
                    for fasta_id in self.fasta_sequences:
                        if fasta_id.startswith(base_id + '|'):
                            seq_str = self.fasta_sequences[fasta_id]
                            self.logger.debug(f"Found fuzzy match: {qid} -> {fasta_id}")
                            break
            
            if not seq_str or q_start <= 0 or cons_len <= 0:
                cache[qid] = ''
                if not seq_str:
                    self.logger.debug(f"No sequence found for {qid}")
                continue
            try:
                cache[qid] = self._extract_ring_sequence(seq_str, q_start, cons_len)
            except Exception as e:
                self.logger.debug(f"Extract failed for {qid}: {e}")
                cache[qid] = ''
        
        # Map back to original dataframe
        eSeq_series = df['query_id'].map(cache).fillna('')
        
        # Validation
        if self.config.validate_sequences and len_col in df.columns:
            lens = pd.to_numeric(df[len_col], errors='coerce')
            mism = (eSeq_series.str.len() != lens).fillna(False)
            mism_count = int(mism.sum())
            if mism_count:
                self.logger.warning(f"eSeq length mismatch in {mism_count} rows (non-fatal)")
        
        return eSeq_series
    
    def _aggregate_uecc_by_ename(self, df: pd.DataFrame) -> pd.DataFrame:
        """Aggregate UeccDNA sequences by eName (matches original step6_menagerie.py)."""
        if 'eName' not in df.columns:
            self.logger.warning("No 'eName' column found, skipping aggregation")
            return df
        
        # Handle NaN/empty eName values - keep them as-is
        valid_ename_mask = df['eName'].notna() & (df['eName'] != '') & (df['eName'] != 'NA')
        
        if not valid_ename_mask.any():
            # No valid eNames to aggregate
            df = df.copy()
            df['cluster_size'] = 1
            df['aggregated'] = False
            return df
        
        # Process aggregatable and non-aggregatable records separately
        aggregatable = df[valid_ename_mask].copy()
        non_aggregatable = df[~valid_ename_mask].copy()
        
        # Group by eName and aggregate
        processed_groups = []
        
        for ename, group in aggregatable.groupby('eName', dropna=False):
            if len(group) == 1:
                # Single record, no aggregation needed
                rep = group.iloc[0].copy()
                rep['cluster_size'] = 1
                rep['aggregated'] = False
            else:
                # Multiple records, aggregate them
                rep = group.iloc[0].copy()  # Use first as representative
                
                # Aggregate eReads (semicolon-separated deduplication)
                if 'eReads' in group.columns:
                    all_reads = []
                    for reads_str in group['eReads'].dropna():
                        if pd.notna(reads_str) and str(reads_str) != 'NA':
                            reads = str(reads_str).split(';')
                            all_reads.extend(r.strip() for r in reads if r.strip())
                    unique_reads = list(dict.fromkeys(all_reads))  # Preserve order while deduplicating
                    rep['eReads'] = ';'.join(unique_reads) if unique_reads else 'NA'
                
                # Aggregate eRepeatNum (sum of unique eccDNA_id combinations)
                if 'eRepeatNum' in group.columns and 'eccDNA_id' in group.columns:
                    unique_combinations = group[['eccDNA_id', 'eRepeatNum']].drop_duplicates()
                    total_repeat_num = 0
                    for _, row in unique_combinations.iterrows():
                        try:
                            repeat_val = pd.to_numeric(row['eRepeatNum'], errors='coerce')
                            if pd.notna(repeat_val):
                                total_repeat_num += repeat_val
                        except (ValueError, TypeError):
                            continue
                    rep['eRepeatNum'] = total_repeat_num
                
                # eStrand consensus (mode)
                if 'eStrand' in group.columns:
                    strand_counts = group['eStrand'].value_counts()
                    rep['eStrand'] = strand_counts.index[0] if not strand_counts.empty else rep['eStrand']
                
                # eccDNA_id selection (first valid one)
                if 'eccDNA_id' in group.columns:
                    valid_ids = group['eccDNA_id'].dropna()
                    valid_ids = valid_ids[valid_ids != 'NA']
                    if not valid_ids.empty:
                        rep['eccDNA_id'] = valid_ids.iloc[0]
                
                # Set aggregation metadata
                rep['cluster_size'] = len(group)
                rep['aggregated'] = True
            
            processed_groups.append(rep)
        
        # Combine processed groups
        if processed_groups:
            aggregated_df = pd.DataFrame(processed_groups)
        else:
            aggregated_df = pd.DataFrame()
        
        # Add non-aggregatable records
        if not non_aggregatable.empty:
            non_aggregatable = non_aggregatable.copy()
            non_aggregatable['cluster_size'] = 1
            non_aggregatable['aggregated'] = False
            
            if not aggregated_df.empty:
                # Ensure column order consistency
                common_cols = list(set(aggregated_df.columns) & set(non_aggregatable.columns))
                aggregated_df = aggregated_df[common_cols]
                non_aggregatable = non_aggregatable[common_cols]
                
                result_df = pd.concat([aggregated_df, non_aggregatable], ignore_index=True)
            else:
                result_df = non_aggregatable
        else:
            result_df = aggregated_df
        
        # Ensure all records have required fields
        if 'cluster_size' not in result_df.columns:
            result_df['cluster_size'] = 1
        if 'aggregated' not in result_df.columns:
            result_df['aggregated'] = False
        
        # Log aggregation summary
        original_count = len(df)
        final_count = len(result_df)
        aggregated_count = len(result_df[result_df.get('aggregated', False)])
        
        self.logger.info(f"UeccDNA aggregation: {original_count} → {final_count} records")
        self.logger.info(f"Aggregated groups: {aggregated_count}, Individual records: {final_count - aggregated_count}")
        
        return result_df
    
    @staticmethod
    def _extract_ring_sequence(seq_str: str, q_start: int, cons_len: int) -> str:
        """Extract circular sequence with proper wrapping."""
        seq_len = len(seq_str)
        start_index = q_start - 1
        
        # Adjust start if needed
        if q_start > cons_len:
            start_index = (q_start - cons_len) - 1
        
        end_index = start_index + cons_len
        
        # Handle wrapping
        if end_index <= seq_len:
            return seq_str[start_index:end_index]
        else:
            return seq_str[start_index:] + seq_str[:(end_index - seq_len)]
    
    def _add_uecc_numbering(self, df: pd.DataFrame) -> pd.DataFrame:
        """Add UeccDNA numbering and create FASTA records."""
        df = df.copy()
        df['eccDNA_id'] = 'NA'
        
        if not self.config.add_numbering or 'eSeq' not in df.columns:
            return df
        
        valid_mask = (df['eSeq'].str.len() > 0)
        valid_df = df[valid_mask].copy()
        
        if valid_df.empty:
            return df
        
        # Create IDs and FASTA records
        for idx in valid_df.index:
            self.uecc_counter += 1
            eccDNA_id = f"U{self.uecc_counter}"
            df.loc[idx, 'eccDNA_id'] = eccDNA_id
            
            if self.config.export_fasta:
                seq = valid_df.loc[idx, 'eSeq']
                query_id = valid_df.loc[idx, 'query_id']
                
                record = SeqRecord(
                    Seq(seq),
                    id=eccDNA_id,
                    description=f"UeccDNA_{query_id}"
                )
                self.fasta_records_u.append(record)
        
        return df
    
    def _add_mecc_numbering(self, df: pd.DataFrame) -> pd.DataFrame:
        """Add MeccDNA numbering and create FASTA records."""
        df = df.copy()
        df['eccDNA_id'] = 'NA'
        
        if not self.config.add_numbering or 'eSeq' not in df.columns:
            return df
        
        # Group by query_id to handle multiple alignments
        for query_id, group in df.groupby('query_id'):
            if group['eSeq'].iloc[0]:  # Check first row for sequence
                self.mecc_counter += 1
                eccDNA_id = f"M{self.mecc_counter}"
                df.loc[group.index, 'eccDNA_id'] = eccDNA_id
                
                if self.config.export_fasta:
                    seq = group['eSeq'].iloc[0]  # Use first sequence
                    record = SeqRecord(
                        Seq(seq),
                        id=eccDNA_id,
                        description=f"MeccDNA_{query_id}"
                    )
                    self.fasta_records_m.append(record)
        
        return df
    
    def _add_cecc_numbering(self, df: pd.DataFrame) -> pd.DataFrame:
        """Add CeccDNA numbering and create FASTA records."""
        df = df.copy()
        df['eccDNA_id'] = 'NA'
        
        if not self.config.add_numbering or 'eSeq' not in df.columns:
            return df
        
        # Group by query_id to handle multiple segments
        for query_id, group in df.groupby('query_id'):
            if group['eSeq'].iloc[0]:  # Check first row for sequence
                self.cecc_counter += 1
                eccDNA_id = f"C{self.cecc_counter}"
                df.loc[group.index, 'eccDNA_id'] = eccDNA_id
                
                if self.config.export_fasta:
                    seq = group['eSeq'].iloc[0]  # Use first sequence
                    record = SeqRecord(
                        Seq(seq),
                        id=eccDNA_id,
                        description=f"CeccDNA_{query_id}"
                    )
                    self.fasta_records_c.append(record)
        
        return df
    
    def process_uecc(self, uecc_files: List[Path], output_dir: Path, 
                     prefix: Optional[str] = None) -> Optional[pd.DataFrame]:
        """Process UeccDNA files."""
        if not uecc_files:
            return None
        
        self.logger.info(f"Processing {len(uecc_files)} UeccDNA files")
        
        dfs = []
        for file_path in uecc_files:
            try:
                df = pd.read_csv(file_path)
                if df.empty:
                    continue
                dfs.append(df)
            except Exception as e:
                self.logger.error(f"Error reading {file_path}: {e}")
        
        if not dfs:
            self.logger.warning("No valid UeccDNA data found")
            return None
        
        # Combine and process
        combined_df = pd.concat(dfs, ignore_index=True)
        combined_df['eSeq'] = self._compute_eSeq_column(combined_df)
        combined_df = self._add_uecc_numbering(combined_df)
        
        # Aggregate by eName if enabled (matches original step6_menagerie.py)
        if self.config.aggregate_uecc:
            combined_df = self._aggregate_uecc_by_ename(combined_df)
        
        # Save results
        output_file = output_dir / f"{prefix}_UeccDNA_processed.csv"
        output_file.parent.mkdir(parents=True, exist_ok=True)
        combined_df.to_csv(output_file, index=False)
        
        # Save FASTA
        if self.config.export_fasta and self.fasta_records_u:
            fasta_file = output_dir / f"{prefix}_UeccDNA_pre.fasta"
            with open(fasta_file, 'w') as handle:
                SeqIO.write(self.fasta_records_u, handle, "fasta")
            self.logger.info(f"UeccDNA FASTA saved: {fasta_file}")
        
        self.logger.info(f"UeccDNA processing complete: {len(combined_df)} records, {len(self.fasta_records_u)} sequences")
        return combined_df
    
    def process_mecc(self, mecc_files: List[Path], output_dir: Path,
                     prefix: Optional[str] = None) -> Optional[pd.DataFrame]:
        """Process MeccDNA files."""
        if not mecc_files:
            return None
        
        self.logger.info(f"Processing {len(mecc_files)} MeccDNA files")
        
        dfs = []
        for file_path in mecc_files:
            try:
                df = pd.read_csv(file_path)
                if df.empty:
                    continue
                dfs.append(df)
            except Exception as e:
                self.logger.error(f"Error reading {file_path}: {e}")
        
        if not dfs:
            self.logger.warning("No valid MeccDNA data found")
            return None
        
        # Combine and process
        combined_df = pd.concat(dfs, ignore_index=True)
        
        # Calculate MatDegree (matches original step6_menagerie.py)
        if self.config.calculate_mat_degree and 'Gap_Percentage' in combined_df.columns:
            gap = pd.to_numeric(combined_df['Gap_Percentage'], errors='coerce')
            combined_df['MatDegree'] = (100 - gap).round(2)
            combined_df.loc[gap.isna(), 'MatDegree'] = np.nan
            self.logger.info(f"Calculated MatDegree for {(combined_df['MatDegree'].notna()).sum()} rows")
        elif 'Gap_Percentage' not in combined_df.columns:
            self.logger.warning("Gap_Percentage column not found, skipping MatDegree calculation")
        
        combined_df['eSeq'] = self._compute_eSeq_column(combined_df)
        combined_df = self._add_mecc_numbering(combined_df)
        
        # Reorder columns to place MatDegree after Gap_Percentage (matches original behavior)
        if 'MatDegree' in combined_df.columns and 'Gap_Percentage' in combined_df.columns:
            cols = combined_df.columns.tolist()
            cols.remove('MatDegree')
            gap_idx = cols.index('Gap_Percentage')
            cols.insert(gap_idx + 1, 'MatDegree')
            combined_df = combined_df[cols]
        
        # Save results
        output_file = output_dir / f"{prefix}_MeccDNA_processed.csv"
        output_file.parent.mkdir(parents=True, exist_ok=True)
        combined_df.to_csv(output_file, index=False)
        
        # Save FASTA
        if self.config.export_fasta and self.fasta_records_m:
            fasta_file = output_dir / f"{prefix}_MeccDNA_pre.fasta"
            with open(fasta_file, 'w') as handle:
                SeqIO.write(self.fasta_records_m, handle, "fasta")
            self.logger.info(f"MeccDNA FASTA saved: {fasta_file}")
        
        self.logger.info(f"MeccDNA processing complete: {len(combined_df)} records, {len(self.fasta_records_m)} sequences")
        return combined_df
    
    def process_cecc(self, cecc_files: List[Path], output_dir: Path,
                     prefix: Optional[str] = None) -> Optional[pd.DataFrame]:
        """Process CeccDNA files."""
        if not cecc_files:
            return None
        
        self.logger.info(f"Processing {len(cecc_files)} CeccDNA files")
        
        dfs = []
        for file_path in cecc_files:
            try:
                df = pd.read_csv(file_path)
                if df.empty:
                    continue
                dfs.append(df)
            except Exception as e:
                self.logger.error(f"Error reading {file_path}: {e}")
        
        if not dfs:
            self.logger.warning("No valid CeccDNA data found")
            return None
        
        # Combine and process
        combined_df = pd.concat(dfs, ignore_index=True)
        combined_df['eSeq'] = self._compute_eSeq_column(combined_df)
        combined_df = self._add_cecc_numbering(combined_df)
        
        # Save results
        output_file = output_dir / f"{prefix}_CeccDNA_processed.csv"
        output_file.parent.mkdir(parents=True, exist_ok=True)
        combined_df.to_csv(output_file, index=False)
        
        # Save FASTA
        if self.config.export_fasta and self.fasta_records_c:
            fasta_file = output_dir / f"{prefix}_CeccDNA_pre.fasta"
            with open(fasta_file, 'w') as handle:
                SeqIO.write(self.fasta_records_c, handle, "fasta")
            self.logger.info(f"CeccDNA FASTA saved: {fasta_file}")
        
        self.logger.info(f"CeccDNA processing complete: {len(combined_df)} records, {len(self.fasta_records_c)} sequences")
        return combined_df
    
    def process_xecc(self, all_classified_ids: Set[str], output_dir: Path,
                     prefix: Optional[str] = None) -> None:
        """Process unclassified sequences (XeccDNA)."""
        if not self.config.process_xecc:
            return
        
        self.logger.info("Processing unclassified sequences (XeccDNA)")
        
        all_fasta_ids = set(self.fasta_sequences.keys())
        xecc_ids = all_fasta_ids - all_classified_ids
        
        if not xecc_ids:
            self.logger.info("No unclassified sequences found")
            return
        
        self.logger.info(f"Found {len(xecc_ids)} unclassified sequences")
        
        # Create XeccDNA FASTA
        xecc_records = []
        for i, seq_id in enumerate(sorted(xecc_ids), 1):
            seq_str = self.fasta_sequences[seq_id]
            # Extract only the first half of the sequence (matches original step6_menagerie.py)
            half_length = len(seq_str) // 2
            
            if self.config.xecc_prefix_queryid:
                new_id = f"{seq_id}{self.config.xecc_id_separator}X{i}"
            else:
                new_id = f"X{i}"
            
            record = SeqRecord(
                Seq(seq_str[:half_length]),  # Only first half
                id=self._sanitize_fasta_id(new_id),
                description=f"XeccDNA_{seq_id}"
            )
            xecc_records.append(record)
        
        # Save XeccDNA FASTA
        if xecc_records:
            xecc_file = output_dir / f"{prefix}_XeccDNA.fasta"
            with open(xecc_file, 'w') as handle:
                SeqIO.write(xecc_records, handle, "fasta")
            self.logger.info(f"XeccDNA FASTA saved: {xecc_file} ({len(xecc_records)} sequences)")
    
    def run_pipeline(
        self,
        fasta_file: Path,
        uecc_files: List[Path],
        mecc_files: List[Path],
        cecc_files: List[Path],
        output_dir: Path,
        prefix: str = "sample",
        process_xecc: bool = False
    ) -> Dict[str, pd.DataFrame]:
        """Run the complete Menagerie processing pipeline."""
        self.logger.info("=" * 60)
        self.logger.info("Menagerie - eccDNA Sequence Processing")
        self.logger.info("=" * 60)
        
        # Set configuration
        self.config.process_xecc = process_xecc
        
        # Load FASTA sequences
        self.load_fasta(fasta_file)
        
        # Process each eccDNA type
        results = {}
        all_classified_ids = set()
        
        # Process UeccDNA
        uecc_df = self.process_uecc(uecc_files, output_dir, prefix)
        if uecc_df is not None:
            results['uecc'] = uecc_df
            all_classified_ids.update(uecc_df['query_id'].unique())
        
        # Process MeccDNA
        mecc_df = self.process_mecc(mecc_files, output_dir, prefix)
        if mecc_df is not None:
            results['mecc'] = mecc_df
            all_classified_ids.update(mecc_df['query_id'].unique())
        
        # Process CeccDNA
        cecc_df = self.process_cecc(cecc_files, output_dir, prefix)
        if cecc_df is not None:
            results['cecc'] = cecc_df
            all_classified_ids.update(cecc_df['query_id'].unique())
        
        # Process XeccDNA (unclassified)
        if process_xecc:
            self.process_xecc(all_classified_ids, output_dir, prefix)
        
        # Final summary
        self.logger.info("=" * 60)
        self.logger.info("Menagerie Processing Summary:")
        self.logger.info(f"  UeccDNA: {len(self.fasta_records_u)} sequences")
        self.logger.info(f"  MeccDNA: {len(self.fasta_records_m)} sequences")
        self.logger.info(f"  CeccDNA: {len(self.fasta_records_c)} sequences")
        if process_xecc:
            unclassified_count = len(self.fasta_sequences) - len(all_classified_ids)
            self.logger.info(f"  XeccDNA: {unclassified_count} sequences")
        self.logger.info("=" * 60)
        
        return results