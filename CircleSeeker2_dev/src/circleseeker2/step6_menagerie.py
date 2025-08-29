#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
menagerie.py
- 处理和分类eccDNA序列（U/M/C三种类型）
- 简化ID命名：U1, M1, C1等
- 分别导出 U / M / C 的 FASTA（*_pre.fasta）

Author: eccDNA Pipeline Team
Date: 2025
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional
import re
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# ==================== CONFIGURATION ====================

class MenagerieConfig:
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


# ==================== EXCEPTIONS ====================

class ProcessingError(Exception):
    pass


# ==================== MENAGERIE PROCESSOR ====================

class MenagerieProcessor:
    def __init__(self, config: Optional[MenagerieConfig] = None):
        self.config = config or MenagerieConfig()
        self.logger = logging.getLogger(__name__)
        self.fasta_sequences: Dict[str, str] = {}

        # 为 U / M / C 分别维护 FASTA 记录
        self.fasta_records_u: List[SeqRecord] = []
        self.fasta_records_m: List[SeqRecord] = []
        self.fasta_records_c: List[SeqRecord] = []

        # 计数器
        self.uecc_counter = 0
        self.mecc_counter = 0
        self.cecc_counter = 0

    # --------- helpers ---------
    @staticmethod
    def _sanitize_fasta_id(id_str: str) -> str:
        id_str = re.sub(r"\s+", "_", str(id_str))
        id_str = re.sub(r"[^A-Za-z0-9._-]", "_", id_str)
        return id_str

    def load_fasta(self, fasta_file: Path):
        self.logger.info(f"Loading FASTA file: {fasta_file}")
        try:
            # 存成纯字符串，切片更快
            self.fasta_sequences = {
                record.id: str(record.seq)
                for record in SeqIO.parse(str(fasta_file), "fasta")
            }
            self.logger.info(f"Loaded {len(self.fasta_sequences):,} sequences")
        except Exception as e:
            raise ProcessingError(f"Failed to load FASTA file: {e}")

    # --------- eSeq O(n) 计算 ---------
    def _compute_eSeq_column(self, df: pd.DataFrame) -> pd.Series:
        if 'query_id' not in df.columns:
            return pd.Series([''] * len(df), index=df.index)

        len_col = next((c for c in ['eLength', 'consLen'] if c in df.columns), None)
        q_start_col = 'q_start' if 'q_start' in df.columns else None
        if not len_col or not q_start_col:
            return pd.Series([''] * len(df), index=df.index)

        first_rows = df[['query_id', q_start_col, len_col]].dropna().drop_duplicates(
            subset=['query_id'], keep='first'
        )

        cache: Dict[str, str] = {}
        get_seq = self.fasta_sequences.get

        qids = first_rows['query_id'].to_numpy()
        q_starts = pd.to_numeric(first_rows[q_start_col], errors='coerce').fillna(0).astype(np.int64).to_numpy()
        cons_lens = pd.to_numeric(first_rows[len_col], errors='coerce').fillna(0).astype(np.int64).to_numpy()

        for qid, q_start, cons_len in zip(qids, q_starts, cons_lens):
            seq_str = get_seq(qid)
            if not seq_str or q_start <= 0 or cons_len <= 0:
                cache[qid] = ''
                continue
            try:
                cache[qid] = self._extract_ring_sequence(seq_str, q_start, cons_len)
            except Exception as e:
                self.logger.debug(f"extract failed for {qid}: {e}")
                cache[qid] = ''

        eSeq_series = df['query_id'].map(cache).fillna('')
        if self.config.validate_sequences and len_col in df.columns:
            lens = pd.to_numeric(df[len_col], errors='coerce')
            mism = (eSeq_series.str.len() != lens).fillna(False)
            mism_count = int(mism.sum())
            if mism_count:
                self.logger.warning(f"eSeq length mismatch in {mism_count} rows (non-fatal)")
        return eSeq_series

    @staticmethod
    def _extract_ring_sequence(seq_str: str, q_start: int, cons_len: int) -> str:
        seq_len = len(seq_str)
        start_index = q_start - 1
        if q_start > cons_len:
            start_index = (q_start - cons_len) - 1
        end_index = start_index + cons_len
        if end_index <= seq_len:
            return seq_str[start_index:end_index]
        return seq_str[start_index:] + seq_str[:(end_index - seq_len)]

    # --------- 编号 & 记录写入（使用简化的ID） ---------
    def _add_uecc_numbering(self, df: pd.DataFrame) -> pd.DataFrame:
        df['eccDNA_id'] = 'NA'
        if not self.config.add_numbering or 'eSeq' not in df.columns:
            return df

        seq_ok = df['eSeq'].astype(bool).to_numpy()
        n_valid = int(seq_ok.sum())
        if n_valid == 0:
            return df

        # 使用简化的ID格式：U1, U2, U3...
        ids = [f"U{i}" for i in range(self.uecc_counter + 1, self.uecc_counter + n_valid + 1)]
        self.uecc_counter += n_valid
        df.loc[seq_ok, 'eccDNA_id'] = ids

        for ecc_id, seq in zip(df.loc[seq_ok, 'eccDNA_id'].to_numpy(), df.loc[seq_ok, 'eSeq'].to_numpy()):
            self.fasta_records_u.append(SeqRecord(Seq(seq), id=ecc_id, description=""))
        return df

    def _add_mecc_numbering(self, df: pd.DataFrame) -> pd.DataFrame:
        df['eccDNA_id'] = 'NA'
        if not self.config.add_numbering or 'query_id' not in df.columns or 'eSeq' not in df.columns:
            return df

        for qid, group in df.groupby('query_id', sort=False):
            mask = group['eSeq'].astype(bool)
            if not mask.any():
                continue
            rep_idx = mask.idxmax()
            rep_seq = df.at[rep_idx, 'eSeq']

            self.mecc_counter += 1
            # 使用简化的ID格式：M1, M2, M3...
            ecc_id = f"M{self.mecc_counter}"
            df.loc[group.index, 'eccDNA_id'] = ecc_id

            self.fasta_records_m.append(SeqRecord(Seq(rep_seq), id=ecc_id, description=""))
        return df

    def _add_cecc_numbering(self, df: pd.DataFrame) -> pd.DataFrame:
        df['eccDNA_id'] = 'NA'
        if not self.config.add_numbering or 'eSeq' not in df.columns:
            return df

        group_col = 'assembly_id' if 'assembly_id' in df.columns else 'query_id'
        if group_col not in df.columns:
            return df

        for gid, group in df.groupby(group_col, sort=False):
            mask = group['eSeq'].astype(bool)
            if not mask.any():
                continue
            rep_idx = mask.idxmax()
            rep_seq = df.at[rep_idx, 'eSeq']

            self.cecc_counter += 1
            # 使用简化的ID格式：C1, C2, C3...
            ecc_id = f"C{self.cecc_counter}"
            df.loc[group.index, 'eccDNA_id'] = ecc_id

            self.fasta_records_c.append(SeqRecord(Seq(rep_seq), id=ecc_id, description=""))
        return df

    # --------- U/M/C processing ---------
    def process_uecc(self, uecc_files: List[Path], output_dir: Path, prefix: Optional[str] = None) -> Optional[pd.DataFrame]:
        if not uecc_files:
            return None
        self.logger.info(f"Processing {len(uecc_files)} UeccDNA files")

        dfs = []
        for file_path in uecc_files:
            try:
                df = pd.read_csv(file_path)
                dfs.append(df)
                self.logger.info(f"Loaded {len(df)} records from {file_path}")
            except Exception as e:
                self.logger.error(f"Failed to load {file_path}: {e}")
        if not dfs:
            return None

        combined_df = pd.concat(dfs, ignore_index=True)
        self.logger.info(f"Combined {len(combined_df)} total UeccDNA records")

        if 'query_id' in combined_df.columns:
            combined_df['eSeq'] = self._compute_eSeq_column(combined_df)
        combined_df = self._add_uecc_numbering(combined_df)

        if self.config.aggregate_uecc and 'eName' in combined_df.columns:
            combined_df = self._aggregate_uecc_by_ename(combined_df)

        output_file = output_dir / (f"{prefix}_UeccDNA_processed.csv" if prefix else "UeccDNA_processed.csv")
        combined_df.to_csv(output_file, index=False)
        self.logger.info(f"Saved processed UeccDNA to {output_file}")
        return combined_df

    def _aggregate_uecc_by_ename(self, df: pd.DataFrame) -> pd.DataFrame:
        self.logger.info("Aggregating UeccDNA by eName")
        if 'eName' not in df.columns:
            self.logger.warning("eName column not found, skipping aggregation")
            return df

        rows = []
        for ename, group in df.groupby('eName', dropna=False):
            if pd.isna(ename) or ename == '':
                rows.extend(group.to_dict('records'))
                continue
            rep = group.iloc[0].to_dict()

            if 'eReads' in group.columns:
                all_reads = []
                for reads in group['eReads'].dropna():
                    s = str(reads)
                    if s:
                        all_reads.extend(s.split(';'))
                rep['eReads'] = ';'.join(dict.fromkeys(all_reads))

            if 'eRepeatNum' in group.columns and 'eccDNA_id' in group.columns:
                uniq = group[['eccDNA_id', 'eRepeatNum']].drop_duplicates()
                repeat_sum = pd.to_numeric(uniq['eRepeatNum'], errors='coerce').sum()
                rep['eRepeatNum'] = int(repeat_sum) if not pd.isna(repeat_sum) else 0

            if 'eStrand' in group.columns:
                strands = group['eStrand'].dropna()
                if len(strands) > 0:
                    rep['eStrand'] = strands.mode(dropna=True).iloc[0]

            if 'eccDNA_id' in group.columns:
                valid_ids = group['eccDNA_id'][(group['eccDNA_id'] != 'NA') & group['eccDNA_id'].notna()]
                if len(valid_ids) > 0:
                    rep['eccDNA_id'] = valid_ids.iloc[0]

            rep['cluster_size'] = len(group)
            rep['aggregated'] = True
            rows.append(rep)

        out = pd.DataFrame(rows)
        if 'aggregated' not in out.columns:
            out['aggregated'] = False
        self.logger.info(f"Aggregated {len(df)} rows to {len(out)} rows")
        return out

    def process_mecc(self, mecc_files: List[Path], output_dir: Path, prefix: Optional[str] = None) -> Optional[pd.DataFrame]:
        if not mecc_files:
            return None
        self.logger.info(f"Processing {len(mecc_files)} MeccDNA files")

        dfs = []
        for file_path in mecc_files:
            try:
                df = pd.read_csv(file_path)
                dfs.append(df)
                self.logger.info(f"Loaded {len(df)} records from {file_path}")
            except Exception as e:
                self.logger.error(f"Failed to load {file_path}: {e}")
        if not dfs:
            return None

        combined_df = pd.concat(dfs, ignore_index=True)
        self.logger.info(f"Combined {len(combined_df)} total MeccDNA records")

        if self.config.calculate_mat_degree and 'Gap_Percentage' in combined_df.columns:
            gap = pd.to_numeric(combined_df['Gap_Percentage'], errors='coerce')
            combined_df['MatDegree'] = (100 - gap).round(2)
            combined_df.loc[gap.isna(), 'MatDegree'] = np.nan
            self.logger.info(f"Calculated MatDegree for {(combined_df['MatDegree'].notna()).sum()} rows")
        elif 'Gap_Percentage' not in combined_df.columns:
            self.logger.warning("Gap_Percentage column not found, skipping MatDegree calculation")

        if 'query_id' in combined_df.columns:
            combined_df['eSeq'] = self._compute_eSeq_column(combined_df)
            combined_df = self._add_mecc_numbering(combined_df)

        if 'MatDegree' in combined_df.columns and 'Gap_Percentage' in combined_df.columns:
            cols = combined_df.columns.tolist()
            cols.remove('MatDegree')
            gap_idx = cols.index('Gap_Percentage')
            cols.insert(gap_idx + 1, 'MatDegree')
            combined_df = combined_df[cols]

        output_file = output_dir / (f"{prefix}_MeccDNA_processed.csv" if prefix else "MeccDNA_processed.csv")
        combined_df.to_csv(output_file, index=False)
        self.logger.info(f"Saved processed MeccDNA to {output_file}")
        return combined_df

    def process_cecc(self, cecc_files: List[Path], output_dir: Path, prefix: Optional[str] = None) -> Optional[pd.DataFrame]:
        if not cecc_files:
            return None
        self.logger.info(f"Processing {len(cecc_files)} CeccDNA files")

        dfs = []
        for file_path in cecc_files:
            try:
                df = pd.read_csv(file_path)
                dfs.append(df)
                self.logger.info(f"Loaded {len(df)} records from {file_path}")
            except Exception as e:
                self.logger.error(f"Failed to load {file_path}: {e}")
        if not dfs:
            return None

        combined_df = pd.concat(dfs, ignore_index=True)
        self.logger.info(f"Combined {len(combined_df)} total CeccDNA records")

        if 'query_id' in combined_df.columns:
            combined_df['eSeq'] = self._compute_eSeq_column(combined_df)
            combined_df = self._add_cecc_numbering(combined_df)

        output_file = output_dir / (f"{prefix}_CeccDNA_processed.csv" if prefix else "CeccDNA_processed.csv")
        combined_df.to_csv(output_file, index=False)
        self.logger.info(f"Saved processed CeccDNA to {output_file}")
        return combined_df

    # --------- 分类别 FASTA 导出（带 _pre 后缀） ---------
    def write_split_fastas_with_pre(self, output_dir: Path, prefix: Optional[str] = None):
        if not self.config.export_fasta:
            self.logger.info("FASTA export disabled by config")
            return

        def _emit(records: List[SeqRecord], tag: str):
            """tag in {'U','M','C'}"""
            name_map = {'U': 'UeccDNA', 'M': 'MeccDNA', 'C': 'CeccDNA'}
            base = name_map[tag]
            if not records:
                self.logger.info(f"No {base} sequences to write")
                return
            output_dir.mkdir(parents=True, exist_ok=True)
            filename = f"{prefix}_{base}_pre.fasta" if prefix else f"{base}_pre.fasta"
            out_path = output_dir / filename
            SeqIO.write(records, str(out_path), "fasta")
            self.logger.info(f"Wrote {len(records)} sequences to {out_path}")

        _emit(self.fasta_records_u, 'U')
        _emit(self.fasta_records_m, 'M')
        _emit(self.fasta_records_c, 'C')

        # Summary
        total = len(self.fasta_records_u) + len(self.fasta_records_m) + len(self.fasta_records_c)
        print("\n" + "=" * 60)
        print("SUMMARY")
        print("=" * 60)
        print(f"UeccDNA sequences: {self.uecc_counter} (records: {len(self.fasta_records_u)})")
        print(f"MeccDNA sequences: {self.mecc_counter} (records: {len(self.fasta_records_m)})")
        print(f"CeccDNA sequences: {self.cecc_counter} (records: {len(self.fasta_records_c)})")
        print(f"Total sequences: {total}")
        print("=" * 60)

    # --------- Xecc (unclassified) ---------
    def process_xecc(self, all_classified_ids: set, output_dir: Path, prefix: Optional[str] = None):
        if not self.config.process_xecc:
            return
        self.logger.info("Processing unclassified sequences (Xecc)")

        all_fasta_ids = set(self.fasta_sequences.keys())
        unclassified_ids = all_fasta_ids - set(map(str, all_classified_ids))

        if not unclassified_ids:
            self.logger.info("No unclassified sequences found")
            return

        self.logger.info(f"Found {len(unclassified_ids)} unclassified sequences")

        unclassified_records: List[SeqRecord] = []
        xecc_counter = 0

        for qid in sorted(unclassified_ids):
            seq_str = self.fasta_sequences[qid]
            half_length = len(seq_str) // 2
            xecc_counter += 1

            if self.config.xecc_prefix_queryid:
                sanitized_qid = self._sanitize_fasta_id(qid)
                # 简化X类的ID格式
                fasta_id = f"{sanitized_qid}{self.config.xecc_id_separator}X{xecc_counter}"
            else:
                fasta_id = f"X{xecc_counter}"

            unclassified_records.append(SeqRecord(Seq(seq_str[:half_length]), id=fasta_id, description=""))

        xecc_filename = f"{prefix}_XeccDNA.fasta" if prefix else "XeccDNA.fasta"
        xecc_path = output_dir / xecc_filename

        try:
            SeqIO.write(unclassified_records, str(xecc_path), "fasta")
            self.logger.info(f"Wrote {len(unclassified_records)} Xecc sequences to {xecc_path}")
        except Exception as e:
            self.logger.error(f"Failed to write XeccDNA file: {e}")

    # --------- logging ---------
    @staticmethod
    def setup_logging(log_level: str):
        numeric_level = getattr(logging, log_level.upper(), None)
        if not isinstance(numeric_level, int):
            raise ValueError(f'Invalid log level: {log_level}')
        logging.basicConfig(
            level=numeric_level,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )


# ==================== MAIN ====================

def main() -> int:
    parser = argparse.ArgumentParser(
        description='Menagerie: Process and classify eccDNA sequences (U/M/C types)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Input arguments
    parser.add_argument('-f', '--fasta', required=True, help='Input FASTA file with doubled circular DNA sequences')

    # Separate input arguments for U/M/C
    parser.add_argument('-u', '--uecc', nargs='+', default=[], help='UeccDNA (Unique-locus) classification CSV files')
    parser.add_argument('-m', '--mecc', nargs='+', default=[], help='MeccDNA (Multi-locus) classification CSV files')
    parser.add_argument('-c', '--cecc', nargs='+', default=[], help='CeccDNA (Chimeric) classification CSV files')

    # Output arguments
    parser.add_argument('-o', '--output-dir', default='menagerie_output', help='Output directory')
    parser.add_argument('-p', '--prefix', default=None, help='Prefix for output files')

    # Options
    parser.add_argument('-X', '--process-xecc', action='store_true', help='Generate XeccDNA.fasta for unclassified sequences')
    parser.add_argument('--no-aggregate', action='store_true', help='Skip UeccDNA aggregation by eName')
    parser.add_argument('--no-mat-degree', action='store_true', help='Skip MatDegree calculation for MeccDNA')
    parser.add_argument('--no-validate', action='store_true', help='Skip sequence validation')
    parser.add_argument('--no-numbering', action='store_true', help='Skip eccDNA numbering & omit FASTA records')
    parser.add_argument('--no-fasta', action='store_true', help='Skip FASTA export')
    parser.add_argument('-l', '--log-level', default='INFO', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'], help='Logging level')

    args = parser.parse_args()

    MenagerieProcessor.setup_logging(args.log_level)
    logger = logging.getLogger(__name__)

    fasta_path = Path(args.fasta)
    if not fasta_path.exists():
        logger.error(f"FASTA file not found: {fasta_path}")
        return 1

    if not (args.uecc or args.mecc or args.cecc):
        logger.error("At least one type of classification file must be provided (-u, -m, or -c)")
        return 1

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    config = MenagerieConfig()
    config.process_xecc = args.process_xecc
    config.validate_sequences = not args.no_validate
    config.add_numbering = not args.no_numbering
    config.export_fasta = not args.no_fasta
    config.aggregate_uecc = not args.no_aggregate
    config.calculate_mat_degree = not args.no_mat_degree

    processor = MenagerieProcessor(config)

    try:
        processor.load_fasta(fasta_path)

        all_classified_ids = set()

        if args.uecc:
            uecc_paths = [Path(p) for p in args.uecc]
            uecc_df = processor.process_uecc(uecc_paths, output_dir, args.prefix)
            if uecc_df is not None and 'query_id' in uecc_df.columns:
                all_classified_ids.update(map(str, uecc_df['query_id'].dropna().unique()))

        if args.mecc:
            mecc_paths = [Path(p) for p in args.mecc]
            mecc_df = processor.process_mecc(mecc_paths, output_dir, args.prefix)
            if mecc_df is not None and 'query_id' in mecc_df.columns:
                all_classified_ids.update(map(str, mecc_df['query_id'].dropna().unique()))

        if args.cecc:
            cecc_paths = [Path(p) for p in args.cecc]
            cecc_df = processor.process_cecc(cecc_paths, output_dir, args.prefix)
            if cecc_df is not None and 'query_id' in cecc_df.columns:
                all_classified_ids.update(map(str, cecc_df['query_id'].dropna().unique()))

        # 分类别导出（U/M/C），文件名带 _pre
        processor.write_split_fastas_with_pre(output_dir, args.prefix)

        # Xecc 可选
        processor.process_xecc(all_classified_ids, output_dir, args.prefix)

        logger.info("Menagerie processing completed successfully")
        return 0

    except Exception as e:
        logger.error(f"Processing failed: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())
