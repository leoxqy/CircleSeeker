#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Astrologer: Advanced Coverage Analysis and High-Coverage Region Detection.

This module provides sophisticated analysis of sequencing coverage patterns to identify
potential circular DNA regions. It integrates mosdepth for coverage calculation and
implements custom algorithms for high-coverage region detection and scoring.

Key features:
- Efficient coverage calculation using mosdepth integration
- Robust statistical filtering of high-coverage regions
- Multi-threaded processing for large datasets
- Comprehensive scoring system for potential circular regions
- Support for various input formats and parameters

Typical usage:
    astrologer = Astrologer(bam_file, prefix)
    astrologer.process()

Copyright (c) 2024 CircleSeeker Team
"""

import os
import sys
import csv
import argparse
import logging
import subprocess
import shutil
import statistics
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

import pysam

class Astrologer:
    """
    A class that integrates the following functionalities:
    1. Uses mosdepth to calculate coverage and decompress the generated per-base.bed file
    2. Performs robust_highcov_filter to filter and score high coverage regions
       - Additionally calculates Avg_Read_Length / Avg_Split_Read_Length
    3. Performs local alignment on these regions and generates local consensus sequences
    """

    def __init__(
        self,
        input_bam=None,
        prefix=None,
        # Parameters related to mosdepth
        mosdepth_threads=4,
        # Parameters related to paths, threads, and cleanup
        temp_dir="./temp",
        keep_temp=False,
        log_file=None,
        log_level="INFO",
        # Parameters related to robust_highcov_filter
        z_threshold=3.0,
        min_length=500,
        gap_threshold=30,
        mapq_threshold=1,
        no_reads_list=False,
        min_score=0.0,
        output_filter="output_filter.tsv",
        # Parameters related to generate_local_consensus
        input_txt="output_filter.tsv",
        ref_fa="ref.fa",
        reads_fa="reads.fa",
        output_consensus="my_output_with_eSeq.csv",
        threads=4,
        minimap2_threads=4,
    ):
        """
        Modified constructor that directly accepts keyword arguments (with default values).
        """
        # Save keyword arguments to self
        self.input_bam = input_bam
        self.prefix = prefix
        self.mosdepth_threads = mosdepth_threads
        self.temp_dir = temp_dir
        self.keep_temp = keep_temp
        self.log_file = log_file
        self.log_level = log_level

        self.z_threshold = z_threshold
        self.min_length = min_length
        self.gap_threshold = gap_threshold
        self.mapq_threshold = mapq_threshold
        self.no_reads_list = no_reads_list
        self.min_score = min_score
        self.output_filter = output_filter

        self.input_txt = input_txt
        self.ref_fa = ref_fa
        self.reads_fa = reads_fa
        self.output_consensus = output_consensus
        self.threads = threads
        self.minimap2_threads = minimap2_threads

        # Prepare some internal attributes
        self.regions = []
        self.coverages = []
        self.final_merged = []

        # These two are constants used in robust_highcov_filter
        self.E_S = 1.0          # Expected split read ratio
        self.penalty_K = 10.0   # Penalty coefficient for weighted_cov/num_reads

        # Initialize logger
        self.logger = self._init_logger()

    def _init_logger(self):
        logger = logging.getLogger("AstrologerPipeline")
        logger.setLevel(getattr(logging, self.log_level.upper(), logging.INFO))
        
        logger.propagate = False

        # Console output
        console_handler = logging.StreamHandler(sys.stderr)
        console_handler.setLevel(getattr(logging, self.log_level.upper(), logging.INFO))
        console_formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')
        console_handler.setFormatter(console_formatter)
        logger.addHandler(console_handler)

        # File output
        if self.log_file:
            file_handler = logging.FileHandler(self.log_file)
            file_handler.setLevel(getattr(logging, self.log_level.upper(), logging.INFO))
            file_formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')
            file_handler.setFormatter(file_formatter)
            logger.addHandler(file_handler)

        return logger

    def _run_cmd(self, cmd, cwd=None):
        """
        A wrapper function to run system commands, raising an exception and logging errors if the command fails
        """
        cmd_str = " ".join(cmd)
        self.logger.info(f"RUN: {cmd_str}")
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            cwd=cwd,
            text=True
        )
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            self.logger.error(f"Command failed: {cmd_str} | {stderr.strip()}")
            raise RuntimeError(f"Command failed: {cmd_str} -> {stderr.strip()}")
        return stdout, stderr

    # ========== Step 1: mosdepth calculates coverage and decompresses the per-base.bed file ==========

    def run_mosdepth(self):
        """
        1) Calls mosdepth to calculate coverage
        2) Decompresses the prefix.per-base.bed.gz file
        """
        if not self.prefix:
            raise ValueError("You must specify the prefix parameter for mosdepth output.")

        cmd_mosdepth = [
            "mosdepth",
            "-t", str(self.mosdepth_threads),
            self.prefix,
            self.input_bam
        ]
        self._run_cmd(cmd_mosdepth)

        per_base_file = f"{self.prefix}.per-base.bed.gz"
        if not os.path.exists(per_base_file):
            self.logger.warning(f"{per_base_file} does not exist. Check mosdepth run.")
        else:
            cmd_gunzip = ["gunzip", "-f", per_base_file]
            self._run_cmd(cmd_gunzip)

    # ========== Step 2: robust_highcov_filter related logic ==========

    def median_absolute_deviation(self, data, scale_factor=1.4826):
        med = statistics.median(data)
        devs = [abs(x - med) for x in data]
        mad = statistics.median(devs)
        return med, mad * scale_factor

    def load_coverage(self, bed_file):
        regions = []
        coverages = []
        with open(bed_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                chrom, start, end, cov_str = line.split()
                start, end = int(start), int(end)
                cov = float(cov_str)
                if end > start:
                    regions.append((chrom, start, end, cov))
                    coverages.append(cov)
        self.regions = regions
        self.coverages = coverages

    def identify_high_coverage_regions_robust(self):
        median_cov, mad_cov = self.median_absolute_deviation(self.coverages)
        candidate_flags = []
        if mad_cov == 0:
            # If MAD is 0, fall back to the simple rule: cov > median_cov * 2
            for (chrom, start, end, cov) in self.regions:
                length = end - start + 1
                if length >= self.min_length and cov > median_cov * 2:
                    candidate_flags.append(True)
                else:
                    candidate_flags.append(False)
        else:
            for (chrom, start, end, cov) in self.regions:
                length = end - start + 1
                if length < self.min_length:
                    candidate_flags.append(False)
                    continue
                robust_z = (cov - median_cov) / mad_cov
                if robust_z > self.z_threshold:
                    candidate_flags.append(True)
                else:
                    candidate_flags.append(False)
        return candidate_flags

    def filter_zero_cov_low(self, regions, candidate_flags):
        filtered = []
        new_flags = []
        for (r, f) in zip(regions, candidate_flags):
            chrom, start, end, cov = r
            if f:
                filtered.append(r)
                new_flags.append(f)
            else:
                if cov > 0:
                    filtered.append(r)
                    new_flags.append(f)
        return filtered, new_flags

    def merge_by_gap(self, regions, candidate_flags):
        combined = list(zip(regions, candidate_flags))
        combined.sort(key=lambda x: (x[0][0], x[0][1]))

        merged = []
        current_cluster = None

        for (r, f) in combined:
            chrom, start, end, cov = r
            if current_cluster is None:
                current_cluster = {
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'high_cov': f
                }
            else:
                # If the gap between the current cluster and the new region is small, merge them
                if (chrom == current_cluster['chrom'] and
                    abs(start - current_cluster['end']) < self.gap_threshold):
                    if end > current_cluster['end']:
                        current_cluster['end'] = end
                    if f:
                        current_cluster['high_cov'] = True
                else:
                    merged.append(current_cluster)
                    current_cluster = {
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'high_cov': f
                    }

        if current_cluster:
            merged.append(current_cluster)

        return merged

    def fetch_reads_in_region(self, chrom, start, end):
        samfile = pysam.AlignmentFile(self.input_bam, "rb")
        reads_set = set()
        split_reads_set = set()
        mapq_sum = 0
        mapq_count = 0
        mapq_ge30_count = 0

        sum_read_length = 0
        sum_split_read_length = 0

        for read in samfile.fetch(chrom, start, end):
            if read.mapping_quality < self.mapq_threshold:
                continue

            read_len = read.infer_read_length() or 0
            reads_set.add(read.query_name)
            mapq_sum += read.mapping_quality
            mapq_count += 1
            sum_read_length += read_len

            if read.mapping_quality >= 30:
                mapq_ge30_count += 1

            # Check for split reads
            has_split_operation = False
            if read.cigartuples:
                for (op, length) in read.cigartuples:
                    if op in [3, 4, 5]:
                        has_split_operation = True
                        break

            if read.has_tag('SA') or read.is_supplementary:
                has_split_operation = True

            if has_split_operation:
                split_reads_set.add(read.query_name)
                sum_split_read_length += read_len

        samfile.close()
        return (
            reads_set,
            split_reads_set,
            mapq_sum,
            mapq_count,
            mapq_ge30_count,
            sum_read_length,
            sum_split_read_length
        )

    def fetch_reads_for_merged(self, merged):
        for m in merged:
            (
                rset, sset, mq_sum, mq_count, mq30_count,
                sum_len, sum_split_len
            ) = self.fetch_reads_in_region(m['chrom'], m['start'], m['end'])

            m['reads'] = rset
            m['split_reads'] = sset
            m['mapq_sum'] = mq_sum
            m['mapq_count'] = mq_count
            m['mapq_ge30_count'] = mq30_count

            # New fields: sum of read lengths and sum of split read lengths in the region
            m['sum_read_length'] = sum_len
            m['sum_split_read_length'] = sum_split_len

        return merged

    def merge_by_reads(self, merged):
        changed = True
        while changed:
            changed = False
            new_merged = []
            while merged:
                current = merged.pop(0)
                merged_in = False
                for nm in new_merged:
                    if nm['chrom'] == current['chrom']:
                        close_enough = (
                            abs(current['start'] - nm['end']) < self.gap_threshold or
                            abs(nm['start'] - current['end']) < self.gap_threshold
                        )
                        # If the reads in the current cluster and the new region overlap, merge them
                        if close_enough and len(nm['reads'].intersection(current['reads'])) > 0:
                            if current['start'] < nm['start']:
                                nm['start'] = current['start']
                            if current['end'] > nm['end']:
                                nm['end'] = current['end']
                            nm['reads'].update(current['reads'])
                            nm['split_reads'].update(current['split_reads'])
                            nm['mapq_sum'] += current['mapq_sum']
                            nm['mapq_count'] += current['mapq_count']
                            nm['mapq_ge30_count'] += current['mapq_ge30_count']
                            nm['sum_read_length'] += current['sum_read_length']
                            nm['sum_split_read_length'] += current['sum_split_read_length']
                            if current['high_cov']:
                                nm['high_cov'] = True
                            merged_in = True
                            changed = True
                            break
                if not merged_in:
                    new_merged.append(current)
            merged = new_merged
        return merged

    def compute_average_coverage(self, final_merged):
        chrom_dict = defaultdict(list)
        for r in self.regions:
            chrom_dict[r[0]].append(r)
        for c in chrom_dict:
            chrom_dict[c].sort(key=lambda x: x[1])

        for m in final_merged:
            c = m['chrom']
            start = m['start']
            end = m['end']
            sum_cov = 0.0
            count_bases = 0
            if c not in chrom_dict:
                m['weighted_cov'] = 0.0
                continue
            for (ch, s, e, cov) in chrom_dict[c]:
                if e < start:
                    continue
                if s > end:
                    break
                overlap_start = max(s, start)
                overlap_end = min(e, end)
                if overlap_end >= overlap_start:
                    overlap_len = overlap_end - overlap_start + 1
                    sum_cov += cov * overlap_len
                    count_bases += overlap_len

            if count_bases > 0:
                m['weighted_cov'] = sum_cov / count_bases
            else:
                m['weighted_cov'] = 0.0

        return final_merged

    def calculate_score_and_write_output(self, merged_list, out_file):
        """
        In the output, add:
          - Avg_Read_Length
          - Avg_Split_Read_Length
        """
        scale_factor = 100.0 / 120.0

        header = [
            "chrom", "start", "end", "length", "num_reads", "high_cov",
            "weighted_cov", "split_read_count", "split_read_ratio",
            "Mean_MAPQ", "High_MAPQ_ratio", "score",
            "Avg_Read_Length", "Avg_Split_Read_Length"
        ]
        if not self.no_reads_list:
            header.append("reads_list")

        self.logger.info(f"Writing filter result to {out_file}")
        out_f = open(out_file, 'w') if out_file else None
        output_line = "\t".join(header) + "\n"
        if out_f:
            out_f.write(output_line)
        else:
            print(output_line.strip())

        for m in merged_list:
            if not m['high_cov']:
                continue

            length = m['end'] - m['start'] + 1
            num_reads = len(m['reads'])
            split_count = len(m['split_reads'])
            split_ratio = float(split_count) / num_reads if num_reads > 0 else 0.0
            weighted_cov = m['weighted_cov']

            if m['mapq_count'] > 0:
                Mean_MAPQ = m['mapq_sum'] / m['mapq_count']
                High_MAPQ_ratio = m['mapq_ge30_count'] / m['mapq_count']
            else:
                Mean_MAPQ = 0.0
                High_MAPQ_ratio = 0.0

            diff = abs(split_ratio - self.E_S)
            original_score = Mean_MAPQ * (1 + High_MAPQ_ratio) * (1 - diff)

            if num_reads > 0:
                ratio = weighted_cov / num_reads
                penalty_factor = 1.0 / (1.0 + (ratio / self.penalty_K))
            else:
                penalty_factor = 1.0

            final_score = original_score * penalty_factor * scale_factor
            if final_score < self.min_score:
                continue

            sum_read_length = m.get('sum_read_length', 0)
            sum_split_read_length = m.get('sum_split_read_length', 0)
            avg_read_len = float(sum_read_length) / num_reads if num_reads > 0 else 0.0
            avg_split_read_len = float(sum_split_read_length) / split_count if split_count > 0 else 0.0

            row_data = [
                m['chrom'],
                str(m['start']),
                str(m['end']),
                str(length),
                str(num_reads),
                "1" if m['high_cov'] else "0",
                f"{weighted_cov:.2f}",
                str(split_count),
                f"{split_ratio:.3f}",
                f"{Mean_MAPQ:.2f}",
                f"{High_MAPQ_ratio:.3f}",
                f"{final_score:.2f}",
                f"{avg_read_len:.2f}",
                f"{avg_split_read_len:.2f}"
            ]

            if not self.no_reads_list:
                reads_list = ",".join(m['reads'])
                row_data.append(reads_list)

            line = "\t".join(row_data) + "\n"
            if out_f:
                out_f.write(line)
            else:
                print(line.strip())

        if out_f:
            out_f.close()

    def run_robust_highcov_filter(self):
        bed_file = f"{self.prefix}.per-base.bed"
        out_file = self.output_filter

        self.logger.info("Start robust_highcov_filter...")
        self.load_coverage(bed_file)
        candidate_flags = self.identify_high_coverage_regions_robust()
        filtered_regions, filtered_flags = self.filter_zero_cov_low(self.regions, candidate_flags)
        pre_merged = self.merge_by_gap(filtered_regions, filtered_flags)
        pre_merged = self.fetch_reads_for_merged(pre_merged)
        merged = self.merge_by_reads(pre_merged)
        self.final_merged = self.compute_average_coverage(merged)
        self.calculate_score_and_write_output(self.final_merged, out_file)
        self.logger.info("robust_highcov_filter done.")

    # ========== Step 3: generate_local_consensus related logic ==========

    def load_all_reads_fa(self, reads_fa):
        reads_dict = {}
        current_id = None
        seq_chunks = []
        with open(reads_fa, 'r') as f:
            for line in f:
                line=line.strip()
                if line.startswith('>'):
                    if current_id and seq_chunks:
                        reads_dict[current_id] = "".join(seq_chunks)
                    current_id = line[1:]
                    seq_chunks = []
                else:
                    seq_chunks.append(line)
        if current_id and seq_chunks:
            reads_dict[current_id] = "".join(seq_chunks)
        return reads_dict

    def write_subreads_fa(self, reads_ids, reads_dict, out_fa):
        count_found = 0
        with open(out_fa, 'w') as f:
            for rid in reads_ids:
                if rid in reads_dict:
                    seq = reads_dict[rid]
                    f.write(f">{rid}\n{seq}\n")
                    count_found += 1
        return count_found

    def rename_subref_header(self, subref_fa, chrom):
        tmp_path = subref_fa + ".tmp"
        with open(subref_fa, 'r') as fin, open(tmp_path, 'w') as fout:
            for i, line in enumerate(fin):
                if i == 0 and line.startswith('>'):
                    fout.write(f">{chrom}\n")
                else:
                    fout.write(line)
        os.replace(tmp_path, subref_fa)

    def run_minimap2_and_write_sam(self, ref_fa, reads_fa, sam_file, cwd):
        cmd = [
            "minimap2",
            "-t", str(self.minimap2_threads),
            "-x", "map-hifi",
            "-a",
            ref_fa,
            reads_fa
        ]
        self.logger.info(" ".join(cmd))
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=cwd, text=True)
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            self.logger.error(f"Minimap2 failed: {' '.join(cmd)} | {stderr.strip()}")
            raise RuntimeError("Minimap2 failed.")
        with open(os.path.join(cwd, sam_file), 'w') as f:
            f.write(stdout)

    def process_region(self, r, reads_dict, ref, temp_dir):
        chrom = r["chrom"]
        start = int(r["start"])
        end = int(r["end"])
        reads_list_str = r.get("reads_list", "")

        start_1 = start + 1
        end_1 = end
        region_str = f"{chrom}:{start_1}-{end_1}"
        self.logger.info(f"{region_str} start")

        region_dir_name = f"region_{chrom}_{start}_{end}"
        region_dir = os.path.join(temp_dir, region_dir_name)
        if not os.path.exists(region_dir):
            os.makedirs(region_dir, exist_ok=True)

        reads_ids = [x.strip() for x in reads_list_str.split(',') if x.strip()]
        reads_ids = list(set(reads_ids))

        if not reads_ids:
            r["eSeq"] = ""
            if not self.keep_temp:
                shutil.rmtree(region_dir)
            self.logger.info(f"{region_str} end")
            return r

        sub_reads_fa = os.path.join(region_dir, "SubReads.fa")
        count_found = self.write_subreads_fa(reads_ids, reads_dict, sub_reads_fa)
        if count_found == 0:
            r["eSeq"] = ""
            if not self.keep_temp:
                shutil.rmtree(region_dir)
            self.logger.info(f"{region_str} end")
            return r

        sub_ref_fa = os.path.join(region_dir, "SubRef.fa")
        # samtools faidx ref region
        try:
            with open(sub_ref_fa, 'w') as out_f:
                cmd = ["samtools", "faidx", ref, region_str]
                p = subprocess.Popen(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)
                _, stderr = p.communicate()
                if p.returncode != 0:
                    r["eSeq"] = ""
                    if not self.keep_temp:
                        shutil.rmtree(region_dir)
                    self.logger.info(f"{region_str} end")
                    return r
        except Exception:
            r["eSeq"] = ""
            if not self.keep_temp:
                shutil.rmtree(region_dir)
            self.logger.info(f"{region_str} end")
            return r

        # Modify header
        self.rename_subref_header(sub_ref_fa, chrom)

        local_sam = "local.sam"
        try:
            self.run_minimap2_and_write_sam("SubRef.fa", "SubReads.fa", local_sam, region_dir)
        except RuntimeError:
            r["eSeq"] = ""
            if not self.keep_temp:
                shutil.rmtree(region_dir)
            self.logger.info(f"{region_str} end")
            return r

        local_bam = "local.bam"
        local_sorted_bam = "local.sorted.bam"

        try:
            self._run_cmd(["samtools", "view", "-bS", "local.sam", "-o", local_bam], cwd=region_dir)
            self._run_cmd(["samtools", "sort", local_bam, "-o", local_sorted_bam], cwd=region_dir)
            self._run_cmd(["samtools", "index", local_sorted_bam], cwd=region_dir)

            local_mpileup = "local_mpileup.vcf.gz"
            self._run_cmd(["bcftools", "mpileup", "-f", "SubRef.fa", local_sorted_bam, "-Oz", "-o", local_mpileup], cwd=region_dir)
            self._run_cmd(["tabix", "-p", "vcf", local_mpileup], cwd=region_dir)

            local_vcf = "local.vcf"
            self._run_cmd(["bcftools", "call", "-c", "-o", local_vcf, local_mpileup], cwd=region_dir)
            self._run_cmd(["bgzip", local_vcf], cwd=region_dir)
            local_vcf_gz = local_vcf + ".gz"
            self._run_cmd(["tabix", "-p", "vcf", local_vcf_gz], cwd=region_dir)

            local_consensus = "local_consensus.fa"
            self._run_cmd(["bcftools", "consensus", "-f", "SubRef.fa", local_vcf_gz, "-o", local_consensus], cwd=region_dir)

            consensus_seq = ""
            with open(os.path.join(region_dir, local_consensus), 'r') as fa:
                seq_lines = []
                for line in fa:
                    line=line.strip()
                    if line.startswith(">"):
                        continue
                    seq_lines.append(line)
                consensus_seq = "".join(seq_lines)

            r["eSeq"] = consensus_seq
        except RuntimeError:
            r["eSeq"] = ""

        if not self.keep_temp:
            shutil.rmtree(region_dir)

        self.logger.info(f"{region_str} end")
        return r

    def run_generate_local_consensus(self):
        if not os.path.exists(self.temp_dir):
            os.makedirs(self.temp_dir, exist_ok=True)

        self.logger.info("Start generate_local_consensus...")

        # Load all reads
        reads_dict = self.load_all_reads_fa(self.reads_fa)

        # Open self.input_txt to read region information
        with open(self.input_txt, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            fieldnames = reader.fieldnames
            if not all(k in fieldnames for k in ["chrom", "start", "end"]):
                self.logger.error("Input txt missing columns: chrom, start, end.")
                sys.exit(1)

            # If there is no reads_list column, handle it according to specific needs
            if "reads_list" not in fieldnames:
                self.logger.warning("No 'reads_list' column in input_txt. We'll proceed but might get empty eSeq.")
            regions = list(reader)

        if "eSeq" not in fieldnames:
            fieldnames.append("eSeq")

        results = []
        with ProcessPoolExecutor(max_workers=self.threads) as executor:
            futures = []
            for r in regions:
                futures.append(executor.submit(self.process_region, r, reads_dict, self.ref_fa, self.temp_dir))
            for fut in futures:
                res = fut.result()
                results.append(res)

        self.logger.info(f"Writing final consensus result to {self.output_consensus}")
        with open(self.output_consensus, 'w', newline='') as out_f:
            writer = csv.DictWriter(out_f, fieldnames=fieldnames)
            writer.writeheader()
            for rr in results:
                writer.writerow(rr)

        self.logger.info("generate_local_consensus done.")

    # ========== Unified execution entry point ==========

    def run_all_steps(self):
        """
        Execute the following steps in order:
        1) mosdepth -> decompress
        2) robust_highcov_filter
        3) generate_local_consensus
        """
        self.run_mosdepth()
        self.run_robust_highcov_filter()
        self.run_generate_local_consensus()


def main():
    """
    Command-line entry point:
    Usage example:
        python astrologer.py \
            --bam input.bam \
            --prefix my_prefix \
            --mosdepth-threads 4 \
            --temp-dir ./temp \
            --keep-temp \
            --z-threshold 3.0 \
            --min-length 500 \
            --gap-threshold 30 \
            --mapq-threshold 1 \
            --no-reads-list \
            --min-score 0.0 \
            --output-filter output_filter.tsv \
            --input-txt output_filter.tsv \
            --ref ref.fa \
            --reads-fa reads.fa \
            --output-consensus final.csv \
            --threads 4 \
            --minimap2-threads 4
    """
    parser = argparse.ArgumentParser(description="Astrologer pipeline (modified to use keyword args).")

    # -- General parameters
    parser.add_argument("--bam", required=True, help="Input BAM file.")
    parser.add_argument("--prefix", required=True, help="Prefix for mosdepth output.")
    parser.add_argument("--mosdepth-threads", type=int, default=4, help="Threads for mosdepth.")
    parser.add_argument("--temp-dir", default="./temp", help="Temporary directory.")
    parser.add_argument("--keep-temp", action="store_true", help="Keep temporary files.")
    parser.add_argument("--log-file", default=None, help="Log file path.")
    parser.add_argument("--log-level", default="INFO", help="Log level (DEBUG, INFO, WARNING, ERROR).")

    # -- robust_highcov_filter parameters
    parser.add_argument("--z-threshold", type=float, default=3.0, help="Z-threshold for robust outlier detection.")
    parser.add_argument("--min-length", type=int, default=500, help="Minimum length for region.")
    parser.add_argument("--gap-threshold", type=int, default=30, help="Gap threshold for merging adjacent regions.")
    parser.add_argument("--mapq-threshold", type=int, default=1, help="Minimum MAPQ to consider a read.")
    parser.add_argument("--no-reads-list", action="store_true", help="If set, do not output the reads list.")
    parser.add_argument("--min-score", type=float, default=0.0, help="Minimum score threshold.")
    parser.add_argument("--output-filter", default="output_filter.tsv", help="robust_highcov_filter output (TSV).")

    # -- generate_local_consensus parameters
    parser.add_argument("--input-txt", default="output_filter.tsv", help="Input TSV with regions info.")
    parser.add_argument("--ref", default="ref.fa", help="Reference FASTA.")
    parser.add_argument("--reads-fa", default="reads.fa", help="FASTA containing all reads.")
    parser.add_argument("--output-consensus", default="my_output_with_eSeq.csv", help="Output CSV with eSeq column.")
    parser.add_argument("--threads", type=int, default=4, help="Number of parallel processes.")
    parser.add_argument("--minimap2-threads", type=int, default=4, help="Threads for minimap2.")

    args = parser.parse_args()

    # Map command-line arguments to Astrologer's keyword arguments
    astro = Astrologer(
        input_bam=args.bam,
        prefix=args.prefix,
        mosdepth_threads=args.mosdepth_threads,
        temp_dir=args.temp_dir,
        keep_temp=args.keep_temp,
        log_file=args.log_file,
        log_level=args.log_level,

        z_threshold=args.z_threshold,
        min_length=args.min_length,
        gap_threshold=args.gap_threshold,
        mapq_threshold=args.mapq_threshold,
        no_reads_list=args.no_reads_list,
        min_score=args.min_score,
        output_filter=args.output_filter,

        input_txt=args.input_txt,
        ref_fa=args.ref,
        reads_fa=args.reads_fa,
        output_consensus=args.output_consensus,
        threads=args.threads,
        minimap2_threads=args.minimap2_threads,
    )

    astro.run_all_steps()
    print(f"Astrologer completed. Results: {args.output_consensus}")

if __name__ == "__main__":
    main()