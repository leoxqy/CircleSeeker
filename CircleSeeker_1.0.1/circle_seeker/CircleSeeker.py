#!/usr/bin/env python3
# coding: utf-8

"""
CircleSeeker: A comprehensive pipeline for detecting and analyzing eccDNA from HiFi sequencing data.

This module serves as the main orchestrator for the CircleSeeker pipeline, coordinating various steps
including TideHunter analysis, BLAST searches, and multiple specialized processors (Carousel, Ringmaster,
Sieve, Juggler, Tamer, Astrologer) to identify and characterize different types of eccDNA:
- UeccDNA (Unique eccDNA)
- MeccDNA (Multiple-copy eccDNA)
- CeccDNA (Composite eccDNA)
- MCeccDNA/Ceccm (Multiple-copy Composite eccDNA)
- Inferred UeccDNA (UeccDNAs identified through additional analysis)

The pipeline manages temporary files, checkpoints, and generates comprehensive reports of the findings.

Version: 1.0.1
License: GNU General Public License v2
Copyright (c) 2024 CircleSeeker Team
"""

import argparse
import os
import subprocess
import time
import shutil
from Bio import SeqIO

from circle_seeker.Carousel import Carousel
from circle_seeker.Ringmaster import Ringmaster
from circle_seeker.Sieve import Sieve
from circle_seeker.Juggler import Juggler
from circle_seeker.Tamer import Tamer
from circle_seeker.Astrologer import Astrologer
from circle_seeker.MergeUecc import MergeUecc
from circle_seeker.MergeMecc import MergeMecc
from circle_seeker.ReportGenerator import ReportGenerator
from circle_seeker.MergeCecc import CeccProcessor
from circle_seeker.FAIProcessor import FAIProcessor
from circle_seeker.MergeUeccInf import UeccAnalyzer

# from Carousel import Carousel
# from Ringmaster import Ringmaster
# from Sieve import Sieve
# from Juggler import Juggler
# from Tamer import Tamer
# from Astrologer import Astrologer
# from MergeUecc import MergeUecc
# from MergeMecc import MergeMecc
# from ReportGenerator import ReportGenerator
# from MergeCecc import CeccProcessor
# from FAIProcessor import FAIProcessor
# from MergeUeccInf import UeccAnalyzer


class CircleSeeker:
    def __init__(self, input_fasta, output_prefix, reference_genome, num_threads=8, keep_tmp=False, enable_X=False):
        self.input_fasta = input_fasta
        self.output_prefix = output_prefix
        self.reference_genome = reference_genome
        self.num_threads = num_threads
        self.keep_tmp = keep_tmp
        self.enable_X = enable_X

        dirname = os.path.dirname(self.output_prefix)
        if dirname == '':
            self.output_dir = os.getcwd()
        else:
            self.output_dir = dirname
        
        os.makedirs(self.output_dir, exist_ok=True)
        self.tmp_dir = os.path.join(self.output_dir, '.tmp')
        os.makedirs(self.tmp_dir, exist_ok=True)

        base_name = os.path.basename(self.output_prefix)

        # Use a fixed temporary directory instead of a random one
        self.temp_dir = os.path.join(self.tmp_dir, f"{base_name}_temp")
        if not os.path.exists(self.temp_dir):
            os.makedirs(self.temp_dir, exist_ok=True)

        self.blast_db = os.path.join(self.temp_dir, f"{base_name}_blast_db")
        self.tidehunter_output_path = os.path.join(self.temp_dir, f"{base_name}.TH.ecc_candidates.txt")
        self.carousel_output_path = os.path.join(self.temp_dir, f"{base_name}.carousel_processed_results.csv")
        self.carousel_read_list_path = os.path.join(self.temp_dir, f"{base_name}.carousel_read_list.csv")
        self.carousel_circular_sequences_path = os.path.join(self.temp_dir, f"{base_name}.carousel_circular_sequences.fasta")
        self.blast_results1_path = os.path.join(self.temp_dir, f"{base_name}.blast1_out_alignment_results.txt")
        self.ringmaster_uecc_output_path = f"{self.temp_dir}/{base_name}_Uecc.csv"
        self.ringmaster_mecc_output_path = f"{self.temp_dir}/{base_name}_Mecc.csv"
        self.ringmaster_cecc_output_path = f"{self.temp_dir}/{base_name}_Cecc.csv"
        self.ringmaster_uecc_fa_path = f"{self.temp_dir}/{base_name}_UeccDNA.fa"
        self.ringmaster_mecc_fa_path = f"{self.temp_dir}/{base_name}_MeccDNA.fa"
        self.ringmaster_cecc_fa_path = f"{self.temp_dir}/{base_name}_CeccDNA.fa"

        # XeccDNA file (will be generated by Ringmaster when enable_X is True)
        self.ringmaster_xecc_fa_path = f"{self.temp_dir}/{base_name}_XeccDNA.fa"

        self.sieve_1_output_unSelected_reads_path = os.path.join(self.temp_dir, f"{base_name}.sieve_1_output.fa")
        self.blast_results2_path = os.path.join(self.temp_dir, f"{base_name}.blast2_out_unSel.results.txt")
        self.samtools_index_path = f"{self.sieve_1_output_unSelected_reads_path}.fai"
        self.juggler_output_tecc_path = os.path.join(self.temp_dir, f"{base_name}.tecc_analysis_results.csv")
        self.juggler_output_read_class_path = os.path.join(self.temp_dir, f"{base_name}.read_classification.csv")
        self.juggler_output_num_linr_path = os.path.join(self.temp_dir, f"{base_name}.Num_LinR.csv")
        self.tamer_output_file = os.path.join(self.temp_dir, f"{base_name}.tecc.enhanced.results.csv")

        self.sieve_2_output_path = os.path.join(self.temp_dir, f"{base_name}.Final_unclassified.fa")
        self.minimap2_output_final_unclassified_sam = os.path.join(self.temp_dir, f"{base_name}_Final_unclassified.sam")
        self.sorted_bam_path = os.path.join(self.temp_dir, f"{base_name}_Final_unclassified_sorted.bam")
        self.astrologer_output_csv = os.path.join(self.temp_dir, f"{base_name}_Final.Inf.process.Uecc.csv")

        self.final_confirmed_uecc_csv = os.path.join(self.temp_dir, f"{base_name}.Final.Confirmed.Uecc.csv")
        self.final_confirmed_uecc_fasta = os.path.join(self.temp_dir, f"{base_name}.Final.Confirmed.Uecc.fasta")
        self.final_confirmed_mecc_csv = os.path.join(self.temp_dir, f"{base_name}.Final.Confirmed.Mecc.csv")
        self.final_confirmed_mecc_fasta = os.path.join(self.temp_dir, f"{base_name}.Final.Confirmed.Mecc.fasta")

        self.final_confirmed_cecc_csv = os.path.join(self.temp_dir, f"{base_name}.Final.Confirmed.Cecc.csv")
        self.final_confirmed_cecc_fasta = os.path.join(self.temp_dir, f"{base_name}.Final.Confirmed.Cecc.fasta")
        self.final_confirmed_mcecc_csv = os.path.join(self.temp_dir, f"{base_name}.Final.Confirmed.MCecc.csv")
        self.final_confirmed_mcecc_fasta = os.path.join(self.temp_dir, f"{base_name}.Final.Confirmed.MCecc.fasta")

        # Final naming for XeccDNA
        self.final_confirmed_xecc_fasta = os.path.join(self.temp_dir, f"{base_name}.Final.Confirmed.Xecc.fasta")

        self.final_inferred_uecc_csv = os.path.join(self.temp_dir, f"{base_name}_Final.Inferred.Uecc.csv")
        self.final_inferred_uecc_fasta = os.path.join(self.temp_dir, f"{base_name}.Final.Inferred.Uecc.fasta")
        self.final_shared_count_csv = os.path.join(self.temp_dir, f"{base_name}.Final.SharedCount.csv")
        self.final_merged_uecc_csv = os.path.join(self.temp_dir, f"{base_name}_Final.Merged.Uecc.csv")

        self.final_report_txt = os.path.join(self.temp_dir, f"{base_name}.report.txt")
        self.final_summary_csv = os.path.join(self.output_dir, f"{base_name}.summary.csv")
        self.final_report_html = os.path.join(self.output_dir, f"{base_name}.report.html")

        self.checkpoint_file = os.path.join(self.output_dir, f"{base_name}.checkpoint")

        self.steps = [
            (1,  "Creating BLAST database",
             [self.blast_db+".nsq", self.blast_db+".nin", self.blast_db+".nhr"], 
             self.create_blast_db),
            (2,  "Running TideHunter",
             [self.tidehunter_output_path], 
             self.run_tidehunter),
            (3,  "Running Carousel",
             [self.carousel_output_path, self.carousel_circular_sequences_path],
             self.run_carousel),
            (4,  "Running BLASTN on circular sequences",
             [self.blast_results1_path],
             lambda: self.run_blastn(self.carousel_circular_sequences_path, self.blast_results1_path)),
            (5,  "Running Ringmaster",
             [self.ringmaster_uecc_output_path, self.ringmaster_mecc_output_path, self.ringmaster_cecc_output_path],
             self.run_ringmaster),
            (6,  "Running Sieve on original input",
             [self.sieve_1_output_unSelected_reads_path],
             self.run_sieve_1),
            (7,  "Running BLASTN on unselected sequences",
             [self.blast_results2_path],
             lambda: self.run_blastn(self.sieve_1_output_unSelected_reads_path, self.blast_results2_path)),
            (8,  "Creating index for unselected sequences",
             [self.samtools_index_path],
             lambda: self.create_fasta_index(self.sieve_1_output_unSelected_reads_path)),
            (9,  "Running Juggler",
             [self.juggler_output_tecc_path, self.juggler_output_read_class_path, self.juggler_output_num_linr_path],
             self.run_juggler),
            (10, "Running Tamer",
             [self.tamer_output_file],
             self.run_tamer),
            (11, "Running Sieve 2 to filter unclassified reads",
             [self.sieve_2_output_path],
             self.run_sieve_2),
            (12, "Running final alignment with minimap2",
             [self.minimap2_output_final_unclassified_sam],
             self.run_minimap2),
            (13, "Sorting and indexing final alignment results",
             [self.sorted_bam_path, self.sorted_bam_path+".bai"],
             self.run_sort_and_index),
            (14, "Running Astrologer for final analysis and scoring",
             [self.astrologer_output_csv],
             self.run_astrologer),
            (15, "Merging UeccDNA results",
             [self.final_confirmed_uecc_csv, self.final_confirmed_uecc_fasta],
             self.run_merge_uecc),
            (16, "Merging MeccDNA results",
             [self.final_confirmed_mecc_csv, self.final_confirmed_mecc_fasta],
             self.run_merge_mecc),
            (17, "Processing inferred UeccDNA",
             [self.final_inferred_uecc_csv, self.final_inferred_uecc_fasta, self.final_shared_count_csv, self.final_merged_uecc_csv],
             self.run_merge_uecc_inferred),
            (18, "Processing Cecc and MCecc",
             [self.final_confirmed_cecc_csv, self.final_confirmed_cecc_fasta, self.final_confirmed_mcecc_csv, self.final_confirmed_mcecc_fasta],
             self.run_merge_cecc),
            (19, "Processing FASTA indices",
             [],
             self.run_FAIProcessor),
            (20, "Generating final analysis report",
             [self.final_report_txt, self.final_summary_csv, self.final_report_html],
             self.run_report_generator),
            (21, "Saving final results and cleaning up",
             [],
             self.save_and_cleanup)
        ]
        self.total_steps = len(self.steps)

    def load_checkpoint(self):
        completed_step = 0
        if os.path.exists(self.checkpoint_file):
            with open(self.checkpoint_file, 'r') as f:
                line = f.readline().strip()
                if line.isdigit():
                    completed_step = int(line)
        return completed_step

    def save_checkpoint(self, step_num):
        with open(self.checkpoint_file, 'w') as f:
            f.write(str(step_num))

    def check_step_done(self, step_outputs):
        if not step_outputs:
            # No output files defined, cannot determine if step is done
            return False
        for output in step_outputs:
            if not os.path.exists(output):
                return False
        return True

    def create_blast_db(self):
        command = f"makeblastdb -in {self.reference_genome} -dbtype nucl -out {self.blast_db}"
        subprocess.run(command, shell=True, check=True)

    def run_tidehunter(self):
        command = f"TideHunter -f 2 -t {self.num_threads} -k 16 -w 1 -p 100 -P 2000000 -e 0.1 {self.input_fasta} > {self.tidehunter_output_path}"
        subprocess.run(command, shell=True, check=True)

    def run_carousel(self):
        if not os.path.exists(self.tidehunter_output_path):
            raise FileNotFoundError(f"File {self.tidehunter_output_path} does not exist")
        carousel = Carousel(self.tidehunter_output_path, self.carousel_output_path, self.carousel_read_list_path, self.carousel_circular_sequences_path)
        carousel.run()

    def run_blastn(self, input_fasta, output_file):
        blastn_command = [
            "blastn",
            "-db", self.blast_db,
            "-query", input_fasta,
            "-out", output_file,
            "-num_threads", str(self.num_threads),
            "-word_size", "100",
            "-evalue", "1e-50",
            "-perc_identity", "99",
            "-outfmt", "6"
        ]
        subprocess.run(blastn_command, check=True)

    def run_ringmaster(self):
        if not os.path.exists(self.blast_results1_path):
            raise FileNotFoundError(f"File {self.blast_results1_path} does not exist")
        if not os.path.exists(self.carousel_circular_sequences_path):
            raise FileNotFoundError(f"File {self.carousel_circular_sequences_path} does not exist")

        base_name = os.path.basename(self.output_prefix)
        ringmaster = Ringmaster(
            blast_results_file=self.blast_results1_path,
            circular_seq_fasta=self.carousel_circular_sequences_path,
            Uecc_output_csv=self.ringmaster_uecc_output_path,
            Mecc_output_csv=self.ringmaster_mecc_output_path,
            Cecc_output_csv=self.ringmaster_cecc_output_path,
            process_xecc=self.enable_X,
            num_threads=self.num_threads,
            prefix=base_name,
            output_dir=self.temp_dir
        )
        ringmaster.run()
        print("Ringmaster completed.")

        # If enable_X is True, check if XeccDNA.fa exists and rename it
        if self.enable_X:
            if os.path.exists(self.ringmaster_xecc_fa_path):
                os.rename(self.ringmaster_xecc_fa_path, self.final_confirmed_xecc_fasta)
                print(f"XeccDNA fasta renamed to: {self.final_confirmed_xecc_fasta}")

    def run_sieve_1(self):
        sieve_1 = Sieve(self.input_fasta, self.tidehunter_output_path, self.sieve_1_output_unSelected_reads_path, 't')
        sieve_1.run()

    def run_sieve_2(self):
        sieve_2 = Sieve(self.sieve_1_output_unSelected_reads_path, self.juggler_output_read_class_path, self.sieve_2_output_path, 'c')
        sieve_2.run()

    def create_fasta_index(self, fasta_file):
        subprocess.run(["samtools", "faidx", fasta_file], check=True)
        print(f"Created index for {fasta_file}")

    def run_juggler(self):
        juggler = Juggler(self.blast_results2_path, self.samtools_index_path, self.juggler_output_num_linr_path, self.juggler_output_tecc_path, self.juggler_output_read_class_path)
        juggler.run()

    def run_tamer(self):
        tamer = Tamer(
            input_tecc_csv=self.juggler_output_tecc_path,
            input_fasta=self.sieve_1_output_unSelected_reads_path,
            output_csv=self.tamer_output_file
        )
        tamer.run()
        print("Tamer completed.")

    def run_minimap2(self):
        command = [
            "minimap2",
            "-t", str(self.num_threads),
            "-x", "map-hifi",
            "-a",
            self.reference_genome,
            self.sieve_2_output_path
        ]
        with open(self.minimap2_output_final_unclassified_sam, 'w') as sam_file:
            subprocess.run(command, check=True, stdout=sam_file)
        print(f"Minimap2 alignment completed. Output: {self.minimap2_output_final_unclassified_sam}")

    def run_sort_and_index(self):
        subprocess.run(["samtools", "sort", "-o", self.sorted_bam_path, self.minimap2_output_final_unclassified_sam], check=True)
        subprocess.run(["samtools", "index", self.sorted_bam_path], check=True)
        print(f"Sorted BAM file: {self.sorted_bam_path}")

    def run_astrologer(self):
        base_name = os.path.basename(self.output_prefix)
        prefix_in_temp = os.path.join(self.temp_dir, base_name)
        
        # Place Astrologer's filter and input_txt in temp directory with consistent naming
        astrologer_filter_csv = os.path.join(self.temp_dir, f"{base_name}.astrologer_filter.csv")
        astrologer_out_csv_in_temp = os.path.join(self.temp_dir, f"{base_name}_Final.Inf.process.Uecc.csv")

        astrologer = Astrologer(
            input_bam=self.sorted_bam_path,
            prefix=prefix_in_temp,
            mosdepth_threads=self.num_threads,
            # Add these two lines to point output_filter and input_txt to temp files
            output_filter=astrologer_filter_csv,
            input_txt=astrologer_filter_csv,
            
            output_consensus=astrologer_out_csv_in_temp,
            ref_fa=self.reference_genome,
            reads_fa=self.sieve_2_output_path,
            keep_temp=self.keep_tmp,
            temp_dir=self.temp_dir,
            threads=self.num_threads,
            minimap2_threads=self.num_threads,
        )
        astrologer.run_all_steps()
        print(f"Astrologer completed. Results: {self.astrologer_output_csv}")
        
    def run_merge_uecc(self):
        merger_uecc = MergeUecc(
            uecc_part1=self.ringmaster_uecc_output_path,
            uecc_part2=self.tamer_output_file,
            output_csv=self.final_confirmed_uecc_csv,
            output_fasta=self.final_confirmed_uecc_fasta
        )
        merger_uecc.run_merge_uecc()
        print("MergeUecc completed.")

    def run_merge_mecc(self):
        merger_mecc = MergeMecc(
            meccdna_part1=self.ringmaster_mecc_output_path,
            meccdna_part2=self.tamer_output_file,
            output_csv=self.final_confirmed_mecc_csv,
            output_fasta=self.final_confirmed_mecc_fasta
        )
        merger_mecc.run_merge_mecc()
        print("MergeMecc completed.")

    def run_merge_uecc_inferred(self):
        processor = UeccAnalyzer(
            inf_process_csv=self.astrologer_output_csv,
            confirmed_csv=self.final_confirmed_uecc_csv,
            inferred_output=self.final_inferred_uecc_csv,
            shared_count_output=self.final_shared_count_csv,
            fasta_output=self.final_inferred_uecc_fasta,
            merged_output=self.final_merged_uecc_csv,
            keep_tmp=self.keep_tmp
        )
        processor.run()
        print("MergeUeccInf (UeccAnalyzer) completed.")

    def run_merge_cecc(self):
        cecc_processor = CeccProcessor(argparse.Namespace(
            input_csv=self.ringmaster_cecc_output_path,
            output_cecc_fasta=self.final_confirmed_cecc_fasta,
            output_mcecc_fasta=self.final_confirmed_mcecc_fasta,
            final_cecc_csv=self.final_confirmed_cecc_csv,
            final_mcecc_csv=self.final_confirmed_mcecc_csv
        ))
        cecc_processor.run()
        print("MergeCecc completed.")

    def run_FAIProcessor(self):
        fai_processor = FAIProcessor(
            ueccdna=self.final_confirmed_uecc_fasta,
            inferred_ueccdna=self.final_inferred_uecc_fasta,
            meccdna=self.final_confirmed_mecc_fasta,
            ceccdna=self.final_confirmed_cecc_fasta,
            mceccdna=self.final_confirmed_mcecc_fasta
        )
        fai_processor.process_all_files()
        print("FAIProcessor completed.")

    def run_report_generator(self):
        input_fai = f"{self.input_fasta}.fai"
        if not os.path.exists(input_fai):
            subprocess.run(["samtools", "faidx", self.input_fasta], check=True)

        base_name = os.path.basename(self.output_prefix)

        report_generator = ReportGenerator(
            fai=input_fai,
            ctcr1=self.carousel_read_list_path,
            ctcr2=self.juggler_output_read_class_path,
            linr=self.juggler_output_num_linr_path,
            uecc_fai=f"{self.final_confirmed_uecc_fasta}.fai",
            mecc_fai=f"{self.final_confirmed_mecc_fasta}.fai",
            cecc_fai=f"{self.final_confirmed_cecc_fasta}.fai",
            mcecc_fai=f"{self.final_confirmed_mcecc_fasta}.fai",
            uecc_inferred_fai=f"{self.final_inferred_uecc_fasta}.fai",
            shared_count_csv=self.final_shared_count_csv,
            output=self.final_report_txt,
            summary_output=self.final_summary_csv,
            html_output=self.final_report_html,
            sample_name=base_name,
            keep_tmp=self.keep_tmp
        )
        report_generator.run()
        print("ReportGenerator completed.")

    def save_and_cleanup(self):
        key_words = ["Confirmed", "Inferred", "report", "Merged", "Xecc"]
        suffix = [".csv", ".html", ".fa", ".fasta", ".txt"]

        for file in os.listdir(self.temp_dir):
            if any(key_word in file for key_word in key_words) and any(file.endswith(s) for s in suffix):
                src = os.path.join(self.temp_dir, file)
                dst = os.path.join(self.output_dir, file)
                shutil.move(src, dst)
                print(f"Saved final result: {dst}")

        if not self.keep_tmp:
            shutil.rmtree(self.temp_dir)
            print(f"Cleaned up temporary directory: {self.temp_dir}")
        else:
            print(f"Temporary files kept in: {self.temp_dir}")

    def run_pipeline(self, start_from=None, force=False):
        start_time = time.time()
        completed_step = self.load_checkpoint() if start_from is None else 0
        actual_start_step = start_from if start_from is not None else (completed_step + 1)

        for step_info in self.steps:
            step_num, step_desc, step_outputs, step_func = step_info

            if step_num < actual_start_step:
                print(f"> Step {step_num}/{self.total_steps} - {step_desc} has been skipped (start_from is set).")
                continue

            print(f"\n> Step {step_num}/{self.total_steps} - {step_desc}")
            # If not force and output files exist, skip this step
            if (not force) and self.check_step_done(step_outputs):
                print(f"Output files for {step_desc} already exist, skipping this step.")
                completed_step = step_num
                self.save_checkpoint(completed_step)
                continue

            # Execute the corresponding step function
            step_func()

            # Update checkpoint
            completed_step = step_num
            self.save_checkpoint(completed_step)

        total_time = time.time() - start_time
        print(f"\n> Pipeline completed in {total_time:.2f} seconds")
        print(f"Final results saved in: {self.output_dir}")

        key_words = ["Confirmed", "Inferred", "report", "Merged", "Xecc"]
        suffix = [".csv", ".html", ".fa", ".fasta", ".txt"]
        all_final_files = [
            f for f in os.listdir(self.output_dir)
            if any(key_word in f for key_word in key_words) and any(f.endswith(s) for s in suffix)
        ]
        print("Saved files:")
        for f in all_final_files:
            print(f" - {f}")

def main():
    parser = argparse.ArgumentParser(description="Wrapper script for eccDNA analysis software package")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file, e.g. demo.10k.fasta")
    parser.add_argument("-p", "--prefix", required=True, help="Prefix for all generated files")
    parser.add_argument("-r", "--reference", required=True, help="Reference genome FASTA file, e.g. T2T-CHM13v2.0_chr.fna")
    parser.add_argument("-o", "--output", default=None, help="Output folder (optional)")
    parser.add_argument("-t", "--num_threads", type=int, default=8, help="Number of threads to use (default: 8)")
    parser.add_argument("-kt", "--keep_tmp", action="store_true", help="Keep temporary files in the .tmp directory")
    parser.add_argument("--enable_X", action="store_true", help="Enable the -X parameter for Ringmaster to process XeccDNA")
    parser.add_argument("--start_from", type=int, default=None, help="Start from this step number (ignore checkpoint and steps before this)")
    parser.add_argument("--force", action="store_true", help="Force re-run steps even if output files exist")

    args = parser.parse_args()

    if args.output:
        os.makedirs(args.output, exist_ok=True)
        output_prefix = os.path.join(args.output, args.prefix)
    else:
        output_prefix = args.prefix

    wrapper = CircleSeeker(args.input, output_prefix, args.reference, num_threads=args.num_threads, keep_tmp=args.keep_tmp, enable_X=args.enable_X)
    wrapper.run_pipeline(start_from=args.start_from, force=args.force)

if __name__ == "__main__":
    main()
