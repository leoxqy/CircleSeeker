import argparse
import os
import subprocess
import time
import multiprocessing
from Bio import SeqIO
from Carousel import Carousel
from Ringmaster import Ringmaster
from Sieve import Sieve
from Juggler import Juggler
from Tamer import Tamer
from Astrologer import Astrologer
from MergeUecc import MergeUecc
from MergeMecc import MergeMecc
from TreatOther import TreatOther
from ReportGenerator import ReportGenerator
from Cleaner import Cleaner

class CircleSeeker:
    def __init__(self, input_fasta, output_prefix, reference_genome, existing_db=None):
        """Initialize the CircleSeeker with input and output parameters."""
        self.input_fasta = input_fasta
        self.output_prefix = output_prefix
        self.reference_genome = reference_genome
        self.existing_db = existing_db
        self.blast_db = existing_db if existing_db else f"{output_prefix}_blast_db"
        
        # Step1: tidehunter output
        self.tidehunter_output_path = None
        
        # Step2: carousel output
        self.carousel_output_path = None
        self.carousel_read_list_path = None
        self.carousel_circular_sequences_path = None
        
        # Step3: blastn_1 output
        self.blast_results1_path = None

        # Step4: ringmaster output
        self.ringmaster_uecc_output_path = None
        self.ringmaster_mecc_output_path = None
        self.ringmaster_cecc_output_path = None
        self.ringmaster_uecc_fa_path = None
        self.ringmaster_mecc_fa_path = None
        self.ringmaster_cecc_fa_path = None

        # Step5: sieve_1 output
        self.sieve_1_output_unSelected_reads_path = None

        # Step6: blastn_2 output
        self.blast_results2_path = None

        # Step7: samtools create index
        self.samtools_index_path = None

        # Step8: juggler output
        self.juggler_output_tecc_path = None
        self.juggler_output_read_class_path = None
        self.juggler_output_num_linr_path = None

        # Step9: tamer output
        self.tamer_output_path = None
        
        # Step10: sieve_2 output
        self.sieve_2_output_path = None

        # Step11: minimap2 output
        self.minimap2_output_final_unclassified_sam = None
        
        # Step12: SortAndIndex output
        self.sorted_bam_path = None
        self.sorted_bam_index_path = None

        # Step13: Astrologer output
        self.astrologer_output_csv = None

        # Step14: MergeUecc output
        self.merged_uecc_csv = None
        self.merged_uecc_fasta = None

        # Step15: MergeMecc output
        self.merged_mecc_csv = None
        self.merged_mecc_fasta = None

        # Step16: Other eccDNA types output
        self.xecc_confirmed_csv = None
        self.cecc_confirmed_csv = None
        self.uecc_inferred_csv = None  # Add this line

        # Step17: ReportGenerator output
        self.final_report = None

        self.total_steps = 18 if existing_db else 19  # Increment total steps
        self.num_threads = max(1, multiprocessing.cpu_count() - 2)  # Use all cores except 2

    def create_blast_db(self):
        """Create a BLAST database from the reference genome."""
        command = f"makeblastdb -in {self.reference_genome} -dbtype nucl -out {self.blast_db}"
        subprocess.run(command, shell=True, check=True)

    def run_tidehunter(self):
        """Execute TideHunter for consensus sequence generation."""
        tidehunter_output_path = f"{self.output_prefix}.TH.ecc_candidates.txt"
        command = f"TideHunter -f 2 -t {self.num_threads} -k 16 -w 1 -p 100 -P 2000000 -e 0.1 {self.input_fasta} > {tidehunter_output_path}"
        self.tidehunter_output_path = tidehunter_output_path
        subprocess.run(command, shell=True, check=True)

    def run_carousel(self):
        """Process TideHunter output and identify circular sequences."""
        carousel_output = f"{self.output_prefix}.carousel_processed_results.csv"
        read_list = f"{self.output_prefix}.carousel_read_list.csv"
        circular_sequences = f"{self.output_prefix}.carousel_circular_sequences.fasta"
        # check if the file exists
        if not os.path.exists(self.tidehunter_output_path):
            raise FileNotFoundError(f"File {self.tidehunter_output_path} does not exist")
        
        carousel = Carousel(self.tidehunter_output_path, carousel_output, read_list, circular_sequences)
        carousel.run()

        self.carousel_output_path = carousel_output
        self.carousel_read_list_path = read_list
        self.carousel_circular_sequences_path = circular_sequences

    def run_blastn(self, input_fasta, output_file):
        """Perform BLASTN search against the created BLAST database using the conda BLASTN tool."""
        # debug check which one is none
        
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
        
        try:
            print(f"... Running BLASTN on {self.num_threads} threads ...")
            subprocess.run(blastn_command, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running BLASTN: {e}")
            raise

    def run_ringmaster(self):
        """Classify eccDNAs based on BLAST results."""
        Uecc_output_csv = f"{self.output_prefix}_Uecc.csv"
        Mecc_output_csv = f"{self.output_prefix}_Mecc.csv"
        Cecc_output_csv = f"{self.output_prefix}_Cecc.csv"
        uecc_fa = f"{self.output_prefix}_UeccDNA.fa"
        mecc_fa = f"{self.output_prefix}_MeccDNA.fa"
        cecc_fa = f"{self.output_prefix}_CeccDNA.fa"
        xecc_fa = f"{self.output_prefix}_XeccDNA.fa"

        # check if the file exists
        if not os.path.exists(self.blast_results1_path):
            raise FileNotFoundError(f"File {self.blast_results1_path} does not exist")
        if not os.path.exists(self.carousel_circular_sequences_path):
            raise FileNotFoundError(f"File {self.carousel_circular_sequences_path} does not exist")

        ringmaster = Ringmaster(self.blast_results1_path, self.carousel_circular_sequences_path, 
                                Uecc_output_csv, Mecc_output_csv, Cecc_output_csv, 
                                uecc_fa, mecc_fa, cecc_fa, xecc_fa, self.num_threads)
        ringmaster.run()

        self.ringmaster_uecc_output_path = Uecc_output_csv
        self.ringmaster_mecc_output_path = Mecc_output_csv
        self.ringmaster_cecc_output_path = Cecc_output_csv
        self.ringmaster_uecc_fa_path = uecc_fa
        self.ringmaster_mecc_fa_path = mecc_fa
        self.ringmaster_cecc_fa_path = cecc_fa
        self.ringmaster_xecc_fa_path = xecc_fa

    def run_sieve_1(self):
        """Extract unselected reads for further analysis."""
        sieve_1_output_unSelected_reads_path = f"{self.output_prefix}.sieve_1_output.txt"
        sieve_1 = Sieve(self.input_fasta, self.tidehunter_output_path, sieve_1_output_unSelected_reads_path,'t')
        sieve_1.run()

        self.sieve_1_output_unSelected_reads_path = sieve_1_output_unSelected_reads_path

    def run_sieve_2(self):
        """Further filter unclassified reads."""
        sieve_2_output_path = f"{self.output_prefix}.Final_unclassified.fa"
        sieve_2 = Sieve(self.sieve_1_output_unSelected_reads_path, 
                        self.juggler_output_read_class_path, 
                        sieve_2_output_path, 
                        'c')
        sieve_2.run()

        self.sieve_2_output_path = sieve_2_output_path

    def create_fasta_index(self, fasta_file):
        """Create an index for the FASTA file using samtools."""
        try:
            subprocess.run(["samtools", "faidx", fasta_file], check=True)
            print(f"Created index for {fasta_file}")
            self.samtools_index_path = f'{fasta_file}.fai'

        except subprocess.CalledProcessError as e:
            print(f"Error creating index for {fasta_file}: {e}")
            raise

    def run_juggler(self):
        """Run Juggler for additional analysis."""
        output_tecc = f"{self.output_prefix}.tecc_analysis_results.csv"
        output_read_class = f"{self.output_prefix}.read_classification.csv"
        output_num_linr = f"{self.output_prefix}.Num_LinR.csv"
        
        blast_results = self.blast_results2_path
        fasta_index = self.samtools_index_path

        juggler = Juggler(blast_results, fasta_index, output_num_linr, output_tecc,output_read_class)
        juggler.run()
        self.juggler_output_tecc_path = output_tecc
        self.juggler_output_read_class_path = output_read_class
        self.juggler_output_num_linr_path = output_num_linr

    def run_tamer(self):
        """Run Tamer to filter out UeccDNA and MeccDNA."""
        output_uecc = f"{self.output_prefix}_uecc_part2.fa"
        output_mecc = f"{self.output_prefix}_mecc_part2.fa"
        input_tecc_csv = self.juggler_output_tecc_path
        input_unSel_reads_fasta = self.sieve_1_output_unSelected_reads_path

        tamer = Tamer(input_tecc_csv, input_unSel_reads_fasta, output_uecc, output_mecc)
        self.tamer_output_uecc_path = output_uecc
        self.tamer_output_mecc_path = output_mecc
        tamer.run()
        

    def run_minimap2(self):
        """Run minimap2 for final alignment of unclassified reads."""
        input_fasta = self.sieve_2_output_path
        output_sam = f"{self.output_prefix}_Final_unclassified.sam"
        
        command = [
            "minimap2",
            "-t", str(self.num_threads),
            "-x", "map-hifi",
            "-a",
            self.reference_genome,
            input_fasta
        ]
        
        try:
            print(f"... Running minimap2 on {self.num_threads} threads ...")
            with open(output_sam, 'w') as sam_file:
                subprocess.run(command, check=True, stdout=sam_file)
            print(f"Minimap2 alignment completed. Output saved to {output_sam}")
            self.minimap2_output_final_unclassified_sam = output_sam
        except subprocess.CalledProcessError as e:
            print(f"Error running minimap2: {e}")
            raise

    def run_sort_and_index(self):
        """Sort and index the final alignment results."""
        sorted_bam = f"{self.output_prefix}_Final_unclassified_sorted.bam"
        
        try:
            print("Sorting SAM file...")
            subprocess.run(["samtools", "sort", "-o", sorted_bam, self.minimap2_output_final_unclassified_sam], check=True)
            print(f"Sorted BAM file saved to {sorted_bam}")
            
            print("Indexing sorted BAM file...")
            subprocess.run(["samtools", "index", sorted_bam], check=True)
            print(f"Index created for {sorted_bam}")
            
            self.sorted_bam_path = sorted_bam
            self.sorted_bam_index_path = f"{sorted_bam}.bai"
        except subprocess.CalledProcessError as e:
            print(f"Error during sort and index: {e}")
            raise

    def run_astrologer(self):
        """Run Astrologer for final analysis and scoring."""
        input_bam = self.sorted_bam_path
        output_csv = f"{self.output_prefix}_circular_dna_results_scored.csv"
        
        astrologer = Astrologer(input_bam, output_csv)
        astrologer.process_data()
        
        self.astrologer_output_csv = output_csv
        print(f"Astrologer analysis completed. Results saved to {output_csv}")

    def run_merge_uecc(self):
        """Merge UeccDNA results from part1 and part2."""
        ueccdna_part1 = self.ringmaster_uecc_output_path
        tecc_analysis = self.juggler_output_tecc_path
        uecc_part1 = self.ringmaster_uecc_fa_path
        uecc_part2 = self.tamer_output_uecc_path
        output_csv = f"{self.output_prefix}_UeccDNA.Confirmed.csv"
        output_fasta = f"{self.output_prefix}_UeccDNA.Confirmed.fa"

        merger = MergeUecc(ueccdna_part1, tecc_analysis, uecc_part1, uecc_part2, output_csv, output_fasta)
        merger.run_merge_uecc()

        self.merged_uecc_csv = output_csv
        self.merged_uecc_fasta = output_fasta
        print(f"UeccDNA results merged. CSV saved to {output_csv}, FASTA saved to {output_fasta}")

    def run_merge_mecc(self):
        """Merge MeccDNA results from part1 and part2."""
        meccdna_part1 = self.ringmaster_mecc_output_path
        tecc_analysis = self.juggler_output_tecc_path
        mecc_part1 = self.ringmaster_mecc_fa_path
        mecc_part2 = self.tamer_output_mecc_path
        output_csv = f"{self.output_prefix}_MeccDNA.Confirmed.csv"
        output_fasta = f"{self.output_prefix}_MeccDNA.Confirmed.fa"

        merger = MergeMecc(meccdna_part1, tecc_analysis, mecc_part1, mecc_part2, output_csv, output_fasta)
        merger.run_merge_mecc()

        self.merged_mecc_csv = output_csv
        self.merged_mecc_fasta = output_fasta
        print(f"MeccDNA results merged. CSV saved to {output_csv}, FASTA saved to {output_fasta}")

    def run_treat_other(self):
        """Process other types of eccDNA, including XeccDNA and CeccDNA."""
        xecc_fa = self.ringmaster_xecc_fa_path
        xecc_fai = f"{self.ringmaster_xecc_fa_path}.fai"
        cecc_csv = self.ringmaster_cecc_output_path
        cecc_fa = self.ringmaster_cecc_fa_path
        uecc_csv = self.astrologer_output_csv
        
        xecc_output = f"{self.output_prefix}_XeccDNA.Confirmed.csv"
        cecc_output = f"{self.output_prefix}_CeccDNA.Confirmed.csv"
        uecc_inferred_output = f"{self.output_prefix}_UeccDNA.Inferred.csv"
        
        # Create index for xecc.fa
        subprocess.run(["samtools", "faidx", xecc_fa], check=True)
        
        # Run TreatOther
        treat_other = TreatOther(
            xecc_fai=xecc_fai,
            xecc_fa=xecc_fa,
            cecc_csv=cecc_csv,
            cecc_fa=cecc_fa,
            uecc_csv=uecc_csv,
            xecc_output=xecc_output,
            cecc_output=cecc_output,
            uecc_inferred_output=uecc_inferred_output
        )
        treat_other.process_eccDNA()
        
        self.xecc_confirmed_csv = xecc_output
        self.cecc_confirmed_csv = cecc_output
        self.uecc_inferred_csv = uecc_inferred_output
        print(f"Other eccDNA types processed. XeccDNA results saved to {xecc_output}, CeccDNA results saved to {cecc_output}, UeccDNA inferred results saved to {uecc_inferred_output}")

    def run_report_generator(self):
        """Generate the final analysis report."""
        fai = f"{self.input_fasta}.fai"
        ctcr1 = self.carousel_read_list_path
        ctcr2 = self.juggler_output_read_class_path
        linr = self.juggler_output_num_linr_path
        uecc = self.merged_uecc_csv
        mecc = self.merged_mecc_csv
        xecc = self.xecc_confirmed_csv
        cecc = self.cecc_confirmed_csv
        uecc_inferred = self.uecc_inferred_csv  # Use the new attribute
        output = f"{self.output_prefix}_eccDNA_analysis_report.txt"

        # Create index for input fasta if it doesn't exist
        if not os.path.exists(fai):
            subprocess.run(["samtools", "faidx", self.input_fasta], check=True)

        report_generator = ReportGenerator(
            fai, ctcr1, ctcr2, linr, uecc, mecc, xecc, cecc, uecc_inferred, output
        )
        report_generator.run()

        self.final_report = output
        print(f"Final analysis report generated and saved to {output}")

    def run_cleaner(self):
        """Run the Cleaner to remove intermediate files."""
        cleaner = Cleaner(self.output_prefix)
        cleaner.run()

    def run_pipeline(self):
        """Execute the complete eccDNA analysis pipeline."""
        start_time = time.time()

        def log_progress(step, description):
            elapsed_time = time.time() - start_time
            print(f"\n> Step {step}/{self.total_steps} - {description} - Elapsed time: {elapsed_time:.2f} seconds")

        step = 1
        if not self.existing_db:
            log_progress(step, "Creating BLAST database")
            self.create_blast_db()
            step += 1

        # Step 2 (or 1 if skipping BLAST DB creation): Run TideHunter
        log_progress(step, "Running TideHunter")
        self.run_tidehunter()
        step += 1

        # Step 3: Run Carousel
        log_progress(step, "Running Carousel")
        self.run_carousel()
        step += 1

        # Step 4: Run BLASTN on circular sequences
        log_progress(step, "Running BLASTN on circular sequences")
        blast_results1_path = f"{self.output_prefix}.blast1_out_alignment_results.txt"  
        self.run_blastn(self.carousel_circular_sequences_path, blast_results1_path)
        self.blast_results1_path = blast_results1_path
        step += 1

        # Step 5: Run Ringmaster
        log_progress(step, "Running Ringmaster")
        self.run_ringmaster()
        step += 1

        # Step 6: CircleSeeker_Sieve_1 Run Sieve on original input
        log_progress(step, "Running Sieve on original input")
        self.run_sieve_1()
        step += 1

        # Step 7: Run BLASTN on unselected sequences
        log_progress(step, "Running BLASTN on unselected sequences")
        blast_results2_path = f"{self.output_prefix}.blast2_out_unSel.results.txt"
        self.run_blastn(self.sieve_1_output_unSelected_reads_path, blast_results2_path)
        self.blast_results2_path = blast_results2_path
        step += 1

        # Step 8: Create index for unselected sequences with samtools
        log_progress(step, "Creating index for unselected sequences")
        self.create_fasta_index(self.sieve_1_output_unSelected_reads_path)
        step += 1

        # Step 9: Run Juggler
        log_progress(step, "Running Juggler")
        self.run_juggler()
        step += 1

        # Step 10: Run Tamer to filter out UeccDNA and MeccDNA
        log_progress(step, "Running Tamer")
        self.run_tamer()
        step += 1
        
        # Step 11: CircleSeeker_Sieve_2 - Further filter unclassified reads
        log_progress(step, "Running Sieve 2 to filter unclassified reads")
        self.run_sieve_2()
        step += 1

        # Step 12: CircleSeeker_FinalAlign - Final alignment of unclassified reads
        log_progress(step, "Running final alignment with minimap2")
        self.run_minimap2()
        step += 1

        # New Step 13: CircleSeeker_SortAndIndex - Sort and index final alignment results
        log_progress(step, "Sorting and indexing final alignment results")
        self.run_sort_and_index()
        step += 1

        # New Step 14: CircleSeeker_Astrologer - Final analysis and scoring
        log_progress(step, "Running Astrologer for final analysis and scoring")
        self.run_astrologer()
        step += 1

        # New Step 15: CircleSeeker_MergeUecc - Merge UeccDNA results
        log_progress(step, "Merging UeccDNA results")
        self.run_merge_uecc()
        step += 1

        # New Step 16: CircleSeeker_MergeMecc - Merge MeccDNA results
        log_progress(step, "Merging MeccDNA results")
        self.run_merge_mecc()
        step += 1

        # New Step 17: CircleSeeker_TreatOther - Process other eccDNA types
        log_progress(step, "Processing other eccDNA types")
        self.run_treat_other()
        step += 1

        # New Step 18: CircleSeeker_ReportGenerator - Generate final analysis report
        log_progress(step, "Generating final analysis report")
        self.run_report_generator()
        step += 1

        # New Step 19: CircleSeeker_Cleaner - Clean up intermediate files
        log_progress(self.total_steps, "Cleaning up intermediate files")
        self.run_cleaner()

        total_time = time.time() - start_time
        print(f"\n> Pipeline completed in {total_time:.2f} seconds")

def main():
    """Parse command-line arguments and run the eccDNA pipeline."""
    parser = argparse.ArgumentParser(description="Wrapper script for eccDNA analysis pipeline")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output_prefix", required=True, help="Output prefix for all generated files")
    parser.add_argument("-r", "--reference", required=True, help="Reference genome FASTA file")
    parser.add_argument("-d", "--db", help="Existing BLAST database name (optional)")
    args = parser.parse_args()

    wrapper = CircleSeeker(args.input, args.output_prefix, args.reference, args.db)
    wrapper.run_pipeline()

if __name__ == "__main__":
    main()