import argparse
import os
import subprocess
import time
import multiprocessing
import tempfile
import shutil
from Bio import SeqIO
# build
from circle_seeker.Carousel import Carousel
from circle_seeker.Ringmaster import Ringmaster
from circle_seeker.Sieve import Sieve
from circle_seeker.Juggler import Juggler
from circle_seeker.Tamer import Tamer
from circle_seeker.Astrologer import Astrologer
from circle_seeker.MergeUecc import MergeUecc
from circle_seeker.MergeMecc import MergeMecc
from circle_seeker.ReportGenerator import ReportGenerator
from circle_seeker.FAIProcessor import FAIProcessor
from circle_seeker.MergeUeccInf import UeccProcessor
from circle_seeker.MergeCecc import CeccProcessor

# Debug
# from Carousel import Carousel
# from Ringmaster import Ringmaster
# from Sieve import Sieve
# from Juggler import Juggler
# from Tamer import Tamer
# from Astrologer import Astrologer
# from MergeUecc import MergeUecc
# from MergeMecc import MergeMecc
# from ReportGenerator import ReportGenerator

# Remove TreatOther import and add new imports
# from circle_seeker.TreatOther import TreatOther
# from MergeUeccInf import UeccProcessor
# from MergeCecc import CeccProcessor
# from MergeMecc import MergeMecc

# Add import at the top with other imports
# from FAIProcessor import FAIProcessor

class CircleSeeker:
    def __init__(self, input_fasta, output_prefix, reference_genome, num_threads=8, keep_tmp=False):
        """Initialize the CircleSeeker with input and output parameters."""
        self.input_fasta = input_fasta
        self.output_prefix = output_prefix
        self.reference_genome = reference_genome
        self.num_threads = num_threads
        self.keep_tmp = keep_tmp

        # Determine the output directory
        self.output_dir = os.path.dirname(self.output_prefix) if self.output_prefix else os.getcwd()

        # Create a hidden '.tmp' directory within the output directory
        self.tmp_dir = os.path.join(self.output_dir, '.tmp')
        os.makedirs(self.tmp_dir, exist_ok=True)

        # Create a temporary directory within the .tmp directory
        self.temp_dir = tempfile.mkdtemp(prefix="circleseeker_", dir=self.tmp_dir)

        # Update all file paths to use the new temp_dir
        self.blast_db = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}_blast_db")
        self.tidehunter_output_path = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}.TH.ecc_candidates.txt")
        self.carousel_output_path = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}.carousel_processed_results.csv")
        self.carousel_read_list_path = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}.carousel_read_list.csv")
        self.carousel_circular_sequences_path = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}.carousel_circular_sequences.fasta")
        self.blast_results1_path = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}.blast1_out_alignment_results.txt")
        self.ringmaster_uecc_output_path = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}_Uecc.csv")
        self.ringmaster_mecc_output_path = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}_Mecc.csv")
        self.ringmaster_cecc_output_path = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}_Cecc.csv")
        self.ringmaster_uecc_fa_path = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}_UeccDNA.fa")
        self.ringmaster_mecc_fa_path = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}_MeccDNA.fa")
        self.ringmaster_cecc_fa_path = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}_CeccDNA.fa")
        # self.ringmaster_xecc_fa_path = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}_XeccDNA.fa")
        self.sieve_1_output_unSelected_reads_path = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}.sieve_1_output.txt")
        self.blast_results2_path = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}.blast2_out_unSel.results.txt")
        self.samtools_index_path = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}.sieve_1_output.txt.fai")
        self.juggler_output_tecc_path = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}.tecc_analysis_results.csv")
        self.juggler_output_read_class_path = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}.read_classification.csv")
        self.juggler_output_num_linr_path = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}.Num_LinR.csv")
        self.tamer_output_uecc_path = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}_uecc_part2.fa")
        self.tamer_output_mecc_path = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}_mecc_part2.fa")
        self.sieve_2_output_path = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}.Final_unclassified.fa")
        self.minimap2_output_final_unclassified_sam = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}_Final_unclassified.sam")
        self.sorted_bam_path = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}_Final_unclassified_sorted.bam")
        self.sorted_bam_index_path = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}_Final_unclassified_sorted.bam.bai")
        self.astrologer_output_csv = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}_circular_dna_results_scored.csv")
        self.merged_uecc_csv = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}_UeccDNA.Confirmed.csv")
        self.merged_uecc_fasta = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}_UeccDNA.Confirmed.fa")
        self.merged_mecc_csv = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}_MeccDNA.Confirmed.csv")
        self.merged_mecc_fasta = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}_MeccDNA.Confirmed.fa")
        self.cecc_confirmed_csv = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}_CeccDNA.Confirmed.csv")
        self.uecc_inferred_csv = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}_UeccDNA.Inferred.csv")
        # self.xecc_confirmed_csv = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}_XeccDNA.Confirmed.csv")
        # self.xecc_fa_output = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}_XeccDNA.Confirmed.fa")
        self.cecc_fa_output = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}_CeccDNA.Confirmed.fa")
        self.final_report = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}_eccDNA_analysis_report.txt")

        # Update file paths for new processors
        self.uecc_inferred_output = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}_UeccDNA.Inferred.csv")
        self.cecc_confirmed_fa = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}_CeccDNA.Confirmed.fa")
        self.mcecc_confirmed_fa = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}_MCeccDNA.Confirmed.fa")
        self.cecc_confirmed_csv = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}_CeccDNA.Confirmed.csv")
        self.mcecc_confirmed_csv = os.path.join(self.temp_dir, f"{os.path.basename(output_prefix)}_MCeccDNA.Confirmed.csv")

        self.total_steps = 21  # Always create a new BLAST database

    def create_blast_db(self):
        """Create a BLAST database from the reference genome."""
        command = f"makeblastdb -in {self.reference_genome} -dbtype nucl -out {self.blast_db}"
        subprocess.run(command, shell=True, check=True)

    def run_tidehunter(self):
        """Execute TideHunter for consensus sequence generation."""
        command = f"TideHunter -f 2 -t {self.num_threads} -k 16 -w 1 -p 100 -P 2000000 -e 0.1 {self.input_fasta} > {self.tidehunter_output_path}"
        subprocess.run(command, shell=True, check=True)

    def run_carousel(self):
        """Process TideHunter output and identify circular sequences."""
        if not os.path.exists(self.tidehunter_output_path):
            raise FileNotFoundError(f"File {self.tidehunter_output_path} does not exist")
        
        carousel = Carousel(self.tidehunter_output_path, self.carousel_output_path, self.carousel_read_list_path, self.carousel_circular_sequences_path)
        carousel.run()

    def run_blastn(self, input_fasta, output_file):
        """Perform BLASTN search against the created BLAST database using the conda BLASTN tool."""
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
        if not os.path.exists(self.blast_results1_path):
            raise FileNotFoundError(f"File {self.blast_results1_path} does not exist")
        if not os.path.exists(self.carousel_circular_sequences_path):
            raise FileNotFoundError(f"File {self.carousel_circular_sequences_path} does not exist")

        ringmaster = Ringmaster(self.blast_results1_path, self.carousel_circular_sequences_path, 
                                self.ringmaster_uecc_output_path, self.ringmaster_mecc_output_path, self.ringmaster_cecc_output_path, 
                                self.ringmaster_uecc_fa_path, self.ringmaster_mecc_fa_path, self.ringmaster_cecc_fa_path,self.num_threads)
        ringmaster.run()

    def run_sieve_1(self):
        """Extract unselected reads for further analysis."""
        sieve_1 = Sieve(self.input_fasta, self.tidehunter_output_path, self.sieve_1_output_unSelected_reads_path, 't')
        sieve_1.run()

    def run_sieve_2(self):
        """Further filter unclassified reads."""
        sieve_2 = Sieve(self.sieve_1_output_unSelected_reads_path, 
                        self.juggler_output_read_class_path, 
                        self.sieve_2_output_path, 
                        'c')
        sieve_2.run()

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
        juggler = Juggler(self.blast_results2_path, self.samtools_index_path, self.juggler_output_num_linr_path, self.juggler_output_tecc_path, self.juggler_output_read_class_path)
        juggler.run()

    def run_tamer(self):
        """Run Tamer to filter out UeccDNA and MeccDNA."""
        tamer = Tamer(self.juggler_output_tecc_path, self.sieve_1_output_unSelected_reads_path, self.tamer_output_uecc_path, self.tamer_output_mecc_path)
        tamer.run()

    def run_minimap2(self):
        """Run minimap2 for final alignment of unclassified reads."""
        command = [
            "minimap2",
            "-t", str(self.num_threads),
            "-x", "map-hifi",
            "-a",
            self.reference_genome,
            self.sieve_2_output_path
        ]
        
        try:
            print(f"... Running minimap2 on {self.num_threads} threads ...")
            with open(self.minimap2_output_final_unclassified_sam, 'w') as sam_file:
                subprocess.run(command, check=True, stdout=sam_file)
            print(f"Minimap2 alignment completed. Output saved to {self.minimap2_output_final_unclassified_sam}")
        except subprocess.CalledProcessError as e:
            print(f"Error running minimap2: {e}")
            raise

    def run_sort_and_index(self):
        """Sort and index the final alignment results."""
        try:
            print("Sorting SAM file...")
            subprocess.run(["samtools", "sort", "-o", self.sorted_bam_path, self.minimap2_output_final_unclassified_sam], check=True)
            print(f"Sorted BAM file saved to {self.sorted_bam_path}")
            
            print("Indexing sorted BAM file...")
            subprocess.run(["samtools", "index", self.sorted_bam_path], check=True)
            print(f"Index created for {self.sorted_bam_path}")
        except subprocess.CalledProcessError as e:
            print(f"Error during sort and index: {e}")
            raise

    def run_astrologer(self):
        """Run Astrologer for final analysis and scoring."""
        astrologer = Astrologer(self.sorted_bam_path, self.astrologer_output_csv)
        astrologer.process_data()
        print(f"Astrologer analysis completed. Results saved to {self.astrologer_output_csv}")

    def run_merge_uecc(self):
        """Merge UeccDNA results from part1 and part2."""
        merger = MergeUecc(self.ringmaster_uecc_output_path, self.juggler_output_tecc_path, self.ringmaster_uecc_fa_path, self.tamer_output_uecc_path, self.merged_uecc_csv, self.merged_uecc_fasta)
        merger.run_merge_uecc()
        print(f"UeccDNA results merged. CSV saved to {self.merged_uecc_csv}, FASTA saved to {self.merged_uecc_fasta}")

    def run_merge_mecc(self):
        """Merge MeccDNA results from part1 and part2."""
        merger = MergeMecc(self.ringmaster_mecc_output_path, self.juggler_output_tecc_path, self.ringmaster_mecc_fa_path, self.tamer_output_mecc_path, self.merged_mecc_csv, self.merged_mecc_fasta)
        merger.run_merge_mecc()
        print(f"MeccDNA results merged. CSV saved to {self.merged_mecc_csv}, FASTA saved to {self.merged_mecc_fasta}")

    def run_merge_uecc_inferred(self):
        """Process inferred UeccDNA results."""
        processor = UeccProcessor(self.astrologer_output_csv, self.uecc_inferred_output)
        processor.process_uecc()
        print(f"UeccDNA inferred results processed. Results saved to {self.uecc_inferred_output}")

    def run_merge_cecc(self):
        """Process Cecc and MCecc results."""
        # Create a namespace object to mimic argparse.Namespace
        class Args:
            def __init__(self, **kwargs):
                self.__dict__.update(kwargs)

        args = Args(
            input_csv=self.ringmaster_cecc_output_path,
            fasta_input=self.ringmaster_cecc_fa_path,
            output_cecc_fasta=self.cecc_confirmed_fa,
            output_mcecc_fasta=self.mcecc_confirmed_fa,
            final_cecc_csv=self.cecc_confirmed_csv,
            final_mcecc_csv=self.mcecc_confirmed_csv
        )
        
        processor = CeccProcessor(args)
        processor.run()
        print(f"Cecc results processed. Results saved to {self.cecc_confirmed_csv} and {self.mcecc_confirmed_csv}")

    def run_report_generator(self):
        """Generate the final analysis report."""
        # Create index for input fasta if it doesn't exist
        input_fai = f"{self.input_fasta}.fai"
        if not os.path.exists(input_fai):
            subprocess.run(["samtools", "faidx", self.input_fasta], check=True)

        # Initialize ReportGenerator with all required parameters
        report_generator = ReportGenerator(
            fai=input_fai,
            ctcr1=self.carousel_read_list_path,
            ctcr2=self.juggler_output_read_class_path,
            linr=self.juggler_output_num_linr_path,
            uecc_fai=f"{self.merged_uecc_fasta}.fai",
            mecc_fai=f"{self.merged_mecc_fasta}.fai",
            cecc_fai=f"{self.cecc_confirmed_fa}.fai",
            mcecc_fai=f"{self.mcecc_confirmed_fa}.fai",
            uecc_inferred=self.uecc_inferred_output,
            output=self.final_report,
            summary_output=os.path.join(self.output_dir, f"{os.path.basename(self.output_prefix)}_summary.csv"),
            html_output=os.path.join(self.output_dir, f"{os.path.basename(self.output_prefix)}_report.html")
        )

        # Run report generation
        report_generator.run()
        print(f"Final analysis report generated and saved to {self.final_report}")
        print(f"Summary CSV saved to {os.path.join(self.output_dir, f'{os.path.basename(self.output_prefix)}_summary.csv')}")
        print(f"HTML report saved to {os.path.join(self.output_dir, f'{os.path.basename(self.output_prefix)}_report.html')}")

    def run_FAIProcessor(self):
        """Process and create indices for all eccDNA FASTA files."""
        fai_processor = FAIProcessor(
            ueccdna=self.merged_uecc_fasta,
            meccdna=self.merged_mecc_fasta,
            ceccdna=self.cecc_confirmed_fa,
            mceccdna=self.mcecc_confirmed_fa
        )
        fai_processor.process_all_files()
        print("FASTA indices processing completed")

    def run_pipeline(self):
        """Execute the complete eccDNA analysis pipeline."""
        start_time = time.time()

        def log_progress(step, description):
            elapsed_time = time.time() - start_time
            print(f"\n> Step {step}/{self.total_steps} - {description} - Elapsed time: {elapsed_time:.2f} seconds")

        step = 1
        log_progress(step, "Creating BLAST database")
        self.create_blast_db()
        step += 1

        log_progress(step, "Running TideHunter")
        self.run_tidehunter()
        step += 1

        log_progress(step, "Running Carousel")
        self.run_carousel()
        step += 1

        log_progress(step, "Running BLASTN on circular sequences")
        self.run_blastn(self.carousel_circular_sequences_path, self.blast_results1_path)
        step += 1

        log_progress(step, "Running Ringmaster")
        self.run_ringmaster()
        step += 1

        log_progress(step, "Running Sieve on original input")
        self.run_sieve_1()
        step += 1

        log_progress(step, "Running BLASTN on unselected sequences")
        self.run_blastn(self.sieve_1_output_unSelected_reads_path, self.blast_results2_path)
        step += 1

        log_progress(step, "Creating index for unselected sequences")
        self.create_fasta_index(self.sieve_1_output_unSelected_reads_path)
        step += 1

        log_progress(step, "Running Juggler")
        self.run_juggler()
        step += 1

        log_progress(step, "Running Tamer")
        self.run_tamer()
        step += 1

        log_progress(step, "Running Sieve 2 to filter unclassified reads")
        self.run_sieve_2()
        step += 1

        log_progress(step, "Running final alignment with minimap2")
        self.run_minimap2()
        step += 1

        log_progress(step, "Sorting and indexing final alignment results")
        self.run_sort_and_index()
        step += 1

        log_progress(step, "Running Astrologer for final analysis and scoring")
        self.run_astrologer()
        step += 1

        log_progress(step, "Merging UeccDNA results")
        self.run_merge_uecc()
        step += 1

        log_progress(step, "Merging MeccDNA results")
        self.run_merge_mecc()
        step += 1

        log_progress(step, "Processing inferred UeccDNA")
        self.run_merge_uecc_inferred()
        step += 1

        log_progress(step, "Processing Cecc and MCecc")
        self.run_merge_cecc()
        step += 1

        log_progress(step, "Processing FASTA indices")
        self.run_FAIProcessor()
        step += 1

        log_progress(step, "Generating final analysis report")
        self.run_report_generator()
        step += 1

        # At the end of the run_pipeline method, replace the existing file saving logic with:
        log_progress(step, "Saving final results and cleaning up")
        
        # Define the files to be saved
        key_words = ["Confirmed", "Inferred", "eccDNA_analysis_report"]
        suffix = [".csv", ".html", ".fa", ".txt"]


        # Move files containing key words to the output directory
        saved_files = []

        for file in os.listdir(self.temp_dir):
            if any(key_word in file for key_word in key_words) and any(file.endswith(suffix) for suffix in suffix):
                src = os.path.join(self.temp_dir, file)
                dst = os.path.join(self.output_dir, file)
                shutil.move(src, dst)
                saved_files.append(file)
                print(f"Saved final result: {dst}")

        if not self.keep_tmp:
            # Clean up the temporary directory
            shutil.rmtree(self.temp_dir)
            print(f"Cleaned up temporary directory: {self.temp_dir}")
        else:
            print(f"Temporary files kept in: {self.temp_dir}")

        total_time = time.time() - start_time
        print(f"\n> Pipeline completed in {total_time:.2f} seconds")
        print(f"Final results saved in: {self.output_dir}")
        print("Saved files:")
        for file in saved_files:
            print(f" - {file}")

def main():
    """Parse command-line arguments and run the eccDNA analysis software package."""
    parser = argparse.ArgumentParser(description="Wrapper script for eccDNA analysis software package")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-p", "--prefix", required=True, help="Prefix for all generated files")
    parser.add_argument("-r", "--reference", required=True, help="Reference genome FASTA file")
    parser.add_argument("-o", "--output", default=None, help="Output folder (optional)")
    parser.add_argument("-t", "--num_threads", type=int, default=8, help="Number of threads to use (default: 8)")
    parser.add_argument("-kt", "--keep_tmp", action="store_true", help="Keep temporary files in the .tmp directory")
    args = parser.parse_args()

    # Modify this part to ensure the output_prefix includes the specified output directory
    if args.output:
        os.makedirs(args.output, exist_ok=True)
        output_prefix = os.path.join(args.output, args.prefix)
    else:
        output_prefix = args.prefix

    wrapper = CircleSeeker(args.input, output_prefix, args.reference, num_threads=args.num_threads, keep_tmp=args.keep_tmp)
    wrapper.run_pipeline()

if __name__ == "__main__":
    main()
