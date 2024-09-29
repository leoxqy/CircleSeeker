import argparse
import os
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from .eccDNA_Carousel import Carousel
from .eccDNA_Ringmaster import Ringmaster
from .eccDNA_Sieve import Sieve
import subprocess

class EccDNAWrapper:
    """Wrapper class for the eccDNA analysis pipeline."""

    def __init__(self, input_fasta, output_prefix):
        """Initialize the EccDNAWrapper with input and output parameters."""
        self.input_fasta = input_fasta
        self.output_prefix = output_prefix
        self.tidehunter_output = f"{output_prefix}.TH.cons.out.0.1.txt"
        self.carousel = Carousel()
        self.ringmaster = Ringmaster()
        self.sieve = Sieve()

    def run_tidehunter(self):
        """Execute TideHunter for consensus sequence generation."""
        command = f"TideHunter -f 2 -t 18 -k 16 -w 1 -p 100 -P 2000000 -e 0.1 {self.input_fasta} > {self.tidehunter_output}"
        subprocess.run(command, shell=True, check=True)

    def run_carousel(self):
        """Process TideHunter output and identify circular sequences."""
        carousel_output = f"{self.output_prefix}.processed_results.csv"
        read_list = f"{self.output_prefix}.read_list.csv"
        circular_sequences = "circular_sequences.fasta"
        self.carousel.run(self.tidehunter_output, carousel_output, read_list, circular_sequences)
        return circular_sequences

    def run_blastn(self, input_fasta, output_file):
        """Perform BLASTN search against hg38 reference."""
        blastn_cline = NcbiblastnCommandline(
            query=input_fasta,
            db="/dssg/home/acct-xyzltt/xyzltt/TGS/ref/hg38/hg38_db",
            out=output_file,
            outfmt=6,
            num_threads=18,
            word_size=100,
            evalue=1e-50,
            perc_identity=99
        )
        stdout, stderr = blastn_cline()

    def run_ringmaster(self, blast_results, circular_sequences):
        """Classify eccDNAs based on BLAST results."""
        mecc_fa = "custom_mecc.fa"
        cecc_fa = "custom_cecc.fa"
        xecc_fa = "custom_xecc.fa"
        self.ringmaster.run(blast_results, circular_sequences, mecc_fa, cecc_fa, xecc_fa)

    def run_sieve(self, input_fasta=None, read_list=None, output_file=None):
        """Extract unselected reads for further analysis."""
        if input_fasta is None:
            input_fasta = self.input_fasta
        if read_list is None:
            read_list = f"{self.output_prefix}.read_list.csv"
        if output_file is None:
            output_file = "unSel.fa"
        self.sieve.run(input_fasta, read_list, output_file)
        return output_file

    def create_fasta_index(self, fasta_file):
        """Create an index for the FASTA file using samtools."""
        subprocess.run(["samtools", "faidx", fasta_file], check=True)

    def run_juggler(self, blast_results, fasta_index):
        """Run eccDNA_Juggler.py for additional analysis."""
        output_file = f"{self.output_prefix}.tecc_analysis_results.csv"
        read_classification = f"{self.output_prefix}.read_classification.csv"
        num_linr = f"{self.output_prefix}.Num_LinR.csv"
        
        subprocess.run([
            "python", "eccDNA_Juggler.py",
            "-i", blast_results,
            "-f", fasta_index,
            "-n", num_linr,
            "-o", output_file,
            "-r", read_classification
        ], check=True)

    def run_pipeline(self):
        """Execute the complete eccDNA analysis pipeline."""
        # Step 1: Run TideHunter
        self.run_tidehunter()

        # Step 2: Run Carousel
        carousel_output, read_list, circular_sequences = self.run_carousel()

        # Step 3: Run BLASTN on circular sequences
        blast_results = self.run_blastn(circular_sequences, f"{self.output_prefix}.results_100_50_99.txt")

        # Step 4: Run Ringmaster
        self.run_ringmaster(blast_results, circular_sequences)

        # Step 5: Run Sieve on original input
        unselected_fasta = self.run_sieve(self.input_fasta, read_list, "unSel.fa")

        # Step 6: Run BLASTN on unselected sequences
        unsel_blast_results = self.run_blastn(unselected_fasta, f"{self.output_prefix}.unSel.results_99.txt")

        # Step 7: Create index for unselected sequences
        self.create_fasta_index(unselected_fasta)

        # Step 8: Run Juggler
        self.run_juggler(unsel_blast_results, f"{unselected_fasta}.fai")

        # Step 9: Run Sieve on unselected sequences
        self.run_sieve(unselected_fasta, "read_classification.csv", "End_unSel.fa")

def main():
    """Parse command-line arguments and run the eccDNA pipeline."""
    parser = argparse.ArgumentParser(description="Wrapper script for eccDNA analysis pipeline")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output_prefix", required=True, help="Output prefix for all generated files")
    args = parser.parse_args()

    wrapper = EccDNAWrapper(args.input, args.output_prefix)
    wrapper.run_pipeline()

if __name__ == "__main__":
    main()
