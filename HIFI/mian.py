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

    def run_sieve(self):
        """Extract unselected reads for further analysis."""
        read_list = f"{self.output_prefix}.read_list.csv"
        unsel_fa = "unSel.fa"
        self.sieve.run(self.input_fasta, read_list, unsel_fa)
        return unsel_fa

    def run_pipeline(self):
        """Execute the complete eccDNA analysis pipeline."""
        # Steps 1-6: Run individual components of the pipeline
        # Step 1: Run TideHunter
        self.run_tidehunter()

        # Step 2: Run Carousel
        circular_sequences = self.run_carousel()

        # Step 3: Run BLASTN on circular sequences
        blast_results = f"{self.output_prefix}.results_100_50_99.txt"
        self.run_blastn(circular_sequences, blast_results)

        # Step 4: Run Ringmaster
        self.run_ringmaster(blast_results, circular_sequences)

        # Step 5: Run Sieve
        unsel_fa = self.run_sieve()

        # Step 6: Run BLASTN on unSel.fa
        unsel_blast_results = f"unSel.{self.output_prefix}.results_99.txt"
        self.run_blastn(unsel_fa, unsel_blast_results)

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
