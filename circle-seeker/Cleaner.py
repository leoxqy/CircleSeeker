import os
import glob
import logging
from datetime import datetime

class Cleaner:
    def __init__(self, output_prefix):
        self.output_prefix = output_prefix
        self.log_filename = f"{output_prefix}_cleanup_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
        self.setup_logging()

    def setup_logging(self):
        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            handlers=[
                                logging.FileHandler(self.log_filename),
                                logging.StreamHandler()
                            ])

    def delete_files(self, file_patterns):
        for pattern in file_patterns:
            matching_files = glob.glob(pattern)
            for file in matching_files:
                if os.path.exists(file):
                    os.remove(file)
                    logging.info(f"Deleted: {file}")
                else:
                    logging.warning(f"File not found: {file}")

    def run(self):
        logging.info("Starting cleanup process for CircleSeeker intermediate files")

        # Intermediate file patterns
        intermediate_file_patterns = [
            f"{self.output_prefix}*blast1_out_alignment_results.txt",
            f"{self.output_prefix}*blast2_out_unSel.results.txt",
            f"{self.output_prefix}*Num_LinR.csv",
            f"{self.output_prefix}*read_classification.csv",
            f"{self.output_prefix}*TH.ecc_candidates.txt",
            f"{self.output_prefix}*carousel_circular_sequences.fasta",
            f"{self.output_prefix}*carousel_processed_results.csv",
            f"{self.output_prefix}*carousel_read_list.csv",
            f"{self.output_prefix}*sieve_1_output.txt",
            f"{self.output_prefix}*tecc_analysis_results.csv",
            f"{self.output_prefix}*uecc_part2.fa",
            f"{self.output_prefix}*mecc_part2.fa",
            f"{self.output_prefix}*Final_unclassified.fa",
            f"{self.output_prefix}*Final_unclassified.sam",
            f"{self.output_prefix}*Final_unclassified_sorted.bam",
            f"{self.output_prefix}*Final_unclassified_sorted.bam.bai",
            f"{self.output_prefix}*circular_dna_results_scored.csv",
            f"{self.output_prefix}*Uecc.csv",
            f"{self.output_prefix}*Mecc.csv",
            f"{self.output_prefix}*Cecc.csv",
            f"{self.output_prefix}*UeccDNA.fa",
            f"{self.output_prefix}*MeccDNA.fa",
            f"{self.output_prefix}*CeccDNA.fa",
            f"{self.output_prefix}*XeccDNA.fa",
            f"{self.output_prefix}*sieve_1_output.txt.fai",
            f"{self.output_prefix}*XeccDNA.fa.fai",
        ]

        # Retained file patterns
        retained_file_patterns = [
            f"{self.output_prefix}*_UeccDNA.Confirmed.csv",
            f"{self.output_prefix}*_UeccDNA.Confirmed.fa",
            f"{self.output_prefix}*_MeccDNA.Confirmed.csv",
            f"{self.output_prefix}*_MeccDNA.Confirmed.fa",
            f"{self.output_prefix}*_XeccDNA.Confirmed.csv",
            f"{self.output_prefix}*_XeccDNA.Confirmed.fa",
            f"{self.output_prefix}*_CeccDNA.Confirmed.csv",
            f"{self.output_prefix}*_CeccDNA.Confirmed.fa",
            f"{self.output_prefix}*_UeccDNA.Inferred.csv",
            f"{self.output_prefix}*_eccDNA_analysis_report.txt"
        ]

        logging.info("Deleting intermediate files...")
        self.delete_files(intermediate_file_patterns)

        logging.info("Cleanup complete. Checking for retained files:")
        for pattern in retained_file_patterns:
            matching_files = glob.glob(pattern)
            if matching_files:
                for file in matching_files:
                    logging.info(f"Retained: {file}")
            else:
                logging.warning(f"No files found matching pattern: {pattern}")

        logging.info("Cleanup process completed")
        logging.info(f"Log file saved as: {self.log_filename}")

if __name__ == "__main__":
    main()
