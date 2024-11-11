import pysam
import csv
import argparse
import subprocess
import tempfile
import os
import logging

class Astrologer:
    def __init__(self, input_bam, output_csv):
        self.input_bam = input_bam
        self.output_csv = output_csv
        # 设置日志格式
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

    def is_soft_clipped(self, read):
        """
        Determine if the read contains soft-clipping CIGAR operations
        """
        return any(cig_op == 4 for cig_op, length in read.cigartuples)

    def check_circular_dna(self, read, region_start, region_end):
        """
        Check if the read is soft-clipped at the 5' or 3' end of the region,
        and if the other part is at the opposite end of the corresponding region
        """
        if not self.is_soft_clipped(read):
            return False

        read_start = read.reference_start
        read_end = read.reference_end

        if (read_start <= region_start + 200 and read_end >= region_end - 200):
            return True

        return False

    def calculate_read_length(self, read):
        """
        Calculate the length of the read based on the CIGAR field, including soft-clipped parts
        """
        length = 0
        if read.cigartuples:
            for cig_op, cig_length in read.cigartuples:
                # Operation codes 0 (M), 1 (I), 4 (S), 7 (=), 8 (X) consume query sequence length
                if cig_op in [0, 1, 4, 7, 8]:
                    length += cig_length
        return length

    def calculate_score(self, coverage, num_reads, avg_read_length, num_split_reads, avg_split_read_length):
        """
        Calculate the score based on coverage, number of reads, average read length,
        number of split reads, and average split read length
        """
        score = (coverage * 0.4) + (num_reads * 0.2) + (avg_read_length * 0.1) + (num_split_reads * 0.2) + (avg_split_read_length * 0.1)
        return round(score, 2)

    def run_bedtools_commands(self):
        """
        Run bedtools commands and return the results
        """
        logging.info(f"Processing BAM file: {self.input_bam}")
        # Create a temporary file
        with tempfile.NamedTemporaryFile(mode='w+t', delete=False) as temp_file:
            temp_filename = temp_file.name

            # Run bedtools commands
            cmd1 = f"bedtools genomecov -ibam {self.input_bam} -bg"
            cmd2 = "bedtools merge -i - -d 50 -c 4 -o sum"
            cmd3 = "awk '$4 > 2'"

            p1 = subprocess.Popen(cmd1.split(), stdout=subprocess.PIPE)
            p2 = subprocess.Popen(cmd2.split(), stdin=p1.stdout, stdout=subprocess.PIPE)
            p3 = subprocess.Popen(cmd3, shell=True, stdin=p2.stdout, stdout=temp_file, stderr=subprocess.PIPE)

            p3.communicate()

        return temp_filename

    def process_data(self):
        # Run bedtools commands
        bed_file = self.run_bedtools_commands()
        logging.info("Starting to analyze regions")

        # Add counters
        regions_processed = 0
        regions_with_splits = 0

        # Open BAM file
        bamfile = pysam.AlignmentFile(self.input_bam, "rb")

        # Create output CSV file
        with open(self.output_csv, "w", newline='') as csvfile:
            fieldnames = ['Chromosome', 'Start', 'End', 'Coverage', 'Num_Reads', 'Avg_Read_Length', 'Num_Split_Reads', 'Avg_Split_Read_Length', 'Score']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()

            # Read BED file and process each region
            with open(bed_file, "r") as bed:
                for line in bed:
                    regions_processed += 1  # Increment counter
                    chrom, start, end, coverage = line.strip().split()
                    start, end = int(start), int(end)
                    coverage = int(float(coverage))

                    reads_in_region = list(bamfile.fetch(chrom, start, end))
                    num_reads = len(reads_in_region)

                    if num_reads > 0:
                        # Calculate total read length
                        total_read_length = sum(self.calculate_read_length(read) for read in reads_in_region)
                        avg_read_length = int(round(total_read_length / num_reads)) if total_read_length > 0 else 0

                        # Calculate split reads
                        split_reads = [read for read in reads_in_region if self.is_soft_clipped(read)]
                        num_split_reads = len(split_reads)

                        if num_split_reads > 0:
                            total_split_read_length = sum(self.calculate_read_length(read) for read in split_reads)
                            avg_split_read_length = int(round(total_split_read_length / num_split_reads)) if total_split_read_length > 0 else 0
                        else:
                            avg_split_read_length = 0

                        # Calculate score
                        score = round(self.calculate_score(coverage, num_reads, avg_read_length, num_split_reads, avg_split_read_length), 2)
                    else:
                        avg_read_length = 0
                        num_split_reads = 0
                        avg_split_read_length = 0
                        score = 0.00

                    # Only write to CSV file when Num_Split_Reads > 1
                    if num_split_reads > 1:
                        regions_with_splits += 1  # Increment counter
                        writer.writerow({
                            'Chromosome': chrom,
                            'Start': start,
                            'End': end,
                            'Coverage': int(coverage),
                            'Num_Reads': int(num_reads),
                            'Avg_Read_Length': int(avg_read_length),
                            'Num_Split_Reads': int(num_split_reads),
                            'Avg_Split_Read_Length': int(avg_split_read_length),
                            'Score': f"{score:.2f}"
                        })

        logging.info(f"Analysis complete. Processed {regions_processed} regions, found {regions_with_splits} regions with split reads")

        # Close BAM file
        bamfile.close()

        # Delete temporary BED file
        os.unlink(bed_file)

def main():
    parser = argparse.ArgumentParser(description="Analyze circular DNA from BAM file")
    parser.add_argument("-i", "--input", required=True, help="Input BAM file")
    parser.add_argument("-o", "--output", required=True, help="Output CSV file")

    args = parser.parse_args()
    logging.info("Starting Astrologer analysis")
    
    astrologer = Astrologer(args.input, args.output)
    astrologer.process_data()
    
    logging.info("Analysis finished successfully")

if __name__ == "__main__":
    main()