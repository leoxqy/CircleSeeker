CircleSeeker: Advanced Circular DNA Analysis Tool

CircleSeeker is a specialized bioinformatics pipeline designed for the identification, classification, and characterization of extrachromosomal circular DNA (eccDNA) using PacBio HiFi long-read sequencing data. By leveraging the high accuracy and long read lengths of PacBio HiFi technology, CircleSeeker streamlines the detection of various eccDNA types (UeccDNA, MeccDNA, CeccDNA, XeccDNA) and integrates a mandatory coverage-based analysis (Astrologer) to further refine and confirm eccDNA. CircleSeeker orchestrates multiple third-party tools (e.g., TideHunter, BLAST, minimap2, samtools, mosdepth, bcftools, etc.) within a single cohesive pipeline, minimizing user overhead for command orchestration and data management.

Key Features
	1.	Robust PacBio HiFi Data Processing
	•	Processes large volumes of long-read sequences with multi-threading for high efficiency.
	•	Specialized modules handle tandem repeats, coverage analyses, and multi-step classification.
	2.	Mandatory Coverage Analysis (Astrologer)
	•	CircleSeeker always runs Astrologer for coverage assessment, read-depth filtering, and local consensus generation using mosdepth and bcftools.
	•	This step refines detection of both confirmed and inferred eccDNA.
	3.	Classification of eccDNA Types
	•	Discovers and categorizes UeccDNA (Unique eccDNA), MeccDNA (Multiple-copy eccDNA), CeccDNA (Chimeric eccDNA), MCeccDNA, and XeccDNA based on alignment patterns and coverage information.
	4.	Modular Architecture
	•	Scripts like Juggler.py, Ringmaster.py, Carousel.py, Astrologer.py, Sieve.py, Tamer.py, and various Merge modules each handle distinct tasks in the pipeline.
	•	Easy to debug or adapt specific modules while preserving the full pipeline’s integrity.
	5.	Comprehensive Analysis Report
	•	Automatically generates HTML and text summaries, final CSV listings, and FASTA files of confirmed or inferred eccDNAs.
	•	Tracks read classifications, coverage profiles, duplication events, and more.
	6.	Checkpoint System
	•	CircleSeeker creates a checkpoint file once each major step completes, enabling the pipeline to resume from the most recent step in case of interruption (e.g., system failure).

Installation

CircleSeeker can be installed through two primary methods:

1. Installation via conda

Using conda ensures that CircleSeeker and its required dependencies (including coverage analysis tools) are installed correctly:

conda install -c bioconda circle_seeker

2. Installation from Source

If you prefer manual installation:

git clone https://github.com/leoxqy/CircleSeeker.git
cd CircleSeeker
pip install .

Important:
	•	Make sure minimap2, samtools, BLAST+, seqkit, mosdepth, bcftools, TideHunter, and other necessary binaries are in your $PATH.
	•	The mandatory coverage-based module (Astrologer) depends on mosdepth, bcftools, tabix, and minimap2.

Dependencies

CircleSeeker relies on:
	•	Python 3.6+
	•	numpy (>=1.22)
	•	pandas
	•	pysam
	•	biopython
	•	minimap2
	•	samtools
	•	BLAST+ (for makeblastdb and blastn)
	•	bedtools (optionally used for certain coverage tasks)
	•	seqkit
	•	TideHunter
	•	pbtk (if applicable in your environment)
	•	mosdepth, bcftools, tabix (all required for Astrologer)

When installed via conda, these tools and libraries typically become available automatically. If installing from source, please confirm they are accessible on your system.

Usage Overview

CircleSeeker executes a multi-step pipeline that always includes coverage analysis (Astrologer). A minimal command might look like:

CircleSeeker \
  -i <input_fasta> \
  -p <output_prefix> \
  -r <reference_genome> \
  [-o <output_folder>] \
  [-t <num_threads>] \
  [--enable_X] \
  [--force] \
  [--start_from <step_num>]

Required Arguments
	•	-i, --input
The input FASTA file containing PacBio HiFi reads.
	•	-p, --prefix
Prefix used in naming output files.
	•	-r, --reference
Reference genome FASTA file.

Optional Arguments
	•	-o, --output
Directory for storing results (defaults to current directory).
	•	-t, --num_threads
Number of threads (default: 8).
	•	--enable_X
Enable detection of XeccDNA (unclassified eccDNA).
	•	--force
Force re-run even if outputs already exist.
	•	--start_from <step_num>
Skip all steps before step_num, overriding any checkpoint.

For a deeper explanation of each internal step, parameter tuning, and example commands, see the CircleSeeker Tutorial.

Pipeline Modules

CircleSeeker subdivides tasks into specific scripts:
	•	Carousel.py: Preprocessing & circularizing raw tandem repeats from TideHunter
	•	Ringmaster.py: Classifying BLAST-aligned reads into Uecc, Mecc, Cecc, etc.
	•	Sieve.py: FASTA filtration based on read lists
	•	Juggler.py: BLAST result parsing, bridging partial alignments, marking single vs. multiple occurrence reads
	•	Tamer.py: Standardizing naming (e.g., eName) and extracting sequences from intermediate FASTA for classification
	•	Astrologer.py: Mandatory coverage-based analysis and local consensus generation
	•	MergeUecc/Mecc/Cecc/UeccInf: Deduplicating & merging final results into consolidated sets
	•	FAIProcessor.py: Validation of FASTA index and cleanup of empty files
	•	ReportGenerator.py: Summarizing final results into HTML/TXT/CSV outputs

Output Files

Once the pipeline completes, several final and intermediate files appear, typically prefixed by your -p argument. A typical output directory might include:
	•	<prefix>.checkpoint
	•	Stores an integer step index, allowing CircleSeeker to resume if interrupted.
	•	<prefix>.Final.Confirmed.Uecc.csv/.fasta
	•	Confirmed UeccDNA records and sequences. Columns commonly include:

eName,eChr,eStart,eEnd,eStrand,eLength,eRepeatNum,MatDegree,eClass,eState,eReadNum,eReads,eSeq


	•	<prefix>.Final.Confirmed.Mecc.csv/.fasta
	•	Confirmed MeccDNA (multiple-copy) data. Example CSV headers:

eName,MapNum,eRepeatNum,eState,eChr,eStart,eEnd,eLength,MatDegree,eClass,eReadNum,eReads,eSeq


	•	<prefix>.Final.Confirmed.Cecc.csv/.fasta
	•	Confirmed CeccDNA (chimeric) data. Example CSV headers:

eName,ChimeraNum,eLength,eRepeatNum,eClass,eState,LocNum,Identity,Alignment_length,Chr,Start,End,Strand,ReadsNum,eReads,eSeq


	•	<prefix>.Final.Confirmed.MCecc.csv/.fasta
	•	Confirmed MCecc (multiple-chimeric) data with slightly different columns.
	•	<prefix>.Final.Confirmed.Xecc.fasta
	•	If --enable_X was used, unclassified circular reads stored as XeccDNA (FASTA only).
	•	<prefix>.Final.Inferred.Uecc.csv/.fasta
	•	EccDNA identified through coverage-based inference (Astrologer). Example CSV columns:

eName,eChr,eStart,eEnd,eLength,Num_Reads,Avg_Read_Length,Num_Split_Reads,Avg_Split_Read_Length,Score,eSeq,eClass,eState


	•	<prefix>.Final.Merged.Uecc.csv
	•	Final merged UeccDNA combining confirmed and inferred sets, typically with columns:

eName,eChr,eStart,eEnd,eLength,Num_Reads,Avg_Read_Length,Num_Split_Reads,Avg_Split_Read_Length,Score,eSeq,eState


	•	<prefix>.report.html / <prefix>.report.txt
	•	Human-readable analysis summaries (HTML & text versions).
	•	<prefix>.summary.csv
	•	Summarized details for all eccDNA classes, typically columns like:

Class,Name,Length,Status



Depending on data scale, these CSV/FASTA files can become large. Intermediate files and logs reside in a .tmp folder (removed by default) unless --keep_tmp is specified.

Example Command and Workflow
	1.	Prepare your inputs:
	•	A FASTA file of reads (e.g., GBM_LB.fasta).
	•	A reference genome (e.g., hg38.fasta or T2T-CHM13v2.0_chr.fna).
	2.	Run:

CircleSeeker \
  -i GBM_LB.fasta \
  -p GBM_LB \
  -r T2T-CHM13v2.0_chr.fna \
  -t 16 \
  --enable_X


	3.	Check your results:
	•	GBM_LB.Final.Confirmed.Uecc.csv/.fasta
	•	GBM_LB.Final.Confirmed.Mecc.csv/.fasta
	•	GBM_LB.Final.Confirmed.Cecc.csv/.fasta
	•	GBM_LB.Final.Confirmed.MCecc.csv/.fasta
	•	GBM_LB.Final.Confirmed.Xecc.fasta
	•	GBM_LB.Final.Inferred.Uecc.csv/.fasta
	•	GBM_LB.Final.Merged.Uecc.csv
	•	GBM_LB.report.html
	•	GBM_LB.report.txt
	•	GBM_LB.summary.csv
	•	GBM_LB.checkpoint
	•	etc.

License

CircleSeeker is licensed under the GNU General Public License v2 (GPLv2). You are free to use, modify, and redistribute under these terms.

Authors and Contact
	•	Yaoxin Zhang (yxzhang@ncgr.com)
	•	Leo Xinqi Yu (leoxqy@hotmail.com)

Please open an issue on GitHub for inquiries, bug reports, or feature requests.

Citation

If you use CircleSeeker in your research, please cite the following (to be updated upon publication):

[Placeholder for citation information]

Contributing

We welcome improvements and new feature ideas. To contribute:
	1.	Fork this repository and create a new branch.
	2.	Make changes and thoroughly test them.
	3.	Submit a Pull Request.

See CONTRIBUTING.md for more details.

Acknowledgments

We would like to acknowledge the developers of the tools and libraries that CircleSeeker depends on, as well as the research community for their valuable suggestions and support.