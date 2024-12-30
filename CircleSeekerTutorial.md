Below is a more in-depth CircleSeeker tutorial that dives into the algorithms driving each module, clarifies why each step is necessary, and explains various abbreviations and terminology used throughout the pipeline (such as Uecc, Mecc, Cecc, Xecc, MCecc, CtcR, LinR, FAI, etc.). This expanded tutorial should help you understand how CircleSeeker reaches its results and why each component is included.

Advanced CircleSeeker Tutorial

1. Overview and Terminology

1.1 What Is eccDNA?
	•	eccDNA stands for extrachromosomal circular DNA, a category of circular DNA molecules found outside the main chromosomal set. These can range from small “microDNAs” to large circular elements contributing to genomic plasticity.

1.2 CircleSeeker’s Mission

CircleSeeker aims to identify, classify, and characterize these eccDNA molecules from PacBio HiFi long-read sequencing data, leveraging:
	•	PacBio HiFi’s high accuracy (phred Q20+),
	•	long read lengths (often 10–25 kb or more),
	•	mandatory coverage-based analysis for confirmation.

1.3 Key Abbreviations
	1.	Uecc – Unique eccDNA: Typically one or few copies, straightforward alignment.
	2.	Mecc – Multiple-copy eccDNA: Repetitive or high-copy sequences forming eccDNA.
	3.	Cecc – Chimeric eccDNA: Derived from multiple distinct genomic loci fused into a single circular molecule.
	4.	MCecc – Multiple-chimeric eccDNA: More complex forms of Cecc with multiple chimeric junctions.
	5.	Xecc – Uncertain eccDNA (also “unclassified eccDNA”): Potential eccDNA sequences that don’t fit easily into other categories. Enabled by --enable_X.
	6.	CtcR – Concatemeric Tandem Copies Reads: Highly repetitive tandem read structures (detected by TideHunter).
	7.	LinR – Linear Reads: Non-circular reads or standard linear reads identified by Juggler or other modules.
	8.	FAI – FASTA Index (e.g., <genome>.fasta.fai), which stores the offsets/lengths for each sequence in a FASTA, used by samtools/mosdepth.
	9.	TideHunter – A specialized algorithm that scans for tandem repeats within long reads.
	10.	mosdepth – A tool for computing coverage information, essential for Astrologer’s coverage-based detection.
	11.	Z-threshold, MAD (Median Absolute Deviation) – Statistical measures used in robust detection of high-coverage or outlier regions in Astrologer.

2. Pipeline Algorithmic Flow

CircleSeeker’s pipeline can be viewed as three major phases:
	1.	Tandem Repeat Discovery (TideHunter & Carousel),
	2.	Alignment & Classification (BLAST, Ringmaster, Juggler, Tamer, Sieve),
	3.	Coverage-based Confirmation (Astrologer) + final merges/report.

Below is an algorithm-focused breakdown of each step.

2.1 Phase 1: Tandem Repeat Discovery
	1.	TideHunter
	•	Algorithm: TideHunter identifies tandem repeats by searching read substrings that self-align in repeating motifs.
	•	Output: A text file listing read coordinates for each tandem region.
	2.	Carousel
	•	Goal: Consolidate and “circularize” tandem repeat segments.
	•	Algorithmic Steps:
	•	Merges overlapping TideHunter hits per read to avoid fragmentation.
	•	Generates 2× read sequences for suspect circular molecules (doubling the read can help detect wrap-around alignments).
	•	Filters out short or poorly matched repeats (e.g., average match < 99%).
	•	Reasoning: Real circular molecules often appear as tandem or near-tandem in long reads. Doubling the sequence helps find a “split boundary” typical of circles.

2.2 Phase 2: Alignment & Classification
	1.	makeblastdb + BLAST
	•	Algorithmic Steps:
	•	Creates a BLAST database of the reference genome.
	•	Aligns the circular candidates (from Carousel) to the reference using blastn.
	•	Metrics: Identity, alignment length, start/end positions, E-value, bit score are collected.
	2.	Ringmaster
	•	Purpose: Categorize each query read’s alignment pattern into:
	•	Uecc if it aligns uniquely or has a single best alignment.
	•	Mecc if multiple strong alignment blocks indicate repeated or multi-copy segments in the genome.
	•	Cecc (Chimeric) if alignments from distinct genomic regions form a single continuous circle.
	•	Xecc if uncertain or unclassified (only if --enable_X is set).
	•	Algorithmic Details:
	•	Sums alignment blocks, checks for partial overlaps and consistent orientation.
	•	Evaluates “gap percentages” or identity gaps to decide if multiple blocks represent repeats (Mecc) or distinct breakpoints (Cecc).
	3.	Sieve (First Pass)
	•	Purpose: Retain or discard reads based on prior classification lists.
	•	Algorithm:
	•	Uses seqkit grep to separate “selected” vs. “unselected” reads.
	•	Why? Because we only want to keep reads that might contain hidden circular structures (or discarding confidently linear reads).
	4.	BLAST (Second Pass)
	•	Align the “unselected” portion once again to ensure we haven’t missed any eccDNA signals the first time.
	5.	Juggler
	•	Algorithmic Steps:
	•	Focuses on multi-occurrence alignments to see if partial fragments chain into a circular arrangement.
	•	Creates alignment “chains” to reconstruct a circle and checks for boundary overlap, repeated start/end patterns, etc.
	6.	Tamer
	•	Purpose: Standardize fields like eName, eSeq, eClass, unify partial data.
	•	Consolidates the partial classification from Juggler (and earlier steps) into a single CSV, including newly recognized circular segments.
	7.	Sieve (Second Pass)
	•	Another round to remove finalized reads from the pool, preparing for coverage-based analysis on the remainder.

2.3 Phase 3: Coverage-Based Confirmation (Astrologer) + Merging
	1.	minimap2
	•	Aligns leftover reads to the reference to generate a new BAM with potentially uncovered circular structures.
	2.	Astrologer
	•	Mandatory coverage-based module.
	•	Algorithmic Steps:
	1.	mosdepth computes per-base coverage.
	2.	Coverage data is read, and “robust outlier detection” methods (e.g., MAD + median, Z-threshold) identify abnormally high-coverage intervals that may indicate circular amplification.
	3.	Local assembly or consensus is generated with bcftools if enough reads overlap a region.
	4.	Newly discovered circles are labeled as inferred eccDNA (sometimes referred to as Inferred Uecc).
	•	Why? Some eccDNAs do not present clear “tandem repeat” signatures but exhibit very high coverage or split-read patterns, only discoverable via coverage stats.
	3.	MergeUecc/Mecc/Cecc/UeccInf
	•	Algorithmic Steps:
	•	Gathers “confirmed” and “inferred” sets, identifies duplicate or overlapping calls.
	•	Resolves naming collisions, merges short/partial circles into a final, deduplicated output.
	•	Example merges:
	•	Uecc + Inferred Uecc => final Uecc set
	•	Mecc => final multiple-copy set
	•	Cecc + MCecc => final chimeric sets
	•	Typically merges are done using maximum coverage, minimal mismatch thresholds, or the best average read identity.
	4.	FAIProcessor
	•	Goal: Validate all final FASTA outputs (Uecc, Mecc, Cecc, etc.), checking for empties and creating .fai indexes.
	5.	ReportGenerator
	•	Generates comprehensive HTML/TXT/CSV summaries:
	•	Summaries of how many reads are CtcR vs. LinR.
	•	Class counts (Uecc, Mecc, Cecc, MCecc, Xecc if enabled).
	•	Overlap analysis for inferred vs. confirmed Uecc.
	•	EccDNA length distribution, mean coverage, other stats.

3. Detailed Data Flow with Example Commands

Let’s assume we have:
	•	myPacBioReads.fasta: The input long-read data (PacBio HiFi).
	•	hg38.fasta: The reference genome (human GRCh38).

	1.	Create Conda Environment (Recommended)

conda create -n circle_env python=3.8
conda activate circle_env
conda install -c bioconda circle_seeker


	2.	Run CircleSeeker

CircleSeeker \
  -i myPacBioReads.fasta \
  -p MyAnalysis \
  -r hg38.fasta \
  -t 16 \
  --enable_X


	3.	Internal Steps:
	1.	makeblastdb -in hg38.fasta ...
	2.	TideHunter -f 2 ... myPacBioReads.fasta > .tmp/myReads.TH.ecc_candidates.txt
	3.	Carousel.py merges tandem hits.
	4.	blastn -db hg38 ... .tmp/myReads.carousel_circular_sequences.fasta
	5.	Ringmaster.py => “UeccDNA.csv,” “MeccDNA.csv,” etc.
	6.	Sieve.py => filter leftover reads => “unselected.fa”
	7.	Second blastn => “unselected.blast.out”
	8.	Juggler.py => refine classification
	9.	Tamer.py => unify CSV fields
	10.	Sieve.py => second pass
	11.	minimap2 -a hg38.fasta .tmp/unselected_sieve2.fa > unselected_sieve2.sam
	12.	Astrologer.py => coverage peaks => local consensus => “Final.Inferred.Uecc.csv/fasta”
	13.	MergeUecc.py, MergeMecc.py => unify partial results
	14.	FAIProcessor.py => validate final FASTA indexes
	15.	ReportGenerator.py => “MyAnalysis.report.html/txt,” “MyAnalysis.summary.csv”

4. Interpreting the Results

After the pipeline finishes:
	1.	MyAnalysis.Final.Confirmed.Uecc.csv/fasta
	•	Algorithmically, these are read segments that:
	•	Have strong single alignments.
	•	Appear to form a circle with no conflicting coverage signals.
	•	Not repeated extensively in the genome.
	2.	MyAnalysis.Final.Confirmed.Mecc.csv/fasta
	•	Multiple-copy eccDNAs (high duplication or multi-mapping patterns).
	•	Possibly repeated tandem arrays or high copy plasmid-like structures.
	3.	MyAnalysis.Final.Confirmed.Cecc.csv/fasta
	•	Chimeric circles combining separate loci. May show breakpoints from different chromosomes or distant regions.
	•	Often partial alignment blocks bridging multiple references.
	4.	MyAnalysis.Final.Inferred.Uecc.csv/fasta
	•	Coverage-based detection. Possibly missed by BLAST if no strong tandem repeat signals.
	•	High or abnormally distributed coverage flagged by Astrologer => local assembly => new circle ID.
	5.	MyAnalysis.Final.Merged.Uecc.csv
	•	Consolidation of confirmed + inferred Uecc.
	•	Avoids duplication via threshold-based merges on alignment or coverage evidence.
	6.	MyAnalysis.Final.Confirmed.Xecc.fasta
	•	“Uncertain/Unclassified eccDNA” if --enable_X was used. Typically partial or incomplete evidence.
	7.	MyAnalysis.report.html/txt + .summary.csv
	•	Summaries of read classification:
	•	CtcR (Concatemeric Tandem Copies Reads) found vs. LinR (Linear Reads) determined.
	•	Count of each eccDNA class.
	•	Overlap analysis: how many Inferred Uecc overlap with Confirmed Uecc.
	•	These results help you see how many circles are “unique,” “multiple,” or “chimeric,” how big they are, and coverage stats.

5. Common Questions
	1.	What if no circles are detected?
	•	Possibly your data lacks true eccDNA or coverage is too low. Check intermediate CSVs for partial alignments or run with a relaxed threshold (--force or lower coverage filters in Astrologer).
	2.	Why do we double the sequence in Carousel?
	•	Doubling helps detect circular boundary overlaps (start->end continuity). If a read is truly circular, aligning its end to its start is crucial.
	3.	Why is Astrologer mandatory?
	•	Some eccDNA appear normal under BLAST but show drastically increased coverage or split-read signals. Requiring coverage analysis ensures you don’t miss those.
	•	E.g., large amplified circular regions typical in certain cancers.
	4.	What about runtime & memory?
	•	BLAST can be memory-intensive for large references. Consider HPC usage or high-memory servers.
	•	Multi-threading improves speed but also increases memory usage.
	5.	Can I skip steps manually?
	•	We strongly recommend letting CircleSeeker handle all steps because coverage-based detection (Astrologer) is integral.
	•	You can override with --start_from <step_num> in advanced debugging, but do so cautiously.

6. Tips for Algorithmic Fine-Tuning
	1.	TideHunter Parameters
	•	Default settings usually suffice, but if your reads have unusual tandem repeats, you might tweak -k (k-mer size) or -p (minimum partial match).
	2.	BLAST
	•	Adjust -word_size or -perc_identity for large or small expected eccDNAs. Larger -word_size speeds alignment but is less sensitive.
	3.	Ringmaster
	•	Check the “gap_percentage” threshold if you see ambiguous calls. A smaller threshold means more strict classification.
	4.	mosdepth (Astrologer)
	•	A larger coverage threshold (e.g., ignoring coverage < 5) might speed up runs but may miss low-frequency circles.
	•	The robust MAD approach is less reliant on user-defined cutoffs, but you can tweak the z_threshold.
	5.	Merging Scripts (MergeUecc, MergeMecc, etc.)
	•	Weights coverage, identity, alignment length to unify or discard overlapping circles.
	•	If you see excessive merging, consider adjusting the tolerance for overlap or coordinate proximity in the scripts.

7. Extending or Combining with Other Tools
	•	IGV or UCSC Genome Browser
	•	Convert CircleSeeker output into BED tracks to visualize circle locations in a genome browser.
	•	Structural Variant Tools
	•	Use dedicated SV callers for cross-verification of breakpoints.
	•	Annotation
	•	If you want to see which genes are encompassed by these circles, integrate with gff/gtf annotation.

8. Conclusion

By dissecting how each sub-module (TideHunter → Carousel → BLAST → Ringmaster → Juggler → Tamer → Astrologer → Merge/Report) interacts, you can fine-tune CircleSeeker for different organisms, coverage depths, or specialized circular DNA analyses. Understanding why a read is labeled Uecc vs. Mecc or how coverage evidence elevates certain circles from “unclear” to “inferred confirmed” can vastly improve the reliability of your eccDNA discoveries.

Core Takeaways:
	•	CircleSeeker combines tandem repeat detection + reference alignment + coverage-based local assembly to capture a broad range of eccDNA forms (unique, multi-copy, chimeric, uncertain).
	•	All major steps rely on well-established algorithms (TideHunter for tandem repeats, BLAST for alignment, robust outlier detection with MAD in coverage analysis, etc.).
	•	Understanding the meanings of abbreviations (Uecc, Mecc, Cecc, MCecc, Xecc, CtcR, LinR, FAI, etc.) helps interpret final CSV and FASTA outputs.

9. Further Reading & Support
	1.	Algorithmic References
	•	TideHunter GitHub for the official tandem repeat discovery approach.
	•	Mosdepth Publication for coverage calculation.
	•	BLAST Documentation for alignment parameters.
	2.	Community
	•	Raise issues or ask questions on CircleSeeker’s GitHub Issues.
	•	We also welcome code contributions and new feature suggestions.
	3.	Citations
	•	If CircleSeeker aids your research, please cite it accordingly (info in the README).

With that, you now have a deeper, algorithm-focused tutorial for CircleSeeker that clarifies the entire process and the significance of each classification type. Good luck, and feel free to ask the community for further guidance or advanced parameter tuning!