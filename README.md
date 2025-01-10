# CircleSeeker

> Advanced Circular DNA Analysis Tool for PacBio HiFi Sequencing Data

[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

## Overview

CircleSeeker is a specialized bioinformatics pipeline designed for the identification, classification, and characterization of extrachromosomal circular DNA (eccDNA) using PacBio HiFi long-read sequencing data. By leveraging the high accuracy and long read lengths of PacBio HiFi technology, CircleSeeker streamlines the detection of various eccDNA types (UeccDNA, MeccDNA, CeccDNA, XeccDNA) and integrates a mandatory coverage-based analysis (Astrologer) to further refine and confirm eccDNA. CircleSeeker orchestrates multiple third-party tools (e.g., TideHunter, BLAST, minimap2, samtools, mosdepth, bcftools, etc.) within a single cohesive pipeline, minimizing user overhead for command orchestration and data management.

## Key Features

### 1. Robust PacBio HiFi Data Processing
- Processes large volumes of long-read sequences with multi-threading for high efficiency
- Specialized modules handle tandem repeats, coverage analyses, and multi-step classification

### 2. Mandatory Coverage Analysis (Astrologer)
- CircleSeeker always runs Astrologer for coverage assessment, read-depth filtering, and local consensus generation using mosdepth and bcftools
- This step refines detection of both confirmed and inferred eccDNA

### 3. Classification of eccDNA Types
- Discovers and categorizes:
  - UeccDNA (Unique eccDNA)
  - MeccDNA (Multiple-copy eccDNA)
  - CeccDNA (Chimeric eccDNA)
  - MCeccDNA
  - XeccDNA
- Based on alignment patterns and coverage information

### 4. Modular Architecture
- Scripts like `Juggler.py`, `Ringmaster.py`, `Carousel.py`, `Astrologer.py`, `Sieve.py`, `Tamer.py`, and various Merge modules each handle distinct tasks
- Easy to debug or adapt specific modules while preserving the full pipeline's integrity

### 5. Comprehensive Analysis Report
- Automatically generates HTML and text summaries, final CSV listings, and FASTA files
- Tracks read classifications, coverage profiles, duplication events, and more

### 6. Checkpoint System
- Creates checkpoint files for each completed major step
- Enables pipeline resumption from the most recent step in case of interruption

## Installation

### Option 1: Via Conda (Recommended)
```bash
conda install -c bioconda circle_seeker
```

### Option 2: From Source
```bash
git clone https://github.com/leoxqy/CircleSeeker.git
cd CircleSeeker
pip install .
```

### Important Notes
- Ensure the following tools are in your `$PATH`:
  - minimap2
  - samtools
  - BLAST+
  - seqkit
  - mosdepth
  - bcftools
  - TideHunter

## Dependencies

### Required Software
- Python 3.6+
- numpy (>=1.22)
- pandas
- pysam
- biopython
- minimap2
- samtools
- BLAST+
- seqkit
- TideHunter
- mosdepth
- bcftools
- tabix
- tqdm

### Optional Tools
- bedtools (for certain coverage tasks)
- pbtk (environment-specific)

## Usage

### Basic Command
```bash
CircleSeeker \
  -i <input_fasta> \
  -p <output_prefix> \
  -r <reference_genome> \
  [-o <output_folder>] \
  [-t <num_threads>] \
  [--enable_X] \
  [--force] \
  [--start_from <step_num>]
```

### Required Arguments
| Argument | Description |
|----------|-------------|
| `-i, --input` | Input FASTA file containing PacBio HiFi reads |
| `-p, --prefix` | Prefix for output files |
| `-r, --reference` | Reference genome FASTA file |

### Optional Arguments
| Argument | Description | Default |
|----------|-------------|---------|
| `-o, --output` | Output directory | Current directory |
| `-t, --num_threads` | Number of threads | 8 |
| `--enable_X` | Enable XeccDNA detection | False |
| `--force` | Force re-run existing outputs | False |
| `--start_from` | Skip steps before step_num | None |

## Pipeline Modules

The pipeline consists of several specialized scripts:
- `Carousel.py`: Preprocessing & circularizing raw tandem repeats
- `Ringmaster.py`: Classifying BLAST-aligned reads
- `Sieve.py`: FASTA filtration
- `Juggler.py`: BLAST result parsing
- `Tamer.py`: Standardizing naming and sequence extraction
- `Astrologer.py`: Coverage-based analysis
- `MergeUecc/Mecc/Cecc/UeccInf`: Result deduplication
- `FAIProcessor.py`: FASTA index validation
- `ReportGenerator.py`: Results summarization

## Output Files

### Main Outputs
- `<prefix>.Final.Confirmed.{Uecc,Mecc,Cecc,MCecc}.{csv,fasta}`
- `<prefix>.Final.Confirmed.Xecc.fasta` (if --enable_X used)
- `<prefix>.Final.Inferred.Uecc.{csv,fasta}`
- `<prefix>.Final.Merged.Uecc.csv`
- `<prefix>.report.{html,txt}`
- `<prefix>.summary.csv`
- `<prefix>.checkpoint`

For detailed file formats and contents, please see our [documentation](https://github.com/leoxqy/CircleSeeker).

## Example Workflow

1. **Prepare Inputs**
   - FASTA file of reads (e.g., `GBM_LB.fasta`)
   - Reference genome (e.g., `T2T-CHM13v2.0_chr.fna`)

2. **Run CircleSeeker**
   ```bash
   CircleSeeker \
     -i GBM_LB.fasta \
     -p GBM_LB \
     -r T2T-CHM13v2.0_chr.fna \
     -t 16 \
     --enable_X
   ```

3. **Check Results**
   - Review output files in your specified directory
   - Examine the HTML report for analysis summary
   - Verify the checkpoint file for completion status

## License

This project is licensed under the GNU General Public License v2 (GPLv2) - see the [LICENSE](LICENSE) file for details.

## Authors

- Yaoxin Zhang (yxzhang@ncgr.com)
- Leo Xinqi Yu (leoxqy@hotmail.com)

## Support

Please [open an issue](https://github.com/leoxqy/CircleSeeker/issues) for support, bug reports, or feature requests.

## Citation

If you use CircleSeeker in your research, please cite the following (to be updated upon publication):

[Placeholder for citation information]

## Contributing

We welcome improvements and new feature ideas. To contribute:

1. Fork this repository and create a new branch
2. Make changes and thoroughly test them
3. Submit a Pull Request

See [CONTRIBUTING.md](CONTRIBUTING.md) for more details.

## Acknowledgments

We would like to acknowledge the developers of the tools and libraries that CircleSeeker depends on, as well as the research community for their valuable suggestions and support.