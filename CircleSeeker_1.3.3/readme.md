# CircleSeeker v1.3.3

CircleSeeker is a specialized bioinformatics pipeline designed for detecting and analyzing extrachromosomal circular DNA (eccDNA) from PacBio HiFi sequencing data. It employs a sophisticated multi-step approach to identify and characterize different types of eccDNA structures with high accuracy and efficiency, leveraging the high accuracy and long read lengths of PacBio HiFi technology.

[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

## Features

- **Comprehensive eccDNA Detection**: Identifies multiple types of eccDNA:
  - UeccDNA (Unique eccDNA)
  - MeccDNA (Multiple-copy eccDNA)
  - CeccDNA (Composite eccDNA)
  - MCeccDNA/Ceccm (Multiple-copy Composite eccDNA)
  - Inferred UeccDNA

- **Advanced Analysis Pipeline**:
  - TideHunter integration for circular sequence detection
  - BLAST-based sequence analysis
  - Minimap2 integration for accurate mapping
  - Multi-threaded processing for improved performance
  - Comprehensive reporting and visualization

- **Robust Architecture**:
  - Modular component design
  - Checkpoint system for process recovery
  - Memory-efficient processing
  - Detailed logging system

## Installation

### Prerequisites

- Python 3.6 or higher
- numpy >=1.22
- pandas
- pysam
- Biopython
- minimap2
- samtools
- BLAST+
- bedtools
- seqkit
- TideHunter
- Biopython
- pandas
- numpy
- tqdm

### Using Conda (Recommended)

```bash
conda install -c bioconda circle_seeker
```

### Manual Installation

```bash
git clone https://github.com/leoxqy/CircleSeeker.git
cd CircleSeeker
pip install .
```

## Usage

### Basic Usage

```bash
python /path/to/CircleSeeker.py -i input.fasta -p prefix -r reference.fa -o output_dir -t num_threads
```

### Example Command

```bash
python /path/to/CircleSeeker.py \
    -i input_data.fasta \
    -p output_prefix \
    -r reference_genome.fa \
    -o output_directory \
    -t 8 \
    --enable_X \
    -kt &
```

### Required Arguments

- `-i, --input`: Input HiFi sequencing data in FASTA format
- `-p, --prefix`: Prefix for output files
- `-r, --reference`: Reference genome in FASTA format

### Optional Arguments

- `-o, --output`: Output directory (default: current directory)
- `-t, --threads`: Number of threads (default: 8)

### Developer Parameters

- `-kt, --keep-tmp`: Keep temporary files
- `--enable_X`: Enable extended analysis features

## Output Files

CircleSeeker generates several output files:

- `{prefix}_uecc.csv`: Unique eccDNA results
- `{prefix}_mecc.csv`: Multiple-copy eccDNA results
- `{prefix}_cecc.csv`: Composite eccDNA results
- `{prefix}_report.html`: Comprehensive analysis report
- `{prefix}_statistics.txt`: Summary statistics

## Pipeline Components

CircleSeeker executes the following steps in sequence:

- **BLAST DB**: Initial database creation for sequence alignment
- **TideHunter**: Circular sequence detection
- **Carousel**: Initial sequence processing and pattern analysis
- **BLASTN**: Alignment of circular sequences
- **Ringmaster**: Classification of UeccDNA, MeccDNA, and CeccDNA
- **Sieve**: FASTA read filtering and processing
- **BLASTN**: Analysis of unselected sequences
- **Indexing**: SAMtools-based sequence indexing
- **Juggler**: Sequence analysis and classification
- **Tamer**: Enhanced data validation
- **Sieve**: Secondary filtering of unclassified reads
- **Minimap2**: Final sequence alignment
- **SAMtools**: Sorting and indexing alignments
- **Astrologer**: Final analysis and confidence scoring
- **Merger**: Consolidation of UeccDNA, MeccDNA, CeccDNA results
- **FAIProcessor**: FASTA index management
- **ReportGenerator**: Comprehensive report generation
- **Cleanup**: Final results saving and temporary file management

## Authors

- Yaoxin Zhang (yxzhang@ncgr.com)
- Leo Xinqi Yu (leoxqy@hotmail.com)

## Contributing

We welcome contributions! Please refer to our [contribution guidelines](CONTRIBUTING.md) for more information.

## Citation

If you use CircleSeeker in your research, please cite:

```
[Citation information to be added]
```

## License

This project is licensed under the GNU General Public License v2 - see the LICENSE file for details.

## Support

For issues, bug reports, or feature requests, please use the [GitHub issue tracker](https://github.com/leoxqy/CircleSeeker/issues).

## Acknowledgments

We thank the developers of the tools and libraries that CircleSeeker depends on, as well as the research community for their valuable feedback and support.
