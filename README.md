<p align="center">
  <img src="docs/images/logo.png" width="300" alt="CircleSeeker">
</p>

# CircleSeeker

[![Version](https://img.shields.io/badge/version-0.9.5-blue.svg)](https://github.com/yaoxinzhang/CircleSeeker)
[![Python](https://img.shields.io/badge/python-≥3.9-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-GPL--2.0-green.svg)](LICENSE)

Comprehensive detection and characterization of extrachromosomal circular DNA (eccDNA) from PacBio HiFi sequencing data.

## Installation

### Quick Install (Recommended)

```bash
# Clone repository
git clone https://github.com/yaoxinzhang/CircleSeeker.git
cd CircleSeeker

# Create and activate environment
conda env create -f environment.yml
conda activate circleseeker

# Install CircleSeeker
pip install -e .

# Verify installation
circleseeker --version
```

### Manual Installation

If you prefer to install dependencies manually:

```bash
# 1. Create environment with Python
conda create -n circleseeker python=3.10
conda activate circleseeker

# 2. Core Python packages (conda-forge)
conda install -c conda-forge -y \
  pip biopython pandas numpy tqdm pyyaml \
  click networkx intervaltree packaging \
  matplotlib python-graphviz pysam

# 3. Bioinformatics tools (bioconda)
conda install -c bioconda -c conda-forge -y \
  minimap2 samtools bcftools \
  tidehunter blast cd-hit \
  varlociraptor cyrcular \
  bedtools mappy

# 4. Install CircleSeeker
cd CircleSeeker
pip install -e .
```

### Install Cresil (Recommended Inference Engine)

Cresil provides faster inference than Cyrcular. Install from source:

```bash
git clone https://github.com/visanuwan/cresil.git
cd cresil
pip install -e .
cresil --help   # Verify installation
```

> **Note**: CircleSeeker will automatically detect available inference engines (Cresil or Cyrcular) and use the best available option.

## Quick Start

```bash
# Basic usage
circleseeker -i reads.fasta -r reference.fa -o output_dir

# With common options
circleseeker \
    -i sample.hifi.fasta \
    -r hg38.fa \
    -o results \
    -p sample_name \
    -t 16
```

## Usage

### Required Arguments

- `-i, --input PATH` - Input HiFi reads (FASTA format)
- `-r, --reference PATH` - Reference genome (FASTA format)

### Optional Arguments

- `-o, --output PATH` - Output directory (default: circleseeker_output)
- `-p, --prefix TEXT` - Output file prefix (default: sample)
- `-t, --threads INT` - Number of threads (default: 8)
- `-c, --config PATH` - Configuration file (YAML format)
- `-n, --noise` - Increase verbosity (-n for INFO, -nn for DEBUG)
- `-h, --help` - Show help message
- `--keep-tmp` - Keep temporary files

## Pipeline Overview

CircleSeeker implements a 16-step analysis pipeline:

### Detection Phase (Steps 1-6)
1. **make_blastdb** - Build BLAST database from reference
2. **tidehunter** - Detect tandem repeats in HiFi reads
3. **tandem_to_ring** - Convert tandem repeats to circular candidates
4. **run_blast** - Map candidates to reference genome
5. **um_classify** - Classify eccDNA as Unique (U) or Multiple (M) origin
6. **cecc_build** - Identify Complex (C) eccDNA

### Processing Phase (Steps 7-10)
7. **umc_process** - Process and cluster U/M/C types
8. **cd_hit** - Remove redundancy (90% identity threshold)
9. **ecc_dedup** - Deduplicate and standardize coordinates
10. **read_filter** - Filter confirmed eccDNA reads

### Inference Phase (Steps 11-13)
11. **minimap2** - Prepare reference index (generates BAM only when Cyrcular is used)
12. **ecc_inference** - Detect eccDNA (Cresil preferred, Cyrcular fallback)
13. **iecc_curator** - Curate inferred eccDNA

### Integration Phase (Steps 14-16)
14. **ecc_unify** - Merge confirmed and inferred results
15. **ecc_summary** - Generate statistics and summaries
16. **ecc_packager** - Package final outputs

## Output Files

```
output_dir/sample_name/
├── sample_name_Confirmed_UeccDNA/
│   ├── sample_name_UeccDNA_C.fasta      # Confirmed unique eccDNA sequences
│   ├── sample_name_UeccDNA.bed          # Genomic coordinates (BED format)
│   └── sample_name_UeccDNA.core.csv     # Detailed information table
├── sample_name_Confirmed_MeccDNA/
│   ├── sample_name_MeccDNA_C.fasta      # Confirmed multiple eccDNA sequences
│   ├── sample_name_MeccSites.bed        # All genomic sites
│   ├── sample_name_MeccBestSite.bed     # Best hit per eccDNA
│   └── sample_name_MeccSites.core.csv   # Detailed information table
├── sample_name_Confirmed_CeccDNA/
│   ├── sample_name_CeccDNA_C.fasta      # Confirmed complex eccDNA sequences
│   ├── sample_name_CeccSegments.bed     # All segments
│   ├── sample_name_CeccJunctions.bedpe  # Junction information
│   └── sample_name_CeccSegments.core.csv # Detailed information table
├── sample_name_Inferred_eccDNA/
│   ├── sample_name_UeccDNA_I.csv        # Inferred simple eccDNA
│   ├── sample_name_UeccDNA_I.fasta      # Inferred simple sequences
│   ├── sample_name_chimeric.csv         # Inferred chimeric eccDNA
│   └── sample_name_CeccDNA_I.fasta      # Inferred complex sequences
├── sample_name_merged_output.csv        # All eccDNA combined
├── sample_name_report.html              # Interactive HTML report
└── sample_name_summary.txt              # Summary statistics
```

### Output Format

The merged output file (`sample_name_merged_output.csv`) contains:

| Column | Description |
|--------|-------------|
| eccDNA_id | Unique identifier (e.g., UeccDNA1, MeccDNA1, CeccDNA1) |
| original | Original eccDNA ID from processing |
| Regions | Genomic coordinates (chr:start-end format) |
| Strand | DNA strand (+/-) |
| Length | eccDNA size in base pairs |
| eccDNA_type | Classification (UeccDNA/MeccDNA/CeccDNA) |
| State | Detection method (Confirmed/Inferred) |
| Seg_total | Number of segments (for complex eccDNA) |
| Hit_count | Number of genomic hits |

## Configuration

Create a custom configuration file in YAML format:

```yaml
# Basic settings
threads: 16

# Tool-specific parameters
tools:
  blast:
    word_size: 100
    evalue: "1e-50"
    perc_identity: 99.0
    soft_masking: false

  tidehunter:
    k: 16
    w: 1
    p: 100
    P: 2000000
    e: 0.1
    f: 2

  minimap2:
    preset: "map-hifi"
```

For a complete list of configuration options, see the [configuration reference](docs/CLI_Reference.md).

## System Requirements

- **Operating System**: Linux or macOS
- **Python**: ≥3.9
- **Memory**: Minimum 16GB (32GB recommended)
- **Storage**: ~50GB free space
- **CPU**: 8+ cores recommended

## Dependencies

### Python Packages
| Package | Version | Description |
|---------|---------|-------------|
| pandas | ≥2.0 | Data manipulation |
| numpy | ≥1.23 | Numerical computing |
| biopython | ≥1.81 | Sequence handling |
| pysam | ≥0.22 | BAM/SAM processing |
| networkx | ≥3.0 | Graph algorithms |
| intervaltree | ≥3.1 | Interval operations |
| click | ≥8.1 | CLI framework |
| pyyaml | ≥6.0 | Configuration parsing |

### External Tools (Required)
| Tool | Version | Description |
|------|---------|-------------|
| TideHunter | ≥1.4.0 | Tandem repeat detection |
| BLAST+ | ≥2.12.0 | Sequence alignment (makeblastdb, blastn) |
| minimap2 | ≥2.24 | Read-to-reference mapping |
| samtools | ≥1.17 | BAM file processing |
| bcftools | ≥1.17 | BCF/VCF processing |
| CD-HIT | ≥4.8.1 | Sequence clustering (cd-hit-est) |
| varlociraptor | ≥5.0.0 | Variant calling (for Cyrcular) |

### Inference Engines (At least one required)
| Tool | Description |
|------|-------------|
| Cresil | **Recommended** - Faster, simpler workflow |
| Cyrcular | Fallback option, requires varlociraptor |

## Troubleshooting

### Resume Interrupted Run

```bash
# Resume from last checkpoint
circleseeker -i reads.fa -r ref.fa -o previous_output_dir
```

### Check Installation

```bash
# Verify CircleSeeker version
circleseeker --version

# Test with minimal data
circleseeker -i test.fa -r ref.fa -o test_output
```

### Common Issues

**Out of Memory**: Reduce thread count or increase available memory
```bash
circleseeker -i input.fa -r ref.fa -t 4
```

**Missing Tools**: Install via conda
```bash
conda install -c bioconda -c conda-forge \
  tidehunter blast minimap2 samtools bcftools cd-hit varlociraptor
```

## Documentation

For detailed documentation, see the `docs/` directory:

- [Pipeline Modules](docs/Pipeline_Modules.md) - Detailed algorithm descriptions
- [CLI Reference](docs/CLI_Reference.md) - Complete command-line options

## Citation

If you use CircleSeeker in your research, please cite:

> Zhang Y, You M, Zhou C, et al. (2024). MMC-seq reveals a vast spatiotemporal eccDNA landscape from single cells to tissues. *Manuscript submitted*.

CircleSeeker is the computational pipeline developed as part of the MMC-seq methodology for comprehensive eccDNA detection and characterization

## License

GPL-2.0. See [LICENSE](LICENSE) for details.

## Contact

- Yaoxin Zhang: yxzhang@ncgr.ac.cn
- Leo Xinqi Yu: leoxqy@hotmail.com

## Acknowledgments

We thank the developers of TideHunter, Cyrcular, and other integrated tools.
