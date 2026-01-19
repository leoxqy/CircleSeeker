<p align="center">
  <img src="docs/images/logo.png" width="300" alt="CircleSeeker">
</p>

# CircleSeeker

[![Version](https://img.shields.io/badge/version-1.0.0-blue.svg)](https://github.com/leoxqy/CircleSeeker)
[![CI](https://github.com/leoxqy/CircleSeeker/actions/workflows/ci.yml/badge.svg)](https://github.com/leoxqy/CircleSeeker/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-≥3.9-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-GPL--3.0-green.svg)](LICENSE)

Comprehensive detection and characterization of extrachromosomal circular DNA (eccDNA) from PacBio HiFi sequencing data.

## eccDNA Classification

CircleSeeker identifies and classifies eccDNA into three categories:

| Type | Full Name | Description |
|------|-----------|-------------|
| **UeccDNA** | Unique eccDNA | Derived from a single genomic locus |
| **MeccDNA** | Multiple eccDNA | Derived from multiple distinct genomic loci |
| **CeccDNA** | Chimeric eccDNA | Formed by fusion of segments from different genomic locations |

## Installation

### Quick Install (Recommended)

```bash
# Clone repository
git clone https://github.com/leoxqy/CircleSeeker.git
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
  tidehunter cd-hit last \
  varlociraptor cyrcular \
  bedtools mappy

# 4. Install CircleSeeker
cd CircleSeeker
pip install -e .
```

### Install Cresil (Recommended Inference Engine)

Cresil provides **higher accuracy** for eccDNA inference compared to Cyrcular. Install from source:

```bash
git clone https://github.com/visanuwan/cresil.git
cd cresil
pip install -e .
cresil --help   # Verify installation
```

> **Note**: CircleSeeker automatically detects available inference engines (Cresil or Cyrcular) and uses the best option. Cresil is preferred for its higher precision.

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
- `-v, --verbose` - Increase verbosity (-v for INFO, -vv for DEBUG)
- `-h, --help` - Show help message
- `--keep-tmp` - Keep temporary files

## Pipeline Overview

CircleSeeker implements a 16-step analysis pipeline with two evidence-driven callers:

- **CtcReads**: Reads containing **Ctc** (**C**oncatemeric **t**andem **c**opies) signals
- **CtcReads-Caller** (Steps 1–10): Produces **Confirmed** U/M/C eccDNA from CtcReads evidence
- **SplitReads-Caller** (Steps 11–13): Produces **Inferred** eccDNA from split-read/junction evidence

### Architecture

```
                              ┌─────────────────────┐
                              │  Input HiFi FASTA   │
                              └──────────┬──────────┘
                                         │
                                         ▼
╔═════════════════════════════════════════════════════════════════════════════╗
║                        CtcReads-Caller (Steps 1-10)                         ║
╠═════════════════════════════════════════════════════════════════════════════╣
║                                                                             ║
║       ┌────────────┐      ┌────────────────┐      ┌────────────┐            ║
║       │ TideHunter │─────▶│ tandem_to_ring │─────▶│  minimap2  │            ║
║       │(Ctc detect)│      │    (double)    │      │(alignment) │            ║
║       └────────────┘      └────────────────┘      └──────┬─────┘            ║
║                                                          │                  ║
║                                                          ▼                  ║
║                                            ┌──────────────────────┐         ║
║                                            │     um_classify      │         ║
║                                            │ (Uecc / Mecc / Cecc) │         ║
║                                            └──────────┬───────────┘         ║
║                                                       │                     ║
║                ┌──────────────────────────────────────┼──────────────────┐  ║
║                │                                      │                  │  ║
║                ▼                                      ▼                  ▼  ║
║        ┌──────────────┐                      ┌──────────────┐    ┌────────┐ ║
║        │     Uecc     │                      │     Mecc     │    │  LAST  │ ║
║        │   (unique)   │                      │  (multiple)  │    │ (Cecc) │ ║
║        └──────┬───────┘                      └──────┬───────┘    └───┬────┘ ║
║               │                                     │                │      ║
║               └─────────────────────────────────────┴────────────────┘      ║
║                                                     │                       ║
║                                                     ▼                       ║
║                                      ┌────────────────────────┐             ║
║                                      │ CD-HIT → dedup → filter│             ║
║                                      └────────────┬───────────┘             ║
║                                                   │                         ║
║                                                   ▼                         ║
║                                     ┌────────────────────────┐              ║
║                                     │ Confirmed eccDNA (UMC) │              ║
║                                     └────────────┬───────────┘              ║
╚══════════════════════════════════════════════════╪══════════════════════════╝
                                                   │
╔══════════════════════════════════════════════════╪══════════════════════════╗
║                                                  │                          ║
║                                                  ▼                          ║
║                               ┌────────────────────────────┐                ║
║                               │ Filtered reads (non-CtcR)  │                ║
║                               └──────────────┬─────────────┘                ║
║                                              │                              ║
║                       SplitReads-Caller (Steps 11-13)                       ║
╠═════════════════════════════════════════════════════════════════════════════╣
║                                              │                              ║
║                                              ▼                              ║
║                 ┌────────────────────────────────────────────────┐          ║
║                 │      Cresil (preferred, higher accuracy)       │          ║
║                 │      Cyrcular (fallback)                       │          ║
║                 │           Split-read junction detection        │          ║
║                 └────────────────────────┬───────────────────────┘          ║
║                                          │                                  ║
║                                          ▼                                  ║
║                               ┌────────────────────┐                        ║
║                               │   Inferred eccDNA  │                        ║
║                               └─────────┬──────────┘                        ║
╚═════════════════════════════════════════╪═══════════════════════════════════╝
                                          │
╔═════════════════════════════════════════╪═══════════════════════════════════╗
║                                         │                                   ║
║                        Integration (Steps 14-16)                            ║
╠═════════════════════════════════════════════════════════════════════════════╣
║                                         │                                   ║
║                                         ▼                                   ║
║          ┌───────────┐       ┌─────────────┐       ┌─────────────┐          ║
║          │ ecc_unify │──────▶│ ecc_summary │──────▶│ ecc_packager│          ║
║          │  (merge)  │       │  (report)   │       │  (output)   │          ║
║          └───────────┘       └─────────────┘       └─────────────┘          ║
║                                                                             ║
║                     Final: merged_output.csv + HTML report                  ║
╚═════════════════════════════════════════════════════════════════════════════╝
```

### CtcReads-Caller (Steps 1-10)

| Step | Module | Description |
|------|--------|-------------|
| 1 | check_dependencies | Validate required tools and dependencies |
| 2 | tidehunter | Detect tandem repeats in HiFi reads |
| 3 | tandem_to_ring | Convert tandem repeats to circular candidates |
| 4 | run_alignment | Map candidates to reference genome (minimap2) |
| 5 | um_classify | Classify eccDNA as Unique (U) or Multiple (M) origin |
| 6 | cecc_build | Identify Chimeric (C) eccDNA using LAST |
| 7 | umc_process | Process and cluster U/M/C types |
| 8 | cd_hit | Remove redundancy (99% identity threshold) |
| 9 | ecc_dedup | Deduplicate and standardize coordinates |
| 10 | read_filter | Filter confirmed eccDNA reads |

### SplitReads-Caller (Steps 11-13)

| Step | Module | Description |
|------|--------|-------------|
| 11 | minimap2 | Prepare reference index |
| 12 | ecc_inference | Detect eccDNA (Cresil preferred, Cyrcular fallback) |
| 13 | curate_inferred_ecc | Curate inferred eccDNA |

### Integration (Steps 14-16)

| Step | Module | Description |
|------|--------|-------------|
| 14 | ecc_unify | Merge confirmed and inferred results |
| 15 | ecc_summary | Generate statistics and summaries |
| 16 | ecc_packager | Package final outputs |

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
│   ├── sample_name_CeccDNA_C.fasta      # Confirmed chimeric eccDNA sequences
│   ├── sample_name_CeccSegments.bed     # All segments
│   ├── sample_name_CeccJunctions.bedpe  # Junction information
│   └── sample_name_CeccSegments.core.csv # Detailed information table
├── sample_name_Inferred_eccDNA/
│   ├── sample_name_UeccDNA_I.csv        # Inferred simple eccDNA
│   ├── sample_name_UeccDNA_I.fasta      # Inferred simple sequences
│   ├── sample_name_chimeric.csv         # Inferred chimeric eccDNA
│   └── sample_name_CeccDNA_I.fasta      # Inferred chimeric sequences
├── sample_name_merged_output.csv        # All eccDNA combined
├── sample_name_report.html              # Interactive HTML report
└── sample_name_summary.txt              # Summary statistics
```

### Output Format

The merged output file (`sample_name_merged_output.csv`) contains:

| Column | Description |
|--------|-------------|
| eccDNA_id | Unique identifier (e.g., UeccDNA1, MeccDNA1, CeccDNA1) |
| original_id | Original eccDNA ID prior to final renumbering |
| Regions | Genomic coordinates (chr:start-end format) |
| Strand | DNA strand (+/-) |
| Length | eccDNA size in base pairs |
| eccDNA_type | Classification (UeccDNA/MeccDNA/CeccDNA) |
| State | Detection method (Confirmed/Inferred) |
| Seg_total | Number of segments (for chimeric eccDNA) |
| Hit_count | Number of genomic hits |

> Note: confidence/evidence fields are recorded in the per-type `*.core.csv` files under `Confirmed_*` directories.

## Configuration

Create a custom configuration file in YAML format:

```yaml
# Basic settings
performance:
  threads: 16

# Tool-specific parameters
tools:
  tidehunter:
    k: 16
    w: 1
    p: 100
    P: 2000000
    e: 0.1
    f: 2

  minimap2_align:  # For candidate alignment (Step 4)
    preset: "map-hifi"
    max_target_seqs: 5
    additional_args: ""
    # Identity filtering with length-based compensation (HiFi optimized)
    min_identity: 99.0           # Base identity threshold (%)
    identity_decay_per_10kb: 0.5 # Identity decay per 10kb of sequence length (%)
    min_identity_floor: 97.0     # Minimum identity floor (%)

  minimap2:  # For read mapping (Step 11)
    preset: "map-hifi"
    additional_args: ""
```

For a complete list of configuration options, see the [Configuration Reference](docs/Configuration_Reference.md).

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
| minimap2 | ≥2.24 | Sequence alignment and read mapping |
| samtools | ≥1.17 | BAM file processing |
| CD-HIT | ≥4.8.1 | Sequence clustering (cd-hit-est) |
| LAST | ≥1250 | High-accuracy CeccDNA detection |

### External Tools (Optional - for Cyrcular inference)

| Tool | Version | Description |
|------|---------|-------------|
| bcftools | ≥1.17 | BCF/VCF processing |
| varlociraptor | ≥5.0.0 | Variant calling |

### Inference Engines (At least one required)

| Tool | Description |
|------|-------------|
| Cresil | **Recommended** - Higher accuracy for eccDNA inference |
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

**FASTQ Input**: CircleSeeker requires FASTA format. Convert FASTQ first:
```bash
seqtk seq -A reads.fastq > reads.fasta
```

**Running in Background**: Use `nohup` instead of `&` to avoid multiprocessing issues:
```bash
# Recommended
nohup circleseeker -i input.fa -r ref.fa -o output -t 16 > run.log 2>&1 &

# Or use screen/tmux
screen -S circleseeker
circleseeker -i input.fa -r ref.fa -o output -t 16
```

**Out of Memory**: Reduce thread count or increase available memory
```bash
circleseeker -i input.fa -r ref.fa -t 4
```

**Missing Tools**: Install via conda
```bash
conda install -c bioconda -c conda-forge \
  tidehunter minimap2 samtools bcftools cd-hit last varlociraptor
```

## Documentation

For detailed documentation, see the `docs/` directory:

- [Pipeline Modules](docs/Pipeline_Modules.md) - Detailed algorithm descriptions
- [CLI Reference](docs/CLI_Reference.md) - Complete command-line options
- [Validation Methodology](docs/Validation_Methodology_en.md) - Synthetic U/M/C validation and recall benchmark

## Citation

If you use CircleSeeker in your research, please cite:

> Zhang Y, You M, Zhou C, et al. (2024). MMC-seq reveals a vast spatiotemporal eccDNA landscape from single cells to tissues. *Manuscript submitted*.

CircleSeeker is the computational pipeline developed as part of the MMC-seq methodology for comprehensive eccDNA detection and characterization.

## License

CircleSeeker is dual-licensed: GNU GPL v3.0 (see [LICENSE](LICENSE)) or a commercial license (see [COMMERCIAL_LICENSE.md](COMMERCIAL_LICENSE.md)).

## Contact

- Yaoxin Zhang: yxzhang@ncgr.ac.cn
- Leo Xinqi Yu: leoxqy@hotmail.com

## Acknowledgments

We thank the developers of TideHunter, Cresil, Cyrcular, LAST, and other integrated tools.
