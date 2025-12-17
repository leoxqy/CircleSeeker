<p align="center">
  <img src="docs/images/logo.png" width="300" alt="CircleSeeker">
</p>

# CircleSeeker

[![Version](https://img.shields.io/badge/version-0.9.3-blue.svg)](https://github.com/yaoxinzhang/CircleSeeker)
[![Python](https://img.shields.io/badge/python-≥3.9-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-GPL--2.0-green.svg)](LICENSE)

Comprehensive detection and characterization of extrachromosomal circular DNA (eccDNA) from PacBio HiFi sequencing data.

## Installation

### Using Conda (Recommended)

```bash
conda create -n circleseeker -c bioconda -c conda-forge circleseeker=0.9.3
conda activate circleseeker
```

### From Source

```bash
git clone https://github.com/yaoxinzhang/CircleSeeker.git
cd CircleSeeker
conda env create -f environment.yml
conda activate circleseeker
pip install -e .
```

### Optional: Install Cresil (preferred inference engine)

Cresil provides the fastest inference workflow. Install it alongside CircleSeeker:

```bash
git clone https://github.com/visanuwan/cresil.git
cd cresil
pip install -e .
cresil --help   # Verify installation
```

Cresil requires a minimap2 `.mmi` index of the reference genome. CircleSeeker will create the index automatically when it is missing.

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
# Minimum configuration
min_ecc_length: 100
max_ecc_length: 1000000
fdr_threshold: 1.0

# Performance
threads: 16
memory_limit: "32G"

# Algorithm parameters
blast_identity: 85
cluster_identity: 90
min_tandem_copies: 2
```

## System Requirements

- **Operating System**: Linux or macOS
- **Python**: ≥3.9
- **Memory**: Minimum 16GB (32GB recommended)
- **Storage**: ~50GB free space
- **CPU**: 8+ cores recommended

## Dependencies

### Python Packages
- numpy ≥1.23
- pandas ≥2.0
- biopython ≥1.81
- pysam ≥0.22
- networkx ≥3.0
- intervaltree ≥3.1

### External Tools
- TideHunter ≥1.4.0
- BLAST+ ≥2.12.0
- minimap2 ≥2.24
- samtools ≥1.17
- bcftools ≥1.17
- CD-HIT ≥4.8.1
- varlociraptor ≥5.0.0
- cyrcular

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
conda install -c bioconda tidehunter blast minimap2
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
