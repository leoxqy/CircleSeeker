# CircleSeeker 2

[![PyPI version](https://badge.fury.io/py/circleseeker.svg)](https://badge.fury.io/py/circleseeker)
[![Conda Version](https://img.shields.io/conda/vn/bioconda/circleseeker.svg)](https://anaconda.org/bioconda/circleseeker)
[![License: GPL v2](https://img.shields.io/badge/License-GPL_v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
[![Python Version](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)

## Overview

CircleSeeker 2 is a comprehensive bioinformatics tool for detecting and characterizing extrachromosomal circular DNA (eccDNA) from PacBio HiFi sequencing data. It provides an end-to-end pipeline for identifying various types of eccDNA with high accuracy and detailed annotations.

### Key Features

- **Comprehensive eccDNA Detection**: Identifies multiple types of circular DNA
  - **UeccDNA** (Unique eccDNA): Single-copy circular DNA elements
  - **MeccDNA** (Multi-region eccDNA): Circles derived from multiple genomic regions
  - **CeccDNA** (Complex eccDNA): Structurally complex circular DNA
  - **XeccDNA** (Extra-chromosomal eccDNA): Optional detection of additional circular elements

- **Advanced Analysis Pipeline**:
  - Tandem repeat detection using TideHunter
  - Homology-based classification with BLAST
  - Coverage analysis with Minimap2 and Mosdepth
  - Redundancy removal with CD-HIT-EST
  - Comprehensive reporting in multiple formats

- **Production-Ready Features**:
  - Checkpoint-based resumption from failures
  - Configurable via YAML files
  - Multi-threaded processing
  - Detailed logging and error handling
  - HTML reports with interactive visualizations

## Installation

### Using Conda (Recommended)

```bash
conda install -c bioconda circleseeker
```

### Using pip

```bash
pip install circleseeker
```

### From Source

```bash
git clone https://github.com/circleseeker/circleseeker2.git
cd circleseeker2
pip install -e .
```

### Dependencies

CircleSeeker requires the following external tools to be installed and accessible in your PATH:

- **TideHunter** ≥ 1.4.3
- **BLAST+** ≥ 2.13.0
- **CD-HIT** ≥ 4.8.1
- **Minimap2** ≥ 2.24
- **Samtools** ≥ 1.16
- **Mosdepth** ≥ 0.3.3

## Quick Start

### Basic Usage

```bash
# Run CircleSeeker on your HiFi reads
circleseeker -i input.fasta -o output_dir -r reference.fasta

# With multiple threads
circleseeker -i input.fasta -o output_dir -r reference.fasta -t 16

# Enable XeccDNA detection
circleseeker -i input.fasta -o output_dir -r reference.fasta --enable-xeccdna
```

### Using Configuration Files

```bash
# Generate a default configuration file
circleseeker --generate-config > my_config.yaml

# Run with custom configuration
circleseeker -i input.fasta -o output_dir -r reference.fasta -c my_config.yaml
```

### Resume from Checkpoint

```bash
# If a run fails, resume from the last successful step
circleseeker -i input.fasta -o output_dir -r reference.fasta --resume
```

## Output Files

CircleSeeker generates organized outputs in the following structure:

```
output_dir/
├── organized_results/
│   ├── {sample}_UeccDNA.bed          # Unique eccDNA coordinates
│   ├── {sample}_UeccDNA.core.csv     # Detailed UeccDNA information
│   ├── {sample}_MeccDNA.sites.bed    # Multi-region eccDNA sites
│   ├── {sample}_MeccDNA.bestsite.bed # Best representative sites
│   ├── {sample}_MeccDNA.core.csv     # Detailed MeccDNA information
│   ├── {sample}_CeccDNA.segments.bed # Complex eccDNA segments
│   ├── {sample}_CeccDNA.junctions.bedpe # Junction coordinates
│   ├── {sample}_InferredUeccDNA.bed  # Inferred unique eccDNA (optional)
│   └── {sample}_file_summary.txt     # Summary of all outputs
├── reports/
│   └── {sample}_eccDNA_report.html   # Interactive HTML report
├── logs/
│   └── circleseeker_{timestamp}.log  # Detailed execution log
└── intermediate/                      # Intermediate files (can be deleted)
```

### Output Format Descriptions

- **BED files**: Standard BED format with eccDNA coordinates
- **CSV files**: Detailed information including:
  - Circle ID and type
  - Genomic coordinates
  - Sequence length
  - Coverage statistics
  - Quality metrics
- **BEDPE files**: Paired-end format for junction coordinates
- **HTML report**: Comprehensive visualization including:
  - Summary statistics
  - Size distributions
  - Chromosomal distributions
  - Coverage analysis
  - Interactive charts

## Advanced Usage

### Pipeline Control

```bash
# Skip specific steps
circleseeker -i input.fasta -o output_dir -r reference.fasta --skip-report --skip-organize

# Run only specific modules
circleseeker run-tidehunter -i input.fasta -o output_dir
circleseeker run-tandem-to-ring -i tidehunter_output.txt -o output_dir
```

### Performance Tuning

```yaml
# config.yaml
performance:
  chunk_size: 10000      # Process reads in chunks
  max_memory: 16         # Maximum memory in GB
  compression: true      # Use compression for intermediate files

quality:
  min_circle_size: 100   # Minimum circle size to report
  min_coverage: 2.0      # Minimum coverage threshold
  max_breakpoints: 10    # Maximum allowed breakpoints
```

## Documentation

- Test guide: `docs/testing/HOW_TO_USE_TESTS.md`
- Technical plan: `docs/Technical_Reconstruction_Plan.md`
- More docs and images: `docs/`

## Project Structure

- `src/circleseeker`: Python package with core, modules, external tool wrappers, reports, and resources.
- `tests`: Pytest suite with unit/integration tests and `test_data/` fixtures.
- `docs`: Long-form docs, guides, and assets (images).
- `scripts`: Dev utilities and test runner scripts.
- `conda-recipe`: Conda packaging recipe.
- `release`: Local build artifacts and distributions.

## Development

- Install dev deps: `make dev`
- Install pre-commit hooks: `pre-commit install`
- Run tests: `make test` (or see `scripts/run_tests.sh`)

## Citation

If you use CircleSeeker in your research, please cite:

```
CircleSeeker 2: Comprehensive detection of extrachromosomal circular DNA from PacBio HiFi sequencing
[Citation details to be added upon publication]
```

## License

CircleSeeker is distributed under the GNU General Public License v2.0. See [LICENSE](LICENSE) file for details.

## Support

- **Documentation**: [https://circleseeker2.readthedocs.io](https://circleseeker2.readthedocs.io)
- **Issues**: [GitHub Issues](https://github.com/circleseeker/circleseeker2/issues)
- **Discussions**: [GitHub Discussions](https://github.com/circleseeker/circleseeker2/discussions)

## Contributing

We welcome contributions! Please see our [Contributing Guide](CONTRIBUTING.md) for details.

## Acknowledgments

CircleSeeker uses a circus-themed naming convention for its modules, making the complex pipeline more memorable and maintainable. Special thanks to all contributors and the bioinformatics community for their support.
