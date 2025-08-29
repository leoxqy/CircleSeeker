# CircleSeeker2

[![Version](https://img.shields.io/badge/version-2.0.0-blue)](pyproject.toml)
[![Python](https://img.shields.io/badge/python-3.9+-green)](pyproject.toml)
[![License](https://img.shields.io/badge/license-GPL--2.0-blue)](LICENSE)

## Overview

CircleSeeker2 is a comprehensive bioinformatics pipeline for detecting extrachromosomal circular DNA (eccDNA) from PacBio HiFi long-read sequencing data. Using a multi-round iterative analysis approach, it achieves high-sensitivity identification of various eccDNA types.

### Key Features

- 🔍 **Multi-round Analysis**: TideHunter → BLAST+ (2x) → Minimap2 for comprehensive detection
- ⚡ **High Performance**: Parallel processing with efficient memory management
- 📊 **Multiple eccDNA Types**: Detects UeccDNA, MeccDNA, CeccDNA, MCeccDNA, and XeccDNA
- 🔧 **Modular Design**: Clean architecture with reusable components
- 📦 **Easy Installation**: Modern Python packaging with pip/conda support
- 📈 **Comprehensive Reports**: HTML, CSV, and FASTA outputs with detailed statistics
- 🔄 **Checkpoint System**: Resume pipeline from any step

## Installation

### Development Installation

```bash
# Clone repository
git clone <repository-url>
cd CircleSeeker2_dev

# Install in development mode
pip install -e .

# Verify installation
circleseeker2 --version
circleseeker2 validate
```

### Dependencies

**Python packages** (automatically installed):
- click>=8.1
- pandas>=2.0
- numpy>=1.23
- biopython>=1.81
- pysam>=0.22
- pyyaml>=6.0
- networkx>=3.0

**External tools** (install separately):
- TideHunter
- BLAST+ (makeblastdb, blastn)
- Minimap2
- Samtools
- Mosdepth
- CD-HIT

## Quick Start

### Basic Usage

```bash
# Run with minimal parameters
circleseeker2 run input.fasta reference.fa

# Run with custom output and threads
circleseeker2 run input.fasta reference.fa \
    -o results \
    -p sample1 \
    -t 16
```

### Using Configuration File

```bash
# Generate configuration template
circleseeker2 init-config -o config.yaml

# Edit config.yaml with your parameters
# Then run with configuration
circleseeker2 run input.fasta reference.fa -c config.yaml
```

### Advanced Usage

```bash
# Resume from checkpoint
circleseeker2 run input.fasta reference.fa -o previous_output

# Start from specific step
circleseeker2 run input.fasta reference.fa --start-from 5

# Stop at specific step
circleseeker2 run input.fasta reference.fa --stop-at 8

# Dry run to see pipeline steps
circleseeker2 run input.fasta reference.fa --dry-run
```

## Pipeline Architecture

### Modern Package Structure

```
src/circleseeker2/
├── __init__.py              # Package initialization
├── __version__.py           # Version information
├── __main__.py             # Package entry point
├── cli.py                  # Click-based CLI
├── config.py               # Configuration management
├── exceptions.py           # Custom exceptions
├── core/
│   └── pipeline.py         # Main pipeline orchestrator
├── external/               # External tool wrappers
│   ├── base.py            # Base wrapper class
│   ├── tidehunter.py      # TideHunter wrapper
│   ├── blast.py           # BLAST+ wrapper
│   ├── minimap2.py        # Minimap2 wrapper
│   ├── samtools.py        # Samtools wrapper
│   └── mosdepth.py        # Mosdepth wrapper
├── modules/                # Analysis modules
│   ├── carousel.py         # Circular DNA analysis
│   ├── gatekeeper.py      # eccDNA classification
│   └── processors/        # Processing utilities
├── algorithms/            # Core algorithms
├── io/                    # Input/Output utilities
├── utils/                 # Utility functions
│   └── validators.py      # Installation validation
├── reports/               # Report generation
└── resources/             # Configuration templates
```

### Pipeline Steps

1. **make_blastdb** - Build BLAST database from reference
2. **tidehunter** - Run TideHunter for tandem repeat detection
3. **carousel** - Process TideHunter output and classify sequences
4. **run_blast** - Run BLAST alignment against reference
5. **gatekeeper** - Classify eccDNA into Uecc/Mecc/Unclassified
6. **trapeze** - Process Cecc candidates (placeholder)
7. **menagerie** - Generate FASTA files (placeholder)
8. **cd_hit** - Remove redundant sequences (placeholder)
9. **harmonizer** - Coordinate results (placeholder)
10. **sieve** - Filter sequences (placeholder)
11. **minimap2** - Map with Minimap2 (placeholder)
12. **mosdepth** - Calculate depth coverage (placeholder)
13. **contortionist** - Split regions (placeholder)
14. **juggler** - Final processing (placeholder)

## Configuration

### YAML Configuration Format

```yaml
# Input files
input_file: "input.fasta"
reference: "reference.fa"
output_dir: "results"
prefix: "sample"

# Feature flags
enable_xecc: false

# Performance settings
performance:
  threads: 8
  max_memory: "16G"

# Quality control parameters
quality:
  min_quality_score: 0.99
  min_coverage: 10
  min_identity: 99.0

# External tool parameters
tools:
  tidehunter:
    k: 16
    w: 1
    p: 100
  blast:
    word_size: 100
    evalue: "1e-50"
    perc_identity: 99.0
```

## Output Files

```
output_dir/
├── config.yaml                    # Effective configuration
├── sample.checkpoint              # Pipeline checkpoint
├── sample_processed.csv           # Carousel results
├── sample_circular.fasta          # Circular sequences
├── sample_blast_results.tsv       # BLAST alignments
├── sample_uecc.csv                # UeccDNA results
├── sample_mecc.csv                # MeccDNA results
└── sample_unclassified.csv        # Unclassified results
```

## Commands

- `circleseeker2 run` - Run the main pipeline
- `circleseeker2 init-config` - Generate configuration template
- `circleseeker2 validate` - Validate installation
- `circleseeker2 benchmark` - Run performance benchmarks
- `circleseeker2 --help` - Show help information

## Development Status

✅ **Completed**:
- Modern package structure and CLI
- Configuration management system
- Pipeline architecture with checkpoints
- External tool wrappers
- Core modules (TideHunter, Carousel, Gatekeeper)
- Validation and testing framework

🚧 **In Progress**:
- Additional pipeline steps (steps 5-13)
- Comprehensive testing
- Documentation completion

## Migration from v1

The original CircleSeeker2 script has been completely restructured:

- **Modular Architecture**: Separated into logical modules
- **Configuration System**: YAML-based configuration with validation
- **Checkpoint System**: Resume from any step
- **Modern CLI**: Click-based interface with multiple commands
- **Error Handling**: Proper exception handling and logging
- **Package Structure**: Follows modern Python packaging standards

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make changes following the existing code style
4. Add tests for new functionality
5. Submit a pull request

## License

GNU General Public License v2.0 - see [LICENSE](LICENSE) file for details.

## Acknowledgments

CircleSeeker2 builds upon excellent bioinformatics tools:
- [TideHunter](https://github.com/yangao07/TideHunter) - Tandem repeat detection
- [BLAST+](https://blast.ncbi.nlm.nih.gov/) - Sequence alignment
- [Minimap2](https://github.com/lh3/minimap2) - Long-read alignment
- [samtools](https://github.com/samtools/samtools) - SAM/BAM manipulation