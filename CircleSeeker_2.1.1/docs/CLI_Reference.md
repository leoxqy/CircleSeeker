# CircleSeeker CLI Reference

## Overview

CircleSeeker v2.1.1 - Comprehensive eccDNA detection from PacBio HiFi sequencing data.

## Basic Usage

```bash
circleseeker -i <reads.fa> -r <ref.fa> [options]
```

Or use the capitalized alias:
```bash
CircleSeeker -i <reads.fa> -r <ref.fa> [options]
```

## Required Parameters

- `-i, --input PATH` - Input FASTA file containing HiFi reads
- `-r, --reference PATH` - Reference genome FASTA file

## Optional Parameters

### Output Configuration
- `-o, --output PATH` - Output directory (default: `circleseeker_output`)
- `-p, --prefix TEXT` - Output file prefix (default: `sample`)

### Runtime Configuration
- `-c, --config PATH` - Configuration file in YAML format
- `-t, --threads INTEGER` - Number of threads to use (default: 8)

### Logging and Debugging
- `-n, --noise` - Increase log verbosity (use `-n` for INFO, `-nn` for DEBUG)
- `-h, --help` - Show help message and exit
- `--keep-tmp` - Retain temporary working directory after completion
- `--debug` - Enable advanced options and debug logging

### Version Information
- `-v, --version` - Show CircleSeeker version and exit

## Advanced Options (--debug mode only)

When `--debug` is enabled, additional options become available:

- `--start-from INT` - Start from a specific step (1-based index)
- `--stop-at INT` - Stop at a specific step (inclusive)
- `--resume` - Resume from last checkpoint
- `--force` - Force re-run all steps (ignore checkpoints)
- `--generate-config` - Print default configuration YAML and exit
- `--show-steps` - Show pipeline steps and exit
- `--dry-run` - Show actions without executing
- `--skip-report` - Skip HTML report generation
- `--log-output PATH` - Path for additional log file

## Pipeline Steps

CircleSeeker executes the following 16 steps:

1. **make_blastdb** - Build BLAST database from reference
2. **tidehunter** - Detect tandem repeats in HiFi reads
3. **tandem_to_ring** - Convert tandem repeats to circular candidates
4. **run_blast** - BLAST candidates against reference
5. **um_classify** - Classify into UeccDNA and MeccDNA
6. **cecc_build** - Build complex eccDNA (CeccDNA)
7. **umc_process** - Process U/M/C types and generate FASTA
8. **cd_hit** - Deduplicate sequences with CD-HIT
9. **organize_coords** - Organize and deduplicate coordinates
10. **read_filter** - Filter confirmed eccDNA reads
11. **minimap2** - Align filtered reads to reference
12. **cyrcular_calling** - Detect circular DNA with Cyrcular
13. **iecc_curator** - Curate inferred eccDNA
14. **ecc_unify** - Unify all eccDNA results
15. **ecc_summary** - Generate summary statistics
16. **ecc_packager** - Package final outputs

## Example Commands

### Basic run with default settings:
```bash
circleseeker -i sample.fasta -r hg38.fa -o results/
```

### Run with more threads and keep temporary files:
```bash
circleseeker -i sample.fasta -r hg38.fa -o results/ -t 16 --keep-tmp
```

### Run with custom configuration:
```bash
circleseeker -i sample.fasta -r hg38.fa -c custom_config.yaml -o results/
```

### Debug mode with verbose logging:
```bash
circleseeker -i sample.fasta -r hg38.fa --debug -nn -o debug_results/
```

### Generate default configuration file:
```bash
circleseeker --debug --generate-config > my_config.yaml
```

### Show pipeline steps without running:
```bash
circleseeker --debug --show-steps
```

## Output Structure

```
output_dir/
├── 00_raw_data/          # Input data links
├── 01_tidehunter/        # Tandem repeat detection
├── 02_tandem_to_ring/    # Circular candidates
├── 03_blast/             # BLAST results
├── 04_um_classify/       # U/M classification
├── 05_cecc_build/        # Complex eccDNA
├── 06_umc_process/       # Processed U/M/C
├── 07_alignment/         # minimap2 alignments
├── 08_cyrcular/          # Cyrcular results
├── 09_varlociraptor/     # Variant calling (FDR=0.05)
├── 10_dedup/             # CD-HIT deduplication
├── 11_final/             # Final results
├── logs/                 # Execution logs
├── report/               # HTML report
├── checkpoint/           # Pipeline checkpoints
└── config.yaml           # Used configuration
```

## Configuration File

CircleSeeker uses YAML configuration files. Key sections include:

```yaml
# Input/Output
input_file: reads.fasta
reference: reference.fa
output_dir: output/
prefix: sample

# Runtime
threads: 8
keep_tmp: false

# Analysis parameters
min_ecc_length: 100
max_ecc_length: 1000000
fdr_threshold: 0.05  # For varlociraptor

# Feature flags
enable_xecc: true  # Detect unknown eccDNA types
```

## Notes

- CircleSeeker automatically saves checkpoints and can resume interrupted runs
- The FDR threshold for varlociraptor is now set to 0.05 by default (previously 1.0)
- All external tools must be available in PATH
- Minimum Python version: 3.9

## See Also

- Main README.md for installation instructions
- pyproject.toml for dependency versions
- GitHub repository for issues and updates