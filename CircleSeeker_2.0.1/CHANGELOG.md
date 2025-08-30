# Changelog

All notable changes to CircleSeeker will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.1] - 2025-08-29

### Added
- Comprehensive README documentation
- HTML report generation with interactive visualizations (Playbill module)
- Automated output organization (Propmaster module)
- Support for XeccDNA detection (optional)
- Configuration file support (YAML)
- Checkpoint-based pipeline resumption
- Multi-threading support for performance optimization

### Changed
- Unified package naming to `circleseeker` (lowercase)
- Improved CLI with both direct run and subcommand options
- Enhanced error handling and logging
- Optimized memory usage for large datasets

### Fixed
- Package import consistency issues
- Command-line entry point compatibility
- Memory leaks in long-running analyses

## [2.0.0] - 2024-12-01

### Added
- Complete pipeline redesign with modular architecture
- Circus-themed module naming convention
- Support for PacBio HiFi sequencing data
- Three-stage analysis pipeline (Ringmaster, Astrologer, Final)
- Comprehensive eccDNA classification (UeccDNA, MeccDNA, CeccDNA)
- Integration with external tools (TideHunter, BLAST, CD-HIT, Minimap2, Mosdepth)
- Multiple output formats (BED, CSV, FASTA, BEDPE)

### Changed
- Migrated from monolithic script to modular package structure
- Improved algorithm efficiency for large-scale data
- Enhanced filtering and quality control mechanisms
- Standardized output formats

### Removed
- Legacy command-line interface
- Deprecated analysis methods

## [1.0.0] - 2023-06-15

### Added
- Initial release of CircleSeeker
- Basic eccDNA detection from sequencing data
- Simple command-line interface
- Basic output in TSV format

---

For more detailed information about each release, please refer to the [GitHub Releases](https://github.com/circleseeker/circleseeker2/releases) page.