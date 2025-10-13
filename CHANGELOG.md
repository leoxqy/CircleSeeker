# Changelog

All notable changes to CircleSeeker will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.9.2] - 2025-01-11

### Fixed
- **Critical Bug Fix**: Fixed `read_filter` module to correctly handle semicolon-separated read IDs in core CSV files
  - Previously, multiple read IDs separated by semicolons were treated as a single ID
  - Now properly splits and processes each individual read ID
  - This ensures all reads from CD-HIT clusters are correctly filtered

- **Cecc Filtering Logic**: Fixed segment count filtering in `cecc_build` module
  - Changed from `len(g) <= min_segments` to `len(g) < min_segments`
  - Now correctly keeps eccDNA with ≥2 segments (previously required ≥3)
  - Aligns with Cecc definition: complex eccDNA must have at least 2 segments

### Added
- Added `logging` import to `read_filter.py` for proper type hints
- Enhanced documentation with inline comments explaining filtering logic

### Changed
- Improved docstring for `min_segments` parameter to be more explicit
- Updated version badge and installation instructions in README.md

### Technical Details
- **Files Modified**:
  - `src/circleseeker/modules/read_filter.py`: Lines 19, 113-117, 146, 216
  - `src/circleseeker/modules/cecc_build.py`: Lines 290, 502
  - `src/circleseeker/__version__.py`: Version updated to 0.9.2
  - `pyproject.toml`: Version updated to 0.9.2
  - `README.md`: Version badge and conda installation command

- **Impact**:
  - More accurate read filtering (no false negatives from clustered reads)
  - Correct Cecc classification (no false negatives from 2-segment eccDNA)
  - Better alignment with biological definitions of eccDNA types

## [0.9.1] - 2025-01-10

### Added
- Initial public release
- Comprehensive 16-step eccDNA detection pipeline
- Support for dual inference engines (Cresil/Cyrcular)
- HTML report generation
- Checkpoint and resume functionality

### Features
- Detection of Unique, Multiple, and Complex eccDNA types
- Integration of confirmed and inferred eccDNA
- Comprehensive output formats (CSV, BED, BEDPE, FASTA)
- Extensive test coverage (420+ unit tests)
- Detailed documentation

---

## Upgrade Instructions

### From 0.9.1 to 0.9.2

No breaking changes. Simply update your installation:

```bash
# If installed via conda
conda update circleseeker

# If installed from source
cd CircleSeeker
git pull
pip install -e .
```

**Important**: If you have existing results from version 0.9.1, we recommend re-running from step 10 (read_filter) onwards to benefit from the bug fixes:

```bash
circleseeker -i input.fa -r reference.fa -o existing_output_dir --start-from read_filter
```

This will reprocess the read filtering and subsequent steps with the corrected logic.
