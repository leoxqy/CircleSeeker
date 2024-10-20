# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.3] - 2024-10-20

### Added
- Implemented skipping of reads with more than 5 repeated regions in the Carousel section. This optimization prevents excessive processing time for large groups that were causing significant delays.
- Integrated progress bars using tqdm:
  1. In the Carousel section for processing repeated readNames.
  2. In the Ringmaster section for analyzing groups.

### Changed
- Modified the Carousel processing to skip combining rows when the number of repeated regions exceeds 5, improving overall performance for large datasets.

### Performance
- Significantly reduced processing time for datasets with large numbers of repeated regions by implementing the skipping mechanism in the Carousel section.

## [1.0.3] - 2024-10-19

### Fixed
- Added conditional checks in `MergeMecc.py` and `MergeUecc.py` to handle empty or non-existent FASTA files, preventing errors when trying to open empty files.
