# Changelog

All notable changes to CircleSeeker will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.10.3] - 2026-01-16

### Added
- **XeccDNA source reads inclusion**: XeccDNA source reads are now included in inference input
  - Extracts original read names from XeccDNA FASTA (ring ID format: `{readName}|{repN}|{consLen}|{copyNum}|circular`)
  - Appends corresponding source reads to filtered FASTA before inference
  - Gives XeccDNA reads another chance to be detected as Inferred eccDNA
- **ResultKeys constant**: Added `XECC_SOURCE_READS_ADDED` to track XeccDNA reads appended to inference input

### Fixed
- **iecc_curator column names**: Fixed column name inconsistency
  - `eccdna_type` → `eccDNA_type` (matches downstream ecc_unify expectations)
  - `state` → `State` (matches confirmed eccDNA column naming)
- **Assert to explicit exception**: Replaced `assert reference_mmi is not None` with proper `PipelineError` raise
- **keep_tmp CLI parameter**: Fixed type annotation (`Optional[bool]` → `bool`) and config priority handling
- **Unused type imports**: Cleaned up unused `Dict`, `List`, `Set`, `Tuple` imports across all modules
- **check_docs_sync.py**: Removed references to deleted documentation files (Simulation_Validation)

### Changed
- **XeccDNA output handling**: XeccDNA.fasta is no longer copied to final output
  - Source reads are now in inference pipeline, detected as Inferred eccDNA if valid
  - XeccDNA.fasta retained in temp directory for debugging only
- **Log files location**: External tool logs now saved to `logs/` subdirectory
  - Previously logs like `minimap2_index.log` were saved to output root
  - Now organized under `{output}/logs/` for cleaner output structure

### Removed
- **Inference fallback to original input**: Removed incorrect fallback logic in minimap2 and Cresil steps
  - Inference now correctly requires filtered reads from earlier pipeline steps
  - Skips inference if filtered FASTA is missing (instead of using unfiltered input)
- **External tool timeout**: Removed default timeout for external bioinformatics tools
  - Cresil, Varlociraptor, and other tools now run without timeout limits
  - Large datasets may require extended runtime; timeouts caused premature termination

## [0.10.2] - 2026-01-15

### Added
- **CI coverage reporting**: Added pytest-cov with 80% threshold and Codecov integration
- **Type checking in CI**: Added mypy typecheck job for Python 3.12

### Fixed
- **Type annotations**: Unified type annotations across 64 Python files
  - Replaced `List[T]` with `list[T]`, `Dict[K,V]` with `dict[K,V]`, etc.
  - Consistent with Python 3.9+ and PEP 585
- **Documentation version numbers**: Updated version references from 0.10.0 to 0.10.2
- **README broken links**: Fixed Simulation_Validation → Validation_Methodology link
- **pyproject.toml license format**: Fixed deprecated license string to `{ text = "GPL-3.0-only" }`

### Changed
- **GitHub repository links**: Updated all links to new owner `leoxqy/CircleSeeker`

## [0.10.1] - 2026-01-14

### Added
- **Adaptive secondary mapping thresholds**: Improved Uecc classification with three-layer adaptive logic
  - Relative threshold (`u_secondary_max_ratio`): Allow small secondary mappings relative to primary coverage
  - High coverage tolerance (`u_high_coverage_threshold`): Relax thresholds when primary coverage ≥ 98%
  - MAPQ weighting (`u_high_mapq_threshold`): Increase tolerance when MAPQ ≥ 50 indicates unique mapping
- **Simulation auto-scaling**: Automatically scale genome size to maintain consistent repeat density
  - `auto_scale_genome`: Enable/disable automatic genome scaling (default: True)
  - `target_repeat_coverage`: Target repeat coverage ratio (default: 1.5%)
  - `max_repeat_coverage`: Maximum repeat coverage when auto-scaling disabled (default: 30%)

### Fixed
- **Large-scale recall regression**: Significantly improved recall on 10,000+ scale datasets
  - Same conditions (5 Mb genome, ~100% repeat coverage): 87.46% → 91.96% (+4.5%)
  - With auto-scaling (1.5% repeat coverage): Achieves 100% recall
- **False negative reduction**: Reduced FN count by 36% in extreme repeat density scenarios

### Changed
- **Secondary mapping veto logic**: `_u_has_significant_secondary_mapping` now considers primary coverage strength and MAPQ when evaluating secondary evidence

## [0.10.0] - 2026-01-13

### Added
- **MAPQ-aware Uecc gating**: Preserve minimap2 PAF `mapq` into alignment TSV and optionally veto Uecc calls below `tools.um_classify.mapq_u_min` (default 0 disables)
- **Simulation validation**: Add synthetic U/M/C validation and recall benchmark scripts (`tests/simulation/`) and documentation
- **Confidence scoring**: Emit `confidence_score` plus evidence fields (MAPQ/identity/coverage/2nd-locus) in confirmed/merged tables for traceable filtering
- **Negative control**: Add low-MAPQ negative-control tests to cap high-confidence false positives

### Fixed
- **Performance regression**: Avoid redundant CECC analysis on already-classified reads and speed up `um_classify` pandas-heavy code paths
- **Temporary directory safety**: Validate `runtime.tmp_dir` and restrict auto-cleanup to subdirectories under the output directory
- **Coordinate conventions**: Standardized subject coordinates to 0-based half-open (`start0/end0`) across minimap2 alignment conversion and downstream modules
- **CD-HIT outputs**: Pipeline now detects representative FASTA whether `cd-hit-est -o` is given with or without a `.fasta` suffix
- **Pipeline logs**: Step numbers in runtime logs are now 1-based (consistent with `--show-steps`, `--start-from`, `--stop-at`)
- **Large FASTA merge**: Stream large combined FASTA writes to avoid high memory use
- **Minimap2 index handling**: Reuse existing `.mmi` next to reference when available; otherwise build in the pipeline temp directory
- **Cresil reference indexing**: Fallback `.fai` creation via output-dir symlink when reference directory is not writable

### Changed
- **Alignment TSV format**: Append a trailing `mapq` column when converting minimap2 PAF output; U/M classifier accepts both 13- and 14-column inputs for compatibility
- **Documentation**: Updated CLI/pipeline docs to reflect minimap2-only alignment (BLAST references removed)
- **Packaging metadata**: Aligned version/license metadata across README/conda recipe and code
- **Pipeline structure**: Split monolithic `pipeline.py` into step executors + IO contracts with schema validation

## [0.9.8] - 2025-12-31

### Changed
- **BLAST removed**: Removed BLAST as alignment option, minimap2 is now the only aligner
  - Deleted `external/blast.py` and related test files
  - Renamed `minimap2_blast.py` → `minimap2_align.py`
  - Removed `--aligner` and `--blast-word-mode` CLI options
  - Removed `makeblastdb` pipeline step
- **Pipeline step renaming**: `run_blast` → `run_alignment` for clarity
- **Step numbering**: Pipeline now has 16 steps (1-16)
  - Step 1: check_dependencies (new)
  - Step 2: tidehunter
  - Step 3-16: remaining steps
- **Output filenames**: `*_blast_results.tsv` → `*_alignment_results.tsv`
- **Generic naming**: Renamed BLAST-specific identifiers to generic alignment names
  - `BLAST_COLUMNS` → `ALIGNMENT_COLUMNS`
  - `read_blast_results()` → `read_alignment_results()`
  - `classify_blast_results()` → `classify_alignment_results()`

### Added
- **Step 1: check_dependencies**: Dependency checking is now a formal pipeline step
  - Runs as the first step before tidehunter
  - Validates all required tools before pipeline execution

### Fixed
- **Best alignment selection**: Fixed `ecc_dedup.py` to use `alignment_length` instead of `bit_score` when bit_score is 0 (minimap2 output)
  - Priority: bit_score (if non-zero) > alignment_length > first row

## [0.9.6] - 2025-12-29

### Added
- **Pipeline stderr coverage tests**: Added unit tests to ensure minimap2 and varlociraptor pipelines do not pipe stderr (prevents deadlocks)
- **Minimap2 alignment option**: Added a minimap2-backed alignment path for candidate mapping (PAF converted to BLAST TSV) with `tools.aligner=minimap2`
- **CLI aligner and BLAST presets**: Added `--aligner` and `--blast-word-mode` to switch aligner and word_size presets without editing config

### Fixed
- **Subprocess pipeline deadlocks**: Avoided piping stderr in minimap2→samtools and varlociraptor→bcftools pipelines
- **Report version string**: ecc_summary now uses the package version instead of a hard-coded v2.1.0

### Changed
- **License**: Switched to GPL-3.0-only with a commercial license option (dual-licensed) and updated metadata/docs accordingly
- **Version metadata**: Updated to 0.9.6 across packaging and documentation
- **BLAST defaults**: Adjusted `word_size` to 50 and `max_target_seqs` to 200 for speed/quality balance
- **Default aligner**: Set `tools.aligner` to `minimap2` for candidate alignment by default

### Removed
- **ReadFilterAdapter**: Removed unused adapter that referenced a non-existent ReadFilter class

## [0.9.5] - 2025-12-26

### Added
- **Inference engine display**: Shows which inference tool (Cresil/Cyrcular) will be used at startup
  - Displays in configuration summary before pipeline execution
  - Shows warning if only Cyrcular is available (Cresil recommended)

### Fixed
- **CLI code refactoring**: Extracted common logic to `_execute_pipeline()` function
  - Eliminates code duplication between `cli()` and `run()` commands
  - Reduces maintenance burden and potential behavior divergence
- **`--show-steps` side effect**: No longer creates output directories when only viewing steps
  - Uses lightweight `_show_pipeline_steps()` that reads class-level STEPS directly
- **Bare except clauses**: Replaced with specific exception types
  - `display.py`: 2 occurrences fixed
  - `ecc_unify.py`: 3 occurrences fixed
  - `external_tools.py`: 4 occurrences fixed
- **Empty f-strings**: Removed 14 unnecessary f-string prefixes across modules
- **Trailing whitespace**: Cleaned up 4 files
- **`import *` usage**: Replaced with explicit imports in `modules/tools/__init__.py`
- **Duplicate import**: Removed redundant `import pandas` in `column_standards.py`
- **TideHunter case sensitivity**: Now correctly detects both `TideHunter` and `tidehunter` executables
  - Bioconda may install the tool with different name casing on different systems
  - Added `_find_tidehunter_executable()` helper function
  - Updated dependency checker and validators to check alternative names
- **Project URLs**: Fixed GitHub repository URLs in `pyproject.toml` and `conda-recipe/meta.yaml`
- **Config file parameter handling**: Now properly reads all parameters from config files
  - CLI arguments take precedence over config values
  - Config values take precedence over hardcoded defaults
  - Fixed: output_dir, prefix, threads now correctly use config file values when not specified via CLI
  - Added `cfg.validate()` call before dependency check (fail fast on missing files)
  - `--keep-tmp/--no-keep-tmp` flag pair now explicitly overrides config file `keep_tmp` setting
  - Error message clarifies config file usage
- **Dependency checker improvements**:
  - Added min_version checking with version comparison
  - Added TideHunter to dependency list with alt_names support
  - bcftools/varlociraptor now conditionally required when only Cyrcular is available (no Cresil)
  - Added cd-hit-est to `validate --full` check
  - Version warnings displayed but don't block execution
- **Removed dead code**: `--skip-organize` flag removed from `run` subcommand (was never wired up)
- **IeccCurator export fix**: Replaced non-existent class with actual functions
  - Exports: `curate_ecc_tables`, `generate_fasta_sequences`, `write_curated_tables`, etc.
  - Fixed import paths in `modules/tools/__init__.py`
- **Removed duplicate config template**: Consolidated to single source in `resources/__init__.py`

### Changed
- **Version metadata**: Updated to 0.9.5 across all files including `conda-recipe/meta.yaml`

## [0.9.4] - 2025-12-25

### Added
- **Chimeric eccDNA Overlap Detection**: New segment-wise reciprocal overlap method for detecting redundant chimeric eccDNA
  - Replaces exact string matching with flexible segment-by-segment comparison
  - Uses 99% reciprocal overlap threshold with ±10bp coordinate tolerance
  - Significantly reduces false negatives caused by minor coordinate differences (2-10bp)
  - Backward compatible: both 'exact' and 'overlap' methods available via `method` parameter

- **Pre-flight Dependency Checker**: Comprehensive tool dependency checking at startup
  - Checks all required tools (minimap2, samtools, blastn, makeblastdb, cd-hit-est) before pipeline execution
  - Verifies at least one inference tool (cresil or cyrcular) is available
  - Clear, actionable error messages with installation hints
  - Prevents wasting computation time on incomplete installations
  - Distinguishes between required and optional dependencies

### Fixed
- **`verbose` parameter**: Fixed always-on verbose output in ecc_packager module (19 occurrences)
- **ResultKeys constant mismatch**: Fixed AttributeError when overview table is missing (TSV vs CSV naming)

### Removed
- **`--skip-report` flag**: Removed from CLI and pipeline (reports are always generated and packaged)

### Changed
- **BLAST soft masking**: Disabled soft masking in BLAST alignment (inherited from 0.9.3)
  - Added `-soft_masking false` to BLAST runner
  - Ensures confirmed candidates can align even when inputs are soft-masked

- **README Configuration**: Updated configuration examples to match actual implementation
  - Removed undocumented parameters (min_ecc_length, fdr_threshold, etc.)
  - Added correct tool-specific configuration structure
  - Now reflects actual config.py structure

### Technical Details
- **Files Modified**:
  - `src/circleseeker/modules/ecc_unify.py`:
    - Added `method` parameter to `find_redundant_chimeric()` (default: 'overlap')
    - New function: `find_redundant_chimeric_overlap()` - segment-wise matching
    - New function: `find_redundant_chimeric_exact()` - legacy exact matching
    - New helper functions: `_segments_match()`, `reciprocal_overlap_ok()`
  - `src/circleseeker/external/blast.py`:
    - Added `soft_masking` parameter to `BlastN.run_blast()` (default: False)
    - Added `soft_masking` to `BlastRunner.__init__()` and `BlastRunner.run()`
  - `src/circleseeker/config.py`:
    - Added `soft_masking: False` to BLAST configuration
  - `tests/unit/test_ecc_unify.py`:
    - New tests for chimeric overlap detection edge cases
  - Version files: `src/circleseeker/__version__.py`, `pyproject.toml`

- **Impact**:
  - **Reduced false negative rate**: Expected to drop from ~20% to ~5% for chimeric eccDNA
  - **Better biological accuracy**: Tolerates natural coordinate variations from different tools
  - **Backward compatible**: Old method still available via `method='exact'`

### Example

Before (0.9.3 and earlier):
```python
# Two biologically identical chimeric eccDNA with 5bp coordinate difference
Cecc1: chr1:1000-2000;chr2:3000-4000  # Confirmed
Cecc2: chr1:1000-2000;chr2:3005-4005  # Inferred (5bp shift)
# Result: Cecc2 NOT detected as redundant (false negative)
```

After (0.9.4):
```python
# Same data with new overlap method
# Result: Cecc2 correctly detected as redundant ✅
# Reason: 99% overlap + within 10bp tolerance
```

---

## [0.9.3] - 2025-12-16

### Changed
- **BLAST nucleotide alignment**: Lowercase (soft-masked) regions are no longer suppressed during seeding/search. Added `-soft_masking false` to the BLAST runner so confirmed candidates can align even when inputs are soft-masked.
- Updated version metadata and README badges to 0.9.3.

### Technical Details
- **Files Modified**:
  - `src/circleseeker/external/blast.py`: Added explicit `-soft_masking false` option to `blastn` invocation via `BlastRunner`.
  - `src/circleseeker/__version__.py`, `pyproject.toml`, `README.md`: Bumped version to 0.9.3.

---

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

### From 0.9.2/0.9.3 to 0.9.4

No breaking changes. The new overlap detection method is automatically enabled by default.

```bash
# If installed via conda
conda activate circleseeker
pip install --upgrade circleseeker

# Or reinstall from source
cd CircleSeeker
pip install -e .
```

### Configuration

The new overlap method is enabled by default. To use the legacy exact matching:

```python
# In your code (if using CircleSeeker as a library)
from circleseeker.modules.ecc_unify import find_redundant_chimeric

redundant_ids = find_redundant_chimeric(
    inferred_df,
    confirmed_df,
    method='exact'  # Use legacy method
)
```

---

## Versioning

CircleSeeker follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html):
- **MAJOR** version for incompatible API changes
- **MINOR** version for backward-compatible functionality additions
- **PATCH** version for backward-compatible bug fixes
