# CircleSeeker Pipeline Modules (v0.9.5)

This document summarizes the 16-step CircleSeeker pipeline for English-speaking users.

---

## 1. Overview

CircleSeeker organizes the workflow into four stages:

| Stage | Steps | Purpose |
|-------|-------|---------|
| **Detection** | 1-6 | Discover circular DNA candidates from HiFi reads |
| **Processing** | 7-10 | Consolidate, de-duplicate, and filter U/M/C outputs |
| **Inference** | 11-13 | Run Cresil (or Cyrcular) and curate inferred eccDNA |
| **Integration** | 14-16 | Merge results, generate reports, and package deliverables |

Intermediates reside in `<output>/.tmp_work/`, and the final artifacts are copied into the packaged directory layout.

---

## 2. Step Summary

| # | Name | Type | Main Input | Main Output |
|---|------|------|------------|-------------|
| 1 | make_blastdb | External (BLAST+) | Reference FASTA | BLAST database |
| 2 | tidehunter | External | HiFi reads FASTA | Tandem repeat consensus |
| 3 | tandem_to_ring | Internal | TideHunter output | Candidate FASTA/CSV |
| 4 | run_blast | External (BLAST+) | Candidates & DB | BLAST TSV |
| 5 | um_classify | Internal | BLAST TSV | `um_classify.uecc.csv`, `um_classify.mecc.csv` |
| 6 | cecc_build | Internal | Unclassified hits | `cecc_build.csv` |
| 7 | umc_process | Internal | U/M/C tables | Harmonized CSV & FASTA |
| 8 | cd_hit | External (CD-HIT) | FASTA | Non-redundant FASTA |
| 9 | ecc_dedup | Internal | CSV/FASTA | Deduplicated coordinates |
|10 | read_filter | Internal | Deduped outputs | Confirmed eccDNA FASTA |
|11 | minimap2 | External (minimap2) | Reference | `.mmi` index |
|12 | ecc_inference | External + Internal | Cresil/Cyrcular inputs | Inferred TSV/FASTA |
|13 | iecc_curator | Internal | Inferred results | Curated CSV/FASTA |
|14 | ecc_unify | Internal | Confirmed & inferred CSV | Unified table, overlap stats |
|15 | ecc_summary | Internal | Unified table, processed CSV | HTML/TXT summary, all-fasta |
|16 | ecc_packager | Internal | Finalized files | End-user directory tree |

---

## 3. Stage Details

### 3.1 Detection (Steps 1-6)

1. **make_blastdb** - Build BLAST index from the reference genome.
2. **tidehunter** - Identify tandem repeats characteristic of rolling-circle amplification.
3. **tandem_to_ring** - Convert repeats into circular candidates via overlap graph analysis.
4. **run_blast** - Align candidates back to the reference to obtain genomic context.
5. **um_classify** - Categorize alignments into UeccDNA or MeccDNA, retaining unclassified records.
6. **cecc_build** - Assemble complex eccDNA (CeccDNA) structures from multi-segment hits.

### 3.2 Processing (Steps 7-10)

7. **umc_process** - Normalize U/M/C results, produce FASTA/CSV inventories.
8. **cd_hit** - Remove redundant sequences at 90% identity.
9. **ecc_dedup** - Harmonize coordinates and collapse duplicates.
10. **read_filter** - Filter confirmed eccDNA reads into an exportable FASTA.

### 3.3 Inference (Steps 11-13)

11. **minimap2** - Ensure the reference `.mmi` index exists (create it if needed).
12. **ecc_inference** - Prefer Cresil; fall back to Cyrcular. Missing `.fai` triggers a `samtools faidx` attempt.
13. **iecc_curator** - Clean and standardize inferred outputs.

### 3.4 Integration (Steps 14-16)

14. **ecc_unify** - Merge confirmed and inferred results, detecting redundant chimeric eccDNA. v0.9.4 introduces segment-wise overlap detection:
    - Uses 99% reciprocal overlap threshold with +/-10bp coordinate tolerance
    - Compares chimeric structures segment by segment, reducing false negatives from minor coordinate differences
    - Generates unified CSV and overlap metrics

15. **ecc_summary** - Aggregate all FASTA into `<prefix>_all.fasta` and create HTML/TXT reports.
16. **ecc_packager** - Copy and rename deliverables into the final customer-facing structure, respecting `--keep-tmp`.

---

## 4. Operational Notes

- **External dependencies** Ensure `tidehunter`, `blastn`, `cd-hit-est`, `minimap2`, `samtools`, `cresil` (or `cyrcular`) are available in `PATH`. v0.9.4 adds a pre-flight dependency check that validates all required tools before pipeline execution, with clear error messages and installation hints.
- **Checkpointing** Each step updates `<prefix>.checkpoint`. Use `--resume` to continue or `--force` to rerun from scratch.
- **Configuration** Parameters are managed through `circleseeker.config.Config`; command-line arguments override YAML.
- **Debug tooling** `circleseeker --debug --show-steps` reveals progress; `show-checkpoint` lists checkpoint details.

---

## 5. Further Reading

- `src/circleseeker/core/pipeline.py` - Pipeline orchestration and state management
- `src/circleseeker/config.py` - Configuration schema and defaults
- `src/circleseeker/modules/` - Implementation of internal modules
- `tests/unit/` - Unit tests demonstrating expected behaviors

If documentation diverges from actual behavior, please open an issue or submit a pull request.
