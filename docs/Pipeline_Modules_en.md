# CircleSeeker Pipeline Modules (v1.0.0)

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

From an evidence perspective, the same 16 steps can be understood as two callers plus a shared integration stage:

- **CtcReads**: reads carrying **Ctc** (**C**oncatemeric **t**andem **c**opies) signals (tracked as CtcR-* classes in `tandem_to_ring.csv`).
- **CtcReads-Caller** (Steps 1–10): produces **Confirmed** U/M/C eccDNA from CtcReads evidence.
- **SplitReads-Caller** (Steps 11–13): infers eccDNA from split-read/junction evidence (Cresil preferred, Cyrcular fallback) and produces **Inferred** eccDNA.
- **Integration** (Steps 14–16): de-redundancy, merging, reporting, and packaging for delivery.

Intermediates reside in `<output>/.tmp_work/`, and the final artifacts are copied into the packaged directory layout.

---

## 2. Step Summary

| # | Name | Type | Main Input | Main Output |
|---|------|------|------------|-------------|
| 1 | check_dependencies | Internal | Config/environment | Dependency report (fails fast) |
| 2 | tidehunter | External | HiFi reads FASTA | Tandem repeat consensus |
| 3 | tandem_to_ring | Internal | TideHunter output | Candidate FASTA/CSV |
| 4 | run_alignment | External (minimap2/LAST) | Candidates & reference | Alignment TSV |
| 5 | um_classify | Internal | Alignment TSV | `um_classify.uecc.csv`, `um_classify.mecc.csv`, `um_classify.unclassified.csv` |
| 6 | cecc_build | Internal | Unclassified alignments | `cecc_build.csv` |
| 7 | umc_process | Internal | U/M/C tables | Harmonized CSV & FASTA |
| 8 | cd_hit | External (CD-HIT) | FASTA | Non-redundant FASTA |
| 9 | ecc_dedup | Internal | CSV/FASTA | Deduplicated coordinates |
|10 | read_filter | Internal | Raw reads, tandem_to_ring.csv | Inference input FASTA |
|11 | minimap2 | External (minimap2) | Reference | `.mmi` index |
|12 | ecc_inference | External + Internal | Cresil/Cyrcular inputs | Inferred TSV/FASTA |
|13 | curate_inferred_ecc | Internal | Inferred results | Curated CSV/FASTA |
|14 | ecc_unify | Internal | Confirmed & inferred CSV | Unified table, overlap stats |
|15 | ecc_summary | Internal | Unified table, processed CSV | HTML/TXT summary, all-fasta |
|16 | ecc_packager | Internal | Finalized files | End-user directory tree |

---

## 3. Stage Details

### 3.1 Detection (CtcReads-Caller, Steps 1-6)

1. **check_dependencies** - Validate required external tools and at least one inference engine (Cresil or Cyrcular).
2. **tidehunter** - Identify tandem repeats characteristic of rolling-circle amplification.
3. **tandem_to_ring** - Convert repeats into circular candidates via overlap graph analysis.
4. **run_alignment** - Align candidates back to the reference using minimap2 (default) or LAST (`tools.alignment.aligner`). Emits a BLAST outfmt 6-like TSV with an additional trailing `mapq` column and supports length-compensated identity thresholds.
5. **um_classify** - Categorize alignments into UeccDNA or MeccDNA using a ring-coverage + locus-clustering model with optional U-call MAPQ gating and secondary-mapping veto. Emits `um_classify.uecc.csv`, `um_classify.mecc.csv`, and `um_classify.unclassified.csv` (see `docs/UMC_Classification_Model.md`).
6. **cecc_build** - Run Cecc detection on `um_classify.unclassified.csv`. LAST is preferred for doubled-sequence repeat detection; if LAST or required inputs are missing, fall back to the graph-based method. Emits `cecc_build.csv`.

### 3.2 Processing (CtcReads-Caller, Steps 7-10)

7. **umc_process** - Normalize U/M/C results and produce FASTA/CSV outputs; when `enable_xecc=true`, XeccDNA FASTA is also emitted.
8. **cd_hit** - Remove redundant sequences at 99% identity.
9. **ecc_dedup** - Harmonize coordinates and collapse duplicates.
10. **read_filter** - Use Sieve with `tandem_to_ring.csv` to filter CtcR reads and build the inference FASTA (keeps all reads if CSV is missing).

### 3.3 Inference (SplitReads-Caller, Steps 11-13)

11. **minimap2** - Ensure the reference `.mmi` index exists; BAM alignment is skipped when Cresil is selected.
12. **ecc_inference** - Prefer Cresil; fall back to Cyrcular on failure. Missing `.fai` triggers `samtools faidx`, and failures write `<prefix>_inference_failed.txt`.
13. **curate_inferred_ecc** - Clean and standardize inferred outputs; generates FASTA when reference + pysam are available (`iecc_curator`).

### 3.4 Integration (Integration, Steps 14-16)

14. **ecc_unify** - Merge confirmed and inferred results, detecting redundant chimeric eccDNA via segment-wise reciprocal overlap (default 0.99 with +/-10bp tolerance); generates unified CSV and overlap metrics.

15. **ecc_summary** - Aggregate all FASTA into `<prefix>_all.fasta` and create HTML/TXT reports.
16. **ecc_packager** - Copy and rename deliverables into the final customer-facing structure, respecting `--keep-tmp`. XeccDNA outputs are not packaged; their source reads are folded into inference input.

---

## 4. Operational Notes

- **External dependencies** Ensure `tidehunter`, `cd-hit-est`, `minimap2`, `samtools`, `cresil` (or `cyrcular`) are available in `PATH`. If you enable LAST (`tools.alignment.aligner=last` or CeccBuild v4), install LAST (`lastdb`/`lastal`).
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
