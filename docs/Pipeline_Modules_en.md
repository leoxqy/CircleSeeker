# CircleSeeker Pipeline Modules (v1.1.0)

This document summarizes the 16-step CircleSeeker pipeline, organized into 5 phases, for English-speaking users.

---

## 1. Overview

CircleSeeker organizes the workflow into five phases:

| Phase | Steps | Purpose |
|-------|-------|---------|
| **Preprocessing** | 1-3 | Dependency checks, tandem repeat detection, circular candidate extraction |
| **CtcReads-Caller** | 4-9 | Alignment, U/M/C classification, deduplication — producing **Confirmed** eccDNA |
| **SplitReads-Caller** | 10-13 | Read filtering, split-read inference — producing **Inferred** eccDNA |
| **Integration** | 14-15 | Merge results and generate summary report |
| **Packaging** | 16 | Organize final output directory |

From an evidence perspective, the same 16 steps can be understood as two callers plus preprocessing/integration/packaging:

- **CtcReads**: reads carrying **Ctc** (**C**oncatemeric **t**andem **c**opies) signals (tracked as CtcR-* classes in `tandem_to_ring.csv`).
- **CtcReads-Caller** (Steps 4–9): produces **Confirmed** U/M/C eccDNA from CtcReads evidence.
- **SplitReads-Caller** (Steps 10–13): infers eccDNA from split-read/junction evidence using built-in SplitReads-Core and produces **Inferred** eccDNA.
- **Preprocessing** (Steps 1–3): prepares inputs for both callers.
- **Integration** (Steps 14–15): de-redundancy merging and reporting.
- **Packaging** (Step 16): assembles the final deliverable directory.

Intermediates reside in `<output>/.tmp_work/`, and the final artifacts are copied into the packaged directory layout.

---

## 2. Step Summary

| # | Name | Phase | Type | Main Input | Main Output |
|---|------|-------|------|------------|-------------|
| 1 | check_dependencies | Preprocessing | Internal | Config/environment | Dependency report (fails fast) |
| 2 | tidehunter | Preprocessing | External | HiFi reads FASTA | Tandem repeat consensus |
| 3 | tandem_to_ring | Preprocessing | Internal | TideHunter output | Candidate FASTA/CSV |
| 4 | run_alignment | CtcReads-Caller | External (minimap2) | Candidates & reference | Alignment TSV |
| 5 | um_classify | CtcReads-Caller | Internal | Alignment TSV | `um_classify.uecc.csv`, `um_classify.mecc.csv`, `um_classify.unclassified.csv` |
| 6 | cecc_build | CtcReads-Caller | Internal | Unclassified alignments | `cecc_build.csv` |
| 7 | umc_process | CtcReads-Caller | Internal | U/M/C tables | Harmonized CSV & FASTA |
| 8 | cd_hit | CtcReads-Caller | External (CD-HIT) | FASTA | Non-redundant FASTA |
| 9 | ecc_dedup | CtcReads-Caller | Internal | CSV/FASTA | Deduplicated coordinates |
|10 | read_filter | SplitReads-Caller | Internal | Raw reads, tandem_to_ring.csv | Inference input FASTA |
|11 | minimap2 | SplitReads-Caller | External (minimap2) | Reference | `.mmi` index |
|12 | ecc_inference | SplitReads-Caller | Internal | SplitReads-Core | Inferred TSV/FASTA |
|13 | curate_inferred_ecc | SplitReads-Caller | Internal | Inferred results | Curated CSV/FASTA |
|14 | ecc_unify | Integration | Internal | Confirmed & inferred CSV | Unified table, overlap stats |
|15 | ecc_summary | Integration | Internal | Unified table, processed CSV | HTML/TXT summary |
|16 | ecc_packager | Packaging | Internal | Finalized files | End-user directory tree |

---

## 3. Phase Details

### 3.1 Preprocessing (Steps 1-3)

1. **check_dependencies** - Validate required external tools (SplitReads-Core is built-in).
2. **tidehunter** - Identify tandem repeats characteristic of rolling-circle amplification.
3. **tandem_to_ring** - Convert repeats into circular candidates via overlap graph analysis.

### 3.2 CtcReads-Caller (Steps 4-9)

4. **run_alignment** - Align candidates back to the reference using minimap2. Emits a BLAST outfmt 6-like TSV with an additional trailing `mapq` column and supports length-compensated identity thresholds.
5. **um_classify** - Categorize alignments into UeccDNA or MeccDNA using a ring-coverage + locus-clustering model with optional U-call MAPQ gating and secondary-mapping veto. Emits `um_classify.uecc.csv`, `um_classify.mecc.csv`, and `um_classify.unclassified.csv` (see `docs/UMC_Classification_Model.md`).
6. **cecc_build** - Run Cecc detection on `um_classify.unclassified.csv`. LAST is preferred for doubled-sequence repeat detection; if LAST or required inputs are missing, fall back to the graph-based method. Emits `cecc_build.csv`.
7. **umc_process** - Normalize U/M/C results and produce FASTA/CSV outputs; when `enable_xecc=true`, XeccDNA FASTA is also emitted.
8. **cd_hit** - Remove redundant sequences at 99% identity.
9. **ecc_dedup** - Harmonize coordinates and collapse duplicates.

### 3.3 SplitReads-Caller (Steps 10-13)

10. **read_filter** - Use Sieve with `tandem_to_ring.csv` to filter CtcR reads and build the inference FASTA (keeps all reads if CSV is missing).
11. **minimap2** - Ensure the reference `.mmi` index exists for SplitReads-Core.
12. **ecc_inference** - Run built-in SplitReads-Core for split-read based eccDNA detection. The algorithm has two phases: Trim (alignment and filtering) + Identify (graph-based circular detection). Missing `.fai` triggers `samtools faidx`, and failures write `<prefix>_inference_failed.txt`. See [`docs/SplitReads_Core_Algorithm_en.md`](SplitReads_Core_Algorithm_en.md) for algorithm details.
13. **curate_inferred_ecc** - Clean and standardize inferred outputs; generates FASTA when reference + pysam are available (`iecc_curator`). Output includes `num_split_reads` (supporting reads), `prob_present` (presence probability), and `hifi_abundance` (coverage estimate).

### 3.4 Integration (Steps 14-15)

14. **ecc_unify** - Merge confirmed and inferred results, detecting redundant chimeric eccDNA via segment-wise reciprocal overlap (default 0.99 with +/-10bp tolerance); generates unified CSV and overlap metrics.
15. **ecc_summary** - Aggregate all FASTA into `<prefix>_all.fasta` and create HTML/TXT reports.

### 3.5 Packaging (Step 16)

16. **ecc_packager** - Copy and rename deliverables into the final customer-facing structure, respecting `--keep-tmp`. XeccDNA outputs are not packaged; their source reads are folded into inference input.

---

## 4. Operational Notes

- **External dependencies** Ensure `tidehunter`, `cd-hit-est`, `minimap2`, `samtools`, and `lastal` are available in `PATH`. SplitReads-Core is built-in and requires only the `mappy` Python package.
- **Checkpointing** Each step updates `<prefix>.checkpoint`. Use `--resume` to continue or `--force` to rerun from scratch.
- **Configuration** Parameters are managed through `circleseeker.config.Config`; command-line arguments override YAML.
- **Debug tooling** `circleseeker --show-steps` reveals the 5-phase progress view; `show-checkpoint` lists checkpoint details.

---

## 5. Further Reading

- [`docs/UMC_Classification_Model_en.md`](UMC_Classification_Model_en.md) - U/M/C classification model and CeccBuild algorithm
- [`docs/SplitReads_Core_Algorithm_en.md`](SplitReads_Core_Algorithm_en.md) - SplitReads-Core inference algorithm details
- [`docs/Output_Format_Reference_en.md`](Output_Format_Reference_en.md) - Output file format reference
- `src/circleseeker/core/pipeline.py` - Pipeline orchestration and state management
- `src/circleseeker/config.py` - Configuration schema and defaults
- `src/circleseeker/modules/` - Implementation of internal modules
- `tests/unit/` - Unit tests demonstrating expected behaviors

If documentation diverges from actual behavior, please open an issue or submit a pull request.
