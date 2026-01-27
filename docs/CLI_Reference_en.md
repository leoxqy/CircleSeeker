# CircleSeeker CLI Reference (v1.1.0)

This document describes the CircleSeeker 1.1.0 command-line interface for English-speaking users.

---

## 1. Quick Start

```bash
circleseeker -i reads.fasta -r reference.fa -o results/
```

Required parameters:

- `-i, --input PATH` HiFi reads FASTA
- `-r, --reference PATH` Reference genome FASTA

Common optional flags:

- `-o, --output PATH` Output directory (`circleseeker_output` by default)
- `-p, --prefix TEXT` Filename prefix (`sample` by default)
- `-t, --threads INT` Number of threads (default `min(8, CPU_count×2)`)
- `-c, --config PATH` Configuration file (YAML format)
- `--keep-tmp` Retain the temporary `.tmp_work` directory (default: remove); overrides `keep_tmp` in config
- `--turbo` Enable turbo mode (RAM-backed temp directory for faster I/O)
- `--preset CHOICE` Sensitivity preset (`relaxed` / `balanced` / `strict`)
- `--show-steps` List all 16 steps without executing them
- `--dry-run` Show planned operations without running

---

## 2. Logging & Debugging

- `-v, --verbose` Increase log verbosity (-v for INFO, -vv for DEBUG; default WARNING)
- `--debug` Unlock advanced options and hidden subcommands (does not affect log level)
- `--help-advanced` Show advanced/debug options and exit
- `-h, --help` Show help and exit
- `-V, --version` Print CircleSeeker version

---

## 3. Advanced Options (Debug Mode Only)

These flags are available only when `--debug` is present:

- `--start-from INT` Begin at a specific step (1-based index)
- `--stop-at INT` Stop after the specified step (inclusive)
- `--resume` Resume from the last checkpoint
- `--force` Ignore checkpoints and rerun the full pipeline
- `--log-file PATH` Write logs to an additional file

---

## 4. Subcommands

| Subcommand | Description |
|------------|-------------|
| `init-config` | Generate default configuration (`--stdout` for terminal output, `--output-file` for file) |
| `show-checkpoint` | Inspect checkpoint information (`-d` to specify output directory) |
| `validate` | Verify installation and dependencies (`--full` for extended checks) |

```bash
circleseeker init-config --stdout
circleseeker init-config --output-file config.yaml
circleseeker show-checkpoint -d results/ -p sample
circleseeker validate
```

---

## 5. Execution Lifecycle

1. **Dependency check** At startup, CircleSeeker verifies that all required external tools (minimap2, samtools, cd-hit-est, lastal) are available. SplitReads-Core is built-in and requires only the `mappy` Python package. If dependencies are missing, clear error messages with installation hints are shown.
2. **Temporary workspace** Intermediates live under `<output>/.tmp_work/` (configurable via `runtime.tmp_dir`, supports relative or absolute paths). Use `--keep-tmp` to retain them.
3. **Config & checkpoints** During execution, `config.yaml` and `<prefix>.checkpoint` are saved in the output directory; they are automatically cleaned up on successful completion. Use `--keep-tmp` to preserve them for debugging or resuming interrupted runs.
4. **Auto indexing** Missing `.mmi` or `.fai` triggers `minimap2 -d` and `samtools faidx` automatically.
5. **Packaging** `ecc_packager` copies/renames deliverables into the final structure documented in the README.

---

## 6. Evidence-Driven Callers

For clarity, CircleSeeker describes the workflow as two evidence-driven callers:

- **CtcReads**: reads carrying **Ctc** (**C**oncatemeric **t**andem **c**opies) signals (tracked as CtcR-* classes in `tandem_to_ring.csv`).
- **CtcReads-Caller** (Steps 1–10): produces **Confirmed** U/M/C eccDNA from CtcReads evidence.
- **SplitReads-Caller** (Steps 11–13): infers eccDNA from split-read/junction evidence using built-in SplitReads-Core and produces **Inferred** eccDNA.
- **Integration** (Steps 14–16): de-redundancy, merging, reporting, and packaging for delivery.

## 7. Step Overview

> Quick mapping: Steps 1–10 are **CtcReads-Caller**; Steps 11–13 are **SplitReads-Caller**; Steps 14–16 are **Integration**.

| # | Name | Purpose |
|---|------|---------|
| 1 | `check_dependencies` | Verify required tools and inference engine |
| 2 | `tidehunter` | Detect tandem repeats |
| 3 | `tandem_to_ring` | Convert repeats to circular candidates |
| 4 | `run_alignment` | Align candidates back to the reference (minimap2 or LAST) |
| 5 | `um_classify` | Coverage + locus-based U/Mecc classification |
| 6 | `cecc_build` | LAST-first complex eccDNA detection (graph fallback) |
| 7 | `umc_process` | Consolidate U/M/C outputs |
| 8 | `cd_hit` | Remove redundant sequences |
| 9 | `ecc_dedup` | Harmonize coordinates |
|10 | `read_filter` | Filter CtcR reads and build inference FASTA |
|11 | `minimap2` | Build reference index for SplitReads-Core |
|12 | `ecc_inference` | SplitReads-Core inference (built-in) |
|13 | `curate_inferred_ecc` | Curate inferred eccDNA tables (`iecc_curator`) |
|14 | `ecc_unify` | Merge confirmed and inferred tables using segment overlap algorithm for chimeric redundancy detection |
|15 | `ecc_summary` | Generate statistics and reports |
|16 | `ecc_packager` | Package final deliverables |

---

## 8. Output Layout

```
<output>/
├── <prefix>_merged_output.csv
├── <prefix>_summary.txt
├── <prefix>_report.html
├── <prefix>_Confirmed_UeccDNA/
├── <prefix>_Confirmed_MeccDNA/
├── <prefix>_Confirmed_CeccDNA/
├── <prefix>_Inferred_eccDNA/
└── logs/ (if additional logging enabled)
```

Subdirectories follow the naming described in the README "Output Files" section.

---

## 9. Example Commands

```bash
# Basic run
circleseeker -i sample.hifi.fasta -r hg38.fa -o results/ -p sample

# Use YAML overrides and keep intermediates
circleseeker -i sample.hifi.fasta -r hg38.fa -c configs/sample.yaml --keep-tmp

# View pipeline steps (no execution)
circleseeker --show-steps

# Dry run (preview)
circleseeker --dry-run -i sample.hifi.fasta -r hg38.fa -o results/

# Sensitivity preset
circleseeker --preset strict -i sample.hifi.fasta -r hg38.fa -o results/

# Debug: resume from checkpoint
circleseeker --debug --resume -i sample.hifi.fasta -r hg38.fa -o results/

# Generate default config
circleseeker init-config --stdout > default_config.yaml
circleseeker init-config --output-file config.yaml
```

---

## 10. FAQ

- **"Requires --debug" errors** Only `--start-from`, `--stop-at`, `--resume`, `--force`, and `--log-file` require `--debug`.
- **Dependency check fails** The error message indicates which tools are missing. Install via conda: `conda install -c bioconda -c conda-forge <tool_name>`.
- **Missing indexes** Automatic generation should cover `.mmi`/`.fai`; otherwise run `minimap2 -d ...` or `samtools faidx ...`.
- **Resuming** Ensure `<prefix>.checkpoint` exists, then rerun with `--debug --resume`.

---

See `docs/Pipeline_Modules_en.md` and the README for deeper context. Contributions are welcome if the docs drift from actual behavior.
