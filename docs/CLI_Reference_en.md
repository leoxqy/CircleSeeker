# CircleSeeker CLI Reference (v0.9.5)

This document describes the CircleSeeker 0.9.5 command-line interface for English-speaking users.

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
- `-t, --threads INT` Number of threads (default 8)
- `-c, --config PATH` Configuration file (YAML format)
- `--keep-tmp / --no-keep-tmp` Retain or remove the temporary `.tmp_work` directory (default: remove); explicitly overrides the `keep_tmp` setting in config file

---

## 2. Logging & Debugging

- `-n, --noise` Increase log verbosity (once for INFO, twice for DEBUG)
- `--debug` Enable debug mode and expose advanced options
- `-h, --help` Show help and exit
- `-v, --version` Print CircleSeeker version

---

## 3. Advanced Options (Debug Mode Only)

These flags are available only when `--debug` is present:

- `--start-from INT` Begin at a specific step (1-based index)
- `--stop-at INT` Stop after the specified step (inclusive)
- `--resume` Resume from the last checkpoint
- `--force` Ignore checkpoints and rerun the full pipeline
- `--generate-config` Print default configuration YAML and exit
- `--show-steps` List all 16 steps without executing them
- `--dry-run` Show planned operations without running
- `--log-output PATH` Write logs to an additional file

---

## 4. Hidden Subcommands

Visible only when `--debug` is set:

| Subcommand | Description |
|------------|-------------|
| `run` | Legacy-compatible entry point |
| `init-config` | Write default configuration to disk |
| `show-checkpoint` | Inspect checkpoint information |
| `validate` | Verify installation and dependencies |

```bash
circleseeker --debug init-config -o config.yaml
circleseeker --debug show-checkpoint -o results/ -p sample
```

---

## 5. Execution Lifecycle

1. **Dependency check** At startup, CircleSeeker verifies that all required external tools (minimap2, samtools, blastn, makeblastdb, cd-hit-est) and at least one inference engine (cresil or cyrcular) are available. If dependencies are missing, clear error messages with installation hints are shown.
2. **Temporary workspace** Intermediates live under `<output>/.tmp_work/` (configurable via `runtime.tmp_dir`, supports relative or absolute paths). Use `--keep-tmp` to retain them.
3. **Config & checkpoints** During execution, `config.yaml` and `<prefix>.checkpoint` are saved in the output directory; they are automatically cleaned up on successful completion. Use `--keep-tmp` to preserve them for debugging or resuming interrupted runs.
4. **Auto indexing** Missing `.mmi` or `.fai` triggers `minimap2 -d` and `samtools faidx` automatically.
5. **Packaging** `ecc_packager` copies/renames deliverables into the final structure documented in the README.

---

## 6. Step Overview

| # | Name | Purpose |
|---|------|---------|
| 1 | `make_blastdb` | Build BLAST database |
| 2 | `tidehunter` | Detect tandem repeats |
| 3 | `tandem_to_ring` | Convert repeats to circular candidates |
| 4 | `run_blast` | Align candidates back to the reference |
| 5 | `um_classify` | Classify into UeccDNA / MeccDNA |
| 6 | `cecc_build` | Assemble complex eccDNA |
| 7 | `umc_process` | Consolidate U/M/C outputs |
| 8 | `cd_hit` | Remove redundant sequences |
| 9 | `ecc_dedup` | Harmonize coordinates |
|10 | `read_filter` | Filter confirmed eccDNA reads |
|11 | `minimap2` | Prepare reference index / alignments |
|12 | `ecc_inference` | Cresil inference (Cyrcular fallback) |
|13 | `iecc_curator` | Curate inferred eccDNA |
|14 | `ecc_unify` | Merge confirmed and inferred tables using segment overlap algorithm for chimeric redundancy detection |
|15 | `ecc_summary` | Generate statistics and reports |
|16 | `ecc_packager` | Package final deliverables |

---

## 7. Output Layout

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

## 8. Example Commands

```bash
# Basic run
circleseeker -i sample.hifi.fasta -r hg38.fa -o results/ -p sample

# Use YAML overrides and keep intermediates
circleseeker -i sample.hifi.fasta -r hg38.fa -c configs/sample.yaml --keep-tmp

# Debug utilities
circleseeker --debug --show-steps
circleseeker --debug --resume -i sample.hifi.fasta -r hg38.fa -o results/

# Generate default config
circleseeker --debug --generate-config > default_config.yaml
```

---

## 9. FAQ

- **"Requires --debug" errors** Add `--debug` to unlock advanced options.
- **Dependency check fails** The error message indicates which tools are missing. Install via conda: `conda install -c bioconda -c conda-forge <tool_name>`.
- **Missing indexes** Automatic generation should cover `.mmi`/`.fai`; otherwise run `minimap2 -d ...` or `samtools faidx ...`.
- **Resuming** Ensure `<prefix>.checkpoint` exists, then rerun with `--resume`.

---

See `docs/Pipeline_Modules_en.md` and the README for deeper context. Contributions are welcome if the docs drift from actual behavior.
