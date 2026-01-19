# CircleSeeker Configuration Reference

This manual provides a complete list of all configuration options for CircleSeeker v1.0.0, including default values, valid ranges, and usage instructions.

---

## Configuration File Format

CircleSeeker uses YAML format configuration files. Specify via `-c` or `--config`:

```bash
circleseeker -i reads.fasta -r reference.fa -c config.yaml
```

Generate default configuration template:

```bash
circleseeker --debug --generate-config > config.yaml
```

> Note: CLI `--preset` (debug only) is applied after loading the config file but before CLI overrides.

---

## 1. Main Configuration Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_file` | Path | None | **Required**, input HiFi reads file path |
| `reference` | Path | None | **Required**, reference genome FASTA file path |
| `output_dir` | Path | `circleseeker_output` | Output directory path |
| `prefix` | string | `sample` | Output file prefix |
| `enable_xecc` | bool | `true` | Generate XeccDNA (unclassified rings) and append source reads to inference input |
| `threads` | int | `8` | Number of threads (same as `performance.threads`) |

### Skip Flags

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `skip_tidehunter` | bool | `false` | Skip TideHunter tandem repeat detection step |
| `skip_carousel` | bool | `false` | Skip tandem-to-ring conversion |
| `skip_organize` | bool | `false` | Skip result organization step |

---

## 2. Runtime Configuration (`runtime`)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `log_level` | string | `WARNING` | Log level: `DEBUG`, `INFO`, `WARNING`, `ERROR` |
| `log_file` | Path | None | Log file path (optional) |
| `tmp_dir` | Path | `.tmp_work` | Temporary directory; relative paths are placed under output directory |
| `keep_tmp` | bool | `false` | Keep temporary files directory |
| `checkpoint_policy` | string | `continue` | Policy on config changes: `continue`, `reset`, `fail` |
| `enable_progress` | bool | `true` | Enable tqdm progress display |

Notes:
- Relative `runtime.tmp_dir` must be a subdirectory (not `.` and no `..`).
- Use absolute paths if you want temp files outside the output directory.

---

## 3. Performance Configuration (`performance`)

| Parameter | Type | Default | Range | Description |
|-----------|------|---------|-------|-------------|
| `threads` | int | `8` | >= 1 | Number of parallel threads |

---

## 4. Tool Configuration (`tools`)

### 4.1 TideHunter (`tools.tidehunter`)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `k` | int | `16` | k-mer size |
| `w` | int | `1` | Window size |
| `p` | int | `100` | Minimum repeat unit length |
| `P` | int | `2000000` | Maximum repeat unit length |
| `e` | float | `0.1` | Maximum error rate |
| `f` | int | `2` | Output format (1=single-line, 2=detailed) |
| `c` | int | `2` | Minimum copy number |

---

### 4.2 tandem_to_ring (`tools.tandem_to_ring`)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `min_ave_match` | float | `99.0` | Minimum average match percentage |

---

### 4.3 U/M Classifier (`tools.um_classify`)

#### Core Threshold Parameters

| Parameter | Type | Default | Range | Description |
|-----------|------|---------|-------|-------------|
| `gap_threshold` | float | `10.0` | 0-100 | Alignment gap percentage threshold |
| `theta_full` | float | `0.95` | 0-1 | Full-length coverage threshold (fraction) |
| `theta_u` | float | `0.95` | 0-1 | UeccDNA coverage threshold |
| `theta_m` | float | `0.95` | 0-1 | MeccDNA coverage threshold |
| `theta_u2_max` | float | `0.05` | 0-1 | Max 2nd-locus coverage for UeccDNA (1.0 disables) |
| `theta_locus` | float | `0.95` | 0-1 | Reciprocal overlap threshold for locus clustering |
| `pos_tol_bp` | int | `50` | >= 0 | Position tolerance (bp) |

#### U Classification Auxiliary Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `mapq_u_min` | int | `0` | Minimum MAPQ for UeccDNA (0 disables) |
| `u_secondary_min_frac` | float | `0.01` | Secondary coverage fraction threshold |
| `u_secondary_min_bp` | int | `50` | Secondary coverage length threshold (bp) |
| `u_contig_gap_bp` | int | `1000` | Gap threshold within same contig (bp) |
| `u_secondary_max_ratio` | float | `0.05` | Secondary/primary coverage ratio threshold |
| `u_high_coverage_threshold` | float | `0.98` | High coverage threshold for relaxed veto |
| `u_high_mapq_threshold` | int | `50` | High MAPQ threshold for relaxed veto |

#### Mecc Filters (Optional)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `mapq_m_ambiguous_threshold` | int | `0` | If all loci MAPQ below this, avoid Mecc call (0 disables) |
| `mecc_identity_gap_threshold` | float | `0` | Reject Mecc when max-2nd identity gap exceeds threshold (0 disables) |

#### Ambiguity Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `delta_uc` | float | `0.05` | U/C boundary ambiguity tolerance |
| `epsilon_mc` | float | `0.05` | M/C boundary ambiguity tolerance |

#### Legacy Compatible Parameters (Deprecated but Accepted)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `min_full_length_coverage` | float | `95.0` | Legacy full-length coverage (percentage) |
| `max_identity_gap_for_mecc` | float | `5.0` | Legacy maximum identity gap for MeccDNA |

---

### 4.4 CeccDNA Builder (`tools.cecc_build`)

CeccBuild v4 prefers LAST; it falls back to a graph-based method if LAST or inputs are missing.

| Parameter | Type | Default | Range | Description |
|-----------|------|---------|-------|-------------|
| `overlap_threshold` | float | `0.95` | 0-1 | Segment overlap threshold |
| `min_segments` | int | `2` | >= 1 | Minimum number of segments |
| `edge_tolerance` | int | `20` | >= 1 | Query edge tolerance (bp) |
| `tau_gap` | int | `20` | >= 1 | Gap tolerance (bp) |
| `position_tolerance` | int | `50` | >= 0 | Position tolerance (bp) |
| `half_query_buffer` | int | `50` | >= 0 | Mid-point buffer for doubled sequences (bp) |
| `locus_overlap_threshold` | float | `0.95` | 0-1 | Reciprocal overlap threshold for loci |
| `theta_chain` | float | `0.95` | 0-1 | Chain coverage threshold (converted to `min_match_degree`) |
| `min_match_degree` | float | `95.0` | 0-100 | Minimum match degree (percent) |
| `max_rotations` | int | `20` | >= 1 | Maximum rotation count |

---

### 4.5 Alignment Config (`tools.alignment`)

Controls the aligner used by `run_alignment`.

> Note: The current version only supports minimap2 for `run_alignment`; LAST is only used for complex eccDNA detection in `cecc_build`.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `aligner` | string | `minimap2` | Currently only `minimap2` supported |
| `min_identity` | float | `99.0` | Minimum identity threshold (percent) |
| `min_alignment_length` | int | `50` | Minimum alignment length (bp) |
| `db_prefix` | Path | None | Prebuilt LAST DB prefix (used by cecc_build only) |

---

### 4.6 Minimap2 Candidate Alignment (`tools.minimap2_align`)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `preset` | string | `map-hifi` | minimap2 preset |
| `max_target_seqs` | int | `5` | Max target sequences |
| `additional_args` | string | `""` | Extra minimap2 args |
| `min_identity` | float | `99.0` | Base identity threshold (percent) |
| `identity_decay_per_10kb` | float | `0.5` | Identity decay per 10kb |
| `min_identity_floor` | float | `97.0` | Identity floor (percent) |
| `split_by_length` | bool | `false` | Split presets by length |
| `split_length` | int | `5000` | Length threshold (bp) |
| `preset_short` | string | `map-hifi` | Preset for short sequences |
| `preset_long` | string | `map-hifi` | Preset for long sequences |

---

### 4.7 Minimap2 Inference Mapping (`tools.minimap2`)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `preset` | string | `map-hifi` | minimap2 preset |
| `additional_args` | string | `""` | Extra minimap2 args |

---

### 4.8 Samtools (`tools.samtools`)

Samtools configuration (reserved).

```yaml
tools:
  samtools: {}
```

---

## 5. Full Configuration Example

```yaml
input_file: /path/to/reads.fasta
reference: /path/to/reference.fa
output_dir: circleseeker_output
prefix: sample

enable_xecc: true
skip_tidehunter: false
skip_carousel: false
skip_organize: false

runtime:
  log_level: INFO
  log_file: null
  tmp_dir: .tmp_work
  keep_tmp: false
  checkpoint_policy: continue
  enable_progress: true

performance:
  threads: 8

tools:
  tidehunter:
    k: 16
    w: 1
    p: 100
    P: 2000000
    e: 0.1
    f: 2
    c: 2

  tandem_to_ring:
    min_ave_match: 99.0

  um_classify:
    gap_threshold: 10.0
    theta_full: 0.95
    theta_u: 0.95
    theta_m: 0.95
    theta_u2_max: 0.05
    mapq_u_min: 0
    u_secondary_min_frac: 0.01
    u_secondary_min_bp: 50
    u_contig_gap_bp: 1000
    u_secondary_max_ratio: 0.05
    u_high_coverage_threshold: 0.98
    u_high_mapq_threshold: 50
    theta_locus: 0.95
    pos_tol_bp: 50
    delta_uc: 0.05
    epsilon_mc: 0.05

  cecc_build:
    overlap_threshold: 0.95
    min_segments: 2
    edge_tolerance: 20
    tau_gap: 20
    position_tolerance: 50
    half_query_buffer: 50
    locus_overlap_threshold: 0.95
    theta_chain: 0.95
    min_match_degree: 95.0
    max_rotations: 20

  alignment:
    aligner: minimap2
    min_identity: 99.0
    min_alignment_length: 50
    db_prefix: null

  minimap2_align:
    preset: map-hifi
    max_target_seqs: 5
    additional_args: ""
    min_identity: 99.0
    identity_decay_per_10kb: 0.5
    min_identity_floor: 97.0
    split_by_length: false
    split_length: 5000
    preset_short: map-hifi
    preset_long: map-hifi

  minimap2:
    preset: map-hifi
    additional_args: ""

  samtools: {}
```

---

## 6. Configuration Priority

Priority order (high to low):

1. **CLI arguments** (e.g., `-t 16` overrides config)
2. **Preset** (`--preset` applied after config load)
3. **Config file** (`-c config.yaml`)
4. **Code defaults**

Notes:
- `tools` uses a **deep-merge** strategy: only provided keys override defaults.
- `--keep-tmp` overrides `runtime.keep_tmp`.
- Relative `runtime.tmp_dir` is placed under `output_dir`.

---

## 7. Related Docs

- [CLI Reference](CLI_Reference_en.md) - Command-line options
- [Pipeline Modules](Pipeline_Modules_en.md) - 16-step workflow
- [Output Format Reference](Output_Format_Reference_en.md) - Output file schemas
