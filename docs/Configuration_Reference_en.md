# CircleSeeker Configuration Reference

This manual provides a complete list of all configuration options for CircleSeeker v0.10.3, including default values, valid ranges, and usage instructions.

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

---

## 1. Main Configuration Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_file` | Path | None | **Required**, input HiFi reads file path |
| `reference` | Path | None | **Required**, reference genome FASTA file path |
| `output_dir` | Path | `circleseeker_output` | Output directory path |
| `prefix` | string | `sample` | Output file prefix |
| `enable_xecc` | bool | `true` | Enable extended eccDNA detection |
| `threads` | int | `8` | Number of threads (same as `performance.threads`) |

### Skip Flags

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `skip_tidehunter` | bool | `false` | Skip TideHunter tandem repeat detection step |
| `skip_carousel` | bool | `false` | Skip tandem-to-ring conversion step (alias: `skip_tandem_to_ring`) |
| `skip_organize` | bool | `false` | Skip result organization step |

---

## 2. Runtime Configuration (`runtime`)

Configuration options that control program runtime behavior.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `log_level` | string | `WARNING` | Log level: `DEBUG`, `INFO`, `WARNING`, `ERROR` |
| `log_file` | Path | None | Log file path (optional) |
| `tmp_dir` | Path | `.tmp_work` | Temporary files directory; relative paths are placed under output directory |
| `keep_tmp` | bool | `false` | Whether to keep temporary files directory |
| `checkpoint_policy` | string | `continue` | Checkpoint policy on config changes: `continue`, `reset`, `fail` |
| `enable_progress` | bool | `true` | Whether to enable tqdm progress bar display |

### Example

```yaml
runtime:
  log_level: INFO
  log_file: /path/to/circleseeker.log
  tmp_dir: .tmp_work
  keep_tmp: false
  checkpoint_policy: continue
  enable_progress: true
```

---

## 3. Performance Configuration (`performance`)

Configuration options that control resource usage.

| Parameter | Type | Default | Range | Description |
|-----------|------|---------|-------|-------------|
| `threads` | int | `8` | >= 1 | Number of parallel threads |

### Example

```yaml
performance:
  threads: 16
```

---

## 4. Tool Configuration (`tools`)

Parameter configuration for external tools and internal modules.

### 4.1 TideHunter (`tools.tidehunter`)

Parameter configuration for the tandem repeat detection tool.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `k` | int | `16` | k-mer size |
| `w` | int | `1` | Window size |
| `p` | int | `100` | Minimum repeat unit length |
| `P` | int | `2000000` | Maximum repeat unit length |
| `e` | float | `0.1` | Maximum error rate |
| `f` | int | `2` | Minimum repeat count |

### Example

```yaml
tools:
  tidehunter:
    k: 16
    w: 1
    p: 100
    P: 2000000
    e: 0.1
    f: 2
```

---

### 4.2 U/M Classifier (`tools.um_classify`)

Parameter configuration for the UeccDNA and MeccDNA classification module.

#### Core Threshold Parameters

| Parameter | Type | Default | Range | Description |
|-----------|------|---------|-------|-------------|
| `gap_threshold` | float | `10.0` | 0-100 | Alignment gap percentage threshold |
| `theta_full` | float | `0.95` | 0-1 | Full-length coverage threshold (fraction) |
| `theta_u` | float | `0.95` | 0-1 | UeccDNA coverage threshold |
| `theta_m` | float | `0.95` | 0-1 | MeccDNA coverage threshold |
| `theta_u2_max` | float | `0.05` | 0-1 | Maximum coverage for UeccDNA secondary locus (set to 1.0 to disable) |
| `theta_locus` | float | `0.95` | 0-1 | Locus coverage threshold |
| `pos_tol_bp` | int | `50` | >= 0 | Position tolerance (bp) |

#### U Classification Auxiliary Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `mapq_u_min` | int | `0` | Minimum MAPQ threshold for UeccDNA (0 disables) |
| `u_secondary_min_frac` | float | `0.01` | Minimum fraction threshold for secondary alignments |
| `u_secondary_min_bp` | int | `50` | Minimum length for secondary alignments (bp) |
| `u_contig_gap_bp` | int | `1000` | Gap threshold within the same contig (bp) |

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

### Example

```yaml
tools:
  um_classify:
    gap_threshold: 10.0
    theta_full: 0.95
    theta_u: 0.95
    theta_m: 0.95
    theta_u2_max: 0.05
    mapq_u_min: 20
    pos_tol_bp: 50
```

---

### 4.3 CeccDNA Builder (`tools.cecc_build`)

Parameter configuration for the complex eccDNA (CeccDNA) detection module.

| Parameter | Type | Default | Range | Description |
|-----------|------|---------|-------|-------------|
| `overlap_threshold` | float | `0.95` | 0-1 | Segment overlap threshold |
| `min_segments` | int | `2` | >= 2 | Minimum number of segments |
| `edge_tolerance` | int | `20` | >= 0 | Query edge tolerance (bp) |
| `tau_gap` | int | `20` | >= 0 | Gap tolerance (bp) |
| `position_tolerance` | int | `50` | >= 0 | Position tolerance (bp); used for legacy closure checks |
| `locus_overlap_threshold` | float | `0.95` | 0-1 | Reciprocal overlap threshold for genomic intervals |
| `theta_chain` | float | `0.95` | 0-1 | Chain coverage threshold |
| `max_rotations` | int | `20` | >= 1 | Maximum rotation count |

#### Legacy Compatible Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `min_match_degree` | float | `95.0` | Legacy minimum match degree (percentage) |

### Example

```yaml
tools:
  cecc_build:
    overlap_threshold: 0.95
    min_segments: 2
    edge_tolerance: 20
    theta_chain: 0.95
    max_rotations: 20
```

---

### 4.4 Minimap2 Short-read Alignment (`tools.minimap2_align`)

Minimap2 configuration for aligning candidate fragments to the reference genome.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `preset` | string | `sr` | minimap2 preset mode, `sr` recommended for short sequences |
| `max_target_seqs` | int | `200` | Maximum number of target sequences |
| `additional_args` | string | `""` | Additional minimap2 arguments |

### Example

```yaml
tools:
  minimap2_align:
    preset: sr
    max_target_seqs: 200
    additional_args: "-x asm5"
```

---

### 4.5 Minimap2 Index/Alignment (`tools.minimap2`)

Minimap2 configuration for building reference index and inference-stage alignment.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `preset` | string | `map-hifi` | minimap2 preset mode, `map-hifi` recommended for HiFi reads |
| `additional_args` | string | `""` | Additional minimap2 arguments |

### Example

```yaml
tools:
  minimap2:
    preset: map-hifi
    additional_args: ""
```

---

### 4.6 Samtools (`tools.samtools`)

Samtools tool configuration (currently reserved, no specific parameters).

```yaml
tools:
  samtools: {}
```

---

## 5. Complete Configuration Example

```yaml
# CircleSeeker Configuration File Example
input_file: /path/to/reads.fasta
reference: /path/to/reference.fa
output_dir: circleseeker_output
prefix: sample

# Feature flags
enable_xecc: true
skip_tidehunter: false
skip_carousel: false
skip_organize: false

# Runtime configuration
runtime:
  log_level: INFO
  log_file: null
  tmp_dir: .tmp_work
  keep_tmp: false
  checkpoint_policy: continue
  enable_progress: true

# Performance configuration
performance:
  threads: 8

# Tool configuration
tools:
  tidehunter:
    k: 16
    w: 1
    p: 100
    P: 2000000
    e: 0.1
    f: 2

  um_classify:
    gap_threshold: 10.0
    theta_full: 0.95
    theta_u: 0.95
    theta_m: 0.95
    theta_u2_max: 0.05
    mapq_u_min: 0
    pos_tol_bp: 50
    delta_uc: 0.05
    epsilon_mc: 0.05

  cecc_build:
    overlap_threshold: 0.95
    min_segments: 2
    edge_tolerance: 20
    tau_gap: 20
    theta_chain: 0.95
    max_rotations: 20

  minimap2_align:
    preset: sr
    max_target_seqs: 200
    additional_args: ""

  minimap2:
    preset: map-hifi
    additional_args: ""

  samtools: {}
```

---

## 6. Configuration Priority

Configuration parameter priority from highest to lowest:

1. **Command-line arguments**: e.g., `-t 16` overrides `threads` in config file
2. **Configuration file**: specified via `-c config.yaml`
3. **Default values**: default values defined in code

### Notes

- `tools` configuration uses **deep merge** strategy: only specified parameters are overridden, unspecified ones retain default values
- `--keep-tmp` command-line argument overrides `runtime.keep_tmp` in config file
- Temporary directory `tmp_dir` with relative path resolves to a subdirectory of `output_dir`

---

## 7. Related Documentation

- [CLI Reference](CLI_Reference_en.md) - Command-line argument details
- [Pipeline Modules](Pipeline_Modules_en.md) - 16-step pipeline details
- [Output Format Reference](Output_Format_Reference_en.md) - Output file format definitions
