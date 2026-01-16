# CircleSeeker 配置参考手册

本手册完整列出 CircleSeeker v0.10.3 所有配置选项，包括默认值、取值范围及使用说明。

---

## 配置文件格式

CircleSeeker 使用 YAML 格式的配置文件。通过 `-c` 或 `--config` 参数指定：

```bash
circleseeker -i reads.fasta -r reference.fa -c config.yaml
```

生成默认配置模板：

```bash
circleseeker --debug --generate-config > config.yaml
```

---

## 1. 主配置参数

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `input_file` | Path | 无 | **必选**，输入 HiFi reads 文件路径 |
| `reference` | Path | 无 | **必选**，参考基因组 FASTA 文件路径 |
| `output_dir` | Path | `circleseeker_output` | 输出目录路径 |
| `prefix` | string | `sample` | 输出文件前缀 |
| `enable_xecc` | bool | `true` | 是否启用扩展 eccDNA 检测功能 |
| `threads` | int | `8` | 线程数（同 `performance.threads`） |

### 跳过标志

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `skip_tidehunter` | bool | `false` | 跳过 TideHunter 串联重复检测步骤 |
| `skip_carousel` | bool | `false` | 跳过串联转环状序列步骤（别名：`skip_tandem_to_ring`） |
| `skip_organize` | bool | `false` | 跳过结果组织步骤 |

---

## 2. 运行时配置 (`runtime`)

控制程序运行时行为的配置项。

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `log_level` | string | `WARNING` | 日志级别，可选：`DEBUG`、`INFO`、`WARNING`、`ERROR` |
| `log_file` | Path | 无 | 日志文件路径（可选） |
| `tmp_dir` | Path | `.tmp_work` | 临时文件目录；相对路径会置于输出目录下 |
| `keep_tmp` | bool | `false` | 是否保留临时文件目录 |
| `checkpoint_policy` | string | `continue` | 配置变更时的检查点策略：`continue`（继续）、`reset`（重置）、`fail`（失败） |
| `enable_progress` | bool | `true` | 是否启用 tqdm 进度条显示 |

### 示例

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

## 3. 性能配置 (`performance`)

控制资源使用的配置项。

| 参数 | 类型 | 默认值 | 取值范围 | 说明 |
|------|------|--------|----------|------|
| `threads` | int | `8` | >= 1 | 并行线程数 |

### 示例

```yaml
performance:
  threads: 16
```

---

## 4. 工具配置 (`tools`)

外部工具和内部模块的参数配置。

### 4.1 TideHunter (`tools.tidehunter`)

串联重复检测工具的参数配置。

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `k` | int | `16` | k-mer 大小 |
| `w` | int | `1` | 窗口大小 |
| `p` | int | `100` | 最小重复单元长度 |
| `P` | int | `2000000` | 最大重复单元长度 |
| `e` | float | `0.1` | 最大错误率 |
| `f` | int | `2` | 最小重复次数 |

### 示例

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

### 4.2 U/M 分类器 (`tools.um_classify`)

UeccDNA 与 MeccDNA 分类模块的参数配置。

#### 核心阈值参数

| 参数 | 类型 | 默认值 | 取值范围 | 说明 |
|------|------|--------|----------|------|
| `gap_threshold` | float | `10.0` | 0-100 | 比对间隙百分比阈值 |
| `theta_full` | float | `0.95` | 0-1 | 全长覆盖率阈值（分数形式） |
| `theta_u` | float | `0.95` | 0-1 | UeccDNA 覆盖率阈值 |
| `theta_m` | float | `0.95` | 0-1 | MeccDNA 覆盖率阈值 |
| `theta_u2_max` | float | `0.05` | 0-1 | UeccDNA 次优位点最大覆盖率（设为 1.0 禁用） |
| `theta_locus` | float | `0.95` | 0-1 | 位点覆盖率阈值 |
| `pos_tol_bp` | int | `50` | >= 0 | 位置容差（bp） |

#### U 分类辅助参数

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `mapq_u_min` | int | `0` | UeccDNA 最小 MAPQ 阈值（0 表示禁用） |
| `u_secondary_min_frac` | float | `0.01` | 次级比对最小比例阈值 |
| `u_secondary_min_bp` | int | `50` | 次级比对最小长度（bp） |
| `u_contig_gap_bp` | int | `1000` | 同一 contig 内的间隙阈值（bp） |

#### 模糊判定参数

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `delta_uc` | float | `0.05` | U/C 边界模糊容差 |
| `epsilon_mc` | float | `0.05` | M/C 边界模糊容差 |

#### 旧版兼容参数（已弃用但仍接受）

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `min_full_length_coverage` | float | `95.0` | 旧版全长覆盖率（百分比形式） |
| `max_identity_gap_for_mecc` | float | `5.0` | 旧版 MeccDNA 最大同一性间隙 |

### 示例

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

### 4.3 CeccDNA 构建器 (`tools.cecc_build`)

复杂 eccDNA (CeccDNA) 检测模块的参数配置。

| 参数 | 类型 | 默认值 | 取值范围 | 说明 |
|------|------|--------|----------|------|
| `overlap_threshold` | float | `0.95` | 0-1 | 片段重叠阈值 |
| `min_segments` | int | `2` | >= 2 | 最小片段数 |
| `edge_tolerance` | int | `20` | >= 0 | query 端边缘容差（bp） |
| `tau_gap` | int | `20` | >= 0 | 间隙容差（bp） |
| `position_tolerance` | int | `50` | >= 0 | 位置容差（bp）；用于旧版闭环检查 |
| `locus_overlap_threshold` | float | `0.95` | 0-1 | 基因组区间互反重叠阈值 |
| `theta_chain` | float | `0.95` | 0-1 | 链覆盖率阈值 |
| `max_rotations` | int | `20` | >= 1 | 最大旋转次数 |

#### 旧版兼容参数

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `min_match_degree` | float | `95.0` | 旧版最小匹配度（百分比形式） |

### 示例

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

### 4.4 Minimap2 短序列比对 (`tools.minimap2_align`)

用于候选片段比对到参考基因组的 minimap2 配置。

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `preset` | string | `sr` | minimap2 预设模式，短序列推荐 `sr` |
| `max_target_seqs` | int | `200` | 最大目标序列数 |
| `additional_args` | string | `""` | 额外的 minimap2 参数 |

### 示例

```yaml
tools:
  minimap2_align:
    preset: sr
    max_target_seqs: 200
    additional_args: "-x asm5"
```

---

### 4.5 Minimap2 索引/比对 (`tools.minimap2`)

用于构建参考索引和推断阶段比对的 minimap2 配置。

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `preset` | string | `map-hifi` | minimap2 预设模式，HiFi reads 推荐 `map-hifi` |
| `additional_args` | string | `""` | 额外的 minimap2 参数 |

### 示例

```yaml
tools:
  minimap2:
    preset: map-hifi
    additional_args: ""
```

---

### 4.6 Samtools (`tools.samtools`)

Samtools 工具配置（目前为预留，暂无特定参数）。

```yaml
tools:
  samtools: {}
```

---

## 5. 完整配置示例

```yaml
# CircleSeeker 配置文件示例
input_file: /path/to/reads.fasta
reference: /path/to/reference.fa
output_dir: circleseeker_output
prefix: sample

# 功能开关
enable_xecc: true
skip_tidehunter: false
skip_carousel: false
skip_organize: false

# 运行时配置
runtime:
  log_level: INFO
  log_file: null
  tmp_dir: .tmp_work
  keep_tmp: false
  checkpoint_policy: continue
  enable_progress: true

# 性能配置
performance:
  threads: 8

# 工具配置
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

## 6. 配置优先级

配置参数的优先级从高到低为：

1. **命令行参数**：如 `-t 16` 会覆盖配置文件中的 `threads`
2. **配置文件**：通过 `-c config.yaml` 指定
3. **默认值**：代码中定义的默认值

### 注意事项

- 配置文件中的 `tools` 配置采用**深度合并**策略：只覆盖指定的参数，未指定的保持默认值
- `--keep-tmp` 命令行参数会覆盖配置文件中的 `runtime.keep_tmp`
- 临时目录 `tmp_dir` 如使用相对路径，会解析为 `output_dir` 的子目录

---

## 7. 相关文档

- [CLI 使用手册](CLI_Reference.md) - 命令行参数详解
- [管道模块说明](Pipeline_Modules.md) - 16 步流程详解
- [输出格式参考](Output_Format_Reference.md) - 输出文件格式定义
