# CircleSeeker 输出格式参考手册

本手册定义 CircleSeeker v1.0.0 的主要输出文件结构与字段说明。

---

## 输出目录结构

`ecc_packager` 打包后的标准目录结构如下：

```
<output>/
└── <prefix>/
    ├── <prefix>_Confirmed_UeccDNA/
    │   ├── <prefix>_UeccDNA_C.fasta
    │   ├── <prefix>_UeccDNA.bed
    │   └── <prefix>_UeccDNA.core.csv
    ├── <prefix>_Confirmed_MeccDNA/
    │   ├── <prefix>_MeccDNA_C.fasta
    │   ├── <prefix>_MeccSites.bed
    │   ├── <prefix>_MeccBestSite.bed
    │   └── <prefix>_MeccSites.core.csv
    ├── <prefix>_Confirmed_CeccDNA/
    │   ├── <prefix>_CeccDNA_C.fasta
    │   ├── <prefix>_CeccSegments.bed
    │   ├── <prefix>_CeccJunctions.bedpe
    │   └── <prefix>_CeccSegments.core.csv
    ├── <prefix>_Inferred_eccDNA/
    │   ├── <prefix>_UeccDNA_I.csv
    │   ├── <prefix>_UeccDNA_I.fasta
    │   ├── <prefix>_chimeric.csv
    │   └── <prefix>_CeccDNA_I.fasta
    ├── <prefix>_merged_output.csv
    ├── <prefix>_report.html
    └── <prefix>_summary.txt
```

> 说明：`<prefix>_UeccDNA_I.csv` 在打包阶段由 `<prefix>_simple.csv` 重命名而来。

---

## 1. 合并输出（`*_merged_output.csv`）

包含所有 Confirmed + Inferred 的最终汇总表。

| 字段 | 类型 | 说明 |
|------|------|------|
| `eccDNA_id` | string | 最终编号（UeccDNA1/MeccDNA1/CeccDNA1） |
| `original_id` | string | 合并/重编号前的原始 ID |
| `Regions` | string | 基因组坐标（CeccDNA 为多段，用 `;` 分隔） |
| `Strand` | string | 链方向（CeccDNA 多段用 `;` 分隔） |
| `Length` | int | eccDNA 长度（bp） |
| `eccDNA_type` | string | `UeccDNA` / `MeccDNA` / `CeccDNA` |
| `State` | string | `Confirmed` / `Inferred` |
| `Seg_total` | int | 片段数（U/M 为 1） |
| `Hit_count` | int | 基因组命中数（Mecc 典型 >1） |

---

## 2. Confirmed 输出

### 2.1 UeccDNA

#### `*_UeccDNA.core.csv`

| 字段 | 说明 |
|------|------|
| `eccDNA_id` | 唯一 ID |
| `chr` / `start0` / `end0` / `strand` | 0-based 坐标与方向 |
| `length` | 序列长度（bp） |
| `match_degree` | 匹配度（0-100） |
| `confidence_score` / `query_cov_best` / `query_cov_2nd` | 证据打分与覆盖度（可能为空） |
| `mapq_best` / `identity_best` | 最佳比对质量（可能为空） |
| `low_mapq` / `low_identity` | 低质量标记（可能为空） |
| `copy_number` / `repeat_number` | 拷贝数/重复数估计 |
| `eccdna_type` | `Uecc` |
| `num_merged` / `merged_from_ids` | 合并信息 |
| `reads_count` / `read_name` | 支持 reads 数与列表 |

#### `*_UeccDNA.bed`

列顺序：`chrom, chromStart, chromEnd, name, score, strand, length, eccdna_type, repeat_number, match_degree`

#### `*_UeccDNA_C.fasta`

FASTA header 格式：
```
>UeccDNA1|chr:start-end(strand)|length=...|repeats=...|reads=...
```

---

### 2.2 MeccDNA

#### `*_MeccSites.core.csv`

在 Uecc 核心字段基础上追加：

| 字段 | 说明 |
|------|------|
| `hit_index` | 当前位点序号（1-based） |
| `hit_count` | 同一 eccDNA 命中总数 |

#### `*_MeccSites.bed`

列顺序：`chrom, chromStart, chromEnd, name, score, strand, length, eccdna_type, copy_number`

#### `*_MeccBestSite.bed`

列顺序（保留 legacy 字段名）：`chrom, chromStart, chromEnd, name, score, strand, eLength, eClass, copyNum`

#### `*_MeccDNA_C.fasta`

FASTA header 格式：
```
>MeccDNA1|multi_loci:{sites}_sites|length=...|copies=...|reads=...
```

---

### 2.3 CeccDNA

#### `*_CeccSegments.core.csv`

在 Uecc 核心字段基础上追加：

| 字段 | 说明 |
|------|------|
| `seg_index` / `seg_total` | 片段序号与总数 |
| `junction_role` | `head` / `middle` / `tail` |

#### `*_CeccSegments.bed`

列顺序：`chrom, chromStart, chromEnd, name, score, strand, length, eccdna_type, copy_number, match_degree`

#### `*_CeccJunctions.bedpe`

列顺序：`chrom1, start1, end1, chrom2, start2, end2, name, score, strand1, strand2`

#### `*_CeccDNA_C.fasta`

FASTA header 格式：
```
>CeccDNA1|segments:{seg_total}|junctions:{chr-chain}|length=...|copies=...|reads=...
```

---

## 3. Inferred 输出

### 3.1 `*_UeccDNA_I.csv`

字段：
`eccDNA_id, chr, start0, end0, strand, length, eccdna_type, state, num_split_reads, prob_present, prob_artifact, hifi_abundance`

### 3.2 `*_chimeric.csv`

字段：
`eccDNA_id, chr, start0, end0, strand, length, eccdna_type, state, seg_index, seg_total, junction_role, read_count, num_split_reads, prob_present, prob_artifact, hifi_abundance`

### 3.3 FASTA

- `*_UeccDNA_I.fasta`：`>ID|chr:start-end|length=...|type=simple`
- `*_CeccDNA_I.fasta`：`>ID|segments=N|length=...|type=chimeric|seg1:chr:start-end;...`

---

## 4. 统计与报告

- `*_summary.txt`：文本统计摘要
- `*_report.html`：HTML 报告
