# CircleSeeker 输出格式参考手册

本手册完整定义 CircleSeeker v0.10.3 所有输出文件的格式和字段说明。

---

## 输出目录结构

完成运行后，输出目录包含以下结构：

```
<output>/
├── <prefix>_merged_output.csv          # 所有 eccDNA 汇总表
├── <prefix>_report.html                # 交互式 HTML 报告
├── <prefix>_summary.txt                # 统计摘要
├── <prefix>_Confirmed_UeccDNA/         # 确认的简单 eccDNA
│   ├── <prefix>_UeccDNA_C.csv
│   └── <prefix>_UeccDNA_C.fasta
├── <prefix>_Confirmed_MeccDNA/         # 确认的多拷贝 eccDNA
│   ├── <prefix>_MeccDNA_C.csv
│   └── <prefix>_MeccDNA_C.fasta
├── <prefix>_Confirmed_CeccDNA/         # 确认的复杂 eccDNA
│   ├── <prefix>_CeccDNA_C.csv
│   ├── <prefix>_CeccDNA_C.fasta
│   ├── <prefix>_CeccJunctions.bedpe
│   └── <prefix>_CeccSegments.core.csv
└── <prefix>_Inferred_eccDNA/           # 推断的 eccDNA
    ├── <prefix>_UeccDNA_I.csv
    ├── <prefix>_UeccDNA_I.fasta
    ├── <prefix>_chimeric.csv
    └── <prefix>_CeccDNA_I.fasta
```

> 说明：Confirmed_* 目录为 **CtcReads-Caller** 的输出；Inferred_eccDNA 目录为 **SplitReads-Caller** 的输出。  
> 其中 **CtcReads** 指含 **Ctc**（**C**oncatemeric **t**andem **c**opies）信号的 reads（在 `tandem_to_ring.csv` 中以 CtcR-* 分类体现）。

---

## 1. 合并输出文件（merged_output.csv）

主输出文件，包含所有确认和推断的 eccDNA 汇总信息。

### 字段定义

| 字段 | 类型 | 说明 |
|------|------|------|
| `eccDNA_id` | string | 唯一标识符（如 UeccDNA1、MeccDNA1、CeccDNA1） |
| `original_id` | string | 重新编号前的原始 ID |
| `Regions` | string | 基因组坐标（格式：`chr:start-end`；CeccDNA 可有多个，用 `;` 分隔） |
| `Strand` | string | DNA 链（`+` 或 `-`；CeccDNA 可有多个） |
| `Length` | int | eccDNA 长度（bp） |
| `eccDNA_type` | string | 分类类型：`UeccDNA`、`MeccDNA`、`CeccDNA` |
| `State` | string | 检测状态：`Confirmed`（CtcReads-Caller）或 `Inferred`（SplitReads-Caller） |
| `Seg_total` | int | 片段数（CeccDNA 有多个，U/M 为 1） |
| `Hit_count` | int | 基因组命中数（MeccDNA 特有） |
| `confidence_score` | float | 置信度评分 [0,1]（越高越可靠） |
| `query_cov_best` | float | 最优位点/链的查询覆盖率 [0,1] |
| `query_cov_2nd` | float | 次优位点/链的查询覆盖率 [0,1] |
| `mapq_best` | int | 最优比对的 MAPQ 值 [0-60] |
| `identity_best` | float | 最优比对的同一性百分比 [0-100] |
| `low_mapq` | bool | 低 MAPQ 标记（mapq_best < 20） |
| `low_identity` | bool | 低同一性标记（identity_best < 95） |

> **注意**：置信度和证据字段主要针对 **Confirmed** 条目填充；**Inferred** 条目这些字段可能为空。

### 示例

```csv
eccDNA_id,original_id,Regions,Strand,Length,eccDNA_type,State,Seg_total,Hit_count,confidence_score,query_cov_best,query_cov_2nd,mapq_best,identity_best,low_mapq,low_identity
UeccDNA1,U001,chr1:10000-10500,+,500,UeccDNA,Confirmed,1,1,0.98,0.99,0.02,60,99.5,False,False
MeccDNA1,M001,chr2:20000-20800,+,800,MeccDNA,Confirmed,1,3,0.85,0.97,0.95,55,98.2,False,False
CeccDNA1,C001,chr1:5000-5200;chr3:8000-8300,+;-,500,CeccDNA,Confirmed,2,2,0.92,0.96,0.03,58,99.0,False,False
```

---

## 2. 确认的 eccDNA 输出

### 2.1 UeccDNA（简单 eccDNA）

#### CSV 文件格式（*_UeccDNA_C.csv）

| 字段 | 类型 | 说明 |
|------|------|------|
| `eccDNA_id` | string | 唯一标识符 |
| `chr` | string | 染色体名称 |
| `start0` | int | 0-based 起始坐标 |
| `end0` | int | 0-based 终止坐标（半开区间） |
| `strand` | string | DNA 链（+/-） |
| `length` | int | 序列长度（bp） |
| `reads` | string | 支持 reads 名称（分号分隔） |
| `copy_number` | float | 拷贝数估计 |
| `match_degree` | float | 匹配度 [0-100] |

#### FASTA 文件格式

```
>UeccDNA1 chr1:10000-10500(+) length=500
ATCGATCG...
```

### 2.2 MeccDNA（多拷贝 eccDNA）

#### CSV 文件格式（*_MeccDNA_C.csv）

| 字段 | 类型 | 说明 |
|------|------|------|
| `eccDNA_id` | string | 唯一标识符 |
| `chr` | string | 染色体名称 |
| `start0` | int | 0-based 起始坐标 |
| `end0` | int | 0-based 终止坐标 |
| `strand` | string | DNA 链 |
| `length` | int | 序列长度 |
| `reads` | string | 支持 reads 名称 |
| `copy_number` | float | 拷贝数估计 |
| `hit_index` | int | 当前命中索引（1-based） |
| `hit_count` | int | 总命中数 |

### 2.3 CeccDNA（复杂 eccDNA）

#### 主 CSV 文件（*_CeccDNA_C.csv）

| 字段 | 类型 | 说明 |
|------|------|------|
| `eccDNA_id` | string | 唯一标识符 |
| `chr` | string | 片段所在染色体 |
| `start0` | int | 片段 0-based 起始 |
| `end0` | int | 片段 0-based 终止 |
| `strand` | string | 片段链方向 |
| `length` | int | 片段长度 |
| `seg_index` | int | 片段索引（1-based） |
| `seg_total` | int | 总片段数 |
| `junction_role` | string | 连接角色：`head`、`body`、`tail` |
| `reads` | string | 支持 reads 名称 |
| `copy_number` | float | 拷贝数估计 |

#### 片段详情文件（*_CeccSegments.core.csv）

包含 CeccDNA 每个片段的详细信息，字段与主 CSV 相同。

#### BEDPE 连接文件（*_CeccJunctions.bedpe）

标准 BEDPE 格式，记录 CeccDNA 片段间的连接信息：

| 列 | 说明 |
|----|------|
| chrom1 | 片段1 染色体 |
| start1 | 片段1 起始（0-based） |
| end1 | 片段1 终止 |
| chrom2 | 片段2 染色体 |
| start2 | 片段2 起始 |
| end2 | 片段2 终止 |
| name | 连接名称（eccDNA_id） |
| score | 评分 |
| strand1 | 片段1 链方向 |
| strand2 | 片段2 链方向 |

---

## 3. 推断的 eccDNA 输出

### 3.1 推断的 UeccDNA（*_UeccDNA_I.csv）

字段与确认的 UeccDNA 基本相同，但置信度字段可能为空。

### 3.2 推断的嵌合体（*_chimeric.csv）

包含推断引擎检测到的潜在嵌合体 eccDNA。

### 3.3 推断的 CeccDNA（*_CeccDNA_I.fasta）

FASTA 格式，包含推断的复杂 eccDNA 序列。

---

## 4. 统计报告

### 摘要文件（*_summary.txt）

纯文本格式，包含运行统计：

```
CircleSeeker Analysis Summary
=============================
Sample: sample_name
Date: 2024-01-15 10:30:00

Total eccDNA detected: 150
  - UeccDNA (Confirmed): 80
  - MeccDNA (Confirmed): 30
  - CeccDNA (Confirmed): 15
  - Inferred: 25

Length distribution:
  - <500 bp: 45
  - 500-1000 bp: 60
  - 1000-5000 bp: 35
  - >5000 bp: 10
```

### HTML 报告（*_report.html）

交互式报告，包含：
- eccDNA 类型分布饼图
- 长度分布直方图
- 染色体分布条形图
- 可交互的数据表格

---

## 5. 坐标系统说明

CircleSeeker 在内部和输出文件中采用 **0-based 半开区间** 坐标系统：

- `start0`：包含的起始位置（0-based）
- `end0`：不包含的终止位置（0-based）
- 区间表示：`[start0, end0)`

### 示例

```
序列：  A T C G A T C G
位置：  0 1 2 3 4 5 6 7

区间 [2, 5) 表示：C G A（位置 2、3、4）
长度 = end0 - start0 = 5 - 2 = 3
```

### 与 1-based 坐标的转换

- BED 格式：直接使用（BED 本身是 0-based）
- FASTA header：显示为 `start0+1` 到 `end0`（1-based 包含区间）
- 与 UCSC/Ensembl 浏览器对接时：`start1 = start0 + 1`

---

## 6. 字段命名标准

CircleSeeker 使用统一的列命名标准（`ColumnStandard`）：

| 标准名称 | 旧版别名 | 说明 |
|----------|----------|------|
| `eccDNA_id` | - | 主标识符 |
| `chr` | `eChr` | 染色体 |
| `start0` | `eStart0`, `eStart` | 0-based 起始 |
| `end0` | `eEnd0`, `eEnd` | 0-based 终止 |
| `strand` | `eStrand` | 链方向 |
| `length` | `eLength`, `Length` | 长度 |
| `reads` | `eReads`, `readName` | reads 名称 |
| `copy_number` | `copyNum` | 拷贝数 |
| `match_degree` | `MatDegree` | 匹配度 |
| `eccdna_type` | `eClass` | eccDNA 类型 |
| `state` | `State` | 检测状态 |
| `hit_count` | `Hit_count` | 命中数 |
| `seg_total` | `Seg_total` | 片段总数 |

---

## 7. 相关文档

- [CLI 使用手册](CLI_Reference.md) - 命令行参数详解
- [配置参考手册](Configuration_Reference.md) - 配置选项详解
- [管道模块说明](Pipeline_Modules.md) - 16 步流程详解
