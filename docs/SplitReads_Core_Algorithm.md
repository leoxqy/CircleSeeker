# SplitReads-Core 算法说明（v1.1.0）

本文档详细描述 CircleSeeker 内置的 SplitReads-Core 模块，用于基于 split-read 证据推断 eccDNA。

> 说明：SplitReads-Core 的算法灵感来源于 [CReSIL](https://github.com/Peppermint-Lab/CReSIL)，并针对 PacBio HiFi 数据进行了优化。

---

## 1. 概述

SplitReads-Core 是一个基于 split-read 比对模式和断点证据的 eccDNA 检测模块。它通过分析 reads 在参考基因组上的比对模式，利用图算法检测环状 DNA 结构。

**主要特点**：
- 针对 HiFi 数据优化（使用 `map-hifi` preset）
- 无需外部依赖工具（使用 Python 的 mappy 进行比对）
- 基于图算法检测成环模式
- 支持简单型（单片段）和嵌合型（多片段）eccDNA 检测

---

## 2. 两阶段流程

### 2.1 Trim 阶段

**目的**：将原始 reads 比对到参考基因组，提取高质量的比对片段。

**流程**：

```
原始 reads (FASTA)
    ↓
minimap2 比对 (map-hifi preset)
    ↓
过滤低质量比对 (mapq >= 20)
    ↓
合并/修剪重叠/间隙区域
    ↓
trim_df (比对片段表)
```

**输出字段**：

| 字段 | 说明 |
|------|------|
| readid | 原始 read ID |
| q_start / q_end | query 上的起止位置 |
| r_start / r_end | 参考基因组上的起止位置 |
| ref | 染色体名 |
| strand | 比对方向 (+/-) |
| mapq | 比对质量 |
| order | 同一 read 内的片段顺序 |

**关键参数**：

| 参数 | 说明 | 默认值 |
|------|------|--------|
| preset | minimap2 preset | map-hifi |
| mapq | 最小比对质量 | 20 |
| allow_gap | 允许的 gap 大小 | 10 bp |
| allow_overlap | 允许的重叠大小 | 10 bp |

### 2.2 Identify 阶段

**目的**：基于 trim 阶段的比对结果，通过图算法识别环状结构。

**流程**：

```
trim_df (比对片段表)
    ↓
计算基因组覆盖度
    ↓
合并邻近区域 → merge regions
    ↓
过滤低覆盖度区域 (depth >= 5)
    ↓
分析断点方向
    ↓
构建断点图
    ↓
检测环状子图
    ↓
eccDNA 候选列表
```

---

## 3. 图算法详解

### 3.1 Merge Region 生成

将基因组上覆盖度足够的连续区域合并为 "merge region"：

```
条件：
- 平均覆盖度 >= average_depth (默认 5)
- 区域长度 >= min_region_size (默认 200 bp)

输出格式：
mergeid = "{chrom}_{start}_{end}"
```

### 3.2 断点方向分析

对于同一个 read 的多个比对片段，分析相邻片段之间的链方向转换：

```python
def check_breakpoint_direction(df_sorted_by_order):
    """
    检查相邻比对片段的链方向。

    返回: (region1, region2, strand_transition, is_valid)

    strand_transition 格式: "{strand1}_{strand2}"
    例如: "+_+", "+_-", "-_+", "-_-"
    """
```

### 3.3 构建断点图

```python
G = nx.MultiDiGraph()

# 节点: merge regions
# 边: 相邻比对片段之间的连接

for each read:
    for adjacent_pair in read.sorted_alignments:
        if is_valid_breakpoint(pair):
            G.add_edge(region1, region2, weight=1)

# 过滤低支持度的边
graph_filtered = {k: v for k, v in graph.items()
                  if v >= breakpoint_depth}  # 默认 5
```

### 3.4 环检测

```python
# 获取连通分量
subgraphs = nx.connected_components(G.to_undirected())

for subgraph in subgraphs:
    nodes = list(subgraph.nodes())

    if len(nodes) == 1:
        # 单节点: 简单型 eccDNA (自环)
        is_cyclic = has_self_loop(subgraph)

    elif len(nodes) == 2:
        # 双节点: 检查双向边
        is_cyclic = has_bidirectional_edge(subgraph)

    else:
        # 多节点: 使用 cycle_basis 检测环
        cycles = nx.cycle_basis(subgraph)
        is_cyclic = len(cycles) > 0
```

### 3.5 链闭合验证

对于检测到的环，验证链方向是否能形成有效闭合：

```python
def check_strand_closure(strand_transitions):
    """
    验证链转换序列是否能形成闭合环。

    规则：
    - +_+ 和 -_- 保持方向
    - +_- 和 -_+ 翻转方向
    - 有效闭合需要偶数次翻转
    """
    flips = count(trans for trans in transitions
                  if trans in ["+_-", "-_+"])
    return flips % 2 == 0
```

---

## 4. 输出格式

### 4.1 eccDNA_final.txt

| 字段 | 说明 |
|------|------|
| id | eccDNA ID (如 ec1, ec2) |
| merge_region | 区域列表 (如 "chr1_100_200_+,chr2_300_400_-") |
| merge_len | 总长度 |
| num_region | 片段数 |
| ctc | 是否为串联重复型 |
| numreads | 支持的 reads 数 |
| totalbase | 总覆盖碱基数 |
| coverage | 估计覆盖度 |

### 4.2 分类逻辑

```
num_region == 1  →  简单型 (IUeccDNA)
num_region >= 2  →  嵌合型 (ICeccDNA)
```

---

## 5. 参数清单

| 参数 | 说明 | 默认值 |
|------|------|--------|
| **比对参数** | | |
| preset | minimap2 preset | map-hifi |
| mapq | 最小比对质量 | 20 |
| exclude_chrs | 排除的染色体 | "" |
| **Trim 参数** | | |
| allow_gap | 允许的 gap 大小 | 10 bp |
| allow_overlap | 允许的重叠大小 | 10 bp |
| **Identify 参数** | | |
| min_region_size | 最小区域大小 | 200 bp |
| overlap_check_size | 边界检查大小 | 50 bp |
| breakpoint_depth | 最小断点支持度 | 5 |
| average_depth | 最小平均覆盖度 | 5.0 |

---

## 6. 与 CtcReads-Caller 的区别

| 特性 | CtcReads-Caller | SplitReads-Core |
|------|-----------------|-----------------|
| 证据类型 | 串联重复 (TideHunter) | Split-read 比对模式 |
| 输出状态 | Confirmed | Inferred |
| 适用场景 | 高覆盖度、明显滚环扩增 | 低覆盖度、无明显重复 |
| 检测方法 | 覆盖度模型 + LAST | 断点图 + NetworkX |

---

## 7. 依赖

- **mappy**: Python minimap2 绑定
- **networkx**: 图算法库
- **pybedtools**: BED 文件操作（需要系统 bedtools）
- **intervaltree**: 区间树操作

---

## 8. 参考

- [CReSIL](https://github.com/Peppermint-Lab/CReSIL) - 原始算法实现
- [minimap2](https://github.com/lh3/minimap2) - 序列比对
- [NetworkX](https://networkx.org/) - 图算法
