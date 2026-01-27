# U/M/C eccDNA 覆盖度模型与 CeccBuild 说明

本文件描述 CircleSeeker 当前实现中的 U/M 分类覆盖度模型，以及 CeccBuild 的检测策略。U/M 判别基于覆盖度与位点聚类；Cecc 默认由 LAST 驱动的重复模式检测完成，图结构链式模型仅作为回退方案。

> 说明：文档聚焦 um_classify/cecc_build 的实际逻辑；如需完整代码细节，请参阅 `src/circleseeker/modules/um_classify.py` 与 `src/circleseeker/modules/cecc_build.py`。

---

## 1. 输入与符号

- `S`：TideHunter 的共识序列（环的一圈），长度 `L`（通常取 `consLen`）。
- `Q = S + S`：双倍化 query，长度 `2L`，用于处理环起点不确定性。
- 比对输出集合 `A = {a_i}`：每条命中 `a_i` 至少包含：
  - query 区间 `q_i = [qs_i, qe_i)`
  - 参考位点 `g_i`（chr/start/end/strand）
  - 可选质量信息（identity/mapq/aln_len）

---

## 2. Locus 聚类（同位点合并）

为避免边界抖动导致的错误拆分，需要将相近的基因组命中合并为同一 locus。

推荐规则（参数化）：

- `chr` 相同且 `strand` 相同
- 互惠重叠率 `overlap_len / min(len_i, len_j) ≥ θ_locus`
- 可选边界容差：`|start_i - start_j| ≤ pos_tol_bp` 且 `|end_i - end_j| ≤ pos_tol_bp`

默认建议：`θ_locus=0.95`，`pos_tol_bp=50`。

---

## 3. 覆盖度与 U/M 判别

### 3.1 覆盖度

对任意命中集合 `X`，覆盖度定义为：

`Cov(X) = | ⋃ proj(q_i) (i∈X) | / L`

其中 `proj(q_i)` 为将 query 区间投影到环坐标 `[0, L)` 后的区间集合。

### 3.2 Uecc 判别

- `U_cov = max_ℓ Cov(ℓ)`
- `U_cov_2nd`：第二高 locus 覆盖度

判定条件：

- `U_cov ≥ θ_U`
- `U_cov_2nd ≤ θ_U2_max`

附加约束（可选）：

- `mapq_u_min`：要求最佳 locus MAPQ ≥ 阈值
- `span_ratio_min`：要求基因组跨度/共识长度比例 ≥ 阈值（默认 0.95）
  - 用于减少 CeccDNA 被误分类为 UeccDNA 的情况
  - 当 eccDNA 实际由多个片段组成时，单个 locus 的跨度通常小于共识长度
  - 计算公式：`span_ratio = best_genomic_span / cons_len`
- 次级比对 veto（防止多染色体/非连续映射）：
  - `u_secondary_min_frac` / `u_secondary_min_bp`
  - `u_contig_gap_bp`
  - `u_secondary_max_ratio`
  - `u_high_coverage_threshold` / `u_high_mapq_threshold`

### 3.3 Mecc 判别

- `M_count = #{ℓ | Cov(ℓ) ≥ θ_M}`
- 判定：`M_count ≥ 2`

可选过滤（默认关闭）：

- `mapq_m_ambiguous_threshold`：若所有 locus MAPQ 低于阈值，则拒绝 Mecc
- `mecc_identity_gap_threshold`：若最高与次高 identity 差超过阈值，则拒绝 Mecc

---

## 4. Cecc 检测

### 4.1 CeccBuild 概述

CeccBuild 使用 LAST 对双倍化序列进行高精度比对，通过图结构分析识别嵌合型 eccDNA。

**核心思路**：
- 输入：双倍化序列 `Q = S + S`（长度 `2L`）
- 比对：使用 LAST（lastal | last-split）获得高精度、最优链比对
- 分析：基于图结构检测环状模式

### 4.2 图结构算法（Graph-based Detection）

#### 4.2.1 Locus 聚类

将邻近的比对合并为 locus（同 U/M 分类）：
- 同一染色体、同一链方向
- 基因组坐标重叠或相邻（容差 `position_tolerance`）

#### 4.2.2 构建有向图

```
G = (V, E)
V = {locus_1, locus_2, ..., locus_n}
E = {(locus_i, locus_j, strand_trans) | 在 query 序列上相邻}
```

- **节点**：每个 locus 作为一个节点，携带 `(chr, start, end, strand)` 信息
- **边**：按 query 坐标排序后，相邻的 locus 之间连边
- **链转换**：记录相邻比对的链方向变化（`+_+`, `+_-`, `-_+`, `-_-`）

#### 4.2.3 环检测（双倍序列模式）

对于双倍序列 `[A-B-C][A-B-C]`：
- 比对模式应为：`A' → B → C → A → B'`（`'` 表示部分比对）
- 检测方法：寻找 locus 序列的重复模式
- 验证条件：
  - 第一半和第二半的 locus 序列应匹配（旋转等价）
  - 链转换需满足闭合条件（偶数次翻转）

#### 4.2.4 链闭合验证

```python
# 链转换闭合检查
flips = count(trans for trans in transitions if trans in ["+_-", "-_+"])
valid = (flips % 2 == 0)  # 偶数次翻转才能形成有效闭合
```

### 4.3 LAST 比对流程

```bash
# 1. 建立索引
lastdb -P threads db_prefix reference.fa

# 2. 比对 + 最优链选择
lastal -P threads db_prefix query.fa | last-split > output.maf
```

**last-split 的作用**：为每个 query 区域选择最优比对链，过滤多重映射噪声。

### 4.4 关键参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `min_match_degree` | 最小匹配度（覆盖度 × 100） | 90.0 |
| `overlap_threshold` | 互惠重叠阈值 | 0.95 |
| `min_segments` | 最少片段数 | 2 |
| `position_tolerance` | locus 合并位置容差 | 100 bp |
| `edge_tolerance` | query 边界容差 | 10 bp |
| `half_query_buffer` | 双倍序列中点缓冲 | 50 bp |

### 4.5 输出

每个检测到的 CeccDNA 包含：
- `eccDNA_id`：唯一标识符
- `segments`：各片段的基因组坐标列表
- `junction_roles`：片段角色（head/middle/tail）
- `match_degree`：匹配度分数
- `CeccClass`：分类（`Cecc-InterChr` / `Cecc-IntraChr`）

---

## 5. 参数清单（默认值）

| 参数 | 含义 | 默认 |
|------|------|------|
| `theta_u` / `theta_m` / `theta_full` | U/M 覆盖度阈值 | 0.95 |
| `theta_u2_max` | U 第二位点覆盖上限 | 0.05 |
| `theta_locus` | locus 互惠重叠阈值 | 0.95 |
| `pos_tol_bp` | locus 边界容差 | 50 |
| `mapq_u_min` | U 的 MAPQ 下限 | 0 |
| `span_ratio_min` | U 的基因组跨度/共识长度比例下限 | 0.95 |
| `u_secondary_min_frac` / `u_secondary_min_bp` | 次级比对阈值 | 0.01 / 50 |
| `u_contig_gap_bp` | U 连续性阈值 | 1000 |
| `u_secondary_max_ratio` | 次级/主位点覆盖比 | 0.05 |
| `u_high_coverage_threshold` | 高覆盖放宽阈值 | 0.98 |
| `u_high_mapq_threshold` | 高 MAPQ 放宽阈值 | 50 |
| `mapq_m_ambiguous_threshold` | Mecc MAPQ 过滤 | 0 (禁用) |
| `mecc_identity_gap_threshold` | Mecc identity 差过滤 | 0 (禁用) |
| `overlap_threshold` | Cecc 互惠重叠阈值 | 0.95 |
| `min_segments` | Cecc 最少片段数 | 2 |
| `tau_gap` / `edge_tolerance` | 链 gap 容差 | 20 |
| `min_match_degree` | Cecc 匹配度阈值 | 95.0 |
| `half_query_buffer` | 双倍序列中点缓冲 | 50 |

---

## 6. 备注

- `delta_uc` / `epsilon_mc` 等歧义参数在当前实现中保留但未单独输出文件；可用于后续扩展或研究用途。
