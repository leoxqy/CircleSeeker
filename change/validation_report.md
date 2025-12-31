# CircleSeeker 模拟数据验证报告

**生成日期**: 2025-12-31
**数据来源**: change/ 文件夹模拟数据
**CircleSeeker版本**: 0.9.8

---

## 1. 数据概览

### 1.1 模拟数据统计

| 项目 | 数值 |
|------|------|
| 总reads数 | 50,000 |
| 独立eccDNA数 | 587 |
| UeccDNA | 492 个 (41,563 reads) |
| MeccDNA | 86 个 (6,075 reads) |
| CeccDNA | 8 个 (2,357 reads) |

### 1.2 检测结果统计

| 类型 | 检测数 | 真实数 |
|------|--------|--------|
| UeccDNA | 489 | 492 |
| MeccDNA | 74 | 86 |
| CeccDNA | 15 | 8 |
| 总计 | 578 | 587 |

---

## 2. 检出率分析

### 2.1 按eccDNA个体计算

| 类型 | 真实数 | 检出数 | 检出率 | 状态 |
|------|--------|--------|--------|------|
| UeccDNA | 492 | 456 | 92.7% | ✓ 良好 |
| MeccDNA | 86 | 34 | 39.5% | ⚠ 需关注 |
| CeccDNA | 8 | 7 | 87.5% | ✓ 良好 |

### 2.2 Reads级别混淆矩阵

```
真实类型    总reads    正确分类    错误分类    未检出      正确率
─────────────────────────────────────────────────────────────
UeccDNA     41,563    35,671      192       5,700      86.0%
MeccDNA      6,075        68      216       5,791       1.1%  ❌
CeccDNA      2,357        30      616       1,711       1.3%  ❌
```

### 2.3 详细分类流向

```
真实类型    ->UeccDNA  ->MeccDNA  ->CeccDNA  ->XeccDNA  未检出
───────────────────────────────────────────────────────────────
UeccDNA      35,671       192          0          0     5,700
MeccDNA         216        68          0          0     5,791
CeccDNA         616         0         30          0     1,711
```

---

## 3. 关键问题诊断

### 3.1 问题1: MeccDNA大量漏检 (严重)

**现象**:
- 6,075条MeccDNA reads中，仅68条被正确分类 (1.1%)
- 5,791条完全未被检测到 (95.3%)
- 216条被误分类为UeccDNA

**漏检MeccDNA列表 (reads数最多的)**:

| eccDNA ID | 总reads | 未检出 | 平均repeat |
|-----------|---------|--------|------------|
| MeccDNA_000575 | 1,109 | 1,108 | 5.2 |
| MeccDNA_000239 | 805 | 805 | 7.0 |
| MeccDNA_000772 | 663 | 621 | 6.6 |
| MeccDNA_000601 | 416 | 416 | 6.5 |
| MeccDNA_000769 | 419 | 399 | 9.2 |

**可能原因**:
- MeccDNA的多拷贝特征未被正确识别
- minimap2比对后的多位点分析逻辑可能存在问题
- 这些reads的repeat_count都>=2，不是TideHunter检测问题

### 3.2 问题2: CeccDNA被错误分类 (严重)

**现象**:
- CeccDNA_000059 的616条reads被100%误分类为UeccDNA
- 仅30条CeccDNA reads被正确识别
- 1,711条CeccDNA reads未被检出

**各CeccDNA检出详情**:

| eccDNA ID | 总reads | 检出数 | 检出率 | 单体长度(bp) |
|-----------|---------|--------|--------|--------------|
| CeccDNA_000059 | 616 | 616 | 100%* | 6,277 |
| CeccDNA_000019 | 374 | 5 | 1.3% | 5,046 |
| CeccDNA_000032 | 378 | 9 | 2.4% | 10,457 |
| CeccDNA_000029 | 337 | 8 | 2.4% | 9,994 |
| CeccDNA_000098 | 315 | 0 | 0% | 5,793 |
| CeccDNA_000043 | 231 | 4 | 1.7% | 6,681 |
| CeccDNA_000030 | 105 | 3 | 2.9% | 6,231 |
| CeccDNA_000064 | 1 | 1 | 100% | 2,245 |

*注: CeccDNA_000059被错误分类为UeccDNA

**可能原因**:
- 嵌合结构特征未被正确识别
- CeccDNA的跨染色体segment检测逻辑可能存在问题

### 3.3 问题3: repeat_count<2 的reads全部漏检 (预期行为)

**现象**:
- 1,802条repeat_count<2的reads全部未检出
- 检出率: 0%

**说明**:
- 这是预期行为，TideHunter需要至少2次重复才能检测串联重复结构
- 不需要修复

### 3.4 问题4: 长单体eccDNA检出率较低

**按单体长度统计检出率**:

| 单体长度范围 | 检出reads | 未检出reads | 检出率 |
|--------------|-----------|-------------|--------|
| 0-500 bp | 480 | 196 | 71.0% |
| 500-1000 bp | 3,449 | 1,066 | 76.4% |
| 1000-2000 bp | 10,973 | 2,074 | 84.1% |
| 2000-3000 bp | 11,505 | 3,579 | 76.3% |
| 3000-5000 bp | 8,417 | 2,891 | 74.4% |
| 5000-10000 bp | 1,839 | 1,221 | 60.1% |
| 10000-50000 bp | 133 | 375 | **26.2%** |

**可能原因**:
- TideHunter对长单体串联重复检测能力有限
- 可能需要调整TideHunter参数

---

## 4. 漏检eccDNA详细列表

### 4.1 漏检UeccDNA (36个)

| eccDNA ID | reads数 | 平均repeat |
|-----------|---------|------------|
| UeccDNA_000152 | 633 | 1.0 |
| UeccDNA_001670 | 518 | 6.5 |
| UeccDNA_003855 | 272 | 7.8 |
| UeccDNA_002418 | 271 | 9.9 |
| UeccDNA_000773 | 255 | 18.3 |
| UeccDNA_004199 | 236 | - |
| UeccDNA_004368 | 154 | - |
| UeccDNA_000145 | 116 | 1.0 |
| ... | ... | ... |

### 4.2 漏检MeccDNA (52个)

| eccDNA ID | reads数 | 平均repeat |
|-----------|---------|------------|
| MeccDNA_000239 | 805 | 7.0 |
| MeccDNA_000601 | 416 | 6.5 |
| MeccDNA_000871 | 212 | 11.3 |
| MeccDNA_000793 | 169 | 13.0 |
| MeccDNA_000786 | 114 | 10.8 |
| MeccDNA_000581 | 111 | - |
| MeccDNA_000939 | 103 | - |
| MeccDNA_000429 | 87 | - |
| ... | ... | ... |

### 4.3 漏检CeccDNA (1个)

| eccDNA ID | reads数 | 平均repeat |
|-----------|---------|------------|
| CeccDNA_000098 | 315 | 3.5 |

---

## 5. XeccDNA分析

- 检测到 **4,385** 条XeccDNA reads
- 这些reads在truth.tsv中**未找到对应记录**
- 可能是背景reads或假阳性
- 需要进一步分析这些reads的来源

---

## 6. 结论与建议

### 6.1 总体评估

| 检测类型 | 评估 | 说明 |
|----------|------|------|
| UeccDNA | ✓ 良好 | 92.7%个体检出，86%reads正确分类 |
| MeccDNA | ✗ 严重问题 | 仅39.5%个体检出，1.1%reads正确分类 |
| CeccDNA | ⚠ 分类问题 | 87.5%个体检出，但分类准确率仅1.3% |

### 6.2 建议排查方向

1. **MeccDNA检测逻辑**
   - 检查minimap2比对后多位点识别逻辑
   - 检查MeccDNA分类条件是否过于严格
   - 分析为什么大量MeccDNA reads完全未进入分类流程

2. **CeccDNA分类逻辑**
   - 检查嵌合结构检测条件
   - 分析为什么CeccDNA_000059被全部误分为UeccDNA
   - 检查跨染色体segment的识别逻辑

3. **TideHunter参数优化**
   - 考虑调整参数以提高长单体eccDNA的检出率
   - 评估是否可以降低最小repeat次数要求

### 6.3 后续分析计划

- [ ] 分析MeccDNA reads在pipeline各步骤的流向
- [ ] 检查CeccDNA_000059与其他CeccDNA的结构差异
- [ ] 分析XeccDNA reads的真实来源
- [ ] 代码审查MeccDNA和CeccDNA的分类逻辑

---

## 附录: 分析脚本

分析使用Python脚本对truth.tsv和检测结果进行比对，主要步骤：

1. 从truth.tsv提取真实eccDNA ID和reads映射
2. 从各检测结果CSV提取检出reads
3. 通过read_id进行交叉比对
4. 计算混淆矩阵和检出率

---

## 7. 根因分析

### 7.1 MeccDNA检测问题 - 代码定位

**问题代码**: `src/circleseeker/modules/um_classify.py` 第436-456行

```python
def _is_full_length_repeat(self, group: pd.DataFrame) -> bool:
    # 要求覆盖度 >= 95% 且至少2个full-length copies
    coverages = (valid_group["Rlength"] / valid_group[ColumnStandard.LENGTH]) * 100
    full_length_count = (coverages >= self.min_full_length_coverage).sum()
    return full_length_count >= 2  # 默认 min_full_length_coverage = 95%
```

**问题分析**:
1. MeccDNA的consensus序列比对到基因组重复区域时
2. 由于重复序列之间存在变异，单个比对覆盖度很难全部达到95%
3. 结果：大量MeccDNA被分类到"unclassified"而非"Mecc"
4. unclassified中的数据没有后续的MeccDNA识别逻辑

**建议修复**:
- 方案A: 降低 `min_full_length_coverage` 阈值到85-90%
- 方案B: 改变判断逻辑为"有多个比对位置"即可，不要求全部95%覆盖

### 7.2 CeccDNA分类问题 - 代码定位

**分类流程**:
```
um_classify.py:
  单比对位置 -> UeccDNA (直接)
  多比对位置且覆盖度>=95% -> MeccDNA
  其他 -> unclassified -> cecc_build.py -> CeccDNA
```

**CeccDNA_000059被误分为UeccDNA的可能原因**:

1. **比对结果只有单位置**: CeccDNA的嵌合结构如果在比对时只输出一个最佳比对（minimap2默认行为），就会在um_classify阶段直接被判定为UeccDNA

2. **Overlap处理问题**: 如果多个segment比对到相近位置，经过overlap处理后可能只剩1个，导致被误判为UeccDNA

**建议修复**:
- 确保minimap2输出所有secondary alignments（-N参数）
- 在um_classify之前保留完整的多位点比对信息
- CeccDNA检测应该在UeccDNA分类之前进行

### 7.3 数据流问题总结

```
当前流程问题:
┌──────────────────────────────────────────────────────────────────┐
│ TideHunter -> minimap2 -> um_classify                            │
│                               ↓                                   │
│              ┌─────────────────┴─────────────────┐               │
│              ↓                 ↓                 ↓               │
│           UeccDNA           MeccDNA         unclassified         │
│           (错误地包含了    (条件太严格      (包含大量            │
│            CeccDNA)        漏掉95%的       MeccDNA和CeccDNA)     │
│                            MeccDNA)                               │
│                                                 ↓                 │
│                                            cecc_build             │
│                                         (只处理嵌合结构,          │
│                                          不识别MeccDNA)           │
└──────────────────────────────────────────────────────────────────┘
```

### 7.4 建议的修复优先级

| 优先级 | 问题 | 修复建议 | 预期效果 |
|--------|------|----------|----------|
| P0 | MeccDNA 95%漏检 | 降低覆盖度阈值或改变判断逻辑 | MeccDNA检出率从5%提升到80%+ |
| P1 | CeccDNA误分类 | 调整minimap2参数和分类顺序 | CeccDNA准确率从1%提升到80%+ |
| P2 | unclassified处理 | 添加MeccDNA fallback逻辑 | 减少漏检 |

---

*报告生成: Claude Code*
