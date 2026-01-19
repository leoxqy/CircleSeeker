# CircleSeeker Pipeline 模块说明（v1.0.0）

本文件概述 CircleSeeker 1.0.0 的 16 个管线步骤、输入输出关系以及关键实现要点，帮助使用者理解整体流程与模块职责。

---

## 1. 全局视图

CircleSeeker 以 4 个阶段串联 16 个步骤，既包含外部工具调用，也包含内建的 Python 模块：

| 阶段 | 步骤范围 | 目标 |
|------|----------|------|
| **检测阶段** | 1-6 | 基于 HiFi reads 和参考基因组，识别潜在的环状 DNA 候选 |
| **处理阶段** | 7-10 | 对 U/M/C 三类候选进行聚合、去重与过滤 |
| **推断阶段** | 11-13 | 使用 Cresil（或 Cyrcular）补充推断结果并进行整理 |
| **整合阶段** | 14-16 | 合并所有信息、生成报告并打包产物 |

从证据来源角度，16 个步骤也可以更直观地理解为“两条 Caller + 一个整合阶段”：

- **CtcReads**：含 **Ctc**（**C**oncatemeric **t**andem **c**opies）信号的 reads（对应 `tandem_to_ring.csv` 中的 CtcR-* 分类）。
- **CtcReads-Caller**（步骤 1-10）：基于 CtcReads 证据产出 **Confirmed** U/M/C eccDNA。
- **SplitReads-Caller**（步骤 11-13）：基于 split-reads/junction 证据进行推断（Cresil 优先，Cyrcular 备用），产出 **Inferred** eccDNA。
- **Integration**（步骤 14-16）：对 Confirmed/Inferred 做去冗余合并、统计与打包交付。

运行过程中，所有中间文件写入 `<output>/.tmp_work/`，最终由 `ecc_packager` 复制到目标目录结构。

---

## 2. 步骤总览

| 序号 | 名称 | 类型 | 主要输入 | 主要输出 |
|------|------|------|----------|----------|
| 1 | check_dependencies | 内部 | 配置/环境 | 依赖检查报告（失败则中止） |
| 2 | tidehunter | 外部 | HiFi reads FASTA | 重复片段列表、共识序列 |
| 3 | tandem_to_ring | 内部 | TideHunter 输出 | 环状候选 FASTA/CSV |
| 4 | run_alignment | 外部（minimap2/LAST） | 候选 FASTA、参考基因组 | alignment TSV |
| 5 | um_classify | 内部 | alignment TSV | `um_classify.uecc.csv` / `mecc.csv` |
| 6 | cecc_build | 内部 | 未分类对齐记录 | `cecc_build.csv` |
| 7 | umc_process | 内部 | U/M/C CSV | 标准化表格、FASTA |
| 8 | cd_hit | 外部（CD-HIT） | FASTA | 去冗余 FASTA |
| 9 | ecc_dedup | 内部 | 各类 CSV/FASTA | 去重后的坐标与序列 |
|10 | read_filter | 内部 | 原始 reads、tandem_to_ring.csv | 推断输入 FASTA |
|11 | minimap2 | 外部（minimap2） | 参考基因组 | `.mmi` 索引（必要时自动生成） |
|12 | ecc_inference | 外部 + 内部 | Cresil/Cyrcular 所需文件 | 推断结果 TSV/FASTA |
|13 | curate_inferred_ecc | 内部 | 推断结果 | 精炼后的表格与 FASTA（内部调用 `iecc_curator`） |
|14 | ecc_unify | 内部 | 确认/推断 CSV | 合并表格、重叠统计 |
|15 | ecc_summary | 内部 | 合并表、原始统计 | HTML 报告、TXT 摘要、all.fasta |
|16 | ecc_packager | 内部 | 各类目录与单文件 | 最终输出目录树 |

---

## 3. 阶段详解

### 3.1 检测阶段（CtcReads-Caller，步骤 1-6）

1. **check_dependencies**
   启动时预检外部工具与推断引擎（至少需要 Cresil 或 Cyrcular 之一）。缺失依赖会中止运行并给出安装提示。

2. **tidehunter**
   调用 TideHunter 检测 reads 中的串联重复，识别潜在滚环扩增事件，生成共识序列与 repeat 信息。

3. **tandem_to_ring**
   将重复序列转换为候选环状 DNA。通过图结构分析 overlaps，识别简单读段、复杂读段等类型，输出 FASTA 与候选元数据。

4. **run_alignment**
   使用 minimap2（默认）或 LAST（`tools.alignment.aligner`）将候选序列比对回参考基因组，输出 BLAST outfmt 6-like 的 TSV，末尾追加 `mapq` 字段供下游解析。minimap2 支持长度补偿的 identity 门控（`min_identity` / `identity_decay_per_10kb` / `min_identity_floor`）。

5. **um_classify**
   对对齐结果进行解析，基于“环坐标覆盖度 + locus 聚类”的模型划分 UeccDNA / MeccDNA；可选通过 `mapq_u_min` 与次级比对 veto 提升 Uecc 可靠性。输出 `um_classify.uecc.csv` / `um_classify.mecc.csv` / `um_classify.unclassified.csv`。判别模型详见 `docs/UMC_Classification_Model.md`。

6. **cecc_build**
   基于 `um_classify.unclassified.csv` 进行 Cecc 检测。优先使用 LAST 对双倍序列进行重复模式检测；若 LAST 或关键输入缺失则回退图方法。输出 `cecc_build.csv`。

### 3.2 处理阶段（CtcReads-Caller，步骤 7-10）

7. **umc_process**
   统一整理 U/M/C 三类结果，生成标准化 CSV 与 FASTA；若 `enable_xecc=true` 会额外导出 XeccDNA FASTA。

8. **cd_hit**
   对 FASTA 进行 99% identity 去冗余，减少重复候选。

9. **ecc_dedup**
   整理并去重坐标信息，生成规范化 BED/CSV，确保重复记录合并到统一 ID。

10. **read_filter**
    使用 Sieve 根据 `tandem_to_ring.csv` 过滤 CtcR reads，生成推断输入 FASTA；若 CSV 缺失则保留全部 reads。

### 3.3 推断阶段（SplitReads-Caller，步骤 11-13）

11. **minimap2**
    检查并生成参考基因组 `.mmi` 索引，用于推断阶段。若选择 Cresil，则通常跳过 BAM 比对，仅保留索引。

12. **ecc_inference**
    - 优先调用 Cresil：对 reads 进行 trim/identify，若缺失 `.fai` 会尝试执行 `samtools faidx` 自动补建。
    - 若 Cresil 不可用或执行失败，则尝试 Cyrcular 作为备选。
    - 推断失败会写入 `<prefix>_inference_failed.txt` 供排查。

13. **curate_inferred_ecc**
    对推断结果进行清洗、格式化，生成简单/嵌合 CSV；若参考基因组可用且 pysam 可用，会额外生成 FASTA（内部调用 `iecc_curator` 模块）。

### 3.4 整合阶段（Integration，步骤 14-16）

14. **ecc_unify**
    汇总确认与推断结果，基于片段互惠重叠检测冗余嵌合体（默认阈值 0.99、坐标容差 10bp），并产出 `*_unified.csv` 与重叠统计文件。

15. **ecc_summary**
    - 合并所有 FASTA，生成 `<prefix>_all.fasta`；
    - 调用 `EccSummary` 计算样本统计并输出 HTML/TXT 报告；
    - 缺失的统计文件会自动填充默认值以保证后续步骤可运行。

16. **ecc_packager**
    按固定结构复制并重命名文件，形成最终产物目录（Confirmed_U/M/C、Inferred_eccDNA、报告等），同时清理临时文件（除非启用 `--keep-tmp`）。XeccDNA 不再拷贝到最终目录，源 reads 已并入推断输入。

---

## 4. 运行时注意事项

- **外部依赖** 请确认 `tidehunter`, `cd-hit-est`, `minimap2`, `samtools`, `cresil`（或 `cyrcular`）均在 `PATH` 内。若启用 LAST（`tools.alignment.aligner=last` 或 CeccBuild v4），需额外安装 LAST（`lastdb`/`lastal`）。
- **检查点机制** 每个步骤开始前都会写入 `<prefix>.checkpoint`，一旦失败可通过 `--resume` 从最近步骤继续；`--force` 可清除状态后重跑。
- **配置继承** 所有参数均由 `circleseeker.config.Config` 驱动，可通过 YAML 或 CLI 覆写。文档中的默认值等同于仓库内置配置。
- **调试工具** `circleseeker --debug --show-steps` 可查看当前步骤状态；`show-checkpoint` 子命令能列出历史执行记录。

---

## 5. 深入阅读

- `src/circleseeker/core/pipeline.py` 主管线实现与状态管理
- `src/circleseeker/config.py` 配置结构及默认值
- `src/circleseeker/modules/` 各内部模块的详细算法
- `tests/unit/` 单元测试示例，有助于理解输入输出契约

如文档与实际行为存在差异，欢迎提交 issue 或 pull request。
