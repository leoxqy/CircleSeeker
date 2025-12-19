# CircleSeeker Pipeline 模块说明（v0.9.3）

本文件概述 CircleSeeker 0.9.3 的 16 个管线步骤、输入输出关系以及关键实现要点，帮助使用者理解整体流程与模块职责。

---

## 1. 全局视图

CircleSeeker 以 4 个阶段串联 16 个步骤，既包含外部工具调用，也包含内建的 Python 模块：

| 阶段 | 步骤范围 | 目标 |
|------|----------|------|
| **检测阶段** | 1–6 | 基于 HiFi reads 和参考基因组，识别潜在的环状 DNA 候选 |
| **处理阶段** | 7–10 | 对 U/M/C 三类候选进行聚合、去重与过滤 |
| **推断阶段** | 11–13 | 使用 Cresil（或 Cyrcular）补充推断结果并进行整理 |
| **整合阶段** | 14–16 | 合并所有信息、生成报告并打包产物 |

运行过程中，所有中间文件写入 `<output>/.tmp_work/`，最终由 `ecc_packager` 复制到目标目录结构。

---

## 2. 步骤总览

| 序号 | 名称 | 类型 | 主要输入 | 主要输出 |
|------|------|------|----------|----------|
| 1 | make_blastdb | 外部（BLAST+） | 参考基因组 FASTA | `*.nsq` 等 BLAST 数据库 |
| 2 | tidehunter | 外部 | HiFi reads FASTA | 重复片段列表、共识序列 |
| 3 | tandem_to_ring | 内部 | TideHunter 输出 | 环状候选 FASTA/CSV |
| 4 | run_blast | 外部（BLAST+） | 候选 FASTA、数据库 | BLAST TSV |
| 5 | um_classify | 内部 | BLAST TSV | `um_classify.uecc.csv` / `mecc.csv` |
| 6 | cecc_build | 内部 | 未分类 BLAST 记录 | `cecc_build.csv` |
| 7 | umc_process | 内部 | U/M/C CSV | 标准化表格、FASTA |
| 8 | cd_hit | 外部（CD-HIT） | FASTA | 去冗余 FASTA |
| 9 | ecc_dedup | 内部 | 各类 CSV/FASTA | 去重后的坐标与序列 |
|10 | read_filter | 内部 | 去重结果、reads | 确认的 eccDNA reads FASTA |
|11 | minimap2 | 外部（minimap2） | 参考基因组 | `.mmi` 索引（必要时自动生成） |
|12 | ecc_inference | 外部 + 内部 | Cresil/Cyrcular 所需文件 | 推断结果 TSV/FASTA |
|13 | iecc_curator | 内部 | 推断结果 | 精炼后的表格与 FASTA |
|14 | ecc_unify | 内部 | 确认/推断 CSV | 合并表格、重叠统计 |
|15 | ecc_summary | 内部 | 合并表、原始统计 | HTML 报告、TXT 摘要、all.fasta |
|16 | ecc_packager | 内部 | 各类目录与单文件 | 最终输出目录树 |

---

## 3. 阶段详解

### 3.1 检测阶段（步骤 1–6）

1. **make_blastdb**  
   使用 `makeblastdb` 为参考基因组建立核酸数据库。若数据库已存在会跳过。

2. **tidehunter**  
   调用 TideHunter 检测 reads 中的串联重复，识别潜在滚环扩增事件，生成共识序列与 repeat 信息。

3. **tandem_to_ring**  
   将重复序列转换为候选环状 DNA。通过图结构分析 overlaps，识别简单读段、复杂读段等类型，输出 FASTA 与候选元数据。

4. **run_blast**  
   使用 BLAST 将候选序列比对回参考基因组，获取 identity、覆盖度与基因组坐标。

5. **um_classify**  
   对 BLAST 结果进行解析，按单一来源或多重来源划分为 UeccDNA / MeccDNA，同时保留未分类记录供后续处理。

6. **cecc_build**  
   针对未分类结果进行多段识别，构建复杂环状结构（CeccDNA）的片段组合及连接信息。

### 3.2 处理阶段（步骤 7–10）

7. **umc_process**  
   统一整理 U/M/C 三类结果，生成标准化 CSV 与 FASTA，并统计对应数量。

8. **cd_hit**  
   对 FASTA 进行 90% identity 去冗余，减少重复候选。

9. **ecc_dedup**  
   整理并去重坐标信息，生成规范化 BED/CSV，确保重复记录合并到统一 ID。

10. **read_filter**  
    依据过滤规则筛选确认的 eccDNA reads，输出与 downstream 兼容的 FASTA。

### 3.3 推断阶段（步骤 11–13）

11. **minimap2**  
    检查并生成参考基因组 `.mmi` 索引，用于 Cresil 或 Cyrcular。若索引缺失会自动运行 `minimap2 -d`。

12. **ecc_inference**  
    - 优先调用 Cresil：对 reads 进行 trim/identify，若缺失 `.fai` 会尝试执行 `samtools faidx` 自动补建。  
    - 若 Cresil 不可用且检测到 Cyrcular，则使用后者作为备选。  
    输出推断的简单/嵌合 eccDNA。

13. **iecc_curator**  
    对推断结果进行清洗、格式化，生成 CSV/FASTA，并与既有标准列名称对齐。

### 3.4 整合阶段（步骤 14–16）

14. **ecc_unify**  
    汇总确认与推断结果，计算交并集指标，产出统一的 `*_unified.csv` 与重叠统计文件。

15. **ecc_summary**  
    - 合并所有 FASTA，生成 `<prefix>_all.fasta`；
    - 调用 `EccSummary` 计算样本统计并输出 HTML/TXT 报告；
    - 缺失的统计文件会自动填充默认值以保证后续步骤可运行。

16. **ecc_packager**  
    按固定结构复制并重命名文件，形成最终产物目录（Confirmed_U/M/C、Inferred_eccDNA、报告等），同时清理临时文件（除非启用 `--keep-tmp`）。

---

## 4. 运行时注意事项

- **外部依赖** 请确认 `tidehunter`, `blastn`, `cd-hit-est`, `minimap2`, `samtools`, `cresil`（或 `cyrcular`）均在 `PATH` 内。启动时会通过 `ExternalTool` 基类检查。
- **检查点机制** 每个步骤开始前都会写入 `<prefix>.checkpoint`，一旦失败可通过 `--resume` 从最近步骤继续；`--force` 可清除状态后重跑。
- **配置继承** 所有参数均由 `circleseeker.config.Config` 驱动，可通过 YAML 或 CLI 覆写。文档中的默认值等同于仓库内置配置。
- **调试工具** `circleseeker --debug --show-steps` 可查看当前步骤状态；`show-checkpoint` 子命令能列出历史执行记录。

---

## 5. 深入阅读

- `src/circleseeker/core/pipeline.py` 主管线实现与状态管理
- `src/circleseeker/config.py` 配置结构及默认值
- `src/circleseeker/modules/` 各内部模块的详细算法
- `tests/unit/` 单元测试示例，有助于理解输入输出契约
