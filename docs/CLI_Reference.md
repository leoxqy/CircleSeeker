# CircleSeeker CLI 使用手册（v1.1.0）

本手册针对 CircleSeeker 1.1.0 的命令行界面，涵盖常用参数、调试选项、运行时行为与输出结构。默认安装后可通过 `circleseeker` 或 `CircleSeeker` 两个入口调用。

---

## 1. 快速开始

```bash
circleseeker -i reads.fasta -r reference.fa -o results/
```

基本参数：

- `-i, --input PATH` 必选；输入 HiFi reads（FASTA）
- `-r, --reference PATH` 必选；参考基因组 FASTA

常用可选项：

- `-o, --output PATH` 输出目录，默认 `circleseeker_output`
- `-p, --prefix TEXT` 输出文件前缀，默认 `sample`
- `-t, --threads INT` 线程数，默认 `min(8, CPU核数×2)`
- `-c, --config PATH` 配置文件路径（YAML 格式）
- `--keep-tmp` 保留临时目录（`.tmp_work`），默认删除；可显式覆盖配置文件中的 `keep_tmp` 设置
- `--turbo` 启用 turbo 模式（使用 RAM-backed 临时目录加速 I/O）
- `--preset CHOICE` 灵敏度预设（`relaxed` / `balanced` / `strict`）
- `--show-steps` 查看 16 个步骤及状态（不执行）
- `--dry-run` 仅展示计划执行的操作，不实际运行

---

## 2. 日志与调试

- `-v, --verbose` 提升日志级别；一次为 INFO，两次为 DEBUG（默认 WARNING）
- `--debug` 解锁高级选项与隐藏子命令（不影响日志级别）
- `--help-advanced` 显示高级/调试选项说明并退出
- `-h, --help` 显示帮助并退出
- `-V, --version` 打印版本信息

---

## 3. 高级选项（需 `--debug`）

在 `--debug` 模式下，以下选项才会被接受：

- `--start-from INT` 从指定步骤开始执行（1 基索引）
- `--stop-at INT` 执行到指定步骤后停止（包含该步骤）
- `--resume` 从上次检查点恢复
- `--force` 忽略检查点，重新执行所有步骤
- `--log-file PATH` 额外写入日志文件

---

## 4. 子命令

| 子命令 | 说明 |
|--------|------|
| `init-config` | 生成默认配置文件（`--stdout` 输出到终端，`--output-file` 写入文件） |
| `show-checkpoint` | 查看现有运行的检查点信息（`-d` 指定输出目录） |
| `validate` | 检查安装与依赖环境（`--full` 可启用更完整检查） |

调用示例：

```bash
circleseeker init-config --stdout
circleseeker init-config --output-file config.yaml
circleseeker show-checkpoint -d results/ -p sample
circleseeker validate
```

---

## 5. 执行生命周期

1. **依赖检查** 启动时自动检测必需的外部工具（minimap2、samtools、cd-hit-est、lastal）。SplitReads-Core 是内置模块，仅需 `mappy` Python 包。若缺少依赖，程序会给出清晰的错误提示和安装建议。
2. **临时目录** 运行时所有中间文件写入 `<output>/.tmp_work/`（可通过 `runtime.tmp_dir` 配置，支持相对路径或绝对路径）。完成后根据 `--keep-tmp` 选择是否保留。
3. **配置与检查点** 运行期间在输出目录保存 `config.yaml` 与 `<prefix>.checkpoint`；成功完成后会自动清理。使用 `--keep-tmp` 可保留这些文件以便调试或恢复中断的运行。
4. **自动索引**
   - 若找不到参考基因组 `.mmi` 索引，会调用 `minimap2 -d` 自动生成；
   - 若 SplitReads-Core 推断阶段缺少 `.fai`，会尝试通过 `samtools faidx` 自动创建。
5. **最终输出** `ecc_packager` 会将确认与推断结果、报告和汇总拷贝到 `<output>/<prefix>_*/` 目录结构，布局与 README 中的"Output Files"一致。

---

## 6. 两条 Caller（证据来源视角）

为便于理解与展示，本项目将 eccDNA 的分析拆分为两条证据链路：

- **CtcReads**：含 **Ctc**（**C**oncatemeric **t**andem **c**opies）信号的 reads（在 `tandem_to_ring.csv` 中以 CtcR-* 分类体现）。
- **CtcReads-Caller**（步骤 1-10）：基于 CtcReads 证据产出 **Confirmed** U/M/C eccDNA。
- **SplitReads-Caller**（步骤 11-13）：基于 split-reads/junction 证据使用内置 SplitReads-Core 推断 eccDNA，产出 **Inferred** eccDNA。
- **Integration**（步骤 14-16）：对 Confirmed/Inferred 做去冗余合并、统计与打包交付。

## 7. 16 个步骤一览

> 快速定位：步骤 1-10 属于 **CtcReads-Caller**；步骤 11-13 属于 **SplitReads-Caller**；步骤 14-16 属于 **Integration**。

| 序号 | 名称 | 作用 |
|-----|------|------|
| 1 | `check_dependencies` | 检查外部工具与推断引擎可用性 |
| 2 | `tidehunter` | 检测 HiFi reads 中的串联重复 |
| 3 | `tandem_to_ring` | 将重复结构转换为候选环状序列 |
| 4 | `run_alignment` | 使用 minimap2（或 LAST）将候选片段比对至参考基因组 |
| 5 | `um_classify` | 基于覆盖度与位点聚类区分 UeccDNA / MeccDNA |
| 6 | `cecc_build` | LAST 优先的复杂 eccDNA 检测（不可用时回退图方法） |
| 7 | `umc_process` | 生成 U/M/C FASTA 及汇总 |
| 8 | `cd_hit` | 去除冗余序列 |
| 9 | `ecc_dedup` | 合并并标准化坐标 |
|10 | `read_filter` | 过滤 CtcR reads，生成推断输入 FASTA |
|11 | `minimap2` | 构建参考索引供 SplitReads-Core 使用 |
|12 | `ecc_inference` | 使用内置 SplitReads-Core 进行推断 |
|13 | `curate_inferred_ecc` | 整理推断结果（内部调用 `iecc_curator`） |
|14 | `ecc_unify` | 合并确认与推断数据，使用片段重叠算法检测冗余嵌合体 eccDNA |
|15 | `ecc_summary` | 统计并生成 HTML/TXT 报告（含合并 FASTA） |
|16 | `ecc_packager` | 复制/重命名产物，形成最终目录 |

运行中可通过 `circleseeker --show-steps` 查看步骤列表。

---

## 8. 输出结构概览

完成后 `<output>/` 目录中典型内容如下：

```
<output>/
├── <prefix>_merged_output.csv
├── <prefix>_summary.txt
├── <prefix>_report.html
├── <prefix>_Confirmed_UeccDNA/
├── <prefix>_Confirmed_MeccDNA/
├── <prefix>_Confirmed_CeccDNA/
├── <prefix>_Inferred_eccDNA/
└── logs/（若启用额外日志）
```

各子目录的文件命名与 README 中的"Output Files"章节保持一致。

---

## 9. 使用示例

```bash
# 默认运行
circleseeker -i sample.hifi.fasta -r hg38.fa -o results/ -p sample

# 使用 YAML 覆写配置，并保留临时目录
circleseeker -i sample.hifi.fasta -r hg38.fa -c configs/sample.yaml --keep-tmp

# 查看步骤（不执行）
circleseeker --show-steps

# 预览运行计划（dry run）
circleseeker --dry-run -i sample.hifi.fasta -r hg38.fa -o results/

# 使用灵敏度预设
circleseeker --preset strict -i sample.hifi.fasta -r hg38.fa -o results/

# 调试模式下恢复执行
circleseeker --debug --resume -i sample.hifi.fasta -r hg38.fa -o results/

# 生成默认配置模板
circleseeker init-config --stdout > default_config.yaml
circleseeker init-config --output-file config.yaml
```

---

## 10. 常见问题

- **命令报错 "需要 --debug 才能使用某选项"**：请加上 `--debug` 再运行。仅 `--start-from`、`--stop-at`、`--resume`、`--force`、`--log-file` 需要 `--debug`。
- **依赖检查失败**：程序会提示缺少哪些工具。请通过 conda 安装缺失的依赖：`conda install -c bioconda -c conda-forge <tool_name>`。
- **参考基因组缺少索引**：程序会尝试自动创建 `.mmi` 和 `.fai`；如失败，请手工运行 `minimap2 -d ref.fa.mmi ref.fa` 或 `samtools faidx ref.fa`。
- **中途终止想要恢复**：确认输出目录下存在 `<prefix>.checkpoint`，加 `--debug --resume` 重新运行即可。

---

如需更多背景，请参阅 `docs/Pipeline_Modules.md` 与主仓库 README。若发现文档或行为不一致，欢迎提交 issue。
