# CircleSeeker CLI 使用手册（v0.9.1）

本手册针对 CircleSeeker 0.9.1 的命令行界面，涵盖常用参数、调试选项、运行时行为与输出结构。默认安装后可通过 `circleseeker` 或 `CircleSeeker` 两个入口调用。

---

## 1. 快速开始

```bash
circleseeker -i reads.fasta -r reference.fa -o results/
```

基本参数：

- `-i, --input PATH` 必选；输入 HiFi reads（FASTA）
- `-r, --reference PATH` 必选；参考基因组 FASTA

常用可选项：

- `-o, --output PATH` 输出目录，默认 `circleseeker_output`
- `-p, --prefix TEXT` 输出文件前缀，默认 `sample`
- `-t, --threads INT` 线程数，默认 8
- `--keep-tmp` 保留运行时的临时目录（`.tmp_work`）

---

## 2. 日志与调试

- `-n, --noise` 提升日志级别；一次为 INFO，两次为 DEBUG
- `--debug` 开启调试模式，同时解锁高级参数与隐藏子命令
- `-h, --help` 显示帮助并退出
- `-v, --version` 打印版本信息

---

## 3. 高级选项（需 `--debug`）

在 `--debug` 模式下，以下选项才会被接受：

- `--start-from INT` 从指定步骤开始执行（1 基索引）
- `--stop-at INT` 执行到指定步骤后停止（包含该步骤）
- `--resume` 从上次检查点恢复
- `--force` 忽略检查点，重新执行所有步骤
- `--generate-config` 输出默认配置 YAML 并退出
- `--show-steps` 查看 16 个步骤及状态
- `--dry-run` 仅展示计划执行的操作，不实际运行
- `--skip-report` 跳过 HTML 报告生成
- `--log-output PATH` 额外写入日志文件

---

## 4. 隐藏子命令（仅在 `--debug` 模式显示）

| 子命令 | 说明 |
|--------|------|
| `run` | 与顶级命令等价的兼容入口，保留旧参数别名 |
| `init-config` | 将默认配置写入指定路径 |
| `show-checkpoint` | 查看现有运行的检查点信息 |
| `validate` | 检查安装与依赖环境 |
| `benchmark` | 执行性能基准测试（需要额外输入数据） |

调用示例：

```bash
circleseeker --debug init-config -o config.yaml
circleseeker --debug show-checkpoint -o results/ -p sample
```

---

## 5. 执行生命周期

1. **临时目录** 运行时所有中间文件写入 `<output>/.tmp_work/`。完成后根据 `--keep-tmp` 选择是否保留。
2. **配置与检查点** 每次运行会保存 `config.yaml` 与 `<prefix>.checkpoint`，方便恢复或复查。
3. **自动索引** 
   - 若找不到参考基因组 `.mmi` 索引，会调用 `minimap2 -d` 自动生成；
   - 若 Cresil 推断阶段缺少 `.fai`，会尝试通过 `samtools faidx` 自动创建。
4. **最终输出** `ecc_packager` 会将确认与推断结果、报告和汇总拷贝到 `<output>/<prefix>_*/` 目录结构，布局与 README 中的“Output Files”一致。

---

## 6. 16 个步骤一览

| 序号 | 名称 | 作用 |
|-----|------|------|
| 1 | `make_blastdb` | 基于参考基因组构建 BLAST 数据库 |
| 2 | `tidehunter` | 检测 HiFi reads 中的串联重复 |
| 3 | `tandem_to_ring` | 将重复结构转换为候选环状序列 |
| 4 | `run_blast` | 候选片段比对至参考基因组 |
| 5 | `um_classify` | 区分 UeccDNA 与 MeccDNA |
| 6 | `cecc_build` | 识别复杂类型（CeccDNA） |
| 7 | `umc_process` | 生成 U/M/C FASTA 及汇总 |
| 8 | `cd_hit` | 去除冗余序列 |
| 9 | `ecc_dedup` | 合并并标准化坐标 |
|10 | `read_filter` | 过滤确认的 eccDNA reads |
|11 | `minimap2` | 构建参考索引 / 生成 Cyrcular 所需比对 |
|12 | `ecc_inference` | Cresil 优先的推断流程（Cyrcular 备用） |
|13 | `iecc_curator` | 整理推断结果，生成 FASTA/CSV |
|14 | `ecc_unify` | 合并确认与推断数据 |
|15 | `ecc_summary` | 统计并生成 HTML/TXT 报告（含合并 FASTA） |
|16 | `ecc_packager` | 复制/重命名产物，形成最终目录 |

运行中可通过 `circleseeker --debug --show-steps` 查看当前状态与历史完成情况。

---

## 7. 输出结构概览

完成后 `<output>/` 目录中典型内容如下：

```
<output>/
├── <prefix>.checkpoint
├── <prefix>_merged_output.csv
├── <prefix>_summary.txt
├── <prefix>_report.html
├── <prefix>_Confirmed_UeccDNA/
├── <prefix>_Confirmed_MeccDNA/
├── <prefix>_Confirmed_CeccDNA/
├── <prefix>_Inferred_eccDNA/
└── logs/（若启用额外日志）
```

各子目录的文件命名与 README 中的“Output Files”章节保持一致。

---

## 8. 使用示例

```bash
# 默认运行
circleseeker -i sample.hifi.fasta -r hg38.fa -o results/ -p sample

# 使用 YAML 覆写配置，并保留临时目录
circleseeker -i sample.hifi.fasta -r hg38.fa -c configs/sample.yaml --keep-tmp

# 调试模式下查看步骤与恢复执行
circleseeker --debug --show-steps
circleseeker --debug --resume -i sample.hifi.fasta -r hg38.fa -o results/

# 生成默认配置模板
circleseeker --debug --generate-config > default_config.yaml
```

---

## 9. 常见问题

- **命令报错 “需要 --debug 才能使用某选项”**：请加上 `--debug` 再运行。
- **参考基因组缺少索引**：程序会尝试自动创建 `.mmi` 和 `.fai`；如失败，请手工运行 `minimap2 -d ref.fa.mmi ref.fa` 或 `samtools faidx ref.fa`。
- **中途终止想要恢复**：确认输出目录下存在 `<prefix>.checkpoint`，加 `--resume` 重新运行即可。

---

如需更多背景，请参阅 `docs/Pipeline_Modules.md` 与主仓库 README。若发现文档或行为不一致，欢迎提交 issue。 
