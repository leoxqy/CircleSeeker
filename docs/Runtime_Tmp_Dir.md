# 使用 RAM 临时目录 (tmpfs) 加速 CircleSeeker

本文档说明如何将 CircleSeeker 的临时文件放置在 RAM 文件系统（如 `/dev/shm`）中，
以减少慢速磁盘上的 I/O 时间。

## 适用场景

- 大量 I/O 操作的步骤（LAST、minimap2 + samtools sort、CD-HIT、大型 FASTA/TSV 文件）
- HDD 或网络存储成为瓶颈时

## 方法一：使用 --turbo 模式（推荐）

最简单的方式是使用 `--turbo` 命令行参数：

```bash
circleseeker -i reads.fasta -r ref.fa -o results/ --turbo -t 16
```

**Turbo 模式特性**：
- 自动检查 `/dev/shm` 空间（默认需要 10GB）
- 在 `/dev/shm` 中创建临时目录
- 在输出目录创建 `.tmp` 软链接，方便访问中间文件
- 空间不足时自动降级到普通模式
- 仅在 Linux 系统上有效（macOS 无 `/dev/shm`）

**清理行为**：
- 默认：运行结束后删除 `/dev/shm` 内容和软链接
- `--keep-tmp`：将文件从 `/dev/shm` 移动到 `output_dir/tmp/` 保留

**配置文件方式**：

```yaml
runtime:
  turbo_mode: true
  turbo_min_space_gb: 10.0  # 最小空间要求 (GB)
  keep_tmp: false
```

## 方法二：手动配置（高级）

1) 在 `/dev/shm` 中创建任务专用的临时目录

```bash
export CS_TMPDIR="/dev/shm/cs_${SLURM_JOB_ID:-$$}"
mkdir -p "$CS_TMPDIR"
```

2) 将 `runtime.tmp_dir` 设置为该绝对路径

```yaml
runtime:
  tmp_dir: /dev/shm/cs_JOBID
  keep_tmp: true
```

注意：YAML 值不会进行 shell 变量展开。请将 `cs_JOBID` 替换为步骤 1 中创建的实际路径。

3) 保持 `output_dir` 在磁盘上，以确保最终结果持久化

```bash
circleseeker -i reads.fasta -r ref.fa -o results/ -c config.yaml -t 16
```

## 清理

`/dev/shm` 中的临时文件在重启或任务结束后不会持久化。
如果任务失败，建议手动清理：

```bash
rm -rf /dev/shm/cs_${SLURM_JOB_ID}
```

## 限制和注意事项

- `/dev/shm` 是易失性存储。除非您准备好丢失最终结果，否则不要将 `output_dir` 放在那里。
- 如果 `/dev/shm` 空间用完，流水线在写入临时文件时会失败。

## LAST 索引处理

LAST 数据库（CeccBuild 使用）现在会在**参考基因组旁边**创建
（例如 `ref.fa` → `ref.lastdb.*`），并在多次运行间重复使用。
这种一次性索引显著加快了对同一参考基因组的重复分析。

LAST 比对的临时文件（查询序列、比对结果）仍然写入 `runtime.tmp_dir`，
可以从 RAM 存储中受益。

## 快速检查

```bash
df -hT /dev/shm
```
