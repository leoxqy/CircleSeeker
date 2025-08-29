# CircleSeeker2 测试指南

## 基本测试命令

你提供的命令是正确的！以下是完整的测试方案：

### 1. 基础测试（推荐首次运行）

```bash
python -m circleseeker2 run \
  -i step1_Test.HeLa_10k.fasta \
  -r chm13v2.0.fa \
  -o circleseeker2_complete_test \
  -p HeLa_10k \
  -t 8 \
  --verbose
```

### 2. 查看流水线步骤

```bash
python -m circleseeker2 run --show-steps
```

输出应该显示16个步骤：
- Step 0-13: 原有的eccDNA分析步骤
- Step 14: playbill (报告生成)
- Step 15: propmaster (文件整理)

### 3. 分步测试

#### a. 只运行前置步骤 (0-5)
```bash
python -m circleseeker2 run \
  -i step1_Test.HeLa_10k.fasta \
  -r chm13v2.0.fa \
  -o test_partial \
  -p HeLa_10k \
  -t 8 \
  --stop-at 5 \
  --verbose
```

#### b. 恢复运行 (从步骤6继续)
```bash
python -m circleseeker2 run \
  -i step1_Test.HeLa_10k.fasta \
  -r chm13v2.0.fa \
  -o test_partial \
  -p HeLa_10k \
  -t 8 \
  --start-from 6 \
  --verbose
```

#### c. 只运行报告和整理步骤
```bash
python -m circleseeker2 run \
  -i step1_Test.HeLa_10k.fasta \
  -r chm13v2.0.fa \
  -o test_report_only \
  -p HeLa_10k \
  -t 8 \
  --start-from 14 \
  --verbose
```

### 4. 高级测试选项

#### a. 包含XeccDNA分析
```bash
python -m circleseeker2 run \
  -i step1_Test.HeLa_10k.fasta \
  -r chm13v2.0.fa \
  -o test_with_xecc \
  -p HeLa_10k \
  -t 8 \
  --enable-xecc \
  --verbose
```

#### b. 强制重新运行（忽略检查点）
```bash
python -m circleseeker2 run \
  -i step1_Test.HeLa_10k.fasta \
  -r chm13v2.0.fa \
  -o test_force \
  -p HeLa_10k \
  -t 8 \
  --force \
  --verbose
```

#### c. 保留临时文件
```bash
python -m circleseeker2 run \
  -i step1_Test.HeLa_10k.fasta \
  -r chm13v2.0.fa \
  -o test_keep_tmp \
  -p HeLa_10k \
  -t 8 \
  --keep-tmp \
  --verbose
```

#### d. 跳过某些步骤
```bash
python -m circleseeker2 run \
  -i step1_Test.HeLa_10k.fasta \
  -r chm13v2.0.fa \
  -o test_skip \
  -p HeLa_10k \
  -t 8 \
  --skip-report \
  --skip-organize \
  --verbose
```

### 5. 使用配置文件

创建配置文件 `config.yaml`:
```yaml
input_file: step1_Test.HeLa_10k.fasta
reference: chm13v2.0.fa
output_dir: circleseeker2_config_test
prefix: HeLa_10k

performance:
  threads: 8
  max_memory: 16G

quality:
  min_eccdna_size: 100
  max_eccdna_size: 1000000
  min_identity: 99.0

runtime:
  keep_tmp: false
  log_level: INFO
```

运行：
```bash
python -m circleseeker2 run -c config.yaml --verbose
```

### 6. 调试模式（最详细日志）

```bash
python -m circleseeker2 run \
  -i step1_Test.HeLa_10k.fasta \
  -r chm13v2.0.fa \
  -o test_debug \
  -p HeLa_10k \
  -t 8 \
  -vv \
  --log-file debug.log
```

## 预期输出

成功运行后，应该在输出目录看到：

```
circleseeker2_complete_test/
├── organized_results/          # 整理后的结果
│   ├── HeLa_10k_UeccDNA.fa
│   ├── HeLa_10k_MeccDNA.fa
│   ├── HeLa_10k_CeccDNA.fa
│   ├── HeLa_10k_eccDNA_report.html
│   └── HeLa_10k_file_summary.txt
├── HeLa_10k.checkpoint        # 检查点文件
├── HeLa_10k_eccDNA_report.html # 分析报告
└── [其他中间文件]
```

## 故障排除

### 1. 如果运行失败
- 检查日志文件了解错误详情
- 使用 `--force` 强制重新运行
- 检查输入文件路径是否正确

### 2. 内存不足
- 减少线程数: `-t 4`
- 调整配置中的 `max_memory` 参数

### 3. 某个步骤失败
- 使用 `--start-from` 从失败步骤重新开始
- 检查该步骤的依赖是否满足

### 4. 查看检查点状态
```bash
python -m circleseeker2 show-checkpoint -o circleseeker2_complete_test
```

## 性能优化建议

1. **线程数**: 根据CPU核心数设置，通常设为核心数-2
2. **内存**: 建议至少16GB RAM
3. **磁盘空间**: 预留输入文件大小的10倍空间

## 快速验证

最简单的验证命令：
```bash
# 只运行前3个步骤快速测试
python -m circleseeker2 run \
  -i step1_Test.HeLa_10k.fasta \
  -r chm13v2.0.fa \
  -o quick_test \
  -p test \
  -t 4 \
  --stop-at 3 \
  --verbose
```

如果前3步成功，说明环境配置正确，可以运行完整流程。