# 模拟验证与召回率基准（U/M/C）

本仓库提供一套“合成数据”端到端验证，用于回归测试 CircleSeeker 在 **candidate 对齐 → U/M 分类 → C 识别** 这条链路上的表现，并输出 read-level 的 recall/precision/F1。

## 相关脚本

- `tests/simulation/eccdna_simulator.py`：生成模拟参考基因组 + U/Mecc/Cecc 三类 reads + ground truth CSV
- `tests/simulation/run_validation.py`：一键运行 *模拟 → minimap2 对齐 → U/M 分类 → Cecc 检测 → 指标计算*
- `tests/simulation/validation_metrics.py`：读取 ground truth 与 `uecc.csv/mecc.csv/cecc.csv`，计算指标（以 read_id 为主）
- `tests/simulation/run_recall_benchmark.py`：多次模拟（不同 seed）重复运行并汇总 recall（均值/标准差）

## Unclassified 是什么？

在 `tests/simulation/run_validation.py` 里：

1. `UMeccClassifier.run()` 会输出 `uecc_df / mecc_df / unclassified_df`。
2. `unclassified_df` 会保存成 `unclassified.csv`，**并作为 CeccBuild 的输入**，用于后续 C 识别。

因此：

- `unclassified.csv` 主要是 **“U/M 阶段没被判定为 U 或 M 的 reads 的集合”**，其中包含：
  - 大部分真实的 Cecc reads（它们本来就需要下一步 `cecc_build` 来识别）
  - 少量 U/Mecc 边界样本（例如 U 的 second-locus 证据偏强导致被 veto，或 Mecc 未达到 multi-locus 覆盖阈值）
  - 一些对齐质量较差的行（`quality_category`/`Gap_Percentage` 等字段可辅助排查）
- **是否“最终漏检”**应以 `uecc.csv/mecc.csv/cecc.csv` 是否包含该 read 为准：只有三者都没有命中的 reads，才算最终未分类/未检出。

## 为什么 recall 会波动？

在相同代码与相同环境下，**同一个 seed 的结果是可复现的**；不同 seed 会改变“模拟数据难度”，从而带来 recall 波动。常见原因包括：

- 不同 seed 生成的 U/Mecc/Cecc 的长度、拷贝数、跨染色体比例等分布略有差异，使更多/更少样本落在阈值边界（边界样本更容易留在 `Unclassified`）。
- U 的“唯一性约束”（second-locus 上限等）偏保守：边界样本数量变化会直接体现在 U/M 的 FN 上。

建议用 **多次重复 + 均值/标准差** 的方式报告模拟 recall，而不是追求单次 100%。

## 如何运行（推荐）

重复跑 3 次（默认 1000U/1000M/1000C）并汇总：

```bash
python tests/simulation/run_recall_benchmark.py --runs 3
```

输出会写到 `.tmp_work/simulation_benchmarks/.../benchmark_summary.json`（默认被 git 忽略）。

只跑一次并生成 `validation_report.json`：

```bash
python tests/simulation/run_validation.py \
  --simulation-dir .tmp_work/sim_one/simulation \
  --results-dir .tmp_work/sim_one/results \
  --num-uecc 1000 --num-mecc 1000 --num-cecc 1000 \
  --seed 42
```

## 参考结果（示例）

在一组 1000U/1000M/1000C、seed=42/43/44 的三次运行中（read-level）：

- `Uecc` recall：均值约 `0.976`
- `Mecc` recall：均值约 `0.987`
- `Cecc` recall：均值约 `0.999`
- overall recall：均值约 `0.987`

注：这是“模拟回归基线”的示例值，具体数值会随参数/版本变化；更重要的是长期趋势与波动范围是否可控。

---

## Pytest 回归基线（可回归）

为了把模拟验证纳入可回归测试，本仓库提供一个**小规模、固定 seed** 的回归基线：

- 基线文件：`tests/simulation/baselines/synthetic_regression.json`
- 测试用例：`tests/simulation/test_simulation_regression_baseline.py`
- 更新脚本：`tests/simulation/update_synthetic_baseline.py`

运行：

```bash
pytest -q tests/simulation/test_simulation_regression_baseline.py
```

> 说明：该用例使用“合成对齐 TSV”（不依赖外部 minimap2）来回归 U/M/C 分类逻辑，断言 overall 的 recall/precision/F1 不发生大幅退化。

## 一键更新基线

当你**有意**调整了算法/参数并希望更新基线期望值时：

```bash
python tests/simulation/update_synthetic_baseline.py
```

脚本会重写 `tests/simulation/baselines/synthetic_regression.json` 中的 `expected.overall` 字段。

## 负控/反例：假阳性上限（重点关注 U）

负控用例用于约束“高重复/低复杂度导致的低 MAPQ 假阳性”：

- 通过构造 `MAPQ=0` 的 full-coverage 对齐，确保不会产生**高置信度**的 U 调用（例如 `confidence_score >= 0.8`）。
- 对应测试：`tests/simulation/test_simulation_regression_baseline.py` 中的 `test_negative_control_low_mapq_has_no_high_confidence_calls`。
