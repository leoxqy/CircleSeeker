# CircleSeeker Batch Validation Report

**Test Date**: 2026-01-14 ~ 2026-01-15
**Version**: v1.1.0
**Test Environment**: Linux (fat1 server), Python 3.12.11, 32 threads
**minimap2**: v2.30-r1287

> 注意：本报告为历史批量验证快照，基于当时的模拟数据与参数。算法/默认参数更新后，结果可能不同。
> 如需当前版本的验证结果，请使用 `tests/simulation/` 重新运行并生成报告。

---

## 0. Supplement: CtcReads-Caller Core Validation (Simulation)

**Test Date**: 2026-01-19  
**Script**: `scripts/run_ctcreads_core_validation.py` (wraps `tests/simulation/batch_validation.py`)  
**Scales**: 100, 500, 1,000, 3,000  
**Runs per scale**: 3  
**Threads**: 8  
**minimap2 preset**: `sr`  
**Output**: `server_validation_results/ctcreads_core_validation_20260119_190214/`

### Summary (U/M/C aggregated)

| Metric | Value |
|--------|-------|
| Runs | 12 |
| Total ground truth | 13,800 |
| Total predictions | 13,800 |
| True Positives | 13,800 |
| False Positives | 0 |
| False Negatives | 0 |
| Recall | 100.00% |
| Precision | 100.00% |
| F1 Score | 100.00% |
| Accuracy (TP/(TP+FP+FN)) | 100.00% |

### By type (aggregated across all runs)

| Type | GT | Pred | TP | FP | FN | Recall | Precision | F1 |
|------|----|------|----|----|----|--------|-----------|----|
| Uecc | 9,660 | 9,660 | 9,660 | 0 | 0 | 100.00% | 100.00% | 100.00% |
| Mecc | 2,070 | 2,070 | 2,070 | 0 | 0 | 100.00% | 100.00% | 100.00% |
| Cecc | 2,070 | 2,070 | 2,070 | 0 | 0 | 100.00% | 100.00% | 100.00% |

---

## 1. Test Overview

This report documents the comprehensive validation of CircleSeeker's eccDNA detection performance across multiple data scales using gradient testing methodology.

> **Methodology**: For detailed information on how simulated data is generated, how the three eccDNA types (U/M/C) are constructed, and how evaluation metrics are calculated, see [Validation_Methodology.md](Validation_Methodology.md).

### Test Configuration

| Parameter | Value |
|-----------|-------|
| Scales tested | 100, 500, 1,000, 3,000, 5,000, 10,000, 50,000, 100,000 |
| Runs per scale | 10 |
| Random seeds | 42-51 |
| eccDNA composition | 70% Unique + 15% Multimap + 15% Cluster |
| Total runtime | 689.9 minutes (~11.5 hours) |

---

## 2. Overall Results Summary

| Metric | Value |
|--------|-------|
| **Successful runs** | **80/80 (100%)** |
| Total ground truth | 1,696,000 |
| Total True Positives | 1,695,965 |
| **Total False Positives** | **22** |
| **Total False Negatives** | **35** |
| **Overall Recall** | **99.9979%** |
| **Overall Precision** | **99.9987%** |
| **Overall F1 Score** | **99.9983%** |
| **Overall Accuracy** | **99.9966%** |

---

## 3. Results by Scale

| Scale | Runs | TP | FP | FN | Recall | Precision | F1 | Status |
|-------|------|------|-----|-----|--------|-----------|-----|--------|
| 100 | 10/10 | 1,000 | 0 | 0 | **100.00%** | **100.00%** | **100.00%** | Perfect |
| 500 | 10/10 | 5,000 | 0 | 0 | **100.00%** | **100.00%** | **100.00%** | Perfect |
| 1,000 | 10/10 | 10,000 | 0 | 0 | **100.00%** | **100.00%** | **100.00%** | Perfect |
| 3,000 | 10/10 | 29,999 | 1 | 1 | 99.997% | 99.997% | 99.997% | Excellent |
| 5,000 | 10/10 | 49,997 | 2 | 3 | 99.994% | 99.996% | 99.995% | Excellent |
| 10,000 | 10/10 | 99,999 | 1 | 1 | 99.999% | 99.999% | 99.999% | Excellent |
| 50,000 | 10/10 | 499,986 | 8 | 14 | 99.997% | 99.998% | 99.998% | Excellent |
| 100,000 | 10/10 | 999,984 | 10 | 16 | 99.998% | 99.999% | 99.999% | Excellent |

---

## 4. Performance by eccDNA Type

### 4.1 Unique eccDNA (U)

- **Ground Truth**: 70% of total
- **Detection Rate**: ~100% across all scales
- **False Negatives**: 0
- **False Positives**: Occasional (primarily from boundary cases)

### 4.2 Multimap eccDNA (M)

- **Ground Truth**: 15% of total
- **Detection Rate**: >99.99%
- **Primary error source**: Occasional FN in large-scale tests
- **Note**: Multi-mapping complexity can lead to rare misclassifications

### 4.3 Cluster eccDNA (C)

- **Ground Truth**: 15% of total
- **Detection Rate**: ~100%
- **False Negatives**: 1 (only at 100k scale, run 9)
- **Stability**: Excellent across all scales

---

## 5. Error Analysis

### 5.1 False Positive/Negative Distribution

| Scale | Avg FP/run | Avg FN/run | Primary FN Type |
|-------|------------|------------|-----------------|
| 100-1,000 | 0.00 | 0.00 | None |
| 3,000 | 0.10 | 0.10 | MeccDNA |
| 5,000 | 0.20 | 0.30 | MeccDNA |
| 10,000 | 0.10 | 0.10 | MeccDNA |
| 50,000 | 0.80 | 1.40 | MeccDNA |
| 100,000 | 1.00 | 1.60 | MeccDNA |

### 5.2 Key Observations

1. **Perfect detection at small scales**: 100-1,000 eccDNAs show 100% accuracy
2. **Minimal errors at large scales**: Even at 100,000 eccDNAs, avg FP=1.0, avg FN=1.6
3. **Error source**: FN primarily from Multimap type due to alignment complexity
4. **Excellent precision**: No scale shows precision below 99.99%

---

## 6. Visualization

Performance visualization figures are available in `images/`:

- `validation_combined.png/pdf` - Combined 4-panel summary figure
- `accuracy_summary.png/pdf` - Accuracy by scale with detailed annotations
- `fp_fn_by_scale.png/pdf` - FP/FN distribution analysis
- `metrics_by_scale.png/pdf` - Recall/Precision/F1 comparison

![Validation Summary](images/validation_combined.png)

---

## 7. Conclusions

### 7.1 Strengths

1. **Exceptional accuracy**: >99.99% accuracy across all tested scales
2. **Scalability**: Performance remains stable from 100 to 100,000 eccDNAs
3. **Low false positive rate**: Average 0.275 FP per run across all tests
4. **Robust to randomization**: Consistent results across 10 different random seeds
5. **Production ready**: 80/80 runs completed successfully

### 7.2 Recommendations

1. **Re-run for current versions**: Use `tests/simulation/batch_validation.py` to regenerate metrics for the current codebase.
2. **Large-scale analysis**: Use multiple seeds and report mean ± std to avoid single-run bias.
3. **MeccDNA attention**: For critical applications, verify Multimap classifications manually.

---

## 8. Appendix: Test Environment

| Component | Version/Specification |
|-----------|----------------------|
| Server | fat1 (Linux 3.10.0-693.el7.x86_64) |
| Python | 3.12.11 (conda-forge) |
| CircleSeeker | 1.0.0 |
| minimap2 | 2.30-r1287 |
| Threads | 32 |
| minimap2 timeout | 7200s |

---

*Report generated from batch validation run: gradient_validation_20260114_141315*
