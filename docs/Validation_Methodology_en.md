# CircleSeeker Validation Methodology

**Version**: v1.1.0
**Last Updated**: 2026-01-15

This document describes the simulation-based validation methodology used to evaluate CircleSeeker's eccDNA detection performance.

---

## 1. Overview

CircleSeeker validation uses a **synthetic data approach** that:

1. Generates a simulated reference genome
2. Creates three types of synthetic eccDNA reads with known ground truth
3. Runs the minimal validation chain (minimap2 alignment → um_classify → cecc_build)
4. Compares predictions against ground truth to calculate performance metrics

This methodology allows precise measurement of recall, precision, and F1 scores because the true eccDNA locations are known exactly.

---

## 2. Simulated eccDNA Types

### 2.1 Unique eccDNA (UeccDNA)

**Definition**: Simple circular DNA from a single, unique genomic location.

**Simulation Method**:
- Random genomic region selected (length: 200-5,000 bp)
- Region must map uniquely to the reference genome
- Creates a single circular molecule

**Proportion**: 70% of total eccDNA

### 2.2 Multimap eccDNA (MeccDNA)

**Definition**: Circular DNA from repetitive genomic regions that map to multiple loci.

**Simulation Method**:
- Create "repeat families" with 2-10+ genomic copies
- Each family has a consensus sequence (100-500 bp unit)
- Introduce ~0.5% divergence between copies to simulate real repeats
- eccDNA reads derive from these multi-copy regions

**Proportion**: 15% of total eccDNA

### 2.3 Chimeric eccDNA (CeccDNA)

**Definition**: Complex circular DNA composed of multiple non-contiguous genomic segments.

**Simulation Method**:
- 2-4 genomic segments joined in a single circle
- Segments can be from the same chromosome (intra-chromosomal) or different chromosomes (inter-chromosomal, ~30%)
- Each segment: 100-1,000 bp

**Proportion**: 15% of total eccDNA

---

## 3. Reference Genome Construction

The simulated reference genome is **auto-scaled** based on the number of eccDNAs to maintain realistic repeat density:

| Parameter | Default Value | Description |
|-----------|---------------|-------------|
| Chromosomes | 5 | Number of chromosomes |
| Target repeat coverage | 1.5% | Proportion of genome covered by repeats |
| Max repeat coverage | 30% | Upper limit when auto-scaling is disabled |

**Formula**:
```
genome_size = (num_mecc × copies_per_family × avg_unit_length) / target_coverage
```

This ensures that larger-scale tests don't artificially inflate repeat density.

---

## 4. HiFi Read Simulation

Simulated reads mimic PacBio HiFi characteristics:

| Parameter | Value |
|-----------|-------|
| Error rate | 0.1% (HiFi typical) |
| Read structure | Tandem repeats of eccDNA sequence |
| Copy number | Variable (2-10×) |

---

## 5. Validation Pipeline

```
┌─────────────────────────────────────────────────────────────────┐
│                    Simulation Phase                              │
├─────────────────────────────────────────────────────────────────┤
│  1. Generate reference genome (FASTA)                           │
│  2. Create UeccDNA reads + ground truth                         │
│  3. Create MeccDNA reads + ground truth                         │
│  4. Create CeccDNA reads + ground truth                         │
│  5. Output: simulated_reads.fa, ground_truth.csv                │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                    Detection Phase                               │
├─────────────────────────────────────────────────────────────────┤
│  1. Run the minimal validation chain                            │
│     - minimap2 → alignment                                      │
│     - U/M classification                                        │
│     - CeccDNA detection (LAST preferred; graph fallback)        │
│  2. Output: uecc.csv, mecc.csv, cecc.csv                        │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                    Evaluation Phase                              │
├─────────────────────────────────────────────────────────────────┤
│  1. Match predictions to ground truth by read_id                │
│  2. Calculate TP, FP, FN for each eccDNA type                   │
│  3. Compute recall, precision, F1, accuracy                     │
│  4. Output: validation_report.json                              │
└─────────────────────────────────────────────────────────────────┘
```

---

## 6. Evaluation Metrics

### 6.1 Matching Strategy

Predictions are matched to ground truth using **read_id** as the primary key:

- **True Positive (TP)**: Predicted read_id exists in ground truth with correct type
- **False Positive (FP)**: Predicted read_id not in ground truth, or wrong type
- **False Negative (FN)**: Ground truth read_id not in predictions

### 6.2 Metric Definitions

| Metric | Formula | Description |
|--------|---------|-------------|
| **Recall** | TP / (TP + FN) | Proportion of true eccDNAs detected |
| **Precision** | TP / (TP + FP) | Proportion of predictions that are correct |
| **F1 Score** | 2 × (P × R) / (P + R) | Harmonic mean of precision and recall |
| **Accuracy** | TP / (TP + FP + FN) | Detection accuracy (TN not applicable) |

**Note**: True Negatives (TN) are not used because the search space (all possible genomic regions) is essentially infinite.

---

## 7. Batch Validation Protocol

For comprehensive validation, we use a **gradient testing** approach:

### 7.1 Test Configuration

| Parameter | Value |
|-----------|-------|
| Scales | 100, 500, 1,000, 3,000, 5,000, 10,000, 50,000, 100,000 |
| Runs per scale | 10 |
| Random seeds | 42-51 (sequential) |
| eccDNA ratio | 70% U : 15% M : 15% C |

### 7.2 Test Execution

```bash
python tests/simulation/batch_validation.py \
  --scales 100,500,1000,3000,5000,10000,50000,100000 \
  --runs 10 \
  --base-seed 42 \
  --threads 32
```

### 7.3 Result Aggregation

Results are aggregated at multiple levels:
1. **Per-run**: Individual TP/FP/FN for each test
2. **Per-scale**: Average metrics across 10 runs at each scale
3. **Overall**: Combined metrics across all 80 runs

---

## 8. Key Scripts

| Script | Purpose |
|--------|---------|
| `tests/simulation/eccdna_simulator.py` | Generate synthetic reference + reads + ground truth |
| `tests/simulation/run_validation.py` | Run single end-to-end validation |
| `tests/simulation/batch_validation.py` | Run gradient tests across multiple scales |
| `tests/simulation/validation_metrics.py` | Calculate TP/FP/FN and metrics |

---

## 9. Reproducibility

All simulations are fully reproducible:

- **Fixed random seeds**: Each run uses a deterministic seed (42-51)
- **Deterministic genome**: Same seed produces identical reference genome
- **Deterministic reads**: Same seed produces identical eccDNA reads

To reproduce a specific test:

```bash
python tests/simulation/run_validation.py \
  --simulation-dir output/simulation \
  --results-dir output/results \
  --num-uecc 700 --num-mecc 150 --num-cecc 150 \
  --seed 42
```

---

## 10. Limitations

1. **Simplified genome**: Simulated reference lacks real genomic complexity (transposons, segmental duplications)
2. **Idealized reads**: Simulated reads have uniform error distribution
3. **No biological noise**: Real samples contain non-eccDNA circular molecules
4. **Type boundaries**: Edge cases at U/M/C classification boundaries may differ from real data

Despite these limitations, synthetic validation provides a rigorous, reproducible benchmark for algorithm development and regression testing.

---

## 11. Understanding "Unclassified" Reads

In the validation pipeline, `UMeccClassifier.run()` outputs `uecc_df`, `mecc_df`, and `unclassified_df`. The `unclassified.csv` file contains reads that were **not classified as U or M** at the U/M stage:

- Most true CeccDNA reads (expected to be handled by subsequent `cecc_build` step)
- Borderline U/M samples (e.g., U uniqueness veto or M multi-locus coverage insufficient)
- Low-quality alignment rows

**Final detection status** should be determined by checking `uecc.csv`, `mecc.csv`, and `cecc.csv`. Only reads absent from all three are truly undetected.

---

## 12. Why Recall Fluctuates Across Seeds

With identical code and environment, results are **reproducible for the same seed**. Different seeds change the "difficulty profile" of simulated data:

- Different seeds produce varying length/copy-number/inter-chromosomal distributions
- More/fewer samples fall near decision boundaries
- U's uniqueness constraints are conservative: boundary sample counts directly affect U/M FN

**Recommendation**: Report recall using **multiple runs (mean ± std)** rather than pursuing single 100% runs.

---

## 13. Pytest Regression Baseline

For CI/local regression testing, the repo includes a **small, fixed-seed** baseline:

| Component | Path |
|-----------|------|
| Baseline file | `tests/simulation/baselines/synthetic_regression.json` |
| Test case | `tests/simulation/test_simulation_regression_baseline.py` |
| Update script | `tests/simulation/update_synthetic_baseline.py` |

**Run**:
```bash
pytest -q tests/simulation/test_simulation_regression_baseline.py
```

**Update baseline** (when intentionally changing algorithm/parameters):
```bash
python tests/simulation/update_synthetic_baseline.py
```

---

## 14. Negative Control: False Positive Constraints

The negative control constrains low-MAPQ false positives in repetitive/low-complexity contexts:

- Construct full-coverage alignments with `MAPQ=0`
- Assert no high-confidence U calls are produced (`confidence_score >= 0.8`)
- Test: `test_negative_control_low_mapq_has_no_high_confidence_calls`

---

## 15. References

- Validation scripts: `tests/simulation/`
- Results report: `docs/Validation_Report.md`
- Performance figures: `docs/images/validation_combined.png`
