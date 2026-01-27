# U/M/C eccDNA Coverage Model and CeccBuild Documentation

This document describes the U/M classification coverage model and CeccBuild detection strategy implemented in CircleSeeker. U/M classification is based on coverage and locus clustering; Cecc detection is primarily driven by LAST-based repeat pattern detection, with a graph-based chain model as fallback.

> Note: This document focuses on the actual logic of um_classify/cecc_build. For complete code details, please refer to `src/circleseeker/modules/um_classify.py` and `src/circleseeker/modules/cecc_build.py`.

---

## 1. Input and Notation

- `S`: TideHunter consensus sequence (one circle), length `L` (typically `consLen`).
- `Q = S + S`: Doubled query, length `2L`, used to handle circular origin uncertainty.
- Alignment output set `A = {a_i}`: Each hit `a_i` contains at least:
  - Query interval `q_i = [qs_i, qe_i)`
  - Reference locus `g_i` (chr/start/end/strand)
  - Optional quality information (identity/mapq/aln_len)

---

## 2. Locus Clustering (Same-locus Merging)

To avoid erroneous splitting due to boundary jitter, nearby genomic hits should be merged into the same locus.

Recommended rules (parameterized):

- Same `chr` and same `strand`
- Reciprocal overlap ratio `overlap_len / min(len_i, len_j) >= theta_locus`
- Optional boundary tolerance: `|start_i - start_j| <= pos_tol_bp` and `|end_i - end_j| <= pos_tol_bp`

Default recommendations: `theta_locus=0.95`, `pos_tol_bp=50`.

---

## 3. Coverage and U/M Classification

### 3.1 Coverage

For any hit set `X`, coverage is defined as:

`Cov(X) = | union proj(q_i) (i in X) | / L`

Where `proj(q_i)` is the interval set after projecting the query interval to ring coordinates `[0, L)`.

### 3.2 Uecc Classification

- `U_cov = max_l Cov(l)`
- `U_cov_2nd`: Second highest locus coverage

Classification conditions:

- `U_cov >= theta_U`
- `U_cov_2nd <= theta_U2_max`

Additional constraints (optional):

- `mapq_u_min`: Require best locus MAPQ >= threshold
- `span_ratio_min`: Require genomic span / consensus length ratio >= threshold (default 0.95)
  - Used to reduce CeccDNA being misclassified as UeccDNA
  - When eccDNA is actually composed of multiple segments, a single locus span is typically smaller than consensus length
  - Formula: `span_ratio = best_genomic_span / cons_len`
- Secondary alignment veto (prevent multi-chromosome/non-contiguous mapping):
  - `u_secondary_min_frac` / `u_secondary_min_bp`
  - `u_contig_gap_bp`
  - `u_secondary_max_ratio`
  - `u_high_coverage_threshold` / `u_high_mapq_threshold`

### 3.3 Mecc Classification

- `M_count = #{l | Cov(l) >= theta_M}`
- Classification: `M_count >= 2`

Optional filters (disabled by default):

- `mapq_m_ambiguous_threshold`: If all loci have MAPQ below threshold, reject Mecc (0 to disable)
- `mecc_identity_gap_threshold`: If gap between highest and second-highest identity exceeds threshold, reject Mecc (0 to disable)

---

## 4. Cecc Detection

### 4.1 CeccBuild Overview

CeccBuild uses LAST for high-precision alignment of doubled sequences, identifying chimeric eccDNA through graph structure analysis.

**Core Concept**:
- Input: Doubled sequence `Q = S + S` (length `2L`)
- Alignment: Use LAST (lastal | last-split) for high-precision, optimal chain alignment
- Analysis: Graph-based circular pattern detection

### 4.2 Graph-based Algorithm

#### 4.2.1 Locus Clustering

Merge adjacent alignments into loci (same as U/M classification):
- Same chromosome, same strand
- Genomic coordinates overlap or adjacent (tolerance `position_tolerance`)

#### 4.2.2 Building Directed Graph

```
G = (V, E)
V = {locus_1, locus_2, ..., locus_n}
E = {(locus_i, locus_j, strand_trans) | adjacent in query sequence}
```

- **Nodes**: Each locus as a node, carrying `(chr, start, end, strand)` information
- **Edges**: Connect adjacent loci after sorting by query coordinates
- **Strand Transitions**: Record strand direction changes between adjacent alignments (`+_+`, `+_-`, `-_+`, `-_-`)

#### 4.2.3 Cycle Detection (Doubled Sequence Pattern)

For doubled sequence `[A-B-C][A-B-C]`:
- Alignment pattern should be: `A' → B → C → A → B'` (`'` indicates partial alignment)
- Detection method: Find repeated locus sequence pattern
- Validation conditions:
  - First half and second half locus sequences should match (rotation equivalent)
  - Strand transitions must satisfy closure condition (even number of flips)

#### 4.2.4 Strand Closure Validation

```python
# Strand transition closure check
flips = count(trans for trans in transitions if trans in ["+_-", "-_+"])
valid = (flips % 2 == 0)  # Even number of flips for valid closure
```

### 4.3 LAST Alignment Pipeline

```bash
# 1. Build index
lastdb -P threads db_prefix reference.fa

# 2. Alignment + optimal chain selection
lastal -P threads db_prefix query.fa | last-split > output.maf
```

**Role of last-split**: Selects optimal alignment chain for each query region, filtering multi-mapping noise.

### 4.4 Key Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `min_match_degree` | Minimum match degree (coverage × 100) | 90.0 |
| `overlap_threshold` | Reciprocal overlap threshold | 0.95 |
| `min_segments` | Minimum segment count | 2 |
| `position_tolerance` | Locus merge position tolerance | 100 bp |
| `edge_tolerance` | Query boundary tolerance | 10 bp |
| `half_query_buffer` | Doubled sequence midpoint buffer | 50 bp |

### 4.5 Output

Each detected CeccDNA includes:
- `eccDNA_id`: Unique identifier
- `segments`: List of genomic coordinates for each segment
- `junction_roles`: Segment roles (head/middle/tail)
- `match_degree`: Match degree score
- `CeccClass`: Classification (`Cecc-InterChr` / `Cecc-IntraChr`)

---

## 5. Parameter Reference (Defaults)

| Parameter | Description | Default |
|-----------|-------------|---------|
| `theta_u` / `theta_m` / `theta_full` | U/M coverage thresholds | 0.95 |
| `theta_u2_max` | U second locus coverage upper limit | 0.05 |
| `theta_locus` | Locus reciprocal overlap threshold | 0.95 |
| `pos_tol_bp` | Locus boundary tolerance | 50 |
| `mapq_u_min` | U minimum MAPQ threshold | 0 |
| `span_ratio_min` | U genomic span / consensus length ratio minimum | 0.95 |
| `u_secondary_min_frac` / `u_secondary_min_bp` | Secondary alignment thresholds | 0.01 / 50 |
| `u_contig_gap_bp` | U contiguity threshold | 1000 |
| `u_secondary_max_ratio` | Secondary/primary coverage ratio | 0.05 |
| `u_high_coverage_threshold` | High coverage relaxation threshold | 0.98 |
| `u_high_mapq_threshold` | High MAPQ relaxation threshold | 50 |
| `mapq_m_ambiguous_threshold` | Mecc MAPQ filter | 0 (disabled) |
| `mecc_identity_gap_threshold` | Mecc identity gap filter | 0 (disabled) |
| `overlap_threshold` | Cecc reciprocal overlap threshold | 0.95 |
| `min_segments` | Cecc minimum segment count | 2 |
| `tau_gap` / `edge_tolerance` | Chain gap tolerance | 20 |
| `min_match_degree` | Cecc match degree threshold | 95.0 |
| `half_query_buffer` | Doubled sequence midpoint buffer | 50 |

---

## 6. Notes

- `delta_uc` / `epsilon_mc` and other ambiguity parameters are retained in the current implementation but do not produce separate output files; they may be used for future extensions or research purposes.
