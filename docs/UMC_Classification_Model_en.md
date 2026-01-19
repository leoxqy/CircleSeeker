# U/M/C eccDNA Coverage Model and CeccBuild Documentation

This document describes the U/M classification coverage model and CeccBuild v4 detection strategy implemented in CircleSeeker. U/M classification is based on coverage and locus clustering; Cecc detection is primarily driven by LAST-based repeat pattern detection, with a graph-based chain model as fallback.

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

### 4.1 CeccBuild v4 (LAST-first)

CeccBuild v4 uses LAST for high-precision alignment of doubled sequences, looking for repeat patterns where "first half and second half align to the same genomic location" to identify complex circles.

Key inputs:
- `reference_fasta`
- `tandem_to_ring.fasta` (doubled sequences)

Key parameters:
- `min_match_degree` (or `theta_chain`)
- `overlap_threshold`
- `min_segments`
- `half_query_buffer`

### 4.2 Graph-based Fallback (when LAST unavailable)

When LAST is missing or key inputs are unavailable, falls back to chain coverage model:

- Build candidate chains based on alignment segments
- Control adjacent gaps via `tau_gap`
- Ensure coverage via `min_match_degree`
- Require at least `min_segments` different loci

---

## 5. Parameter Reference (Defaults)

| Parameter | Description | Default |
|-----------|-------------|---------|
| `theta_u` / `theta_m` / `theta_full` | U/M coverage thresholds | 0.95 |
| `theta_u2_max` | U second locus coverage upper limit | 0.05 |
| `theta_locus` | Locus reciprocal overlap threshold | 0.95 |
| `pos_tol_bp` | Locus boundary tolerance | 50 |
| `mapq_u_min` | U minimum MAPQ threshold | 0 |
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
