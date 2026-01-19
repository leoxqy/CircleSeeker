# CircleSeeker Output Format Reference

This document describes the primary output layout and file schemas for CircleSeeker v1.0.0.

---

## Output Directory Layout

The packaged layout produced by `ecc_packager` is:

```
<output>/
└── <prefix>/
    ├── <prefix>_Confirmed_UeccDNA/
    │   ├── <prefix>_UeccDNA_C.fasta
    │   ├── <prefix>_UeccDNA.bed
    │   └── <prefix>_UeccDNA.core.csv
    ├── <prefix>_Confirmed_MeccDNA/
    │   ├── <prefix>_MeccDNA_C.fasta
    │   ├── <prefix>_MeccSites.bed
    │   ├── <prefix>_MeccBestSite.bed
    │   └── <prefix>_MeccSites.core.csv
    ├── <prefix>_Confirmed_CeccDNA/
    │   ├── <prefix>_CeccDNA_C.fasta
    │   ├── <prefix>_CeccSegments.bed
    │   ├── <prefix>_CeccJunctions.bedpe
    │   └── <prefix>_CeccSegments.core.csv
    ├── <prefix>_Inferred_eccDNA/
    │   ├── <prefix>_UeccDNA_I.csv
    │   ├── <prefix>_UeccDNA_I.fasta
    │   ├── <prefix>_chimeric.csv
    │   └── <prefix>_CeccDNA_I.fasta
    ├── <prefix>_merged_output.csv
    ├── <prefix>_report.html
    └── <prefix>_summary.txt
```

> Note: `<prefix>_UeccDNA_I.csv` is renamed from `<prefix>_simple.csv` during packaging.

---

## 1. Merged Output (`*_merged_output.csv`)

Final summary table combining Confirmed + Inferred eccDNA.

| Field | Type | Description |
|-------|------|-------------|
| `eccDNA_id` | string | Final ID (UeccDNA1/MeccDNA1/CeccDNA1) |
| `original_id` | string | Original ID before final renumbering |
| `Regions` | string | Genomic coordinates (CeccDNA uses `;` for multiple segments) |
| `Strand` | string | Strand (CeccDNA uses `;` for multiple segments) |
| `Length` | int | eccDNA length (bp) |
| `eccDNA_type` | string | `UeccDNA` / `MeccDNA` / `CeccDNA` |
| `State` | string | `Confirmed` / `Inferred` |
| `Seg_total` | int | Segment count (U/M = 1) |
| `Hit_count` | int | Genomic hit count |

---

## 2. Confirmed Outputs

### 2.1 UeccDNA

#### `*_UeccDNA.core.csv`

| Field | Description |
|-------|-------------|
| `eccDNA_id` | Unique ID |
| `chr` / `start0` / `end0` / `strand` | 0-based coordinates and strand |
| `length` | Length (bp) |
| `match_degree` | Match degree (0-100) |
| `confidence_score` / `query_cov_best` / `query_cov_2nd` | Evidence scores (may be empty) |
| `mapq_best` / `identity_best` | Best alignment quality (may be empty) |
| `low_mapq` / `low_identity` | Low-quality flags (may be empty) |
| `copy_number` / `repeat_number` | Copy/repeat estimates |
| `eccdna_type` | `Uecc` |
| `num_merged` / `merged_from_ids` | Merge bookkeeping |
| `reads_count` / `read_name` | Read count and list |

#### `*_UeccDNA.bed`

Columns: `chrom, chromStart, chromEnd, name, score, strand, length, eccdna_type, repeat_number, match_degree`

#### `*_UeccDNA_C.fasta`

Header format:
```
>UeccDNA1|chr:start-end(strand)|length=...|repeats=...|reads=...
```

---

### 2.2 MeccDNA

#### `*_MeccSites.core.csv`

Adds the following fields on top of Uecc core fields:

| Field | Description |
|-------|-------------|
| `hit_index` | Site index (1-based) |
| `hit_count` | Total sites for the eccDNA |

#### `*_MeccSites.bed`

Columns: `chrom, chromStart, chromEnd, name, score, strand, length, eccdna_type, copy_number`

#### `*_MeccBestSite.bed`

Columns (legacy names): `chrom, chromStart, chromEnd, name, score, strand, eLength, eClass, copyNum`

#### `*_MeccDNA_C.fasta`

Header format:
```
>MeccDNA1|multi_loci:{sites}_sites|length=...|copies=...|reads=...
```

---

### 2.3 CeccDNA

#### `*_CeccSegments.core.csv`

Adds the following fields on top of Uecc core fields:

| Field | Description |
|-------|-------------|
| `seg_index` / `seg_total` | Segment index and total |
| `junction_role` | `head` / `middle` / `tail` |

#### `*_CeccSegments.bed`

Columns: `chrom, chromStart, chromEnd, name, score, strand, length, eccdna_type, copy_number, match_degree`

#### `*_CeccJunctions.bedpe`

Columns: `chrom1, start1, end1, chrom2, start2, end2, name, score, strand1, strand2`

#### `*_CeccDNA_C.fasta`

Header format:
```
>CeccDNA1|segments:{seg_total}|junctions:{chr-chain}|length=...|copies=...|reads=...
```

---

## 3. Inferred Outputs

### 3.1 `*_UeccDNA_I.csv`

Fields:
`eccDNA_id, chr, start0, end0, strand, length, eccdna_type, state, num_split_reads, prob_present, prob_artifact, hifi_abundance`

### 3.2 `*_chimeric.csv`

Fields:
`eccDNA_id, chr, start0, end0, strand, length, eccdna_type, state, seg_index, seg_total, junction_role, read_count, num_split_reads, prob_present, prob_artifact, hifi_abundance`

### 3.3 FASTA

- `*_UeccDNA_I.fasta`: `>ID|chr:start-end|length=...|type=simple`
- `*_CeccDNA_I.fasta`: `>ID|segments=N|length=...|type=chimeric|seg1:chr:start-end;...`

---

## 4. Reports

- `*_summary.txt`: plain-text summary
- `*_report.html`: HTML report
