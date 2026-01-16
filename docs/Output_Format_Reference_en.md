# CircleSeeker Output Format Reference

This manual provides complete format and field definitions for all CircleSeeker v0.10.3 output files.

---

## Output Directory Structure

After a completed run, the output directory contains the following structure:

```
<output>/
├── <prefix>_merged_output.csv          # All eccDNA summary table
├── <prefix>_report.html                # Interactive HTML report
├── <prefix>_summary.txt                # Summary statistics
├── <prefix>_Confirmed_UeccDNA/         # Confirmed simple eccDNA
│   ├── <prefix>_UeccDNA_C.csv
│   └── <prefix>_UeccDNA_C.fasta
├── <prefix>_Confirmed_MeccDNA/         # Confirmed multi-copy eccDNA
│   ├── <prefix>_MeccDNA_C.csv
│   └── <prefix>_MeccDNA_C.fasta
├── <prefix>_Confirmed_CeccDNA/         # Confirmed complex eccDNA
│   ├── <prefix>_CeccDNA_C.csv
│   ├── <prefix>_CeccDNA_C.fasta
│   ├── <prefix>_CeccJunctions.bedpe
│   └── <prefix>_CeccSegments.core.csv
└── <prefix>_Inferred_eccDNA/           # Inferred eccDNA
    ├── <prefix>_UeccDNA_I.csv
    ├── <prefix>_UeccDNA_I.fasta
    ├── <prefix>_chimeric.csv
    └── <prefix>_CeccDNA_I.fasta
```

> Note: The Confirmed_* directories are produced by **CtcReads-Caller**, while Inferred_eccDNA is produced by **SplitReads-Caller**.  
> **CtcReads** refers to reads carrying **Ctc** (**C**oncatemeric **t**andem **c**opies) signals (tracked as CtcR-* classes in `tandem_to_ring.csv`).

---

## 1. Merged Output File (merged_output.csv)

The main output file containing summary information for all confirmed and inferred eccDNA.

### Field Definitions

| Field | Type | Description |
|-------|------|-------------|
| `eccDNA_id` | string | Unique identifier (e.g., UeccDNA1, MeccDNA1, CeccDNA1) |
| `original_id` | string | Original ID before renumbering |
| `Regions` | string | Genomic coordinates (format: `chr:start-end`; CeccDNA may have multiple, separated by `;`) |
| `Strand` | string | DNA strand (`+` or `-`; CeccDNA may have multiple) |
| `Length` | int | eccDNA length (bp) |
| `eccDNA_type` | string | Classification type: `UeccDNA`, `MeccDNA`, `CeccDNA` |
| `State` | string | Detection state: `Confirmed` (CtcReads-Caller) or `Inferred` (SplitReads-Caller) |
| `Seg_total` | int | Number of segments (multiple for CeccDNA, 1 for U/M) |
| `Hit_count` | int | Number of genomic hits (MeccDNA specific) |
| `confidence_score` | float | Confidence score [0,1] (higher is more reliable) |
| `query_cov_best` | float | Best locus/chain query coverage [0,1] |
| `query_cov_2nd` | float | Second-best locus/chain query coverage [0,1] |
| `mapq_best` | int | Best alignment MAPQ value [0-60] |
| `identity_best` | float | Best alignment identity percentage [0-100] |
| `low_mapq` | bool | Low MAPQ flag (mapq_best < 20) |
| `low_identity` | bool | Low identity flag (identity_best < 95) |

> **Note**: Confidence and evidence fields are primarily populated for **Confirmed** entries; **Inferred** entries may have empty values for these fields.

### Example

```csv
eccDNA_id,original_id,Regions,Strand,Length,eccDNA_type,State,Seg_total,Hit_count,confidence_score,query_cov_best,query_cov_2nd,mapq_best,identity_best,low_mapq,low_identity
UeccDNA1,U001,chr1:10000-10500,+,500,UeccDNA,Confirmed,1,1,0.98,0.99,0.02,60,99.5,False,False
MeccDNA1,M001,chr2:20000-20800,+,800,MeccDNA,Confirmed,1,3,0.85,0.97,0.95,55,98.2,False,False
CeccDNA1,C001,chr1:5000-5200;chr3:8000-8300,+;-,500,CeccDNA,Confirmed,2,2,0.92,0.96,0.03,58,99.0,False,False
```

---

## 2. Confirmed eccDNA Output

### 2.1 UeccDNA (Simple eccDNA)

#### CSV File Format (*_UeccDNA_C.csv)

| Field | Type | Description |
|-------|------|-------------|
| `eccDNA_id` | string | Unique identifier |
| `chr` | string | Chromosome name |
| `start0` | int | 0-based start coordinate |
| `end0` | int | 0-based end coordinate (half-open interval) |
| `strand` | string | DNA strand (+/-) |
| `length` | int | Sequence length (bp) |
| `reads` | string | Supporting read names (semicolon-separated) |
| `copy_number` | float | Copy number estimate |
| `match_degree` | float | Match degree [0-100] |

#### FASTA File Format

```
>UeccDNA1 chr1:10000-10500(+) length=500
ATCGATCG...
```

### 2.2 MeccDNA (Multi-copy eccDNA)

#### CSV File Format (*_MeccDNA_C.csv)

| Field | Type | Description |
|-------|------|-------------|
| `eccDNA_id` | string | Unique identifier |
| `chr` | string | Chromosome name |
| `start0` | int | 0-based start coordinate |
| `end0` | int | 0-based end coordinate |
| `strand` | string | DNA strand |
| `length` | int | Sequence length |
| `reads` | string | Supporting read names |
| `copy_number` | float | Copy number estimate |
| `hit_index` | int | Current hit index (1-based) |
| `hit_count` | int | Total number of hits |

### 2.3 CeccDNA (Complex eccDNA)

#### Main CSV File (*_CeccDNA_C.csv)

| Field | Type | Description |
|-------|------|-------------|
| `eccDNA_id` | string | Unique identifier |
| `chr` | string | Segment chromosome |
| `start0` | int | Segment 0-based start |
| `end0` | int | Segment 0-based end |
| `strand` | string | Segment strand direction |
| `length` | int | Segment length |
| `seg_index` | int | Segment index (1-based) |
| `seg_total` | int | Total number of segments |
| `junction_role` | string | Junction role: `head`, `body`, `tail` |
| `reads` | string | Supporting read names |
| `copy_number` | float | Copy number estimate |

#### Segment Details File (*_CeccSegments.core.csv)

Contains detailed information for each CeccDNA segment, with the same fields as the main CSV.

#### BEDPE Junction File (*_CeccJunctions.bedpe)

Standard BEDPE format recording junction information between CeccDNA segments:

| Column | Description |
|--------|-------------|
| chrom1 | Segment 1 chromosome |
| start1 | Segment 1 start (0-based) |
| end1 | Segment 1 end |
| chrom2 | Segment 2 chromosome |
| start2 | Segment 2 start |
| end2 | Segment 2 end |
| name | Junction name (eccDNA_id) |
| score | Score |
| strand1 | Segment 1 strand direction |
| strand2 | Segment 2 strand direction |

---

## 3. Inferred eccDNA Output

### 3.1 Inferred UeccDNA (*_UeccDNA_I.csv)

Fields are similar to confirmed UeccDNA, but confidence fields may be empty.

### 3.2 Inferred Chimeric (*_chimeric.csv)

Contains potential chimeric eccDNA detected by the inference engine.

### 3.3 Inferred CeccDNA (*_CeccDNA_I.fasta)

FASTA format containing inferred complex eccDNA sequences.

---

## 4. Statistical Reports

### Summary File (*_summary.txt)

Plain text format containing run statistics:

```
CircleSeeker Analysis Summary
=============================
Sample: sample_name
Date: 2024-01-15 10:30:00

Total eccDNA detected: 150
  - UeccDNA (Confirmed): 80
  - MeccDNA (Confirmed): 30
  - CeccDNA (Confirmed): 15
  - Inferred: 25

Length distribution:
  - <500 bp: 45
  - 500-1000 bp: 60
  - 1000-5000 bp: 35
  - >5000 bp: 10
```

### HTML Report (*_report.html)

Interactive report containing:
- eccDNA type distribution pie chart
- Length distribution histogram
- Chromosome distribution bar chart
- Interactive data tables

---

## 5. Coordinate System

CircleSeeker uses **0-based half-open interval** coordinate system internally and in output files:

- `start0`: Inclusive start position (0-based)
- `end0`: Exclusive end position (0-based)
- Interval notation: `[start0, end0)`

### Example

```
Sequence:  A T C G A T C G
Position:  0 1 2 3 4 5 6 7

Interval [2, 5) represents: C G A (positions 2, 3, 4)
Length = end0 - start0 = 5 - 2 = 3
```

### Conversion with 1-based Coordinates

- BED format: Use directly (BED is natively 0-based)
- FASTA header: Displayed as `start0+1` to `end0` (1-based inclusive interval)
- Integration with UCSC/Ensembl browsers: `start1 = start0 + 1`

---

## 6. Column Naming Standards

CircleSeeker uses a unified column naming standard (`ColumnStandard`):

| Standard Name | Legacy Aliases | Description |
|---------------|----------------|-------------|
| `eccDNA_id` | - | Main identifier |
| `chr` | `eChr` | Chromosome |
| `start0` | `eStart0`, `eStart` | 0-based start |
| `end0` | `eEnd0`, `eEnd` | 0-based end |
| `strand` | `eStrand` | Strand direction |
| `length` | `eLength`, `Length` | Length |
| `reads` | `eReads`, `readName` | Read names |
| `copy_number` | `copyNum` | Copy number |
| `match_degree` | `MatDegree` | Match degree |
| `eccdna_type` | `eClass` | eccDNA type |
| `state` | `State` | Detection state |
| `hit_count` | `Hit_count` | Hit count |
| `seg_total` | `Seg_total` | Total segments |

---

## 7. Related Documentation

- [CLI Reference](CLI_Reference_en.md) - Command-line argument details
- [Configuration Reference](Configuration_Reference_en.md) - Configuration options details
- [Pipeline Modules](Pipeline_Modules_en.md) - 16-step pipeline details
