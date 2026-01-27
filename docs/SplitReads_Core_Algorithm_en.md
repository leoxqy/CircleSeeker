# SplitReads-Core Algorithm (v1.1.0)

This document describes the built-in SplitReads-Core module for split-read based eccDNA inference.

> Note: The SplitReads-Core algorithm is inspired by [CReSIL](https://github.com/Peppermint-Lab/CReSIL), with optimizations for PacBio HiFi data.

---

## 1. Overview

SplitReads-Core is an eccDNA detection module based on split-read alignment patterns and breakpoint evidence. It analyzes read alignment patterns against the reference genome and uses graph algorithms to detect circular DNA structures.

**Key Features**:
- Optimized for HiFi data (uses `map-hifi` preset)
- No external tool dependencies (uses Python's mappy for alignment)
- Graph-based circular detection algorithm
- Supports both simple (single-segment) and chimeric (multi-segment) eccDNA detection

---

## 2. Two-Phase Pipeline

### 2.1 Trim Phase

**Purpose**: Align raw reads to reference genome and extract high-quality alignment segments.

**Workflow**:

```
Raw reads (FASTA)
    ↓
minimap2 alignment (map-hifi preset)
    ↓
Filter low-quality alignments (mapq >= 20)
    ↓
Merge/trim overlapping/gap regions
    ↓
trim_df (alignment segment table)
```

**Output Fields**:

| Field | Description |
|-------|-------------|
| readid | Original read ID |
| q_start / q_end | Query start/end positions |
| r_start / r_end | Reference start/end positions |
| ref | Chromosome name |
| strand | Alignment direction (+/-) |
| mapq | Mapping quality |
| order | Segment order within the same read |

**Key Parameters**:

| Parameter | Description | Default |
|-----------|-------------|---------|
| preset | minimap2 preset | map-hifi |
| mapq | Minimum mapping quality | 20 |
| allow_gap | Allowed gap size | 10 bp |
| allow_overlap | Allowed overlap size | 10 bp |

### 2.2 Identify Phase

**Purpose**: Identify circular structures using graph algorithms based on trim phase results.

**Workflow**:

```
trim_df (alignment segment table)
    ↓
Calculate genome coverage
    ↓
Merge adjacent regions → merge regions
    ↓
Filter low-coverage regions (depth >= 5)
    ↓
Analyze breakpoint directions
    ↓
Build breakpoint graph
    ↓
Detect circular subgraphs
    ↓
eccDNA candidate list
```

---

## 3. Graph Algorithm Details

### 3.1 Merge Region Generation

Merge continuous genomic regions with sufficient coverage into "merge regions":

```
Conditions:
- Average coverage >= average_depth (default 5)
- Region length >= min_region_size (default 200 bp)

Output format:
mergeid = "{chrom}_{start}_{end}"
```

### 3.2 Breakpoint Direction Analysis

For multiple alignment segments from the same read, analyze strand transitions between adjacent segments:

```python
def check_breakpoint_direction(df_sorted_by_order):
    """
    Check strand directions of adjacent alignment segments.

    Returns: (region1, region2, strand_transition, is_valid)

    strand_transition format: "{strand1}_{strand2}"
    Examples: "+_+", "+_-", "-_+", "-_-"
    """
```

### 3.3 Building Breakpoint Graph

```python
G = nx.MultiDiGraph()

# Nodes: merge regions
# Edges: connections between adjacent alignment segments

for each read:
    for adjacent_pair in read.sorted_alignments:
        if is_valid_breakpoint(pair):
            G.add_edge(region1, region2, weight=1)

# Filter low-support edges
graph_filtered = {k: v for k, v in graph.items()
                  if v >= breakpoint_depth}  # default 5
```

### 3.4 Cycle Detection

```python
# Get connected components
subgraphs = nx.connected_components(G.to_undirected())

for subgraph in subgraphs:
    nodes = list(subgraph.nodes())

    if len(nodes) == 1:
        # Single node: simple eccDNA (self-loop)
        is_cyclic = has_self_loop(subgraph)

    elif len(nodes) == 2:
        # Two nodes: check bidirectional edge
        is_cyclic = has_bidirectional_edge(subgraph)

    else:
        # Multiple nodes: use cycle_basis for detection
        cycles = nx.cycle_basis(subgraph)
        is_cyclic = len(cycles) > 0
```

### 3.5 Strand Closure Validation

For detected cycles, verify if strand directions can form a valid closure:

```python
def check_strand_closure(strand_transitions):
    """
    Verify if strand transition sequence can form a closed circle.

    Rules:
    - +_+ and -_- maintain direction
    - +_- and -_+ flip direction
    - Valid closure requires even number of flips
    """
    flips = count(trans for trans in transitions
                  if trans in ["+_-", "-_+"])
    return flips % 2 == 0
```

---

## 4. Output Format

### 4.1 eccDNA_final.txt

| Field | Description |
|-------|-------------|
| id | eccDNA ID (e.g., ec1, ec2) |
| merge_region | Region list (e.g., "chr1_100_200_+,chr2_300_400_-") |
| merge_len | Total length |
| num_region | Number of segments |
| ctc | Is tandem repeat type |
| numreads | Supporting read count |
| totalbase | Total covered bases |
| coverage | Estimated coverage |

### 4.2 Classification Logic

```
num_region == 1  →  Simple (IUeccDNA)
num_region >= 2  →  Chimeric (ICeccDNA)
```

---

## 5. Parameter List

| Parameter | Description | Default |
|-----------|-------------|---------|
| **Alignment** | | |
| preset | minimap2 preset | map-hifi |
| mapq | Minimum mapping quality | 20 |
| exclude_chrs | Excluded chromosomes | "" |
| **Trim** | | |
| allow_gap | Allowed gap size | 10 bp |
| allow_overlap | Allowed overlap size | 10 bp |
| **Identify** | | |
| min_region_size | Minimum region size | 200 bp |
| overlap_check_size | Edge check size | 50 bp |
| breakpoint_depth | Minimum breakpoint support | 5 |
| average_depth | Minimum average coverage | 5.0 |

---

## 6. Comparison with CtcReads-Caller

| Feature | CtcReads-Caller | SplitReads-Core |
|---------|-----------------|-----------------|
| Evidence Type | Tandem repeats (TideHunter) | Split-read patterns |
| Output State | Confirmed | Inferred |
| Use Case | High coverage, clear rolling circle | Low coverage, no clear repeats |
| Detection Method | Coverage model + LAST | Breakpoint graph + NetworkX |

---

## 7. Dependencies

- **mappy**: Python minimap2 bindings
- **networkx**: Graph algorithm library
- **pybedtools**: BED file operations (requires system bedtools)
- **intervaltree**: Interval tree operations

---

## 8. References

- [CReSIL](https://github.com/Peppermint-Lab/CReSIL) - Original algorithm implementation
- [minimap2](https://github.com/lh3/minimap2) - Sequence alignment
- [NetworkX](https://networkx.org/) - Graph algorithms
