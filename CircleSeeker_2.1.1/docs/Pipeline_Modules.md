# CircleSeeker Pipeline Modules - Detailed Algorithmic Documentation

## Overview

CircleSeeker v2.1.1 consists of 16 interconnected modules that work together to detect and characterize eccDNA from HiFi sequencing data. This document provides in-depth algorithmic details for internal modules.

## Core Detection Modules

### 1. TideHunter Module (`tidehunter`)
**Purpose**: Detect tandem repeats in HiFi reads that may represent rolling circle amplification products.

**Type**: External tool

**Key Parameters**:
- Minimum period: 60bp
- Maximum period: 100,000bp
- Consensus calling enabled

**Output**: Tandem repeat sequences with consensus

### 2. Tandem-to-Ring Module (`tandem_to_ring`)
**Purpose**: Convert tandem repeat structures into putative circular DNA sequences using graph-based overlap detection.

**Core Algorithm - Graph-based Overlap Detection**:

1. **Overlap Graph Construction**:
   - Builds a directed graph where nodes represent reads
   - Edges represent overlaps between reads (minimum 3bp overlap)
   - Uses NetworkX for graph representation
   - Applies IntervalTree for efficient overlap detection

2. **Read Classification System**:
   - **CtcR (Circle-to-Circle Recombination)**: Reads with overlaps indicating chimeric circles
     - CtcR-perfect: Clean overlap patterns
     - CtcR-inversion: Contains inverted segments
     - CtcR-hybrid: Mixed overlap characteristics
   - **SimpleRead**: Single circular unit without complex overlaps
   - **ComplexRead**: Contains multiple overlapping segments

3. **Processing Pipeline**:
   ```python
   # Key algorithmic steps:
   - Parse TideHunter consensus sequences
   - Build overlap graph using suffix-prefix matching
   - Classify reads based on graph topology
   - Extract circular consensus sequences
   - Apply length filters (100bp - 1Mbp)
   ```

4. **Quality Control**:
   - Filters sequences by copy number (min 2 copies)
   - Validates circular topology through overlap verification
   - Removes low-quality consensus sequences

**Key Data Structures**:
- `ConsensusInfo`: Stores repeat unit, copy number, and sequence metadata
- `ReadNode`: Graph node containing read ID, sequence, and overlap information
- `OverlapGraph`: NetworkX DiGraph with weighted edges for overlap lengths

**Output**: FASTA file with validated circular sequences and classification CSV

### 3. BLAST Module (`run_blast`)
**Purpose**: Map circular candidates back to the reference genome.

**Type**: External tool (BLAST+)

**Key Features**:
- Uses BLAST+ for high-sensitivity alignment
- Identifies genomic origin of eccDNA
- Calculates coverage and identity metrics

**Output**: BLAST alignment results with genomic coordinates

## Classification Modules

### 4. UM Classify Module (`um_classify`)
**Purpose**: Classify eccDNA into Unique (U) or Multiple (M) categories using sophisticated BLAST hit analysis.

**Core Classification Algorithm**:

1. **BLAST Processing Pipeline**:
   ```python
   # Three-stage filtering:
   - Stage 1: Identity filter (>85% identity)
   - Stage 2: Coverage filter (>80% query coverage)
   - Stage 3: Length ratio filter (0.95 < ratio < 1.05)
   ```

2. **Overlap Detection Algorithm**:
   - **Query Overlap**: Detects overlapping regions in query coordinates
   - **Subject Overlap**: Identifies overlapping genomic regions
   - Uses interval arithmetic for precise overlap calculation:
     ```python
     overlap = min(end1, end2) - max(start1, start2)
     overlap_fraction = overlap / min(length1, length2)
     ```

3. **Classification Decision Tree**:
   ```
   if single_blast_hit:
       → UeccDNA (unique origin)
   elif multiple_hits_same_location:
       → UeccDNA (redundant mapping)
   elif hits_on_different_chromosomes:
       → MeccDNA (multi-origin)
   elif overlapping_query_regions > 5%:
       → MeccDNA (chimeric structure)
   else:
       → Analyze gap patterns
   ```

4. **Gap Analysis for Complex Cases**:
   - Calculates gap percentage in query coverage
   - Identifies potential rearrangement breakpoints
   - Gap threshold: 5% for M classification

**Key Metrics Calculated**:
- `match_degree`: Percentage of query covered by BLAST hits
- `gap_percentage`: Percentage of query not covered
- `copy_number`: Estimated from tandem repeat count
- `hit_distribution`: Spatial distribution of BLAST hits

**Output**:
- `um_classify.uecc.csv`: Single-origin eccDNA with genomic coordinates
- `um_classify.mecc.csv`: Multi-origin eccDNA with all hit locations

### 5. CECC Build Module (`cecc_build`)
**Purpose**: Identify and reconstruct Complex eccDNA (CeccDNA) using sweep-line algorithm and pattern recognition.

**Core Algorithm - Sweep-Line Overlap Detection**:

1. **Preprocessing Phase**:
   ```python
   # Data preparation:
   - Parse BLAST results with multiple hits per query
   - Filter hits by identity (>85%) and coverage
   - Sort hits by query position for sweep-line
   ```

2. **Sweep-Line Algorithm**:
   ```python
   # Sweep through query coordinates:
   events = [(start, 'open', hit), (end, 'close', hit)]
   active_set = IntervalTree()

   for position, event_type, hit in sorted_events:
       if event_type == 'open':
           # Check overlaps with active set
           overlaps = active_set.overlap(hit.range)
           if overlaps:
               mark_as_complex_junction(hit)
       else:
           active_set.remove(hit)
   ```

3. **Circular Pattern Recognition**:
   - **Head-Tail Junction Detection**:
     ```python
     if (first_segment.q_start < 100 and
         last_segment.q_end > query_length - 100):
         # Potential circular junction
         validate_circular_continuity()
     ```

   - **Multi-Segment Validation**:
     - Checks if segments form a valid circular path
     - Validates junction points between segments
     - Ensures no gaps > 10% between segments

4. **Complex Structure Assembly**:
   ```python
   # Assembly algorithm:
   1. Group hits by query_id
   2. Order segments by query position
   3. Identify junction types:
      - 'head': First segment in circle
      - 'tail': Last segment in circle
      - 'middle': Internal segments
   4. Calculate total circle length
   5. Validate structural integrity
   ```

5. **Quality Scoring**:
   - Junction quality score based on overlap precision
   - Segment coverage score
   - Structural complexity score (number of segments)

**Key Data Structures**:
- `BlastHit`: Stores alignment details with genomic coordinates
- `CircularSegment`: Represents one segment of complex eccDNA
- `JunctionMap`: Maps segment connections in circular structure

**Output**: CSV with multi-segment eccDNA, including:
- Segment coordinates and order
- Junction types and quality scores
- Total circle length and complexity metrics

### 6. UMC Process Module (`umc_process`)
**Purpose**: Integrate and process all three eccDNA types (U/M/C) with sophisticated clustering algorithms.

**Core Algorithms**:

1. **U-Type Clustering (Single Location)**:
   ```python
   # Location-based clustering:
   signature = f"{chr}:{start}-{end}"
   clusters = group_by_signature(sequences)

   # Aggregate metrics:
   - Sum copy numbers across cluster
   - Average gap percentages
   - Merge read lists
   ```

2. **M-Type Clustering (Multiple Locations)**:
   ```python
   # Multi-location signature generation:
   locations = [(chr1, start1, end1), (chr2, start2, end2), ...]
   signature = ";".join(sorted_locations)

   # Complex matching:
   - All locations must match for clustering
   - Order-independent comparison
   - Tolerance: ±10bp per coordinate
   ```

3. **C-Type Clustering (Complex Segments)**:
   ```python
   # Ordered segment signature:
   segments = order_by_circle_position(segments)
   signature = ";".join([f"{s.chr}:{s.start}-{s.end}"
                         for s in segments])

   # Structural validation:
   - Preserve segment order
   - Validate junction continuity
   - Check circular topology
   ```

4. **Sequence Extraction Algorithm**:
   ```python
   def extract_ring_sequence(seq, start, length):
       # Handle circular wrapping:
       if start > length:
           start = (start - length) - 1

       end = start + length
       if end <= seq_length:
           return seq[start:end]
       else:
           # Wrap around for circular sequence
           part1 = seq[start:]
           part2 = seq[:(end - seq_length)]
           return part1 + part2
   ```

5. **X-Type (Unclassified) Detection**:
   - Identifies sequences not classified as U/M/C
   - Extracts first half of sequence for further analysis
   - Maintains original read IDs for traceability

**Key Features**:
- Memory-efficient sequence library with fallback strategies
- Parallel processing for large datasets
- Comprehensive FASTA ID sanitization
- Cluster size tracking and member aggregation

**Output**:
- Processed CSVs with cluster information
- FASTA files for each eccDNA type
- Cluster statistics and metadata

## Quality Control Modules

### 7. CD-HIT Module (`cd_hit`)
**Purpose**: Remove redundant sequences and cluster similar eccDNA using external CD-HIT tool.

**Type**: External tool

**Parameters**:
- Identity threshold: 90%
- Coverage threshold: 80%
- Memory efficient algorithm
- Word size: 5 (for DNA sequences)

**Output**: Non-redundant eccDNA set with cluster representatives

### 8. ECC Deduplication Module (`ecc_dedup`)
**Purpose**: Perform type-specific deduplication using CD-HIT clusters and organize genomic coordinates.

**Core Deduplication Algorithm**:

1. **Cluster-based Deduplication**:
   ```python
   # Type-specific strategies:
   - Uecc: Keep 1 representative per cluster
   - Mecc: Keep all segments of representative
   - Cecc: Keep all segments of representative

   # Metadata aggregation:
   - Sum copy_numbers from cluster members
   - Average match_degrees
   - Merge read lists (unique only)
   ```

2. **Natural Sorting Algorithm**:
   ```python
   # Extract prefix and numeric parts:
   pattern = r'([A-Za-z]+)(\d+)'
   sort_key = (type_prefix, int(number))
   # Results in: U1, U2, U10 (not U1, U10, U2)
   ```

3. **Coordinate Standardization**:
   - Convert to 0-based coordinates
   - Normalize strand notation (+/-)
   - Validate coordinate integrity
   - Format regions as chr:start-end

**Output**:
- Deduplicated coordinate tables
- Standardized eccDNA IDs
- Aggregated cluster metadata

### 9. Read Filter Module (`read_filter`)
**Purpose**: Separate confirmed eccDNA reads from remaining reads using classification-based filtering.

**Core Filtering Algorithm**:

1. **Classification-based Filtering**:
   ```python
   # Filter criteria:
   CTCR_CLASSES = {
       "CtcR-perfect",
       "CtcR-inversion",
       "CtcR-hybrid"
   }

   # Process:
   1. Load TandemToRing classifications
   2. Identify CtcR variant reads
   3. Filter from FASTA files
   4. Maintain read statistics
   ```

2. **Memory-Efficient Processing**:
   ```python
   # Stream processing for large files:
   for record in fasta_iterator:
       if record.id not in filtered_set:
           write_to_output(record)
       update_statistics()
   ```

3. **Multi-File Handling**:
   - Combines multiple input FASTA files
   - Removes duplicates across files
   - Preserves sequence quality information

**Statistics Tracking**:
- Total reads processed
- CtcR reads filtered
- Retained read percentage
- Classification distribution

**Output**:
- Filtered FASTA file
- Optional samtools index
- Filtering statistics report

## Inference Modules

### 10. Minimap2 Alignment (`minimap2`)
**Purpose**: High-accuracy alignment of filtered reads.

**Type**: External tool

**Settings**:
- HiFi-optimized parameters
- Splice-aware alignment
- Secondary alignment reporting

**Output**: BAM alignment files

### 11. Cyrcular Calling Module (`cyrcular_calling`)
**Purpose**: Detect circular DNA using split-read analysis with statistical validation.

**Type**: External tool integration

**Components**:
- **Cyrcular**: External tool for split-read circular DNA detection
- **Varlociraptor**: Statistical validation (FDR=0.05)
- **BCFtools**: Variant processing and filtering

**Statistical Framework**:
- False Discovery Rate control: 0.05
- Minimum split reads: 3
- Probability present threshold: 0.9

**Output**: Inferred eccDNA with statistical support (BCF/VCF format)

### 12. Inferred ECC Curator (`iecc_curator`)
**Purpose**: Curate and validate inferred eccDNA with sophisticated scoring algorithms.

**Core Curation Algorithm**:

1. **Duplicate Resolution**:
   ```python
   def select_best_duplicate(group):
       # Scoring weights:
       score = (prob_present * 0.5 +
               (num_split_reads / max_splits) * 0.3 +
               (1 - prob_artifact) * 0.2)
       return highest_scoring_entry
   ```

2. **Simple vs Chimeric Classification**:
   ```python
   if segment_count == 1:
       # Simple circle: single continuous region
       type = "IUeccDNA"
       extract_single_region()
   else:
       # Chimeric circle: multiple segments
       type = "ICeccDNA"
       for segment in segments:
           assign_junction_role(segment)
   ```

3. **Junction Role Assignment**:
   - **head**: First segment (connects to tail)
   - **middle**: Internal segments
   - **tail**: Last segment (connects to head)

4. **Quality Metrics**:
   - `prob_present`: Statistical confidence
   - `prob_artifact`: Likelihood of false positive
   - `num_split_reads`: Supporting evidence count
   - `hifi_abundance`: HiFi read abundance

**Strand Determination**:
```python
if start_pos <= end_pos:
    strand = "+"  # Forward
else:
    strand = "-"  # Reverse
    swap(start_pos, end_pos)
```

**Output**:
- Curated simple eccDNA table
- Curated chimeric eccDNA table
- Optional FASTA sequences from reference

## Integration Modules

### 13. ECC Unify Module (`ecc_unify`)
**Purpose**: Merge confirmed and inferred eccDNA results with redundancy detection.

**Core Unification Algorithm**:

1. **Reciprocal Overlap Detection**:
   ```python
   def reciprocal_overlap_ok(a_start, a_end, b_start, b_end,
                            thr=0.99, tol=10):
       overlap = min(a_end, b_end) - max(a_start, b_start)
       if overlap <= 0:
           return False

       # Calculate reciprocal overlap fractions
       overlap_frac_a = overlap / (a_end - a_start)
       overlap_frac_b = overlap / (b_end - b_start)

       # With tolerance for boundary matching
       if abs(a_start - b_start) <= tol and
          abs(a_end - b_end) <= tol:
           return True

       return overlap_frac_a >= thr and overlap_frac_b >= thr
   ```

2. **Chromosome-Indexed Search**:
   ```python
   # Build index for efficient queries:
   chr_index = {
       'chr1': [(start1, end1, id1), (start2, end2, id2), ...],
       'chr2': [...]
   }
   # Sort by start position for sweep-line search
   ```

3. **Redundancy Detection Pipeline**:
   ```python
   # For simple eccDNA:
   1. Check inferred vs confirmed UeccDNA
   2. Apply reciprocal overlap test
   3. Mark redundant with 99% overlap

   # For chimeric eccDNA:
   1. Compare segment signatures
   2. All segments must match
   3. Preserve segment order
   ```

4. **Conflict Resolution**:
   - Confirmed eccDNA takes precedence
   - Higher confidence score wins
   - Preserve maximum information

5. **Final Renumbering**:
   ```python
   # Systematic renumbering:
   - Confirmed: U1, U2, M1, M2, C1, C2...
   - Inferred: IU1, IU2, IM1, IM2, IC1, IC2...
   - Maintain type consistency
   ```

**Integration Strategy**:
- Remove redundant inferred entries
- Merge non-redundant results
- Standardize column formats
- Add source annotations (Confirmed/Inferred)

**Output**:
- Unified eccDNA catalog with all types
- Redundancy report
- Final eccDNA statistics

### 14. ECC Summary Module (`ecc_summary`)
**Purpose**: Generate comprehensive statistics and summaries.

**Analysis Components**:

1. **Statistical Aggregation**:
   - eccDNA count by type (U/M/C)
   - Size distribution analysis
   - Chromosomal distribution
   - Quality score distributions

2. **Genomic Feature Analysis**:
   - Gene overlap detection
   - Repeat element enrichment
   - Hotspot identification

3. **Quality Metrics**:
   - Average match degree
   - Copy number distributions
   - Read support statistics

**Output**:
- Summary tables (CSV/TSV)
- Statistical reports
- Distribution plots data

### 15. ECC Packager Module (`ecc_packager`)
**Purpose**: Organize and package final results.

**Organization Strategy**:

1. **Directory Structure**:
   ```
   final_results/
   ├── sequences/
   │   ├── confirmed_eccDNA.fasta
   │   └── inferred_eccDNA.fasta
   ├── tables/
   │   ├── all_eccDNA.csv
   │   ├── UeccDNA.csv
   │   ├── MeccDNA.csv
   │   └── CeccDNA.csv
   ├── statistics/
   │   └── summary_stats.txt
   └── README.txt
   ```

2. **File Compression**:
   - Large files compressed with gzip
   - Index files for quick access
   - Checksum generation

3. **Documentation Generation**:
   - README with file descriptions
   - Method summary
   - Parameter documentation

**Output**: Well-organized final results package

### 16. Report Generator
**Purpose**: Create interactive HTML reports.

**Type**: External reporting module

**Report Sections**:
- Executive summary
- eccDNA statistics
- Quality metrics
- Visualization plots
- Detailed tables

**Output**: Interactive HTML report

## Module Dependencies

```
TideHunter → Tandem-to-Ring → BLAST
                                ↓
                           UM Classify → CECC Build
                                ↓            ↓
                           UMC Process ←─────┘
                                ↓
                            CD-HIT → ECC Dedup
                                ↓
                           Read Filter
                                ↓
                      Minimap2 → Cyrcular Calling
                                ↓
                          IECC Curator
                                ↓
                           ECC Unify
                                ↓
                      ECC Summary → Report
                                ↓
                         ECC Packager
```

## Algorithm Complexity Analysis

### Time Complexity
- **Tandem-to-Ring**: O(n² log n) for overlap graph construction
- **UM Classify**: O(n × m) where n = queries, m = BLAST hits
- **CECC Build**: O(n log n) for sweep-line algorithm
- **UMC Process**: O(n log n) for clustering operations
- **ECC Unify**: O(n × m) for overlap detection

### Space Complexity
- **Sequence Library**: O(n) for n sequences
- **Overlap Graph**: O(n²) worst case, O(n) typical
- **CD-HIT Clusters**: O(n) for n sequences
- **Chromosome Index**: O(n) for n regions

## Performance Optimization

- **Parallel Processing**: Multiple modules support multi-threading
- **Memory Management**: Stream processing for large files
- **Indexing**: Interval trees and hash maps for fast lookups
- **Caching**: Intermediate results cached for resumption
- **Compression**: Gzip compression for large outputs

## Error Handling

All internal modules include:
- Input validation with informative error messages
- Checkpoint saving for pipeline resumption
- Graceful degradation for missing optional inputs
- Comprehensive logging with debug levels
- Recovery mechanisms for transient failures