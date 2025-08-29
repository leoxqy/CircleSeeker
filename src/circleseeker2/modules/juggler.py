"""
Juggler - Circular DNA Validation Pipeline

This module performs the final validation of circular DNA candidates by:
1. Extracting sequences and reads for each candidate region
2. Creating single and double reference sequences 
3. Aligning reads to both references using Minimap2
4. Analyzing cross-junction alignments for circular validation
5. Generating consensus sequences and final validation metrics

Key features:
- Validates circular DNA through cross-junction analysis
- Creates consensus sequences using pileup or bcftools
- Analyzes alignment patterns for circularity confirmation
- Supports parallel processing for multiple candidates
- Comprehensive validation metrics and statistics

Migrated from step13_Juggler_v1.py to the new CircleSeeker2 architecture.
"""

from __future__ import annotations

import csv
import logging
import os
import shutil
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
import pysam

from circleseeker2.exceptions import PipelineError


# ==================== CONFIGURATION ====================

@dataclass
class ValidationConfig:
    """Configuration for validation parameters."""
    # Thresholds
    min_cross_reads: int = 1
    min_overhang: int = 20
    mapq_threshold: int = 20
    
    # Consensus calling
    consensus_method: str = "pileup"  # "pileup" or "bcftools"
    consensus_mapq: int = 0
    consensus_bq: int = 10
    
    # Alignment
    preset: str = "map-hifi"
    allow_secondary: bool = False
    
    # Processing
    threads: int = 4
    max_workers: int = 1
    keep_temp: bool = False
    
    # Filtering
    only_circular: bool = False
    reads_column: str = "all_read_ids"
    
    # Logging
    log_level: str = "INFO"


@dataclass
class RegionInfo:
    """Information about a genomic region."""
    idx: int
    chrom: str
    start: int
    end: int
    read_ids: List[str]
    is_circular: bool = False
    raw_data: Dict[str, str] = field(default_factory=dict)
    
    @property
    def length(self) -> int:
        return self.end - self.start
    
    @property
    def name(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}"


# ==================== UTILITIES ====================

def check_tool(name: str) -> bool:
    """Check if a tool is available."""
    return shutil.which(name) is not None


def run_command(cmd: List[str], cwd: Optional[Path] = None, 
                capture_output: bool = False) -> Optional[str]:
    """Run external command safely."""
    try:
        if capture_output:
            result = subprocess.run(
                cmd, cwd=str(cwd) if cwd else None,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                check=True, text=True
            )
            return result.stdout
        else:
            subprocess.run(
                cmd, cwd=str(cwd) if cwd else None,
                stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
                check=True
            )
            return None
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Command failed: {' '.join(cmd)}\n{e.stderr}")


def ensure_fasta_index(fasta_path: Path) -> None:
    """Ensure FASTA file has an index."""
    fai_path = Path(str(fasta_path) + ".fai")
    if not fai_path.exists():
        if not check_tool("samtools"):
            raise RuntimeError(f"Missing {fai_path} and samtools not found")
        run_command(["samtools", "faidx", str(fasta_path)])


def parse_id_list(id_string: str, delimiter: str = ";") -> List[str]:
    """Parse a delimited string of IDs."""
    if not id_string or id_string.strip() in ["", ".", "NA", "None"]:
        return []
    
    # Handle both semicolon and comma delimiters
    id_string = id_string.replace(",", delimiter)
    ids = [id.strip() for id in id_string.split(delimiter) if id.strip()]
    
    # Remove potential suffixes like /ccs
    cleaned_ids = []
    for id in ids:
        cleaned_ids.append(id.split()[0])  # Take first part if space-separated
    
    return cleaned_ids


def write_fasta(output_path: Path, records: List[Tuple[str, str]]) -> None:
    """Write FASTA records to file."""
    with open(output_path, "w") as f:
        for header, sequence in records:
            f.write(f">{header}\n")
            # Write sequence in 80-character lines
            for i in range(0, len(sequence), 80):
                f.write(sequence[i:i+80] + "\n")


# ==================== CORE FUNCTIONS ====================

def extract_reads_pysam(read_ids: Set[str], reads_fasta: pysam.FastaFile,
                        output_path: Path, logger: logging.Logger) -> Tuple[int, int]:
    """
    Extract reads using pysam (fast indexed access).
    
    Returns:
        (found_count, missing_count)
    """
    found = 0
    missing = 0
    
    with open(output_path, "w") as out:
        for read_id in read_ids:
            try:
                # Try exact match first
                sequence = reads_fasta.fetch(read_id)
                out.write(f">{read_id}\n{sequence}\n")
                found += 1
            except KeyError:
                # Try without suffix
                base_id = read_id.split("/")[0] if "/" in read_id else read_id
                try:
                    sequence = reads_fasta.fetch(base_id)
                    out.write(f">{base_id}\n{sequence}\n")
                    found += 1
                except KeyError:
                    missing += 1
                    logger.debug(f"Read not found: {read_id}")
    
    return found, missing


def prepare_references(chrom: str, start: int, end: int,
                      ref_fasta: pysam.FastaFile, work_dir: Path) -> Tuple[Path, Path, int]:
    """
    Prepare single and double reference sequences.
    
    Returns:
        (single_ref_path, double_ref_path, midpoint)
    """
    # Extract region sequence
    region_seq = ref_fasta.fetch(chrom, start, end).upper()
    region_len = len(region_seq)
    
    # Single reference
    single_ref = work_dir / "ref_1x.fa"
    single_name = f"ref1x|{chrom}:{start}-{end}"
    write_fasta(single_ref, [(single_name, region_seq)])
    
    # Double reference (concatenated)
    double_ref = work_dir / "ref_2x.fa"
    double_name = f"ref2x|{chrom}:{start}-{end}"
    double_seq = region_seq + region_seq
    write_fasta(double_ref, [(double_name, double_seq)])
    
    midpoint = region_len
    
    return single_ref, double_ref, midpoint


def align_minimap2(ref_path: Path, reads_path: Path, output_path: Path,
                  preset: str, threads: int, to_bam: bool = False,
                  allow_secondary: bool = False) -> None:
    """
    Align reads using minimap2.
    
    Args:
        to_bam: If True, output sorted BAM with index; if False, output SAM
    """
    if not check_tool("minimap2"):
        raise RuntimeError("minimap2 not found")
    
    # Build minimap2 command
    mm2_cmd = [
        "minimap2", "-ax", preset,
        "-t", str(threads),
        "--secondary=no" if not allow_secondary else "--secondary=yes",
        "-Y",  # Use soft clipping for supplementary alignments
        str(ref_path), str(reads_path)
    ]
    
    if to_bam:
        # Pipe through samtools sort
        if not check_tool("samtools"):
            raise RuntimeError("samtools not found")
        
        # Use shell pipeline
        pipeline = " | ".join([
            " ".join(mm2_cmd),
            f"samtools sort -@ {threads} -o {output_path}"
        ])
        
        subprocess.run(
            ["bash", "-c", pipeline],
            check=True, stderr=subprocess.DEVNULL
        )
        
        # Index the BAM
        run_command(["samtools", "index", str(output_path)])
    else:
        # Direct SAM output
        with open(output_path, "w") as out:
            subprocess.run(mm2_cmd, stdout=out, check=True, stderr=subprocess.DEVNULL)


def count_cross_reads(alignment_path: Path, ref_name: str, midpoint: int,
                     min_mapq: int, min_overhang: int) -> Tuple[int, Set[str], float]:
    """
    Count reads that cross the midpoint in double reference.
    
    Returns:
        (total_cross_reads, unique_read_ids, average_span)
    """
    mode = "rb" if alignment_path.suffix.lower() in (".bam", ".cram") else "r"
    cross_reads = {}  # read_id -> span_length
    
    with pysam.AlignmentFile(str(alignment_path), mode) as af:
        for aln in af:
            if aln.is_unmapped or aln.is_secondary:
                continue
            
            if aln.reference_name != ref_name:
                continue
            
            if aln.mapping_quality < min_mapq:
                continue
            
            # Check if alignment crosses midpoint
            start = aln.reference_start
            end = aln.reference_end
            
            if start < midpoint < end:
                left_overhang = midpoint - start
                right_overhang = end - midpoint
                
                if left_overhang >= min_overhang and right_overhang >= min_overhang:
                    span = end - start
                    read_id = aln.query_name.split()[0]
                    
                    # Keep the longest span for each read
                    if read_id not in cross_reads or span > cross_reads[read_id]:
                        cross_reads[read_id] = span
    
    # Calculate statistics
    unique_ids = set(cross_reads.keys())
    total_cross = len(unique_ids)
    
    if cross_reads:
        avg_span = round(np.mean(list(cross_reads.values())), 2)
    else:
        avg_span = 0.0
    
    return total_cross, unique_ids, avg_span


def consensus_pileup(bam_path: Path, ref_name: str, ref_seq: str,
                    min_base_quality: int = 10, min_depth: int = 1) -> str:
    """Generate consensus using pysam pileup (fast method)."""
    consensus = list(ref_seq.upper())
    ref_len = len(ref_seq)
    
    with pysam.AlignmentFile(str(bam_path), "rb") as af:
        for pileup_column in af.pileup(
            ref_name, 0, ref_len,
            truncate=True,
            stepper="all",
            min_base_quality=min_base_quality
        ):
            pos = pileup_column.reference_pos
            if pos >= ref_len:
                continue
            
            # Count bases at this position
            base_counts = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0}
            
            for pileup_read in pileup_column.pileups:
                if pileup_read.is_del or pileup_read.is_refskip:
                    continue
                
                query_pos = pileup_read.query_position
                if query_pos is None:
                    continue
                
                base = pileup_read.alignment.query_sequence[query_pos].upper()
                if base in base_counts:
                    base_counts[base] += 1
                else:
                    base_counts["N"] += 1
            
            # Determine consensus base
            total_depth = sum(base_counts.values())
            if total_depth >= min_depth:
                consensus_base = max(base_counts, key=base_counts.get)
                if base_counts[consensus_base] > 0:
                    consensus[pos] = consensus_base
    
    return "".join(consensus)


def consensus_bcftools(bam_path: Path, ref_path: Path, work_dir: Path,
                      min_mapq: int = 0, min_bq: int = 0) -> str:
    """Generate consensus using bcftools (standard method)."""
    if not check_tool("bcftools"):
        raise RuntimeError("bcftools not found")
    
    # Filter BAM if needed
    filtered_bam = bam_path
    if min_mapq > 0:
        filtered_bam = work_dir / "filtered.bam"
        run_command([
            "samtools", "view",
            "-F", "0x900",  # Exclude secondary and supplementary
            "-q", str(min_mapq),
            "-b", str(bam_path),
            "-o", str(filtered_bam)
        ])
        run_command(["samtools", "index", str(filtered_bam)])
    
    # Generate pileup
    vcf_path = work_dir / "variants.vcf.gz"
    mpileup_cmd = ["bcftools", "mpileup", "-f", str(ref_path)]
    if min_bq > 0:
        mpileup_cmd.extend(["-Q", str(min_bq)])
    mpileup_cmd.extend([str(filtered_bam), "-Oz", "-o", str(vcf_path)])
    run_command(mpileup_cmd)
    
    # Call variants
    called_vcf = work_dir / "called.vcf.gz"
    run_command([
        "bcftools", "call", "-c", "-Oz",
        "-o", str(called_vcf), str(vcf_path)
    ])
    
    # Index VCF
    run_command(["tabix", "-p", "vcf", str(called_vcf)])
    
    # Generate consensus
    consensus_fa = work_dir / "consensus.fa"
    run_command([
        "bcftools", "consensus",
        "-f", str(ref_path),
        str(called_vcf),
        "-o", str(consensus_fa)
    ])
    
    # Read consensus
    with pysam.FastaFile(str(consensus_fa)) as f:
        return f.fetch(f.references[0])


# ==================== REGION PROCESSOR ====================

class RegionProcessor:
    """Process a single genomic region."""
    
    def __init__(self, config: ValidationConfig, ref_fasta: pysam.FastaFile,
                 reads_fasta: pysam.FastaFile, output_dir: Path,
                 logger: logging.Logger):
        self.config = config
        self.ref_fasta = ref_fasta
        self.reads_fasta = reads_fasta
        self.output_dir = output_dir
        self.logger = logger
    
    def process(self, region: RegionInfo) -> Dict[str, any]:
        """Process a single region and return results."""
        self.logger.info(f"Processing region {region.idx}: {region.name}")
        
        # Create working directory
        work_dir = self.output_dir / f"region_{region.chrom}_{region.start}_{region.end}"
        work_dir.mkdir(parents=True, exist_ok=True)
        
        try:
            # Step 1: Extract reads
            if not region.read_ids:
                self.logger.warning(f"{region.name}: No read IDs provided")
                return self._create_result(region, validated=False, reason="no_reads")
            
            reads_file = work_dir / "reads.fa"
            found, missing = extract_reads_pysam(
                set(region.read_ids), self.reads_fasta, reads_file, self.logger
            )
            
            if found == 0:
                self.logger.warning(f"{region.name}: No reads found in FASTA")
                return self._create_result(
                    region, validated=False, reason="no_reads_found",
                    found_reads=found, missing_reads=missing
                )
            
            self.logger.debug(f"{region.name}: Extracted {found} reads, {missing} missing")
            
            # Step 2: Prepare references
            ref_1x, ref_2x, midpoint = prepare_references(
                region.chrom, region.start, region.end,
                self.ref_fasta, work_dir
            )
            
            # Get reference names from FASTA files
            with pysam.FastaFile(str(ref_2x)) as f2x:
                ref_2x_name = f2x.references[0]
            with pysam.FastaFile(str(ref_1x)) as f1x:
                ref_1x_name = f1x.references[0]
                ref_1x_seq = f1x.fetch(ref_1x_name)
            
            # Step 3: Align to double reference and count cross reads
            sam_2x = work_dir / "align_2x.sam"
            align_minimap2(
                ref_2x, reads_file, sam_2x,
                self.config.preset, self.config.threads,
                to_bam=False, allow_secondary=self.config.allow_secondary
            )
            
            cross_count, cross_ids, avg_span = count_cross_reads(
                sam_2x, ref_2x_name, midpoint,
                self.config.mapq_threshold, self.config.min_overhang
            )
            
            self.logger.info(
                f"{region.name}: {cross_count} cross reads, avg_span={avg_span:.2f}"
            )
            
            # Check if passes validation
            if cross_count < self.config.min_cross_reads:
                return self._create_result(
                    region, validated=False, reason="insufficient_cross_reads",
                    found_reads=found, missing_reads=missing,
                    cross_reads=cross_count, unique_cross_ids=cross_ids,
                    avg_span=avg_span
                )
            
            # Step 4: Generate consensus for validated regions
            bam_1x = work_dir / "align_1x.bam"
            align_minimap2(
                ref_1x, reads_file, bam_1x,
                self.config.preset, self.config.threads,
                to_bam=True, allow_secondary=self.config.allow_secondary
            )
            
            # Generate consensus
            if self.config.consensus_method == "pileup":
                consensus_seq = consensus_pileup(
                    bam_1x, ref_1x_name, ref_1x_seq,
                    min_base_quality=self.config.consensus_bq
                )
            else:
                consensus_seq = consensus_bcftools(
                    bam_1x, ref_1x, work_dir,
                    min_mapq=self.config.consensus_mapq,
                    min_bq=self.config.consensus_bq
                )
            
            # Save consensus
            consensus_file = work_dir / "consensus.fa"
            write_fasta(consensus_file, [(f"consensus|{region.name}", consensus_seq)])
            
            return self._create_result(
                region, validated=True,
                found_reads=found, missing_reads=missing,
                cross_reads=cross_count, unique_cross_ids=cross_ids,
                avg_span=avg_span, consensus_sequence=consensus_seq,
                consensus_file=str(consensus_file)
            )
            
        except Exception as e:
            self.logger.error(f"Error processing {region.name}: {e}")
            return self._create_result(
                region, validated=False, reason=f"error: {str(e)}"
            )
        finally:
            # Clean up temporary files
            if not self.config.keep_temp:
                try:
                    for file in work_dir.glob("*.sam"):
                        file.unlink()
                    for file in work_dir.glob("*.fa"):
                        if "consensus" not in file.name:
                            file.unlink()
                    for file in work_dir.glob("*.fai"):
                        file.unlink()
                except Exception:
                    pass
    
    def _create_result(self, region: RegionInfo, validated: bool,
                      reason: str = "", **kwargs) -> Dict[str, any]:
        """Create a result dictionary."""
        result = {
            "chrom": region.chrom,
            "start": region.start,
            "end": region.end,
            "region_length": region.length,
            "validated": validated,
            "reason": reason,
            "found_reads": kwargs.get("found_reads", 0),
            "missing_reads": kwargs.get("missing_reads", 0),
            "cross_reads": kwargs.get("cross_reads", 0),
            "unique_cross_reads": len(kwargs.get("unique_cross_ids", [])),
            "unique_cross_read_ids": ";".join(sorted(kwargs.get("unique_cross_ids", []))),
            "avg_span": round(kwargs.get("avg_span", 0.0), 2),
            "span_ratio": round(kwargs.get("avg_span", 0.0) / region.length, 4) if region.length > 0 else 0.0,
            "consensus_sequence": kwargs.get("consensus_sequence", ""),
            "consensus_file": kwargs.get("consensus_file", ""),
            "is_circular_pattern": region.is_circular
        }
        
        # Add original data
        for key, value in region.raw_data.items():
            if key not in result:
                result[f"orig_{key}"] = value
        
        return result


# ==================== MAIN JUGGLER CLASS ====================

class Juggler:
    """Circular DNA validation processor."""
    
    def __init__(self, config: Optional[ValidationConfig] = None,
                 logger: Optional[logging.Logger] = None):
        """Initialize Juggler validator."""
        self.config = config or ValidationConfig()
        self.logger = logger or logging.getLogger(self.__class__.__name__)
        
        # Check required tools
        self._check_tools()
    
    def _check_tools(self) -> None:
        """Check if required external tools are available."""
        required_tools = ["minimap2", "samtools"]
        
        if self.config.consensus_method == "bcftools":
            required_tools.extend(["bcftools", "tabix"])
        
        missing_tools = [t for t in required_tools if not check_tool(t)]
        if missing_tools:
            raise PipelineError(f"Missing required tools: {', '.join(missing_tools)}")
        
        self.logger.debug("All required tools are available")
    
    def load_regions(self, tsv_path: Path) -> List[RegionInfo]:
        """Load regions from TSV file."""
        regions = []
        
        with open(tsv_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            
            # Check required columns
            required_cols = {"chrom", "start", "end"}
            missing = required_cols - set(reader.fieldnames)
            if missing:
                raise ValueError(f"Missing required columns: {missing}")
            
            if self.config.reads_column not in reader.fieldnames:
                self.logger.warning(f"Reads column '{self.config.reads_column}' not found, will try alternatives")
            
            for idx, row in enumerate(reader, 1):
                # Parse circular pattern
                is_circular = False
                if "is_circular_pattern" in row:
                    is_circ_val = row["is_circular_pattern"].lower()
                    is_circular = is_circ_val in ("yes", "true", "1")
                
                # Skip non-circular if requested
                if self.config.only_circular and not is_circular:
                    continue
                
                # Parse coordinates
                chrom = row["chrom"]
                start = int(row["start"])
                end = int(row["end"])
                
                # Extract read IDs from various possible columns
                read_ids = []
                for col in [self.config.reads_column, "all_read_ids", "support_read_ids", "split_read_ids"]:
                    if col in row and row[col]:
                        read_ids = parse_id_list(row[col])
                        if read_ids:
                            break
                
                regions.append(RegionInfo(
                    idx=idx,
                    chrom=chrom,
                    start=start,
                    end=end,
                    read_ids=read_ids,
                    is_circular=is_circular,
                    raw_data=dict(row)
                ))
        
        self.logger.info(f"Loaded {len(regions)} regions from TSV")
        return regions
    
    def process_region_parallel(self, region: RegionInfo, ref_path: Path, 
                               reads_path: Path, output_dir: Path) -> Dict:
        """Wrapper for parallel processing."""
        # Create new logger for subprocess
        logger = logging.getLogger(f"worker_{os.getpid()}")
        logger.setLevel(getattr(logging, self.config.log_level.upper()))
        
        # Open FASTA files
        ref_fasta = pysam.FastaFile(str(ref_path))
        reads_fasta = pysam.FastaFile(str(reads_path))
        
        try:
            processor = RegionProcessor(
                self.config, ref_fasta, reads_fasta, output_dir, logger
            )
            return processor.process(region)
        finally:
            ref_fasta.close()
            reads_fasta.close()
    
    def process_all_regions(self, regions: List[RegionInfo], ref_path: Path,
                          reads_path: Path, output_dir: Path) -> pd.DataFrame:
        """Process all regions with optional parallelization."""
        
        # Ensure FASTA files are indexed
        ensure_fasta_index(ref_path)
        ensure_fasta_index(reads_path)
        
        # Open FASTA files
        ref_fasta = pysam.FastaFile(str(ref_path))
        reads_fasta = pysam.FastaFile(str(reads_path))
        
        try:
            results = []
            
            if self.config.max_workers > 1:
                # Parallel processing
                self.logger.info(f"Processing {len(regions)} regions with {self.config.max_workers} workers")
                
                with ProcessPoolExecutor(max_workers=self.config.max_workers) as executor:
                    # Submit all tasks
                    futures = {}
                    for region in regions:
                        future = executor.submit(
                            self.process_region_parallel,
                            region, ref_path, reads_path, output_dir
                        )
                        futures[future] = region
                    
                    # Collect results
                    for future in as_completed(futures):
                        region = futures[future]
                        try:
                            result = future.result()
                            results.append(result)
                        except Exception as e:
                            self.logger.error(f"Failed to process {region.name}: {e}")
                            results.append({
                                "chrom": region.chrom,
                                "start": region.start,
                                "end": region.end,
                                "validated": False,
                                "reason": f"error: {str(e)}"
                            })
            else:
                # Sequential processing
                self.logger.info(f"Processing {len(regions)} regions sequentially")
                processor = RegionProcessor(
                    self.config, ref_fasta, reads_fasta, output_dir, self.logger
                )
                
                for region in regions:
                    result = processor.process(region)
                    results.append(result)
        
        finally:
            ref_fasta.close()
            reads_fasta.close()
        
        return pd.DataFrame(results)
    
    def generate_publication_outputs(self, validated_df: pd.DataFrame, 
                                   consensus_sequences: Dict[str, str],
                                   output_prefix: Path) -> None:
        """Generate clean formatted outputs for publication."""
        
        self.logger.info("=" * 60)
        self.logger.info("Generating publication-ready outputs...")
        
        # 1. Clean CSV
        clean_data = []
        for idx, row in validated_df.iterrows():
            # Generate simple ID (IUeccDNA1, IUeccDNA2, etc.)
            circ_id = f"IUeccDNA{len(clean_data) + 1}"
            
            # Calculate length
            length = row.get("region_length", row["end"] - row["start"])
            
            # Create name
            name = f"{row['chrom']}:{row['start']}-{row['end']}"
            
            # Extract fields with fallbacks
            unique_split = 0
            for col in ["unique_split_count", "orig_unique_split_count", "orig_split_read_count"]:
                if col in row and pd.notna(row[col]):
                    unique_split = int(row[col])
                    break
            
            mean_mapq = 0.0
            for col in ["mean_mapq_region", "orig_mean_mapq_region", "orig_mean_mapq_support"]:
                if col in row and pd.notna(row[col]):
                    mean_mapq = float(row[col])
                    break
            
            depth_sum = 0.0
            for col in ["depth_sum", "orig_depth_sum"]:
                if col in row and pd.notna(row[col]):
                    depth_sum = float(row[col])
                    break
            
            clean_record = {
                "ID": circ_id,
                "Chrom": row["chrom"],
                "Start": int(row["start"]),
                "End": int(row["end"]),
                "Length": int(length),
                "Name": name,
                "State": "Inferred",
                "Class": "UeccDNA",  # Using IUeccDNA
                "Found_Reads": int(row.get("found_reads", 0)),
                "Unique_Cross_Reads": int(row.get("unique_cross_reads", 0)),
                "Unique_Split_Reads": unique_split,
                "Mean_MAPQ": round(mean_mapq, 2),
                "Depth_Sum": round(depth_sum, 1),
                "Total_Reads": int(row.get("found_reads", 0))
            }
            clean_data.append(clean_record)
        
        clean_df = pd.DataFrame(clean_data)
        clean_df = clean_df.sort_values(["Chrom", "Start", "End"])
        
        # Reset IDs after sorting
        for i in range(len(clean_df)):
            clean_df.iloc[i, clean_df.columns.get_loc("ID")] = f"IUeccDNA{i + 1}"
        
        # Save clean CSV
        clean_csv = Path(f"{output_prefix}.clean.csv")
        clean_df.to_csv(clean_csv, index=False)
        self.logger.info(f"Saved clean CSV to: {clean_csv}")
        
        # 2. Simplified BED
        bed_file = Path(f"{output_prefix}.simplified.bed")
        with open(bed_file, "w") as f:
            f.write("# Circular DNA regions detected by validation pipeline\n")
            f.write("# chrom\tstart\tend\tname\tscore\tstrand\n")
            
            for _, row in clean_df.iterrows():
                # Use Total_Reads as score (more intuitive than normalized score)
                # Cap at 1000 for visualization
                score = min(1000, int(row["Total_Reads"]))
                bed_line = f"{row['Chrom']}\t{row['Start']}\t{row['End']}\t{row['ID']}\t{score}\t.\n"
                f.write(bed_line)
        
        self.logger.info(f"Saved BED file to: {bed_file}")
        
        # 3. Renamed FASTA with simple IDs
        renamed_fasta = Path(f"{output_prefix}.renamed.fasta")
        with open(renamed_fasta, "w") as f:
            for idx, (_, row) in enumerate(clean_df.iterrows()):
                # Find corresponding consensus sequence
                key = f"{row['Chrom']}_{row['Start']}_{row['End']}"
                if key in consensus_sequences:
                    # Use simple ID format: IUeccDNA1, IUeccDNA2, etc.
                    simple_id = f"IUeccDNA{idx + 1}"
                    header = (f"{simple_id}|{row['Name']}|"
                             f"Length={row['Length']}|"
                             f"Cross_Reads={row['Unique_Cross_Reads']}")
                    
                    seq = consensus_sequences[key]
                    f.write(f">{header}\n")
                    for i in range(0, len(seq), 80):
                        f.write(seq[i:i+80] + "\n")
        
        self.logger.info(f"Saved renamed FASTA to: {renamed_fasta}")
        self.logger.info("Publication-ready outputs generated successfully!")
    
    def save_results(self, df: pd.DataFrame, output_prefix: Path) -> None:
        """Save final results - only publication-ready files."""
        
        # Filter validated regions
        validated_df = df[df["validated"] == True].copy()
        
        if len(validated_df) > 0:
            # Collect consensus sequences
            consensus_sequences = {}
            for _, row in validated_df.iterrows():
                if row["consensus_sequence"]:
                    key = f"{row['chrom']}_{row['start']}_{row['end']}"
                    consensus_sequences[key] = row["consensus_sequence"]
            
            # Generate publication outputs only
            self.generate_publication_outputs(validated_df, consensus_sequences, output_prefix)
        
        # Print summary
        self.logger.info("=" * 60)
        self.logger.info(f"Total regions processed: {len(df)}")
        self.logger.info(f"Validated regions: {len(validated_df)}")
        self.logger.info(f"Failed regions: {len(df) - len(validated_df)}")
        
        if len(df) > len(validated_df):
            failure_reasons = df[df["validated"] == False]["reason"].value_counts()
            self.logger.info("Failure reasons:")
            for reason, count in failure_reasons.items():
                self.logger.info(f"  {reason}: {count}")
        
        # Print validation statistics
        if len(validated_df) > 0:
            self.logger.info("Validation statistics:")
            self.logger.info(f"  Average cross reads: {validated_df['cross_reads'].mean():.1f}")
            self.logger.info(f"  Average span: {validated_df['avg_span'].mean():.1f}")
            self.logger.info(f"  Average span ratio: {validated_df['span_ratio'].mean():.3f}")
    
    def run_validation(self, candidates_tsv: Path, reads_fasta: Path,
                      reference_fasta: Path, output_prefix: Path,
                      max_regions: Optional[int] = None) -> pd.DataFrame:
        """Run the complete validation pipeline."""
        self.logger.info("=" * 60)
        self.logger.info("Juggler - Circular DNA Validation Pipeline")
        self.logger.info("=" * 60)
        
        # Load regions
        regions = self.load_regions(candidates_tsv)
        
        if max_regions:
            regions = regions[:max_regions]
            self.logger.info(f"Limited to first {max_regions} regions")
        
        if not regions:
            self.logger.warning("No regions to process")
            return pd.DataFrame()
        
        # Create output directory
        output_dir = output_prefix.parent / f"{output_prefix.name}_workdir"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Process regions
        results_df = self.process_all_regions(
            regions, reference_fasta, reads_fasta, output_dir
        )
        
        # Save only publication-ready results
        self.save_results(results_df, output_prefix)
        
        # Clean up working directory if not keeping temp files
        if not self.config.keep_temp:
            try:
                shutil.rmtree(output_dir)
                self.logger.info("Cleaned up working directory")
            except Exception:
                pass
        
        self.logger.info("=" * 60)
        self.logger.info("Pipeline completed successfully")
        
        return results_df