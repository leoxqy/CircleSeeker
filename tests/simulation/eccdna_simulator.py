#!/usr/bin/env python3
"""
eccDNA Simulation Generator for CircleSeeker Validation

Generates simulated:
- Reference genome (FASTA)
- UeccDNA reads (100): Simple circular DNA from single genomic location
- MeccDNA reads (20): Multi-copy repeat circular DNA from multiple genomic loci
- CeccDNA reads (20): Chimeric circular DNA from multiple genomic segments

Output:
- Simulated reference genome
- Simulated HiFi reads (FASTA)
- Ground truth annotations (CSV)
"""

import random
import csv
import json
from pathlib import Path
from dataclasses import dataclass, field, asdict
from typing import List, Tuple, Optional, Dict


# =============================================================================
# Configuration
# =============================================================================

@dataclass
class SimulationConfig:
    """Configuration for eccDNA simulation."""
    # Reference genome
    num_chromosomes: int = 5
    chr_length: int = 1_000_000  # 1 Mb per chromosome

    # eccDNA counts
    num_uecc: int = 100
    num_mecc: int = 20
    num_cecc: int = 20

    # UeccDNA parameters
    uecc_length_range: Tuple[int, int] = (200, 5000)

    # MeccDNA parameters (repeat copies across multiple loci)
    mecc_unit_range: Tuple[int, int] = (100, 500)  # repeat unit length
    mecc_copy_range: Tuple[float, float] = (2.0, 10.0)  # copy number
    mecc_repeat_mutation_rate: float = 0.005  # divergence among genomic copies (lowered for better alignment)

    # CeccDNA parameters
    cecc_segments_range: Tuple[int, int] = (2, 4)  # number of segments
    cecc_segment_length_range: Tuple[int, int] = (100, 1000)
    cecc_inter_chr_ratio: float = 0.3  # 30% inter-chromosomal

    # Sequencing simulation
    error_rate: float = 0.001  # HiFi error rate ~0.1%

    # Random seed
    seed: int = 42

    # Output directory
    output_dir: str = "simulation_data"


# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class GenomicRegion:
    """Represents a genomic region."""
    chrom: str
    start: int  # 0-based
    end: int    # 0-based, exclusive
    strand: str = "+"

    @property
    def length(self) -> int:
        return self.end - self.start

    def to_region_string(self) -> str:
        """Convert to region string format: chr:start-end_strand"""
        return f"{self.chrom}:{self.start}-{self.end}_{self.strand}"


@dataclass
class EccDNARecord:
    """Ground truth record for an eccDNA."""
    eccdna_id: str
    eccdna_type: str  # Uecc, Mecc, Cecc
    read_id: str
    regions: List[GenomicRegion]
    length: int
    copy_number: float = 1.0
    sequence: str = ""

    # Additional metadata
    is_inter_chr: bool = False  # For CeccDNA
    num_segments: int = 1

    def to_dict(self) -> dict:
        """Convert to dictionary for CSV output."""
        region_strs = [r.to_region_string() for r in self.regions]
        return {
            "eccdna_id": self.eccdna_id,
            "eccdna_type": self.eccdna_type,
            "read_id": self.read_id,
            "regions": ";".join(region_strs),
            "chr": self.regions[0].chrom if self.regions else "",
            "start0": self.regions[0].start if self.regions else 0,
            "end0": self.regions[0].end if self.regions else 0,
            "strand": self.regions[0].strand if self.regions else "+",
            "length": self.length,
            "copy_number": self.copy_number,
            "num_segments": self.num_segments,
            "is_inter_chr": self.is_inter_chr,
        }


@dataclass
class RepeatFamily:
    """A repeat family embedded in the reference genome for MeccDNA simulation."""

    family_id: str
    unit_length: int
    copies: List[GenomicRegion]


# =============================================================================
# Genome Generator
# =============================================================================

class ReferenceGenomeGenerator:
    """Generate a simulated reference genome."""

    BASES = ['A', 'C', 'G', 'T']

    def __init__(self, config: SimulationConfig):
        self.config = config
        self.genome = {}
        self.repeat_families: List[RepeatFamily] = []
        self.repeat_intervals_by_chrom: Dict[str, List[Tuple[int, int]]] = {}
        random.seed(config.seed)

    def generate(self) -> dict:
        """Generate random reference genome."""
        genome_lists: Dict[str, List[str]] = {}
        for i in range(1, self.config.num_chromosomes + 1):
            chrom_name = f"chr{i}"
            genome_lists[chrom_name] = random.choices(self.BASES, k=self.config.chr_length)

        # Embed repeat families for MeccDNA simulation so that Mecc reads can
        # align to multiple distinct genomic loci.
        self._embed_mecc_repeat_families(genome_lists)

        # Finalize genome as strings
        self.genome = {chrom: ''.join(seq_list) for chrom, seq_list in genome_lists.items()}
        return self.genome

    def _embed_mecc_repeat_families(self, genome_lists: Dict[str, List[str]]) -> None:
        """Embed repeat families into the simulated genome for MeccDNA."""
        num_families = max(0, int(self.config.num_mecc))
        if num_families == 0:
            self.repeat_families = []
            self.repeat_intervals_by_chrom = {chrom: [] for chrom in genome_lists}
            return

        min_unit, max_unit = self.config.mecc_unit_range
        max_copy = int(self.config.mecc_copy_range[1])
        # Ensure enough distinct loci per family to support assembling Mecc reads
        copies_per_family = max(2, max_copy + 2)

        occupied: Dict[str, List[Tuple[int, int]]] = {chrom: [] for chrom in genome_lists}
        repeat_families: List[RepeatFamily] = []

        for fam_idx in range(1, num_families + 1):
            unit_len = random.randint(min_unit, max_unit)
            base_seq = ''.join(random.choices(self.BASES, k=unit_len))

            family_copies: List[GenomicRegion] = []
            for _ in range(copies_per_family):
                chrom = random.choice(list(genome_lists.keys()))
                start, end = self._pick_non_overlapping_interval(
                    chrom=chrom,
                    length=unit_len,
                    chr_len=len(genome_lists[chrom]),
                    occupied=occupied,
                )

                variant_seq = self._mutate_sequence(base_seq, self.config.mecc_repeat_mutation_rate)
                genome_lists[chrom][start:end] = list(variant_seq)

                family_copies.append(GenomicRegion(chrom=chrom, start=start, end=end, strand="+"))

            repeat_families.append(
                RepeatFamily(
                    family_id=f"RF{fam_idx}",
                    unit_length=unit_len,
                    copies=family_copies,
                )
            )

        self.repeat_families = repeat_families
        self.repeat_intervals_by_chrom = occupied

    @staticmethod
    def _mutate_sequence(seq: str, mutation_rate: float) -> str:
        """Introduce substitutions to create similar-but-not-identical repeat copies."""
        if mutation_rate <= 0:
            return seq
        seq_list = list(seq)
        for i, base in enumerate(seq_list):
            if random.random() < mutation_rate:
                alternatives = [b for b in "ACGT" if b != base]
                seq_list[i] = random.choice(alternatives)
        return ''.join(seq_list)

    @staticmethod
    def _pick_non_overlapping_interval(
        chrom: str,
        length: int,
        chr_len: int,
        occupied: Dict[str, List[Tuple[int, int]]],
        max_attempts: int = 2000,
    ) -> Tuple[int, int]:
        """Pick a non-overlapping interval on a chromosome."""
        max_start = chr_len - length
        if max_start < 0:
            return 0, chr_len

        intervals = occupied.setdefault(chrom, [])
        for _ in range(max_attempts):
            start = random.randint(0, max_start)
            end = start + length
            if all(end <= s or start >= e for s, e in intervals):
                intervals.append((start, end))
                return start, end

        # Fallback: allow overlap if placement fails (should be rare with defaults)
        start = random.randint(0, max_start)
        end = start + length
        intervals.append((start, end))
        return start, end

    def get_sequence(self, chrom: str, start: int, end: int, strand: str = "+") -> str:
        """Extract sequence from genome."""
        if chrom not in self.genome:
            raise ValueError(f"Chromosome {chrom} not found")

        seq = self.genome[chrom][start:end]

        if strand == "-":
            seq = self._reverse_complement(seq)

        return seq

    @staticmethod
    def _reverse_complement(seq: str) -> str:
        """Get reverse complement of sequence."""
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
                      'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
        return ''.join(complement.get(b, b) for b in reversed(seq))

    def write_fasta(self, output_path: Path):
        """Write genome to FASTA file."""
        with open(output_path, 'w') as f:
            for chrom, seq in self.genome.items():
                f.write(f">{chrom}\n")
                # Write sequence in 80-character lines
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i+80] + "\n")


# =============================================================================
# eccDNA Simulators
# =============================================================================

class EccDNASimulator:
    """Base class for eccDNA simulation."""

    def __init__(
        self,
        genome: dict,
        config: SimulationConfig,
        excluded_intervals_by_chrom: Optional[Dict[str, List[Tuple[int, int]]]] = None,
    ):
        self.genome = genome
        self.config = config
        self.chromosomes = list(genome.keys())
        self.excluded_intervals_by_chrom = excluded_intervals_by_chrom or {}

    def _random_region(self, min_len: int, max_len: int,
                       chrom: Optional[str] = None) -> GenomicRegion:
        """Generate a random genomic region."""
        if chrom is None:
            chrom = random.choice(self.chromosomes)

        length = random.randint(min_len, max_len)
        chr_len = len(self.genome[chrom])

        # Ensure we have room for the region
        max_start = chr_len - length
        if max_start < 0:
            max_start = 0
            length = chr_len

        start = random.randint(0, max_start)
        end = start + length
        strand = random.choice(["+", "-"])

        return GenomicRegion(chrom, start, end, strand)

    def _random_region_avoiding_excluded(
        self,
        min_len: int,
        max_len: int,
        chrom: Optional[str] = None,
        max_attempts: int = 2000,
    ) -> GenomicRegion:
        """Generate a random region that does not overlap excluded intervals."""
        for _ in range(max_attempts):
            region = self._random_region(min_len, max_len, chrom)
            if not self._overlaps_excluded(region):
                return region
        # Fallback if unable to find a clean region
        return self._random_region(min_len, max_len, chrom)

    def _overlaps_excluded(self, region: GenomicRegion) -> bool:
        """Check if a region overlaps any excluded intervals for its chromosome."""
        intervals = self.excluded_intervals_by_chrom.get(region.chrom, [])
        if not intervals:
            return False
        return any(region.end > s and region.start < e for s, e in intervals)

    def _get_sequence(self, region: GenomicRegion) -> str:
        """Extract sequence for a region."""
        seq = self.genome[region.chrom][region.start:region.end]
        if region.strand == "-":
            seq = ReferenceGenomeGenerator._reverse_complement(seq)
        return seq

    def _add_errors(self, seq: str) -> str:
        """Add sequencing errors to simulate HiFi reads."""
        if self.config.error_rate == 0:
            return seq

        seq_list = list(seq)
        for i in range(len(seq_list)):
            if random.random() < self.config.error_rate:
                # Substitution error
                original = seq_list[i]
                alternatives = [b for b in 'ACGT' if b != original.upper()]
                seq_list[i] = random.choice(alternatives)

        return ''.join(seq_list)


class UeccDNASimulator(EccDNASimulator):
    """Simulate UeccDNA (Unique eccDNA from single location)."""

    def generate(self, count: int) -> List[EccDNARecord]:
        """Generate UeccDNA records."""
        records = []

        for i in range(1, count + 1):
            # Random region
            min_len, max_len = self.config.uecc_length_range
            region = self._random_region_avoiding_excluded(min_len, max_len)

            # Get sequence (circular: no tandem repeat for Uecc)
            sequence = self._get_sequence(region)
            sequence = self._add_errors(sequence)

            # Create read ID with metadata
            read_id = f"sim_uecc_{i:04d}|unit_1|{region.length}|1.0"

            record = EccDNARecord(
                eccdna_id=f"UeccDNA{i}",
                eccdna_type="Uecc",
                read_id=read_id,
                regions=[region],
                length=region.length,
                copy_number=1.0,
                sequence=sequence,
                num_segments=1,
            )
            records.append(record)

        return records


class MeccDNASimulator(EccDNASimulator):
    """Simulate MeccDNA (Multi-copy eccDNA from multiple genomic loci)."""

    def __init__(
        self,
        genome: dict,
        config: SimulationConfig,
        repeat_families: Optional[List[RepeatFamily]] = None,
    ):
        super().__init__(genome, config)
        self.repeat_families = repeat_families or []

    def generate(self, count: int) -> List[EccDNARecord]:
        """Generate MeccDNA records."""
        records = []

        for i in range(1, count + 1):
            # "Mecc" in this project means: the *whole eccDNA sequence* has multiple
            # near-full-length matches on the reference genome (multi-locus mapping).
            # We therefore generate one locus as the eccDNA sequence, and record
            # additional genomic loci that carry the same repeat family.
            if not self.repeat_families:
                # Fallback: generate tandem repeat from single locus
                min_unit, max_unit = self.config.mecc_unit_range
                source_region = self._random_region(min_unit, max_unit)
                unit_seq = self._get_sequence(source_region)
                regions_used = [source_region]
                unit_length = source_region.length

                # Generate tandem repeat even in fallback mode
                min_copy, max_copy = self.config.mecc_copy_range
                copy_number = random.uniform(min_copy, max_copy)
                full_copies = int(copy_number)
                partial_fraction = copy_number - full_copies
                partial_len = int(len(unit_seq) * partial_fraction)
                sequence = unit_seq * full_copies + unit_seq[:partial_len]
            else:
                family = self.repeat_families[(i - 1) % len(self.repeat_families)]
                unit_length = family.unit_length

                # Choose how many distinct loci to expose (>=2).
                min_copy, max_copy = self.config.mecc_copy_range
                copy_number = random.uniform(min_copy, max_copy)
                n_loci = max(2, int(round(copy_number)))
                n_loci = min(n_loci, len(family.copies))

                source_region = random.choice(family.copies)
                other_choices = [r for r in family.copies if r != source_region]
                additional = random.sample(other_choices, k=max(0, n_loci - 1)) if other_choices else []
                regions_used = [source_region] + additional

                # MeccDNA is a tandem repeat: repeat the unit sequence copy_number times
                unit_seq = self._get_sequence(source_region)
                full_copies = int(copy_number)
                partial_fraction = copy_number - full_copies
                partial_len = int(len(unit_seq) * partial_fraction)
                sequence = unit_seq * full_copies + unit_seq[:partial_len]

            sequence = self._add_errors(sequence)

            # Create read ID with metadata
            read_id = f"sim_mecc_{i:04d}|unit_1|{unit_length}|{copy_number:.1f}"

            record = EccDNARecord(
                eccdna_id=f"MeccDNA{i}",
                eccdna_type="Mecc",
                read_id=read_id,
                regions=regions_used,
                length=len(sequence),
                copy_number=round(copy_number, 2),
                sequence=sequence,
                num_segments=len(regions_used),
            )
            records.append(record)

        return records


class CeccDNASimulator(EccDNASimulator):
    """Simulate CeccDNA (Chimeric eccDNA from multiple segments)."""

    def generate(self, count: int) -> List[EccDNARecord]:
        """Generate CeccDNA records."""
        records = []

        # Determine how many should be inter-chromosomal
        num_inter = int(count * self.config.cecc_inter_chr_ratio)

        for i in range(1, count + 1):
            is_inter = i <= num_inter

            # Random number of segments
            min_seg, max_seg = self.config.cecc_segments_range
            num_segments = random.randint(min_seg, max_seg)

            # Generate segments
            regions = []
            sequences = []

            if is_inter:
                # Inter-chromosomal: use different chromosomes
                chroms = random.sample(self.chromosomes,
                                       min(num_segments, len(self.chromosomes)))
            else:
                # Intra-chromosomal: use same chromosome
                chrom = random.choice(self.chromosomes)
                chroms = [chrom] * num_segments

            for j, chrom in enumerate(chroms):
                min_len, max_len = self.config.cecc_segment_length_range
                region = self._random_region_avoiding_excluded(min_len, max_len, chrom)
                regions.append(region)
                sequences.append(self._get_sequence(region))

            # Combine sequences
            full_sequence = ''.join(sequences)
            full_sequence = self._add_errors(full_sequence)

            # Create read ID with metadata
            read_id = f"sim_cecc_{i:04d}|unit_1|{len(full_sequence)}|1.0"

            record = EccDNARecord(
                eccdna_id=f"CeccDNA{i}",
                eccdna_type="Cecc",
                read_id=read_id,
                regions=regions,
                length=len(full_sequence),
                copy_number=1.0,
                sequence=full_sequence,
                is_inter_chr=is_inter,
                num_segments=num_segments,
            )
            records.append(record)

        return records


# =============================================================================
# Output Writers
# =============================================================================

class SimulationOutputWriter:
    """Write simulation outputs."""

    def __init__(self, output_dir: Path):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def write_reads_fasta(self, records: List[EccDNARecord], filename: str):
        """Write reads to FASTA file."""
        output_path = self.output_dir / filename
        with open(output_path, 'w') as f:
            for record in records:
                # Write read with metadata in header
                f.write(f">{record.read_id}\n")
                # Write sequence in 80-character lines
                seq = record.sequence
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i+80] + "\n")
        return output_path

    def write_ground_truth(self, records: List[EccDNARecord], filename: str):
        """Write ground truth annotations to CSV."""
        output_path = self.output_dir / filename

        if not records:
            return output_path

        fieldnames = list(records[0].to_dict().keys())

        with open(output_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for record in records:
                writer.writerow(record.to_dict())

        return output_path

    def write_summary(self, config: SimulationConfig,
                      all_records: List[EccDNARecord], filename: str):
        """Write simulation summary."""
        output_path = self.output_dir / filename

        # Count by type
        type_counts = {}
        for record in all_records:
            t = record.eccdna_type
            type_counts[t] = type_counts.get(t, 0) + 1

        # Calculate stats
        summary = {
            "config": asdict(config),
            "statistics": {
                "total_records": len(all_records),
                "by_type": type_counts,
                "total_bases": sum(r.length for r in all_records),
                "avg_length": sum(r.length for r in all_records) / len(all_records) if all_records else 0,
            },
            "files": {
                "reference": "reference.fa",
                "reads": "simulated_reads.fa",
                "ground_truth": "ground_truth.csv",
                "ground_truth_uecc": "ground_truth_uecc.csv",
                "ground_truth_mecc": "ground_truth_mecc.csv",
                "ground_truth_cecc": "ground_truth_cecc.csv",
            }
        }

        with open(output_path, 'w') as f:
            json.dump(summary, f, indent=2, default=str)

        return output_path


# =============================================================================
# Main Simulation Pipeline
# =============================================================================

def run_simulation(config: Optional[SimulationConfig] = None) -> Path:
    """Run the complete simulation pipeline."""

    if config is None:
        config = SimulationConfig()

    # Set random seed
    random.seed(config.seed)

    output_dir = Path(config.output_dir)
    writer = SimulationOutputWriter(output_dir)

    print("=" * 60)
    print("CircleSeeker eccDNA Simulation")
    print("=" * 60)

    # Step 1: Generate reference genome
    print("\n[1/5] Generating reference genome...")
    genome_gen = ReferenceGenomeGenerator(config)
    genome = genome_gen.generate()
    ref_path = output_dir / "reference.fa"
    genome_gen.write_fasta(ref_path)
    print(f"  - Generated {config.num_chromosomes} chromosomes")
    print(f"  - Total size: {config.num_chromosomes * config.chr_length:,} bp")
    if genome_gen.repeat_families:
        total_repeat_copies = sum(len(f.copies) for f in genome_gen.repeat_families)
        print(
            f"  - Embedded {len(genome_gen.repeat_families)} repeat families "
            f"({total_repeat_copies} total copies) for MeccDNA"
        )
    print(f"  - Output: {ref_path}")

    # Step 2: Simulate UeccDNA
    print(f"\n[2/5] Simulating {config.num_uecc} UeccDNA...")
    uecc_sim = UeccDNASimulator(genome, config, excluded_intervals_by_chrom=genome_gen.repeat_intervals_by_chrom)
    uecc_records = uecc_sim.generate(config.num_uecc)
    print(f"  - Length range: {config.uecc_length_range}")

    # Step 3: Simulate MeccDNA
    print(f"\n[3/5] Simulating {config.num_mecc} MeccDNA...")
    mecc_sim = MeccDNASimulator(genome, config, repeat_families=genome_gen.repeat_families)
    mecc_records = mecc_sim.generate(config.num_mecc)
    print(f"  - Unit length range: {config.mecc_unit_range}")
    print(f"  - Copy number range: {config.mecc_copy_range}")

    # Step 4: Simulate CeccDNA
    print(f"\n[4/5] Simulating {config.num_cecc} CeccDNA...")
    cecc_sim = CeccDNASimulator(genome, config, excluded_intervals_by_chrom=genome_gen.repeat_intervals_by_chrom)
    cecc_records = cecc_sim.generate(config.num_cecc)
    num_inter = sum(1 for r in cecc_records if r.is_inter_chr)
    print(f"  - Segments range: {config.cecc_segments_range}")
    print(f"  - Inter-chromosomal: {num_inter}/{config.num_cecc}")

    # Step 5: Write outputs
    print("\n[5/5] Writing outputs...")

    all_records = uecc_records + mecc_records + cecc_records

    # Combined reads file
    reads_path = writer.write_reads_fasta(all_records, "simulated_reads.fa")
    print(f"  - Reads: {reads_path}")

    # Ground truth files
    gt_all = writer.write_ground_truth(all_records, "ground_truth.csv")
    gt_uecc = writer.write_ground_truth(uecc_records, "ground_truth_uecc.csv")
    gt_mecc = writer.write_ground_truth(mecc_records, "ground_truth_mecc.csv")
    gt_cecc = writer.write_ground_truth(cecc_records, "ground_truth_cecc.csv")
    print(f"  - Ground truth: {gt_all}")

    # Summary
    summary_path = writer.write_summary(config, all_records, "simulation_summary.json")
    print(f"  - Summary: {summary_path}")

    # Print final statistics
    print("\n" + "=" * 60)
    print("Simulation Complete!")
    print("=" * 60)
    print(f"\nGenerated eccDNA:")
    print(f"  - UeccDNA: {len(uecc_records)}")
    print(f"  - MeccDNA: {len(mecc_records)}")
    print(f"  - CeccDNA: {len(cecc_records)}")
    print(f"  - Total: {len(all_records)}")
    print(f"\nOutput directory: {output_dir.absolute()}")

    return output_dir


# =============================================================================
# CLI Interface
# =============================================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate simulated eccDNA data for CircleSeeker validation"
    )
    parser.add_argument(
        "-o", "--output",
        default="tests/simulation/simulation_data",
        help="Output directory (default: tests/simulation/simulation_data)"
    )
    parser.add_argument(
        "--num-uecc", type=int, default=100,
        help="Number of UeccDNA to simulate (default: 100)"
    )
    parser.add_argument(
        "--num-mecc", type=int, default=20,
        help="Number of MeccDNA to simulate (default: 20)"
    )
    parser.add_argument(
        "--num-cecc", type=int, default=20,
        help="Number of CeccDNA to simulate (default: 20)"
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for reproducibility (default: 42)"
    )
    parser.add_argument(
        "--chr-length", type=int, default=1_000_000,
        help="Length of each chromosome in bp (default: 1000000)"
    )
    parser.add_argument(
        "--num-chr", type=int, default=5,
        help="Number of chromosomes (default: 5)"
    )

    args = parser.parse_args()

    config = SimulationConfig(
        output_dir=args.output,
        num_uecc=args.num_uecc,
        num_mecc=args.num_mecc,
        num_cecc=args.num_cecc,
        seed=args.seed,
        chr_length=args.chr_length,
        num_chromosomes=args.num_chr,
    )

    run_simulation(config)
