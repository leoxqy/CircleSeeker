from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
TESTS_DIR = ROOT / "tests"
if str(TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(TESTS_DIR))


from simulation.eccdna_simulator import (  # noqa: E402
    SimulationConfig,
    ReferenceGenomeGenerator,
    MeccDNASimulator,
    UeccDNASimulator,
)


def test_mecc_simulation_uses_multiple_loci():
    config = SimulationConfig(
        num_chromosomes=2,
        chr_length=2000,
        num_uecc=0,
        num_mecc=3,
        num_cecc=0,
        mecc_unit_range=(30, 30),
        mecc_copy_range=(2.0, 4.0),
        # Make copies identical so the expectation is unambiguous.
        mecc_repeat_mutation_rate=0.0,
        error_rate=0.0,
        seed=123,
    )

    genome_gen = ReferenceGenomeGenerator(config)
    genome = genome_gen.generate()

    simulator = MeccDNASimulator(genome, config, repeat_families=genome_gen.repeat_families)
    records = simulator.generate(config.num_mecc)

    assert len(records) == config.num_mecc
    for record in records:
        loci = {(r.chrom, r.start, r.end, r.strand) for r in record.regions}
        assert record.num_segments == len(record.regions)
        assert len(loci) == record.num_segments
        assert len(loci) >= 2

        # MeccDNA sequence is a single repeat unit that maps to multiple loci
        source_region = record.regions[0]
        unit_seq = (
            genome[source_region.chrom][source_region.start : source_region.end]
            if source_region.strand == "+"
            else ReferenceGenomeGenerator._reverse_complement(
                genome[source_region.chrom][source_region.start : source_region.end]
            )
        )
        unit_length = source_region.length

        assert record.length == len(record.sequence)
        assert record.length == unit_length
        assert record.sequence == unit_seq


def test_uecc_simulation_avoids_mecc_repeat_intervals():
    config = SimulationConfig(
        num_chromosomes=2,
        chr_length=2000,
        num_uecc=50,
        num_mecc=5,
        num_cecc=0,
        uecc_length_range=(60, 80),
        mecc_unit_range=(30, 30),
        mecc_copy_range=(2.0, 3.0),
        error_rate=0.0,
        seed=7,
    )

    genome_gen = ReferenceGenomeGenerator(config)
    genome = genome_gen.generate()

    simulator = UeccDNASimulator(
        genome,
        config,
        excluded_intervals_by_chrom=genome_gen.repeat_intervals_by_chrom,
    )
    records = simulator.generate(config.num_uecc)

    for record in records:
        region = record.regions[0]
        excluded = genome_gen.repeat_intervals_by_chrom.get(region.chrom, [])
        assert all(region.end <= start or region.start >= end for start, end in excluded)
