"""Tests for read_filter module - Sieve FASTA filtering."""

from pathlib import Path
import importlib.util
import sys
import types

import pytest

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


def _load_read_filter_module():
    package_name = "circleseeker.modules"
    module_name = f"{package_name}.read_filter"

    if module_name in sys.modules:
        return sys.modules[module_name]

    import circleseeker  # noqa: F401  # Ensure base package is initialized

    if package_name not in sys.modules:
        package_module = types.ModuleType(package_name)
        package_module.__path__ = [str(SRC / "circleseeker" / "modules")]
        sys.modules[package_name] = package_module

    spec = importlib.util.spec_from_file_location(
        module_name,
        SRC / "circleseeker" / "modules" / "read_filter.py",
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


Sieve = _load_read_filter_module().Sieve
FilterStats = _load_read_filter_module().FilterStats


@pytest.fixture(autouse=True)
def _disable_samtools(monkeypatch):
    module = _load_read_filter_module()
    monkeypatch.setattr(module.shutil, "which", lambda _: None)


def test_filter_stats_properties():
    """Test FilterStats percentage calculations."""
    stats = FilterStats(total_reads=100, filtered_reads=30, retained_reads=70)
    assert stats.filtered_percentage == 30.0
    assert stats.retained_percentage == 70.0

    # Test zero division handling
    empty_stats = FilterStats()
    assert empty_stats.filtered_percentage == 0.0
    assert empty_stats.retained_percentage == 0.0


def test_sieve_initialization():
    """Test Sieve class initialization."""
    sieve = Sieve()
    assert sieve.reads_to_filter == set()
    assert sieve.stats.total_reads == 0


def test_load_tandem_to_ring_classification(tmp_path):
    """Test loading classification from TandemToRing CSV."""
    # Create a classification CSV
    csv_file = tmp_path / "tandem_to_ring.csv"
    csv_file.write_text("\n".join([
        "readName,readClass",
        "read1,CtcR-perfect",
        "read2,Normal",
        "read3,CtcR-inversion",
        "read4,Normal",
    ]))

    sieve = Sieve()
    sieve.load_tandem_to_ring_classification(csv_file)

    # Should have loaded 2 reads to filter (CtcR variants)
    assert len(sieve.reads_to_filter) == 2
    assert "read1" in sieve.reads_to_filter
    assert "read3" in sieve.reads_to_filter
    assert "read2" not in sieve.reads_to_filter
    assert sieve.stats.csv_total_reads == 4
    assert sieve.stats.csv_ctcr_reads == 2


def test_filter_fasta_files_without_classification(tmp_path):
    """Test filtering FASTA files when no classification is loaded (keeps all)."""
    fasta_path = tmp_path / "input.fasta"
    fasta_path.write_text("\n".join([
        ">read1",
        "ACGT",
        ">read2",
        "TTTT",
    ]))

    output_fasta = tmp_path / "filtered.fasta"

    sieve = Sieve()
    # No classification loaded, so nothing should be filtered
    stats = sieve.filter_fasta_files([fasta_path], output_fasta, generate_index=False)

    assert stats.total_reads == 2
    assert stats.retained_reads == 2
    assert stats.filtered_reads == 0

    retained_headers = [
        line
        for line in output_fasta.read_text().splitlines()
        if line.startswith(">")
    ]
    assert retained_headers == [">read1", ">read2"]


def test_filter_fasta_files_merges_inputs(tmp_path):
    """Test that multiple FASTA files are merged correctly."""
    fasta_a = tmp_path / "input_a.fasta"
    fasta_a.write_text("\n".join([
        ">keep1",
        "AAAA",
    ]))

    fasta_b = tmp_path / "input_b.fasta"
    fasta_b.write_text("\n".join([
        ">keep2",
        "GGGG",
    ]))

    missing_fasta = tmp_path / "missing.fasta"
    output_fasta = tmp_path / "filtered.fasta"

    sieve = Sieve()
    stats = sieve.filter_fasta_files(
        [fasta_a, missing_fasta, fasta_b],
        output_fasta,
        generate_index=False,
    )

    assert stats.total_reads == 2
    assert stats.retained_reads == 2
    # Check both reads are in output
    content = output_fasta.read_text()
    assert ">keep1" in content
    assert ">keep2" in content


def test_run_sieve_complete_pipeline(tmp_path):
    """Test the complete sieve pipeline with classification and filtering."""
    # Create classification CSV
    tandem_to_ring_csv = tmp_path / "tandem_to_ring.csv"
    tandem_to_ring_csv.write_text("\n".join([
        "readName,readClass",
        "keep_read,Normal",
        "filter_read,CtcR-perfect",
    ]))

    # Create input FASTA
    fasta_path = tmp_path / "input.fasta"
    fasta_path.write_text("\n".join([
        ">keep_read",
        "ACGT",
        ">filter_read",
        "TTTT",
    ]))

    output_fasta = tmp_path / "filtered.fasta"

    sieve = Sieve()
    stats = sieve.run_sieve(
        tandem_to_ring_csv=tandem_to_ring_csv,
        input_fastas=[fasta_path],
        output_fasta=output_fasta,
    )

    assert stats.csv_ctcr_reads == 1
    assert stats.filtered_reads == 1
    assert stats.retained_reads == 1
    assert output_fasta.read_text().splitlines() == [
        ">keep_read",
        "ACGT",
    ]


def test_run_sieve_with_custom_ctcr_classes(tmp_path):
    """Test filtering with custom CtcR class set."""
    tandem_to_ring_csv = tmp_path / "tandem_to_ring.csv"
    tandem_to_ring_csv.write_text("\n".join([
        "readName,readClass",
        "read1,CustomBad",
        "read2,Normal",
    ]))

    fasta_path = tmp_path / "input.fasta"
    fasta_path.write_text("\n".join([
        ">read1",
        "AAAA",
        ">read2",
        "CCCC",
    ]))

    output_fasta = tmp_path / "filtered.fasta"

    sieve = Sieve()
    stats = sieve.run_sieve(
        tandem_to_ring_csv=tandem_to_ring_csv,
        input_fastas=[fasta_path],
        output_fasta=output_fasta,
        ctcr_classes={"CustomBad"},  # Custom filter class
    )

    assert stats.filtered_reads == 1
    assert stats.retained_reads == 1
    # Only read2 should remain
    content = output_fasta.read_text()
    assert ">read2" in content
    assert ">read1" not in content


def test_load_tandem_to_ring_classification_file_not_found(tmp_path):
    """Test error handling when classification file doesn't exist."""
    sieve = Sieve()
    with pytest.raises(FileNotFoundError):
        sieve.load_tandem_to_ring_classification(tmp_path / "nonexistent.csv")


def test_load_tandem_to_ring_classification_tab_delimited(tmp_path):
    """Test loading tab-delimited classification file."""
    csv_file = tmp_path / "tandem_to_ring.tsv"
    csv_file.write_text("\t".join(["readName", "readClass"]) + "\n" +
                        "\t".join(["read1", "CtcR-hybrid"]))

    sieve = Sieve()
    sieve.load_tandem_to_ring_classification(csv_file)

    assert "read1" in sieve.reads_to_filter
    assert sieve.stats.csv_ctcr_reads == 1


def test_generate_faidx_skips_empty_fasta(tmp_path, monkeypatch):
    """Empty FASTA should skip faidx without calling samtools."""
    module = _load_read_filter_module()
    monkeypatch.setattr(module.shutil, "which", lambda _: "/usr/bin/samtools")
    sieve = module.Sieve()

    empty_fasta = tmp_path / "empty.fasta"
    empty_fasta.write_text("", encoding="utf-8")

    called = {"value": False}

    def fake_run(*_args, **_kwargs):
        called["value"] = True
        raise AssertionError("samtools should not be invoked for empty FASTA")

    monkeypatch.setattr(module.subprocess, "run", fake_run)

    assert sieve.generate_faidx_index(empty_fasta) is False
    assert called["value"] is False
