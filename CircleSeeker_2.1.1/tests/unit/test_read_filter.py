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


@pytest.fixture(autouse=True)
def _disable_samtools(monkeypatch):
    module = _load_read_filter_module()
    monkeypatch.setattr(module.shutil, "which", lambda _: None)


def test_run_sieve_filters_ctcr_reads(tmp_path):
    classification_path = tmp_path / "carousel.csv"
    classification_path.write_text(
        "\n".join([
            "readName,readClass",
            "keep_read,Other",
            "drop_read,CtcR-perfect",
        ])
    )

    fasta_path = tmp_path / "input.fasta"
    fasta_path.write_text(
        "\n".join([
            ">keep_read",
            "ACGT",
            ">drop_read",
            "TTTT",
        ])
    )

    output_fasta = tmp_path / "filtered.fasta"

    sieve = Sieve()
    stats = sieve.run_sieve(classification_path, [fasta_path], output_fasta)

    assert stats.total_reads == 2
    assert stats.retained_reads == 1
    assert stats.filtered_reads == 1

    retained_headers = [
        line
        for line in output_fasta.read_text().splitlines()
        if line.startswith(">")
    ]
    assert retained_headers == [">keep_read"]


def test_load_carousel_classification_tab_delimiter_and_custom_classes(tmp_path):
    classification_path = tmp_path / "carousel.tsv"
    classification_path.write_text(
        "\n".join([
            "ReadName\tReadClass",
            "keep_read\tOther",
            "custom_drop\tCustom",
            "default_drop\tCtcR-perfect",
        ])
    )

    sieve = Sieve()
    sieve.load_carousel_classification(classification_path, {"Custom"})

    assert sieve.reads_to_filter == {"custom_drop"}
    assert sieve.stats.csv_total_reads == 3
    assert sieve.stats.csv_ctcr_reads == 1


def test_filter_fasta_files_merges_inputs_and_updates_stats(tmp_path):
    classification_path = tmp_path / "carousel.csv"
    classification_path.write_text(
        "\n".join([
            "readName,readClass",
            "keep1,Other",
            "drop1,CtcR-hybrid",
        ])
    )

    sieve = Sieve()
    sieve.load_carousel_classification(classification_path)

    fasta_a = tmp_path / "input_a.fasta"
    fasta_a.write_text(
        "\n".join([
            ">keep1",
            "AAAA",
            ">drop1",
            "CCCC",
        ])
    )

    fasta_b = tmp_path / "input_b.fasta"
    fasta_b.write_text(
        "\n".join([
            ">keep_extra",
            "GGGG",
        ])
    )

    missing_fasta = tmp_path / "missing.fasta"
    output_fasta = tmp_path / "filtered.fasta"

    stats = sieve.filter_fasta_files(
        [fasta_a, missing_fasta, fasta_b],
        output_fasta,
        generate_index=False,
    )

    assert stats is sieve.stats
    assert stats.total_reads == 3
    assert stats.retained_reads == 2
    assert stats.filtered_reads == 1
    assert stats.filtered_percentage == pytest.approx(33.3333, rel=1e-3)
    assert stats.retained_percentage == pytest.approx(66.6667, rel=1e-3)
    assert output_fasta.read_text().splitlines() == [
        ">keep1",
        "AAAA",
        ">keep_extra",
        "GGGG",
    ]
