from pathlib import Path
import sys

import pandas as pd
import pytest
from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.modules.tandem_to_ring import TandemToRing
import numpy as np
from circleseeker.utils.column_standards import ColumnStandard


def test_tandem_to_ring_process_subset(tmp_path):
    input_tsv = tmp_path / "step1_subset.tsv"

    # Input format matches TideHunter TSV (11 columns, no header)
    rows = [
        # Effective_Length = (1000 / 1000) * 100 = 100 -> CtcR-perfect
        ["readA", 1, 1.0, 1000, 1, 1000, 4, 99.0, 0, "0,1", "ATCG"],
        # Effective_Length = (800 / 1000) * 100 = 80 -> CtcR-hybrid
        ["readB", 1, 1.0, 1000, 1, 800, 4, 99.0, 0, "0,1", "GGCC"],
    ]
    pd.DataFrame(rows).to_csv(input_tsv, sep="\t", header=False, index=False)

    output_csv = tmp_path / "step2_processed.csv"
    output_fasta = tmp_path / "step2_circular.fasta"

    module = TandemToRing(
        input_file=input_tsv,
        output_file=output_csv,
        circular_fasta=output_fasta,
    )

    df_main, df_classification = module.process()

    produced_csv = pd.read_csv(output_csv).sort_values("readName").reset_index(drop=True)
    expected_df = pd.DataFrame(
        {"readName": ["readA", "readB"], "readClass": ["CtcR-perfect", "CtcR-hybrid"]}
    ).sort_values("readName").reset_index(drop=True)

    pd.testing.assert_frame_equal(produced_csv, expected_df)
    pd.testing.assert_frame_equal(
        df_classification.sort_values("readName").reset_index(drop=True), expected_df
    )

    produced_records = list(SeqIO.parse(output_fasta, "fasta"))
    produced_records.sort(key=lambda rec: rec.id)

    expected_fasta = {
        "readA|1|4|1|circular": "ATCGATCG",
        "readB|1|4|1|circular": "GGCCGGCC",
    }
    assert [rec.id for rec in produced_records] == sorted(expected_fasta.keys())
    for rec in produced_records:
        assert str(rec.seq) == expected_fasta[rec.id]

    assert not df_main.empty
    unique_ids = {f"{uid}|circular" for uid in df_main["unique_id"].tolist()}
    assert unique_ids == {rec.id for rec in produced_records}


class TestReadAndPreprocessData:
    def _write_tsv(self, tmp_path, rows):
        f = tmp_path / "input.tsv"
        pd.DataFrame(rows).to_csv(f, sep="\t", header=False, index=False)
        return f

    def test_basic(self, tmp_path):
        f = self._write_tsv(tmp_path, [
            ["readA", 1, 1.0, 1000, 1, 1000, 100, 99.5, 0, "0,1", "ATCG"],
        ])
        module = TandemToRing(f, tmp_path / "out.csv", tmp_path / "out.fasta")
        df = module.read_and_preprocess_data()
        assert len(df) == 1
        assert "Effective_Length" in df.columns

    def test_filter_low_quality(self, tmp_path):
        f = self._write_tsv(tmp_path, [
            ["readA", 1, 1.0, 1000, 1, 1000, 100, 99.5, 0, "0,1", "ATCG"],
            ["readB", 1, 1.0, 1000, 1, 1000, 100, 50.0, 0, "0,1", "GGCC"],  # Below default 99.0
        ])
        module = TandemToRing(f, tmp_path / "out.csv", tmp_path / "out.fasta")
        df = module.read_and_preprocess_data()
        assert len(df) == 1
        assert df.iloc[0][ColumnStandard.READS] == "readA"

    def test_effective_length_calculation(self, tmp_path):
        f = self._write_tsv(tmp_path, [
            ["readA", 1, 1.0, 1000, 1, 500, 100, 99.5, 0, "0,1", "ATCG"],
        ])
        module = TandemToRing(f, tmp_path / "out.csv", tmp_path / "out.fasta")
        df = module.read_and_preprocess_data()
        expected = (500 / 1000) * 100  # 50.0
        assert abs(df.iloc[0]["Effective_Length"] - expected) < 0.01

    def test_subpos_conversion(self, tmp_path):
        f = self._write_tsv(tmp_path, [
            ["readA", 1, 1.0, 1000, 1, 1000, 100, 99.5, 0, "0,1,2", "ATCG"],
        ])
        module = TandemToRing(f, tmp_path / "out.csv", tmp_path / "out.fasta")
        df = module.read_and_preprocess_data()
        assert "," not in str(df.iloc[0]["subPos"])
        assert "|" in str(df.iloc[0]["subPos"])


class TestSeparateSimpleComplex:
    def _make_df(self, reads_list):
        rows = []
        for r in reads_list:
            rows.append({
                ColumnStandard.READS: r, "repN": 1, "copyNum": 1.0,
                "readLen": 1000, ColumnStandard.START0: 1, ColumnStandard.END0: 1000,
                "consLen": 100, "aveMatch": 99.5, "subPos": "0|1",
                "consSeq": "ATCG", "Effective_Length": 100.0,
                "readName": r, "start": 1, "end": 1000,
            })
        return pd.DataFrame(rows)

    def test_all_simple(self, tmp_path):
        df = self._make_df(["readA", "readB"])
        module = TandemToRing(tmp_path / "in.tsv", tmp_path / "out.csv", tmp_path / "out.fasta")
        simple, complex_ = module.separate_simple_complex(df)
        assert len(simple) == 2
        assert len(complex_) == 0

    def test_all_complex(self, tmp_path):
        df = self._make_df(["readA", "readA", "readB", "readB"])
        module = TandemToRing(tmp_path / "in.tsv", tmp_path / "out.csv", tmp_path / "out.fasta")
        simple, complex_ = module.separate_simple_complex(df)
        assert len(simple) == 0
        assert len(complex_) == 4

    def test_mixed(self, tmp_path):
        df = self._make_df(["readA", "readB", "readB"])
        module = TandemToRing(tmp_path / "in.tsv", tmp_path / "out.csv", tmp_path / "out.fasta")
        simple, complex_ = module.separate_simple_complex(df)
        assert len(simple) == 1
        assert len(complex_) == 2


class TestClassifySimpleReads:
    def _make_simple(self, eff_length):
        return pd.DataFrame({
            ColumnStandard.READS: ["readA"], "Effective_Length": [eff_length],
        })

    def test_perfect(self, tmp_path):
        module = TandemToRing(tmp_path / "in.tsv", tmp_path / "out.csv", tmp_path / "out.fasta")
        df = module.classify_simple_reads(self._make_simple(100.0))
        assert df.iloc[0]["classification"] == "CtcR-perfect"

    def test_hybrid(self, tmp_path):
        module = TandemToRing(tmp_path / "in.tsv", tmp_path / "out.csv", tmp_path / "out.fasta")
        df = module.classify_simple_reads(self._make_simple(85.0))
        assert df.iloc[0]["classification"] == "CtcR-hybrid"

    def test_other(self, tmp_path):
        module = TandemToRing(tmp_path / "in.tsv", tmp_path / "out.csv", tmp_path / "out.fasta")
        df = module.classify_simple_reads(self._make_simple(50.0))
        assert df.iloc[0]["classification"] == "Other"

    def test_boundary_99(self, tmp_path):
        module = TandemToRing(tmp_path / "in.tsv", tmp_path / "out.csv", tmp_path / "out.fasta")
        df = module.classify_simple_reads(self._make_simple(99.0))
        assert df.iloc[0]["classification"] == "CtcR-perfect"


class TestProcessComplexGroupGraphOptimized:
    def _make_module(self, tmp_path):
        return TandemToRing(tmp_path / "in.tsv", tmp_path / "out.csv", tmp_path / "out.fasta")

    def test_no_overlap(self, tmp_path):
        module = self._make_module(tmp_path)
        df = pd.DataFrame({
            ColumnStandard.START0: [1, 500], ColumnStandard.END0: [100, 600],
        })
        result = module.process_complex_group_with_graph_optimized(df)
        assert len(result) == 2

    def test_overlap_keeps_longest(self, tmp_path):
        module = self._make_module(tmp_path)
        df = pd.DataFrame({
            ColumnStandard.START0: [1, 10], ColumnStandard.END0: [200, 150],
        })
        result = module.process_complex_group_with_graph_optimized(df)
        assert len(result) == 1
        assert result.iloc[0][ColumnStandard.END0] == 200

    def test_single_region(self, tmp_path):
        module = self._make_module(tmp_path)
        df = pd.DataFrame({
            ColumnStandard.START0: [1], ColumnStandard.END0: [100],
        })
        result = module.process_complex_group_with_graph_optimized(df)
        assert len(result) == 1

    def test_multiple_components(self, tmp_path):
        module = self._make_module(tmp_path)
        df = pd.DataFrame({
            ColumnStandard.START0: [1, 10, 500, 510],
            ColumnStandard.END0: [200, 150, 700, 650],
        })
        result = module.process_complex_group_with_graph_optimized(df)
        assert len(result) == 2  # One from each component

    def test_transitive_overlap(self, tmp_path):
        module = self._make_module(tmp_path)
        # A overlaps B, B overlaps C -> all in one component
        df = pd.DataFrame({
            ColumnStandard.START0: [1, 100, 200],
            ColumnStandard.END0: [200, 300, 400],
        })
        result = module.process_complex_group_with_graph_optimized(df)
        assert len(result) == 1


class TestAnalyzeConsLenConsistency:
    def _make_module(self, tmp_path):
        return TandemToRing(tmp_path / "in.tsv", tmp_path / "out.csv", tmp_path / "out.fasta")

    def _make_df(self, read_name, consLen_values):
        rows = []
        for cl in consLen_values:
            rows.append({
                ColumnStandard.READS: read_name, "consLen": cl,
                "Effective_Length": 100.0, "group_total_effective_length": 100.0,
            })
        return pd.DataFrame(rows)

    def test_highly_consistent(self, tmp_path):
        module = self._make_module(tmp_path)
        df = self._make_df("readA", [100, 100, 100])
        result = module.analyze_consLen_consistency(df)
        assert result.iloc[0]["consistency_type"] == "highly_consistent"

    def test_consistent(self, tmp_path):
        module = self._make_module(tmp_path)
        df = self._make_df("readA", [100, 102, 99])
        result = module.analyze_consLen_consistency(df)
        assert result.iloc[0]["consistency_type"] == "consistent"

    def test_integer_multiple(self, tmp_path):
        module = self._make_module(tmp_path)
        df = self._make_df("readA", [100, 200, 300])
        result = module.analyze_consLen_consistency(df)
        assert result.iloc[0]["consistency_type"] == "integer_multiple"

    def test_variable(self, tmp_path):
        module = self._make_module(tmp_path)
        df = self._make_df("readA", [100, 500, 37])
        result = module.analyze_consLen_consistency(df)
        assert result.iloc[0]["consistency_type"] == "variable"

    def test_single_value(self, tmp_path):
        module = self._make_module(tmp_path)
        df = self._make_df("readA", [100])
        result = module.analyze_consLen_consistency(df)
        assert result.iloc[0]["consistency_type"] == "highly_consistent"


class TestProcessAndClassifyComplexReads:
    def _make_module(self, tmp_path):
        return TandemToRing(tmp_path / "in.tsv", tmp_path / "out.csv", tmp_path / "out.fasta")

    def _make_complex_df(self, read_name, eff_lengths, consLens=None):
        rows = []
        if consLens is None:
            consLens = [100] * len(eff_lengths)
        for i, (el, cl) in enumerate(zip(eff_lengths, consLens)):
            rows.append({
                ColumnStandard.READS: read_name,
                "Effective_Length": el,
                "consLen": cl, "copyNum": 1.0, "repN": 1,
                ColumnStandard.START0: i * 100, ColumnStandard.END0: (i + 1) * 100,
                "group_total_effective_length": sum(eff_lengths),
            })
        return pd.DataFrame(rows)

    def _make_consistency(self, read_name, consistency_type, num_regions):
        return pd.DataFrame({
            "readName": [read_name],
            "consistency_type": [consistency_type],
            "num_regions": [num_regions],
        })

    def test_inversion(self, tmp_path):
        module = self._make_module(tmp_path)
        df = self._make_complex_df("readA", [40.0, 40.0])
        cons = self._make_consistency("readA", "highly_consistent", 2)
        result = module.process_and_classify_complex_reads(df, cons)
        assert result.iloc[0]["classification"] == "CtcR-inversion"

    def test_below_threshold(self, tmp_path):
        module = self._make_module(tmp_path)
        df = self._make_complex_df("readA", [30.0, 30.0])
        cons = self._make_consistency("readA", "highly_consistent", 2)
        result = module.process_and_classify_complex_reads(df, cons)
        assert result.iloc[0]["classification"] == "Other"

    def test_single_record_perfect(self, tmp_path):
        module = self._make_module(tmp_path)
        df = self._make_complex_df("readA", [100.0])
        cons = self._make_consistency("readA", "variable", 1)
        result = module.process_and_classify_complex_reads(df, cons)
        assert result.iloc[0]["classification"] == "CtcR-perfect"

    def test_single_record_hybrid(self, tmp_path):
        module = self._make_module(tmp_path)
        df = self._make_complex_df("readA", [85.0])
        cons = self._make_consistency("readA", "variable", 1)
        result = module.process_and_classify_complex_reads(df, cons)
        assert result.iloc[0]["classification"] == "CtcR-hybrid"

    def test_not_consistent_multi(self, tmp_path):
        module = self._make_module(tmp_path)
        df = self._make_complex_df("readA", [80.0, 80.0])
        cons = self._make_consistency("readA", "variable", 2)
        result = module.process_and_classify_complex_reads(df, cons)
        assert all(result["classification"] == "CtcR-hybrid")


class TestCreateMainResultsTable:
    def _make_module(self, tmp_path):
        return TandemToRing(tmp_path / "in.tsv", tmp_path / "out.csv", tmp_path / "out.fasta")

    def test_merge(self, tmp_path):
        module = self._make_module(tmp_path)
        simple = pd.DataFrame({
            ColumnStandard.READS: ["readA"], "repN": [1], "copyNum": [1.5],
            "consLen": [100], "consSeq": ["ATCG"], "classification": ["CtcR-perfect"],
        })
        complex_ = pd.DataFrame({
            ColumnStandard.READS: ["readB"], "repN": [2], "copyNum": [2.3],
            "consLen": [200], "consSeq": ["GGCC"], "classification": ["CtcR-hybrid"],
        })
        result = module.create_main_results_table(simple, complex_)
        assert len(result) == 2
        assert "unique_id" in result.columns

    def test_unique_id_format(self, tmp_path):
        module = self._make_module(tmp_path)
        df = pd.DataFrame({
            ColumnStandard.READS: ["readA"], "repN": [1], "copyNum": [1.5],
            "consLen": [100], "consSeq": ["ATCG"], "classification": ["CtcR-perfect"],
        })
        result = module.create_main_results_table(df, pd.DataFrame())
        uid = result.iloc[0]["unique_id"]
        assert "readA" in uid
        assert "|" in uid

    def test_copynum_rounded(self, tmp_path):
        module = self._make_module(tmp_path)
        df = pd.DataFrame({
            ColumnStandard.READS: ["readA"], "repN": [1], "copyNum": [2.7],
            "consLen": [100], "consSeq": ["ATCG"], "classification": ["CtcR-perfect"],
        })
        result = module.create_main_results_table(df, pd.DataFrame())
        assert result.iloc[0]["copyNum"] == 3


class TestCircularizeSequences:
    def _make_module(self, tmp_path):
        return TandemToRing(tmp_path / "in.tsv", tmp_path / "out.csv", tmp_path / "out.fasta")

    def test_doubled_sequence(self, tmp_path):
        module = self._make_module(tmp_path)
        df = pd.DataFrame({
            "unique_id": ["readA|1|100|1"], "consSeq": ["ATCG"],
        })
        records = module.circularize_sequences(df)
        assert len(records) == 1
        assert str(records[0].seq) == "ATCGATCG"

    def test_multiple_records(self, tmp_path):
        module = self._make_module(tmp_path)
        df = pd.DataFrame({
            "unique_id": ["r1|1|100|1", "r2|1|200|1"],
            "consSeq": ["AT", "GC"],
        })
        records = module.circularize_sequences(df)
        assert len(records) == 2
        assert str(records[0].seq) == "ATAT"
        assert str(records[1].seq) == "GCGC"


class TestCreateReadnameClassification:
    def _make_module(self, tmp_path):
        return TandemToRing(tmp_path / "in.tsv", tmp_path / "out.csv", tmp_path / "out.fasta")

    def test_single_read(self, tmp_path):
        module = self._make_module(tmp_path)
        df = pd.DataFrame({
            ColumnStandard.READS: ["readA"], "classification": ["CtcR-perfect"],
        })
        result = module.create_readname_classification(df)
        assert len(result) == 1
        assert result.iloc[0]["readClass"] == "CtcR-perfect"

    def test_perfect_over_hybrid(self, tmp_path):
        module = self._make_module(tmp_path)
        df = pd.DataFrame({
            ColumnStandard.READS: ["readA", "readA"],
            "classification": ["CtcR-hybrid", "CtcR-perfect"],
        })
        result = module.create_readname_classification(df)
        assert len(result) == 1
        assert result.iloc[0]["readClass"] == "CtcR-perfect"

    def test_inversion_over_hybrid(self, tmp_path):
        module = self._make_module(tmp_path)
        df = pd.DataFrame({
            ColumnStandard.READS: ["readA", "readA"],
            "classification": ["CtcR-hybrid", "CtcR-inversion"],
        })
        result = module.create_readname_classification(df)
        assert len(result) == 1
        assert result.iloc[0]["readClass"] == "CtcR-inversion"

    def test_all_other(self, tmp_path):
        module = self._make_module(tmp_path)
        df = pd.DataFrame({
            ColumnStandard.READS: ["readA", "readB"],
            "classification": ["Other", "Other"],
        })
        result = module.create_readname_classification(df)
        assert len(result) == 2
        assert all(result["readClass"] == "Other")
