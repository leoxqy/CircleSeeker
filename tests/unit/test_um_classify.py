from pathlib import Path
import sys

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.modules.um_classify import UMeccClassifier

def test_um_classify_subset(tmp_path):
    blast_tsv = tmp_path / "step3_subset.tsv"

    q_uecc = "readU|rep1|100|1"
    q_mecc = "readM|rep1|100|2"
    q_unclassified = "readX|rep1|100|1"

    rows = [
        # Single high-quality alignment -> Uecc
        [q_uecc, "chr1", 99.0, 100, 0, 0, 1, 100, 100, 199, 0.0, 200.0, "plus"],
        # Two non-overlapping, full-length alignments -> Mecc
        [q_mecc, "chr1", 99.0, 100, 0, 0, 1, 100, 100, 199, 0.0, 200.0, "plus"],
        [q_mecc, "chr1", 98.5, 100, 0, 0, 1, 100, 400, 499, 0.0, 190.0, "plus"],
        # High-quality but not full-length coverage -> unclassified
        [q_unclassified, "chr2", 99.0, 90, 0, 0, 1, 90, 100, 189, 0.0, 180.0, "plus"],
        [q_unclassified, "chr2", 98.0, 90, 0, 0, 1, 90, 300, 389, 0.0, 170.0, "plus"],
    ]

    pd.DataFrame(rows).to_csv(blast_tsv, sep="\t", header=False, index=False)

    classifier = UMeccClassifier()
    uecc_df, mecc_df, unclassified_df = classifier.run(blast_tsv)

    assert set(uecc_df["query_id"].unique()) == {q_uecc}
    assert set(mecc_df["query_id"].unique()) == {q_mecc}
    assert set(unclassified_df["query_id"].unique()) == {q_unclassified}

    assert set(uecc_df["eccdna_type"].unique()) == {"Uecc"}
    assert set(mecc_df["eccdna_type"].unique()) == {"Mecc"}

    assert len(unclassified_df) == 2
    assert {"unclass_reason", "quality_category"}.issubset(unclassified_df.columns)
    assert set(unclassified_df["quality_category"].unique()) == {"High_quality"}


def test_um_classify_invalid_column_count(tmp_path):
    classifier = UMeccClassifier()
    alignment_df = pd.DataFrame([["query1", "subject1"]])

    with pytest.raises(ValueError):
        classifier.classify_alignment_results(alignment_df, tmp_path / "invalid")


def test_um_classify_overlap_components_mecc(tmp_path):
    classifier = UMeccClassifier()
    query_id = "readA|rep1|100|2"

    alignment_df = pd.DataFrame(
        [
            [query_id, "chr1", 99.0, 100, 1, 0, 1, 100, 100, 199, 0, 0, "plus"],
            [query_id, "chr1", 99.0, 100, 1, 0, 1, 100, 150, 249, 0, 0, "plus"],
            [query_id, "chr1", 99.0, 100, 1, 0, 1, 100, 400, 499, 0, 0, "plus"],
            [query_id, "chr1", 99.0, 100, 1, 0, 1, 100, 450, 549, 0, 0, "plus"],
        ]
    )

    uecc_df, mecc_df, unclassified_df = classifier.classify_alignment_results(
        alignment_df, tmp_path / "um_classify"
    )

    assert uecc_df.empty
    assert mecc_df["query_id"].nunique() == 1
    assert len(mecc_df) == 2
    assert unclassified_df.empty


def test_um_classify_uecc_requires_low_second_locus(tmp_path):
    """If a strong 2nd locus exists, the query should not be labeled as Uecc."""
    blast_tsv = tmp_path / "step3_second_locus.tsv"

    query_id = "readU2|rep1|100|1"

    rows = [
        # Best locus: full ring coverage
        [query_id, "chr1", 99.0, 100, 0, 0, 1, 100, 100, 199, 0.0, 200.0, "plus"],
        # Second locus: substantial (90%) ring coverage, but below theta_m=0.95
        [query_id, "chr1", 99.0, 90, 0, 0, 1, 90, 500, 589, 0.0, 190.0, "plus"],
    ]

    pd.DataFrame(rows).to_csv(blast_tsv, sep="\t", header=False, index=False)

    classifier = UMeccClassifier(theta_u=0.95, theta_m=0.95, theta_u2_max=0.05)
    uecc_df, mecc_df, unclassified_df = classifier.run(blast_tsv)

    assert uecc_df.empty
    assert mecc_df.empty
    assert set(unclassified_df["query_id"].unique()) == {query_id}


def test_um_classify_mecc_can_use_lower_theta_m(tmp_path):
    """Mecc can be recovered by lowering theta_m while keeping U strict."""
    blast_tsv = tmp_path / "step3_theta_m.tsv"

    query_id = "readM2|rep1|100|1"

    rows = [
        # Locus A: full ring coverage
        [query_id, "chr1", 99.0, 100, 0, 0, 1, 100, 100, 199, 0.0, 200.0, "plus"],
        # Locus B: 90% ring coverage (should count toward Mecc when theta_m<=0.9)
        [query_id, "chr1", 99.0, 90, 0, 0, 1, 90, 500, 589, 0.0, 190.0, "plus"],
    ]

    pd.DataFrame(rows).to_csv(blast_tsv, sep="\t", header=False, index=False)

    classifier = UMeccClassifier(theta_u=0.95, theta_m=0.85, theta_u2_max=0.05)
    uecc_df, mecc_df, unclassified_df = classifier.run(blast_tsv)

    assert uecc_df.empty
    assert set(mecc_df["query_id"].unique()) == {query_id}
    assert mecc_df["query_id"].nunique() == 1
    assert len(mecc_df) == 2
    assert unclassified_df.empty


def test_um_classify_uecc_vetoes_significant_secondary_mapping(tmp_path):
    """A strong single-locus coverage can still be vetoed by significant secondary evidence."""
    blast_tsv = tmp_path / "step3_secondary_evidence.tsv"

    query_id = "readC|rep1|10000|1"

    rows = [
        # Primary locus: ~99% ring coverage
        [query_id, "chr1", 99.0, 9900, 0, 0, 1, 9900, 1000, 10899, 0.0, 200.0, "plus"],
        # Secondary locus: 1% ring coverage on another chromosome (low-quality by Gap_Percentage)
        [query_id, "chr5", 99.0, 100, 0, 0, 9901, 10000, 5000, 5099, 0.0, 50.0, "plus"],
    ]

    pd.DataFrame(rows).to_csv(blast_tsv, sep="\t", header=False, index=False)

    classifier = UMeccClassifier(theta_u=0.95, theta_m=0.95, theta_u2_max=0.05)
    uecc_df, mecc_df, unclassified_df = classifier.run(blast_tsv)

    assert uecc_df.empty
    assert mecc_df.empty
    assert set(unclassified_df["query_id"].unique()) == {query_id}


def test_um_classify_uecc_can_be_gated_by_mapq(tmp_path):
    """When enabled, mapq_u_min should veto Uecc calls in low-uniqueness mappings."""
    blast_tsv = tmp_path / "step3_mapq_gate.tsv"

    query_id = "readU_mapq|rep1|100|1"

    rows = [
        # Full ring coverage, but mapq is 0 (ambiguous / repetitive mapping)
        [query_id, "chr1", 99.0, 100, 0, 0, 1, 100, 100, 199, 0.0, 200.0, "plus", 0],
    ]

    pd.DataFrame(rows).to_csv(blast_tsv, sep="\t", header=False, index=False)

    classifier = UMeccClassifier(theta_u=0.95, theta_m=0.95, theta_u2_max=0.05, mapq_u_min=20)
    uecc_df, mecc_df, unclassified_df = classifier.run(blast_tsv)

    assert uecc_df.empty
    assert mecc_df.empty
    assert set(unclassified_df["query_id"].unique()) == {query_id}


def test_um_classify_confidence_score_is_monotonic_in_mapq(tmp_path):
    """For the same read, confidence_score should increase as MAPQ increases."""
    query_id = "readU_conf|rep1|100|1"

    def run_with_mapq(mapq: int) -> float:
        blast_tsv = tmp_path / f"step3_confidence_mapq_{mapq}.tsv"
        rows = [
            [query_id, "chr1", 99.0, 100, 0, 0, 1, 100, 100, 199, 0.0, 200.0, "plus", mapq],
        ]
        pd.DataFrame(rows).to_csv(blast_tsv, sep="\t", header=False, index=False)
        classifier = UMeccClassifier(theta_u=0.95, theta_m=0.95, theta_u2_max=0.05, mapq_u_min=0)
        uecc_df, mecc_df, unclassified_df = classifier.run(blast_tsv)
        assert set(uecc_df["query_id"].unique()) == {query_id}
        assert mecc_df.empty
        assert unclassified_df.empty
        assert "confidence_score" in uecc_df.columns
        return float(uecc_df["confidence_score"].iloc[0])

    score_0 = run_with_mapq(0)
    score_20 = run_with_mapq(20)
    score_60 = run_with_mapq(60)
    assert score_0 <= score_20 <= score_60
