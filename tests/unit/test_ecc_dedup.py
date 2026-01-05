from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.modules.ecc_dedup import CDHitClusters, eccDedup

def test_ecc_dedup_subset(tmp_path):
    input_dir = tmp_path / "inputs"
    input_dir.mkdir()

    uecc_csv = input_dir / "uecc_input.csv"
    uecc_clstr = input_dir / "uecc.id2cluster.csv"
    mecc_csv = input_dir / "mecc_input.csv"
    mecc_clstr = input_dir / "mecc.id2cluster.csv"
    cecc_csv = input_dir / "cecc_input.csv"
    cecc_clstr = input_dir / "cecc.id2cluster.csv"

    pd.DataFrame(
        [
            {
                "eccDNA_id": "U1",
                "chr": "chr1",
                "start0": 0,
                "end0": 100,
                "strand": "+",
                "length": 100,
                "match_degree": 99.5,
                "copy_number": 1,
                "reads": "readA",
                "eSeq": "A" * 100,
            },
            {
                "eccDNA_id": "U2",
                "chr": "chr1",
                "start0": 0,
                "end0": 100,
                "strand": "+",
                "length": 100,
                "match_degree": 98.0,
                "copy_number": 2,
                "reads": "readB",
                "eSeq": "C" * 100,
            },
            {
                "eccDNA_id": "U3",
                "chr": "chr2",
                "start0": 200,
                "end0": 300,
                "strand": "-",
                "length": 100,
                "match_degree": 97.0,
                "copy_number": 1,
                "reads": "readC",
                "eSeq": "G" * 100,
            },
        ]
    ).to_csv(uecc_csv, index=False)
    uecc_clstr.write_text(
        "id,cluster,is_representative\n"
        "U1,grp1,false\n"
        "U2,grp1,true\n"
        "U3,grp2,true\n"
    )

    pd.DataFrame(
        [
            {
                "eccDNA_id": "M1",
                "chr": "chr1",
                "start0": 1000,
                "end0": 1100,
                "strand": "+",
                "length": 100,
                "match_degree": 98.0,
                "copy_number": 2,
                "reads": "mread1",
                "alignment_length": 100,
            },
            {
                "eccDNA_id": "M2",
                "chr": "chr1",
                "start0": 2000,
                "end0": 2100,
                "strand": "+",
                "length": 100,
                "match_degree": 99.0,
                "copy_number": 2,
                "reads": "mread2",
                "alignment_length": 100,
            },
            {
                "eccDNA_id": "M2",
                "chr": "chr2",
                "start0": 3000,
                "end0": 3100,
                "strand": "+",
                "length": 100,
                "match_degree": 97.0,
                "copy_number": 2,
                "reads": "mread3",
                "alignment_length": 90,
            },
        ]
    ).to_csv(mecc_csv, index=False)
    mecc_clstr.write_text(
        "id,cluster,is_representative\n"
        "M1,mgrp1,false\n"
        "M2,mgrp1,true\n"
    )

    pd.DataFrame(
        [
            {
                "eccDNA_id": "C1",
                "chr": "chr1",
                "start0": 0,
                "end0": 100,
                "strand": "+",
                "length": 200,
                "match_degree": 90.0,
                "copy_number": 1,
                "reads": "cread1",
            },
            {
                "eccDNA_id": "C2",
                "chr": "chr1",
                "start0": 100,
                "end0": 200,
                "strand": "+",
                "length": 200,
                "match_degree": 95.0,
                "copy_number": 1,
                "reads": "cread2",
            },
            {
                "eccDNA_id": "C2",
                "chr": "chr1",
                "start0": 300,
                "end0": 400,
                "strand": "+",
                "length": 200,
                "match_degree": 95.0,
                "copy_number": 1,
                "reads": "cread3",
            },
        ]
    ).to_csv(cecc_csv, index=False)
    cecc_clstr.write_text(
        "id,cluster,is_representative\n"
        "C1,cgrp1,false\n"
        "C2,cgrp1,true\n"
    )

    processor = eccDedup()
    results = processor.process_all_types(
        uecc_csv=uecc_csv,
        uecc_clstr=uecc_clstr,
        mecc_csv=mecc_csv,
        mecc_clstr=mecc_clstr,
        cecc_csv=cecc_csv,
        cecc_clstr=cecc_clstr,
        output_dir=tmp_path,
        prefix="step8",
        drop_seq=True,
        generate_confirmed_tables=False,
        organize_output_files=False,
    )

    assert set(results.keys()) == {"Uecc", "Mecc", "Cecc"}

    uecc = results["Uecc"]
    assert uecc["cluster_id"].nunique() == 2
    assert uecc["eccDNA_id"].str.startswith("UeccDNA").all()
    grp1 = uecc.loc[uecc["cluster_id"] == "grp1"].iloc[0]
    assert grp1["num_merged"] == 2
    assert grp1["merged_from_ids"] == "U1;U2"

    mecc = results["Mecc"]
    assert mecc["cluster_id"].nunique() == 1
    assert mecc["eccDNA_id"].str.startswith("MeccDNA").all()
    assert len(mecc) == 2  # two loci retained for the representative ID
    assert mecc["num_merged"].dropna().astype(int).eq(2).all()
    assert mecc["merged_from_ids"].iloc[0] == "M1;M2"

    cecc = results["Cecc"]
    assert cecc["cluster_id"].nunique() == 1
    assert cecc["eccDNA_id"].str.startswith("CeccDNA").all()
    assert len(cecc) == 2  # two segments retained for the representative ID
    assert cecc["num_merged"].dropna().astype(int).eq(2).all()
    assert cecc["merged_from_ids"].iloc[0] == "C1;C2"

    expected_outputs = [
        tmp_path / "step8_UeccDNA.core.csv",
        tmp_path / "step8_UeccDNA.bed",
        tmp_path / "step8_UeccDNA_C.fasta",
        tmp_path / "step8_MeccSites.core.csv",
        tmp_path / "step8_MeccSites.bed",
        tmp_path / "step8_MeccBestSite.bed",
        tmp_path / "step8_MeccDNA_C.fasta",
        tmp_path / "step8_CeccSegments.core.csv",
        tmp_path / "step8_CeccSegments.bed",
        tmp_path / "step8_CeccJunctions.bedpe",
        tmp_path / "step8_CeccDNA_C.fasta",
    ]
    for path in expected_outputs:
        assert path.exists(), f"Expected output missing: {path}"


def test_ecc_dedup_csv_cluster_merging(tmp_path):
    """Ensure CSV cluster mappings are honoured when selecting representatives."""
    cluster_csv = tmp_path / "uecc.id2cluster.csv"
    cluster_csv.write_text(
        "id,cluster,is_representative\n"
        "U1,grp1,false\n"
        "U2,grp1,true\n"
        "U3,grp2,\n"
    )

    clusters = CDHitClusters()
    clusters.parse(cluster_csv)

    assert clusters.rep_of["grp1"] == "U2"
    assert clusters.rep_of["grp2"] == "U3"

    uecc_df = pd.DataFrame(
        [
            {
                "eccDNA_id": "U1",
                "eChr": "chr1",
                "eStart": 10,
                "eEnd": 70,
                "eLength": 60,
                "readName": "readA",
                "eRepeatNum": 1,
                "eStrand": "+",
            },
            {
                "eccDNA_id": "U2",
                "eChr": "chr1",
                "eStart": 20,
                "eEnd": 90,
                "eLength": 70,
                "readName": "readB",
                "eRepeatNum": 2,
                "eStrand": "+",
            },
            {
                "eccDNA_id": "U3",
                "eChr": "chr2",
                "eStart": 5,
                "eEnd": 45,
                "eLength": 40,
                "readName": "readC",
                "eRepeatNum": 3,
                "eStrand": "-",
            },
        ]
    )

    processor = eccDedup()
    processed = processor.process_uecc(uecc_df, clusters, "Uecc")
    processed = processed.sort_values("cluster_id").reset_index(drop=True)

    assert len(processed) == 2

    assert {"eccdna_id", "chr", "start_0based", "end_0based", "strand", "copy_number", "repeat_number", "match_degree"}.issubset(processed.columns)

    cluster1 = processed.loc[processed["cluster_id"] == "grp1"].iloc[0]
    assert cluster1["eccdna_id"] == "U2"
    assert cluster1["reads"] == "readA;readB"
    assert cluster1["repeat_number"] == 3
    assert cluster1["copy_number"] == 3
    assert cluster1["num_merged"] == 2
    assert cluster1["merged_from_ids"] == "U1;U2"
    assert cluster1["orig_eccdna_id"] == "U2"
    assert cluster1["eccdna_type"] == "Uecc"
    assert cluster1["strand"] == "+"

    cluster2 = processed.loc[processed["cluster_id"] == "grp2"].iloc[0]
    assert cluster2["eccdna_id"] == "U3"
    assert cluster2["num_merged"] == 1
    assert cluster2["merged_from_ids"] == "U3"
    assert cluster2["strand"] == "-"
    assert cluster2["copy_number"] == 3
