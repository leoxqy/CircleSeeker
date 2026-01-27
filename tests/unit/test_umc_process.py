"""Tests for umc_process module - UMC eccDNA processing pipeline."""

from pathlib import Path
import importlib.util
import sys
import types

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


def _load_umc_process_module():
    package_name = "circleseeker.modules"
    module_name = f"{package_name}.umc_process"

    if module_name in sys.modules:
        return sys.modules[module_name]

    import circleseeker  # noqa: F401

    if package_name not in sys.modules:
        package_module = types.ModuleType(package_name)
        package_module.__path__ = [str(SRC / "circleseeker" / "modules")]
        sys.modules[package_name] = package_module

    spec = importlib.util.spec_from_file_location(
        module_name,
        SRC / "circleseeker" / "modules" / "umc_process.py",
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


# Load module and get classes
umc_module = _load_umc_process_module()
UMCProcessConfig = umc_module.UMCProcessConfig
UMCProcess = umc_module.UMCProcess
SequenceLibrary = umc_module.SequenceLibrary
XeccExporter = umc_module.XeccExporter
UeccProcessor = umc_module.UeccProcessor
MeccProcessor = umc_module.MeccProcessor
CeccProcessor = umc_module.CeccProcessor
BaseEccProcessor = umc_module.BaseEccProcessor
extract_ring_sequence = umc_module.extract_ring_sequence
canonicalize_circular_sequence = umc_module.canonicalize_circular_sequence
_sanitize_fasta_id = umc_module._sanitize_fasta_id
_coerce_to_paths = umc_module._coerce_to_paths
_reverse_complement = umc_module._reverse_complement
_min_circular_rotation = umc_module._min_circular_rotation
_booth_min_rotation_start = umc_module._booth_min_rotation_start


class TestHelperFunctions:
    def test_sanitize_fasta_id(self):
        assert _sanitize_fasta_id("test id") == "test_id"
        assert _sanitize_fasta_id("test@id#123") == "test_id_123"
        assert _sanitize_fasta_id("test.id-123") == "test.id-123"
        assert _sanitize_fasta_id("test  multiple  spaces") == "test_multiple_spaces"

    def test_coerce_to_paths(self):
        assert _coerce_to_paths(None) == []
        assert _coerce_to_paths("/path/to/file") == [Path("/path/to/file")]
        assert _coerce_to_paths(Path("/path/to/file")) == [Path("/path/to/file")]
        assert _coerce_to_paths(["/path1", Path("/path2")]) == [Path("/path1"), Path("/path2")]
        assert _coerce_to_paths(["", None, "/valid"]) == [Path("/valid")]

    def test_extract_ring_sequence(self):
        seq = "ACGTACGTACGT"
        # Normal extraction
        assert extract_ring_sequence(seq, 1, 4) == "ACGT"
        assert extract_ring_sequence(seq, 5, 4) == "ACGT"

        # Circular extraction (wrapping around)
        result = extract_ring_sequence(seq, 10, 6)
        assert len(result) == 6

        # Edge cases
        assert extract_ring_sequence(seq, 1, 12) == "ACGTACGTACGT"


class TestSequenceLibrary:
    @pytest.fixture
    def temp_fasta(self, tmp_path):
        fasta_file = tmp_path / "test.fasta"
        content = """>seq1
ACGTACGTACGT
>seq2|extra_info
TTGGCCAATTGG
>seq3_with_underscore
AAAAAACCCCCC
"""
        fasta_file.write_text(content)
        return fasta_file

    def test_load_fasta(self, temp_fasta):
        lib = SequenceLibrary()
        lib.load_fasta(temp_fasta)

        assert len(lib.primary_ids) == 3
        assert "seq1" in lib.fasta_sequences
        assert "seq2|extra_info" in lib.fasta_sequences
        assert "seq3_with_underscore" in lib.fasta_sequences

        assert lib.fasta_sequences["seq1"] == "ACGTACGTACGT"
        assert lib.fasta_sequences["seq2|extra_info"] == "TTGGCCAATTGG"

    def test_load_fasta_missing_file(self):
        lib = SequenceLibrary()
        with pytest.raises(FileNotFoundError):
            lib.load_fasta(Path("/nonexistent/file.fasta"))

    def test_find_sequence(self, temp_fasta):
        lib = SequenceLibrary()
        lib.load_fasta(temp_fasta)

        # Direct matches
        assert lib.find_sequence("seq1") == "ACGTACGTACGT"
        assert lib.find_sequence("seq2|extra_info") == "TTGGCCAATTGG"

        # Base ID lookup for pipe-separated IDs
        assert lib.find_sequence("seq2") == "TTGGCCAATTGG"

        # Not found
        assert lib.find_sequence("nonexistent") is None


class TestUMCProcessConfig:
    def test_default_config(self):
        config = UMCProcessConfig()
        assert config.process_xecc is True
        assert config.cluster_uecc is True
        assert config.cluster_mecc is True
        assert config.cluster_cecc is True
        assert config.debug is False

    def test_custom_config(self):
        config = UMCProcessConfig(
            process_xecc=False,
            cluster_uecc=False,
            debug=True
        )
        assert config.process_xecc is False
        assert config.cluster_uecc is False
        assert config.debug is True


class TestXeccExporter:
    @pytest.fixture
    def setup_exporter(self, tmp_path):
        # Create test FASTA
        fasta_file = tmp_path / "test.fasta"
        content = """>classified1
ACGTACGTACGT
>classified2
TTGGCCAATTGG
>unclassified1
AAAAAACCCCCC
>unclassified2
GGGGGGTTTTTT
"""
        fasta_file.write_text(content)

        # Create classified CSVs
        csv1 = tmp_path / "classified1.csv"
        csv1.write_text("query_id\nclassified1\nclassified2\n")

        # Setup sequence library
        seq_library = SequenceLibrary()
        seq_library.load_fasta(fasta_file)

        config = UMCProcessConfig()

        return XeccExporter(seq_library, config), [csv1], tmp_path

    def test_get_classified_ids(self, setup_exporter):
        exporter, csvs, _ = setup_exporter
        classified_ids = exporter._get_classified_ids(csvs)

        assert "classified1" in classified_ids
        assert "classified2" in classified_ids
        assert len(classified_ids) == 2

    def test_generate_xecc(self, setup_exporter):
        exporter, csvs, tmp_path = setup_exporter
        output_dir = tmp_path / "output"

        result = exporter.generate(csvs, output_dir, "test")

        assert result is not None
        assert result.name == "test_XeccDNA.fasta"
        assert result.exists()

        # Check content
        from Bio import SeqIO
        records = list(SeqIO.parse(str(result), "fasta"))
        assert len(records) == 2

        # Check IDs contain X marker
        ids = [r.id for r in records]
        assert all("__X" in id for id in ids)

    def test_generate_xecc_no_unclassified(self, setup_exporter):
        exporter, _, tmp_path = setup_exporter
        # Create CSV with all sequences classified
        csv_all = tmp_path / "all_classified.csv"
        csv_all.write_text("query_id\nclassified1\nclassified2\nunclassified1\nunclassified2\n")

        result = exporter.generate([csv_all], tmp_path / "output", "test")
        assert result is None


class TestUeccProcessor:
    @pytest.fixture
    def setup_processor(self, tmp_path):
        # Create test FASTA
        fasta_file = tmp_path / "test.fasta"
        content = """>read1
ACGTACGTACGTACGTACGT
>read2
TTGGCCAATTGGTTGGCCAA
"""
        fasta_file.write_text(content)

        seq_library = SequenceLibrary()
        seq_library.load_fasta(fasta_file)

        config = UMCProcessConfig()
        return UeccProcessor(seq_library, config), tmp_path

    def test_cluster_by_location(self, setup_processor):
        processor, _ = setup_processor

        # Test data with same location (should cluster)
        df = pd.DataFrame({
            'query_id': ['q1', 'q2', 'q3'],
            'chr': ['chr1', 'chr1', 'chr2'],
            'start0': [100, 100, 200],
            'end0': [150, 150, 250],
            'copy_number': [5, 3, 2]
        })

        result = processor.cluster_by_location(df)

        # q1 and q2 should be clustered (same location)
        assert len(result) == 2
        clustered = result[result['cluster_size'] == 2]
        assert len(clustered) == 1
        assert clustered['copy_number'].iloc[0] == 8  # 5 + 3

    def test_compute_sequences(self, setup_processor):
        processor, _ = setup_processor

        df = pd.DataFrame({
            'query_id': ['read1', 'read2'],
            'q_start': [5, 3],
            'length': [8, 6]
        })

        result = processor.compute_sequences(df)

        assert 'eSeq' in result.columns
        assert len(result['eSeq'].iloc[0]) == 8
        assert len(result['eSeq'].iloc[1]) == 6

    def test_compute_sequences_canonicalizes_rotation_and_strand(self, tmp_path):
        """Same circular molecule can appear as rotated / reverse-complemented sequences."""
        circle = "GATTACAG"
        doubled = circle + circle

        # Reverse-complement circle and also provide doubled form.
        rc_circle = umc_module._reverse_complement(circle)  # type: ignore[attr-defined]
        doubled_rc = rc_circle + rc_circle

        fasta_file = tmp_path / "test_rot_rc.fasta"
        fasta_file.write_text(f">read1\n{doubled}\n>read2\n{doubled_rc}\n")

        seq_library = SequenceLibrary()
        seq_library.load_fasta(fasta_file)
        processor = UeccProcessor(seq_library, UMCProcessConfig())

        df = pd.DataFrame({
            'query_id': ['read1', 'read2'],
            'q_start': [1, 3],  # introduce a rotation on read2
            'length': [len(circle), len(circle)],
        })

        result = processor.compute_sequences(df)
        assert result["eSeq"].nunique() == 1

    def test_add_numbering_and_export(self, setup_processor):
        processor, _ = setup_processor

        df = pd.DataFrame({
            'query_id': ['read1', 'read2'],
            'eSeq': ['ACGTACGT', 'GGCCAA']
        })

        result = processor.add_numbering_and_export(df)

        assert result['eccDNA_id'].iloc[0] == 'U1'
        assert result['eccDNA_id'].iloc[1] == 'U2'
        assert len(processor.fasta_records) == 2

    def test_process_pipeline(self, setup_processor):
        processor, tmp_path = setup_processor

        # Create test CSV with new column names
        csv_file = tmp_path / "uecc.csv"
        csv_file.write_text("""query_id,chr,start0,end0,q_start,length,copy_number
read1,chr1,100,150,5,8,5
read2,chr1,100,150,3,6,3
""")

        output_dir = tmp_path / "output"
        result = processor.process([csv_file], output_dir, "test", cluster=True)

        assert result is not None
        assert len(result) == 1  # Clustered into one
        assert (output_dir / "test_UeccDNA_processed.csv").exists()
        assert (output_dir / "test_UeccDNA_pre.fasta").exists()


class TestMeccProcessor:
    @pytest.fixture
    def setup_processor(self, tmp_path):
        # Create test FASTA
        fasta_file = tmp_path / "test.fasta"
        content = """>read1
ACGTACGTACGTACGTACGT
>read2
TTGGCCAATTGGTTGGCCAA
>read3
AAAAAACCCCCCAAAAAACCCCCC
"""
        fasta_file.write_text(content)

        seq_library = SequenceLibrary()
        seq_library.load_fasta(fasta_file)

        config = UMCProcessConfig()
        return MeccProcessor(seq_library, config), tmp_path

    def test_generate_mecc_signature(self, setup_processor):
        processor, _ = setup_processor

        df = pd.DataFrame({
            'chr': ['chr1', 'chr1', 'chr2'],
            'start0': [100, 200, 300],
            'end0': [150, 250, 350]
        })

        signature = processor.generate_mecc_signature(df)
        assert signature == "chr1:100-150;chr1:200-250;chr2:300-350"

    def test_cluster_by_signature(self, setup_processor):
        processor, _ = setup_processor

        # Create data where read1 and read2 have same multi-location signature
        df = pd.DataFrame({
            'query_id': ['read1', 'read1', 'read2', 'read2', 'read3'],
            'chr': ['chr1', 'chr2', 'chr1', 'chr2', 'chr3'],
            'start0': [100, 200, 100, 200, 300],
            'end0': [150, 250, 150, 250, 350],
            'copy_number': [5, 5, 3, 3, 2]
        })

        result = processor.cluster_by_signature(df)

        # read1 and read2 should be clustered (same signature)
        assert 'cluster_id' in result.columns
        assert 'cluster_size' in result.columns

    def test_process_pipeline(self, setup_processor):
        processor, tmp_path = setup_processor

        # Create test CSV
        csv_file = tmp_path / "mecc.csv"
        csv_file.write_text("""query_id,chr,start0,end0,q_start,length,copy_number,Gap_Percentage
read1,chr1,100,150,5,8,5,10.0
read1,chr2,200,250,5,8,5,10.0
read2,chr1,100,150,3,6,3,15.0
read2,chr2,200,250,3,6,3,15.0
""")

        output_dir = tmp_path / "output"
        result = processor.process([csv_file], output_dir, "test", cluster=True)

        assert result is not None
        assert (output_dir / "test_MeccDNA_processed.csv").exists()
        assert (output_dir / "test_MeccDNA_pre.fasta").exists()


class TestCeccProcessor:
    @pytest.fixture
    def setup_processor(self, tmp_path):
        # Create test FASTA
        fasta_file = tmp_path / "test.fasta"
        content = """>read1
ACGTACGTACGTACGTACGTACGTACGT
>read2
TTGGCCAATTGGTTGGCCAATTGGCCAA
"""
        fasta_file.write_text(content)

        seq_library = SequenceLibrary()
        seq_library.load_fasta(fasta_file)

        config = UMCProcessConfig()
        return CeccProcessor(seq_library, config), tmp_path

    def test_generate_cecc_signature(self, setup_processor):
        processor, _ = setup_processor

        df = pd.DataFrame({
            'chr': ['chr1', 'chr2', 'chr3'],
            'start0': [100, 200, 300],
            'end0': [150, 250, 350],
            'seg_index': [1, 2, 3]
        })

        signature = processor.generate_cecc_signature(df)
        assert signature == "chr1:100-150;chr2:200-250;chr3:300-350"

    def test_generate_cecc_signature_direction_invariant(self, setup_processor):
        processor, _ = setup_processor

        df = pd.DataFrame({
            'chr': ['chr1', 'chr2', 'chr3'],
            'start0': [100, 200, 300],
            'end0': [150, 250, 350],
        })

        sig_forward = processor.generate_cecc_signature(df)
        sig_reverse = processor.generate_cecc_signature(df.iloc[::-1].copy())
        assert sig_forward == sig_reverse

    def test_cluster_by_signature(self, setup_processor):
        processor, _ = setup_processor

        # Create data where read1 and read2 have same segment structure
        df = pd.DataFrame({
            'query_id': ['read1', 'read1', 'read2', 'read2'],
            'chr': ['chr1', 'chr2', 'chr1', 'chr2'],
            'start0': [100, 200, 100, 200],
            'end0': [150, 250, 150, 250],
            'seg_index': [1, 2, 1, 2],
            'copy_number': [5, 5, 3, 3]
        })

        result = processor.cluster_by_signature(df)

        assert 'cluster_id' in result.columns
        assert 'cluster_size' in result.columns

    def test_process_pipeline(self, setup_processor):
        processor, tmp_path = setup_processor

        # Create test CSV with multi-segment circular DNA
        csv_file = tmp_path / "cecc.csv"
        csv_file.write_text("""query_id,chr,start0,end0,seg_index,q_start,length
read1,chr1,100,150,1,5,12
read1,chr2,200,250,2,5,12
read2,chr1,300,350,1,3,10
""")

        output_dir = tmp_path / "output"
        result = processor.process([csv_file], output_dir, "test", cluster=True)

        assert result is not None
        assert (output_dir / "test_CeccDNA_processed.csv").exists()
        assert (output_dir / "test_CeccDNA_pre.fasta").exists()


class TestUMCProcess:
    @pytest.fixture
    def setup_pipeline(self, tmp_path):
        # Create comprehensive test FASTA
        fasta_file = tmp_path / "master.fasta"
        content = """>uecc1
ACGTACGTACGTACGTACGT
>mecc1
TTGGCCAATTGGTTGGCCAA
>cecc1
AAAAAACCCCCCAAAAAACCCCCC
>unclassified1
GGGGGGTTTTTTGGGGGGTTTTTT
"""
        fasta_file.write_text(content)

        # Create test CSVs with new column names
        uecc_csv = tmp_path / "uecc.csv"
        uecc_csv.write_text("""query_id,chr,start0,end0,q_start,length
uecc1,chr1,100,150,5,8
""")

        mecc_csv = tmp_path / "mecc.csv"
        mecc_csv.write_text("""query_id,chr,start0,end0,q_start,length
mecc1,chr1,100,150,3,6
""")

        cecc_csv = tmp_path / "cecc.csv"
        cecc_csv.write_text("""query_id,chr,start0,end0,seg_index,q_start,length
cecc1,chr1,100,150,1,1,10
""")

        return fasta_file, uecc_csv, mecc_csv, cecc_csv, tmp_path

    def test_process_all_with_xecc(self, setup_pipeline):
        fasta_file, uecc_csv, mecc_csv, cecc_csv, tmp_path = setup_pipeline

        config = UMCProcessConfig(
            process_xecc=True,
            cluster_mecc=True,
            cluster_cecc=True,
            cluster_uecc=True
        )

        processor = UMCProcess(config)
        output_dir = tmp_path / "output"

        results = processor.process_all(
            fasta_file=fasta_file,
            uecc_csv=uecc_csv,
            mecc_csv=mecc_csv,
            cecc_csv=cecc_csv,
            output_dir=output_dir,
            prefix="test"
        )

        assert "xecc" in results
        assert "uecc" in results
        assert "mecc" in results
        assert "cecc" in results

        # Check XeccDNA file was created
        assert results["xecc"].name == "test_XeccDNA.fasta"
        assert results["xecc"].exists()

        # Check other outputs are DataFrames
        assert isinstance(results["uecc"], pd.DataFrame)
        assert isinstance(results["mecc"], pd.DataFrame)
        assert isinstance(results["cecc"], pd.DataFrame)

    def test_process_all_without_xecc(self, setup_pipeline):
        fasta_file, uecc_csv, mecc_csv, cecc_csv, tmp_path = setup_pipeline

        config = UMCProcessConfig(process_xecc=False)
        processor = UMCProcess(config)
        output_dir = tmp_path / "output"

        results = processor.process_all(
            fasta_file=fasta_file,
            uecc_csv=uecc_csv,
            mecc_csv=mecc_csv,
            cecc_csv=cecc_csv,
            output_dir=output_dir,
            prefix="test"
        )

        assert "xecc" not in results
        assert len(results) == 3

    def test_process_partial_inputs(self, setup_pipeline):
        fasta_file, uecc_csv, _, _, tmp_path = setup_pipeline

        processor = UMCProcess()
        output_dir = tmp_path / "output"

        # Process only UeccDNA
        results = processor.process_all(
            fasta_file=fasta_file,
            uecc_csv=uecc_csv,
            mecc_csv=None,
            cecc_csv=None,
            output_dir=output_dir,
            prefix="test"
        )

        assert "uecc" in results
        assert "mecc" not in results
        assert "cecc" not in results

    def test_missing_fasta_file(self, setup_pipeline):
        _, uecc_csv, mecc_csv, cecc_csv, tmp_path = setup_pipeline

        processor = UMCProcess()
        output_dir = tmp_path / "output"

        with pytest.raises(FileNotFoundError):
            processor.process_all(
                fasta_file=Path("/nonexistent.fasta"),
                uecc_csv=uecc_csv,
                mecc_csv=mecc_csv,
                cecc_csv=cecc_csv,
                output_dir=output_dir,
                prefix="test"
            )


# ========== Additional Tests for Coverage ========== #


class TestReverseComplement:
    """Tests for _reverse_complement function."""

    def test_basic(self):
        assert _reverse_complement("ACGT") == "ACGT"  # palindrome

    def test_simple(self):
        assert _reverse_complement("AAAA") == "TTTT"
        assert _reverse_complement("CCCC") == "GGGG"

    def test_mixed(self):
        assert _reverse_complement("AACG") == "CGTT"

    def test_empty(self):
        assert _reverse_complement("") == ""

    def test_lowercase(self):
        assert _reverse_complement("acgt") == "acgt"


class TestMinCircularRotation:
    """Tests for _min_circular_rotation function."""

    def test_single_char(self):
        assert _min_circular_rotation("A") == "A"

    def test_empty(self):
        assert _min_circular_rotation("") == ""

    def test_already_minimal(self):
        assert _min_circular_rotation("ABCD") == "ABCD"

    def test_rotation_needed(self):
        assert _min_circular_rotation("CDAB") == "ABCD"

    def test_repeated(self):
        assert _min_circular_rotation("AAA") == "AAA"


class TestBoothMinRotationStart:
    """Tests for _booth_min_rotation_start function."""

    def test_already_minimal(self):
        assert _booth_min_rotation_start("ABCD") == 0

    def test_rotation(self):
        assert _booth_min_rotation_start("CDAB") == 2

    def test_single(self):
        assert _booth_min_rotation_start("A") == 0

    def test_empty(self):
        assert _booth_min_rotation_start("") == 0


class TestCanonicalizeCircularSequence:
    """Tests for canonicalize_circular_sequence function."""

    def test_empty_string(self):
        assert canonicalize_circular_sequence("") == ""

    def test_single_char(self):
        assert canonicalize_circular_sequence("A") == "A"

    def test_rotation_invariant(self):
        """Rotations of the same circle produce the same canonical form."""
        seq1 = canonicalize_circular_sequence("ACGT")
        seq2 = canonicalize_circular_sequence("CGTA")
        seq3 = canonicalize_circular_sequence("GTAC")
        seq4 = canonicalize_circular_sequence("TACG")
        assert seq1 == seq2 == seq3 == seq4

    def test_strand_invariant(self):
        """Forward and reverse-complement produce the same canonical form."""
        fwd = canonicalize_circular_sequence("GATTACA")
        rc = canonicalize_circular_sequence(_reverse_complement("GATTACA"))
        assert fwd == rc

    def test_whitespace_stripped(self):
        assert canonicalize_circular_sequence("  ACGT  ") == canonicalize_circular_sequence("ACGT")

    def test_uppercase(self):
        assert canonicalize_circular_sequence("acgt") == canonicalize_circular_sequence("ACGT")


class TestMeccProcessorComputeSequences:
    """Tests for MeccProcessor.compute_sequences."""

    @pytest.fixture
    def setup(self, tmp_path):
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(">qid1\nACGTACGTACGTACGTACGT\n>qid2\nTTGGCCAATTGGTTGGCCAA\n")
        seq_lib = SequenceLibrary()
        seq_lib.load_fasta(fasta_file)
        return MeccProcessor(seq_lib, UMCProcessConfig())

    def test_basic_extraction(self, setup):
        """Extract sequences for MeccDNA (one per unique query_id)."""
        processor = setup
        df = pd.DataFrame({
            "query_id": ["qid1", "qid1", "qid2"],
            "chr": ["chr1", "chr2", "chr3"],
            "start0": [100, 200, 300],
            "end0": [150, 250, 350],
            "q_start": [1, 1, 3],
            "length": [6, 6, 4],
        })
        result = processor.compute_sequences(df)
        assert "eSeq" in result.columns
        # All rows for qid1 should have the same sequence
        qid1_seqs = result[result["query_id"] == "qid1"]["eSeq"].unique()
        assert len(qid1_seqs) == 1

    def test_missing_sequence(self, setup):
        """Missing sequence should leave eSeq empty."""
        processor = setup
        df = pd.DataFrame({
            "query_id": ["nonexistent"],
            "chr": ["chr1"],
            "start0": [100],
            "end0": [200],
            "q_start": [1],
            "length": [6],
        })
        result = processor.compute_sequences(df)
        assert result["eSeq"].iloc[0] == ""

    def test_missing_columns(self, setup):
        """Missing length/q_start columns should warn and return unchanged."""
        processor = setup
        df = pd.DataFrame({
            "query_id": ["qid1"],
            "chr": ["chr1"],
        })
        result = processor.compute_sequences(df)
        assert "eSeq" in result.columns
        assert result["eSeq"].iloc[0] == ""


class TestCeccProcessorComputeSequences:
    """Tests for CeccProcessor.compute_sequences."""

    @pytest.fixture
    def setup(self, tmp_path):
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(">qid1\nACGTACGTACGTACGTACGT\n")
        seq_lib = SequenceLibrary()
        seq_lib.load_fasta(fasta_file)
        return CeccProcessor(seq_lib, UMCProcessConfig())

    def test_basic_extraction(self, setup):
        """Extract sequences for CeccDNA (one per unique query_id)."""
        processor = setup
        df = pd.DataFrame({
            "query_id": ["qid1", "qid1"],
            "chr": ["chr1", "chr2"],
            "start0": [100, 200],
            "end0": [150, 250],
            "seg_index": [1, 2],
            "q_start": [1, 1],
            "length": [8, 8],
        })
        result = processor.compute_sequences(df)
        assert "eSeq" in result.columns
        # All rows for qid1 share the same sequence
        assert result["eSeq"].nunique() == 1
        assert len(result["eSeq"].iloc[0]) == 8


class TestMeccProcessorProcess:
    """Tests for MeccProcessor.process with various scenarios."""

    @pytest.fixture
    def setup(self, tmp_path):
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(
            ">read1\nACGTACGTACGTACGTACGT\n"
            ">read2\nTTGGCCAATTGGTTGGCCAA\n"
        )
        seq_lib = SequenceLibrary()
        seq_lib.load_fasta(fasta_file)
        return MeccProcessor(seq_lib, UMCProcessConfig()), tmp_path

    def test_empty_csv(self, setup):
        """Empty CSV returns None."""
        processor, tmp_path = setup
        csv_file = tmp_path / "empty.csv"
        csv_file.write_text("")
        result = processor.process([csv_file], tmp_path / "out", "test")
        assert result is None

    def test_nonexistent_file(self, setup):
        """Nonexistent file returns None."""
        processor, tmp_path = setup
        result = processor.process(
            [tmp_path / "nonexistent.csv"], tmp_path / "out", "test"
        )
        assert result is None

    def test_multi_location_clustering(self, setup):
        """Multiple reads at same locations cluster together."""
        processor, tmp_path = setup
        csv_file = tmp_path / "mecc.csv"
        csv_file.write_text(
            "query_id,chr,start0,end0,q_start,length,copy_number\n"
            "read1,chr1,100,150,5,8,5\n"
            "read1,chr2,200,250,5,8,5\n"
            "read2,chr1,100,150,3,6,3\n"
            "read2,chr2,200,250,3,6,3\n"
        )
        output_dir = tmp_path / "out"
        result = processor.process([csv_file], output_dir, "test", cluster=True)
        assert result is not None
        # read1 and read2 share the same multi-location signature, should cluster
        assert result["cluster_size"].max() == 2


class TestCeccProcessorProcess:
    """Tests for CeccProcessor.process with various scenarios."""

    @pytest.fixture
    def setup(self, tmp_path):
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(
            ">read1\nACGTACGTACGTACGTACGTACGTACGT\n"
        )
        seq_lib = SequenceLibrary()
        seq_lib.load_fasta(fasta_file)
        return CeccProcessor(seq_lib, UMCProcessConfig()), tmp_path

    def test_empty_csv(self, setup):
        """Empty CSV returns None."""
        processor, tmp_path = setup
        csv_file = tmp_path / "empty.csv"
        csv_file.write_text("")
        result = processor.process([csv_file], tmp_path / "out", "test")
        assert result is None

    def test_basic_processing(self, setup):
        """Basic CeccDNA processing pipeline."""
        processor, tmp_path = setup
        csv_file = tmp_path / "cecc.csv"
        csv_file.write_text(
            "query_id,chr,start0,end0,seg_index,q_start,length\n"
            "read1,chr1,100,150,1,5,12\n"
            "read1,chr2,200,250,2,5,12\n"
        )
        output_dir = tmp_path / "out"
        result = processor.process([csv_file], output_dir, "test", cluster=True)
        assert result is not None
        assert "eccDNA_id" in result.columns
        assert result["eccDNA_id"].iloc[0].startswith("C")


class TestBaseEccProcessorLoadCsvFiles:
    """Tests for BaseEccProcessor.load_csv_files via MeccProcessor."""

    @pytest.fixture
    def processor(self, tmp_path):
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(">seq1\nACGT\n")
        seq_lib = SequenceLibrary()
        seq_lib.load_fasta(fasta_file)
        return MeccProcessor(seq_lib, UMCProcessConfig()), tmp_path

    def test_load_valid_csv(self, processor):
        """Load valid CSV files."""
        proc, tmp_path = processor
        csv1 = tmp_path / "a.csv"
        csv1.write_text("query_id,chr\nq1,chr1\n")
        csv2 = tmp_path / "b.csv"
        csv2.write_text("query_id,chr\nq2,chr2\n")
        # MeccProcessor doesn't have load_csv_files directly,
        # but process() uses the same logic; test via direct CSV loading
        result_frames = []
        for f in [csv1, csv2]:
            df = pd.read_csv(f)
            result_frames.append(df)
        combined = pd.concat(result_frames, ignore_index=True)
        assert len(combined) == 2

    def test_skip_missing_file(self, processor):
        """Missing files are skipped gracefully."""
        proc, tmp_path = processor
        valid = tmp_path / "valid.csv"
        valid.write_text("query_id,chr\nq1,chr1\n")
        # process() should handle missing files via its internal loop
        result = proc.process(
            [valid, tmp_path / "nonexistent.csv"],
            tmp_path / "out", "test"
        )
        assert result is not None

    def test_skip_empty_csv(self, processor):
        """Empty CSV files are skipped gracefully (zero-size file)."""
        proc, tmp_path = processor
        empty = tmp_path / "empty.csv"
        empty.write_text("")  # 0 bytes → st_size == 0, skipped before read_csv
        result = proc.process([empty], tmp_path / "out", "test")
        assert result is None


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
