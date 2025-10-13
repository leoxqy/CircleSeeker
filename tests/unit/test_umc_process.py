from pathlib import Path
import importlib.util
import sys
import types
import tempfile
import pandas as pd
from unittest.mock import Mock, patch
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

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

    import circleseeker  # noqa: F401  # Ensure base package is initialized

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
extract_ring_sequence = umc_module.extract_ring_sequence
_sanitize_fasta_id = umc_module._sanitize_fasta_id
_coerce_to_paths = umc_module._coerce_to_paths


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
        
        # Circular extraction
        assert extract_ring_sequence(seq, 10, 6) == "TACGTA"
        assert extract_ring_sequence(seq, 11, 4) == "GTAC"
        
        # Edge cases
        assert extract_ring_sequence(seq, 1, 12) == "ACGTACGTACGT"
        assert extract_ring_sequence(seq, 7, 8) == "GTACGTAC"


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
        
        assert len(lib.fasta_sequences) == 4  # 3 primary + 1 base ID for seq2
        assert "seq1" in lib.fasta_sequences
        assert "seq2|extra_info" in lib.fasta_sequences
        assert "seq2" in lib.fasta_sequences  # Base ID also stored
        assert "seq3_with_underscore" in lib.fasta_sequences
        
        assert lib.fasta_sequences["seq1"] == "ACGTACGTACGT"
        assert lib.fasta_sequences["seq2|extra_info"] == "TTGGCCAATTGG"
        assert lib.fasta_sequences["seq2"] == "TTGGCCAATTGG"
    
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
        assert lib.find_sequence("seq2") == "TTGGCCAATTGG"
        
        # Partial matches
        assert lib.find_sequence("seq3") == "AAAAAACCCCCC"
        
        # Not found
        assert lib.find_sequence("nonexistent") is None


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
        
        config = UMCProcessConfig(
            xecc_prefix_queryid=True,
            xecc_id_separator="__"
        )
        
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
        
        # Check IDs and sequences
        ids = [r.id for r in records]
        assert all("__X" in id for id in ids)
        assert all(len(r.seq) == 6 for r in records)  # Half length
    
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
    
    def test_aggregate_by_ename(self, setup_processor):
        processor, _ = setup_processor
        
        # Test data with duplicates
        df = pd.DataFrame({
            'query_id': ['q1', 'q2', 'q3', 'q4', 'q5'],
            'eName': ['repeat1', 'repeat1', 'repeat2', None, 'repeat3'],
            'readName': ['r1;r2', 'r3', 'r4', 'r5', 'r6'],
            'eRepeatNum': [5, 3, 2, 1, 4]
        })
        
        result = processor.aggregate_by_ename(df)
        
        assert len(result) == 4  # 3 unique + 1 non-aggregatable
        assert result[result['eName'] == 'repeat1']['eRepeatNum'].iloc[0] == 8
        assert 'r1' in result[result['eName'] == 'repeat1']['readName'].iloc[0]
        assert 'r3' in result[result['eName'] == 'repeat1']['readName'].iloc[0]
    
    def test_compute_sequences(self, setup_processor):
        processor, _ = setup_processor
        
        df = pd.DataFrame({
            'query_id': ['read1', 'read2'],
            'q_start': [5, 3],
            'eLength': [8, 6]
        })
        
        result = processor.compute_sequences(df)
        
        assert result['eSeq'].iloc[0] == 'ACGTACGT'
        assert result['eSeq'].iloc[1] == 'GGCCAA'
    
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
    
    def test_process_full_pipeline(self, setup_processor):
        processor, tmp_path = setup_processor
        
        # Create test CSV
        csv_file = tmp_path / "uecc.csv"
        csv_file.write_text("""query_id,eName,q_start,eLength,readName,eRepeatNum
read1,repeat1,5,8,r1;r2,5
read2,repeat1,3,6,r3,3
""")
        
        output_dir = tmp_path / "output"
        result = processor.process([csv_file], output_dir, "test", aggregate=True)
        
        assert result is not None
        assert len(result) == 1  # Aggregated
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
            'eChr': ['chr1', 'chr1', 'chr2'],
            'eStart': [100, 200, 300],
            'eEnd': [150, 250, 350]
        })
        
        signature = processor.generate_mecc_signature(df)
        assert signature == "chr1:100-150;chr1:200-250;chr2:300-350"
    
    def test_cluster_by_signature(self, setup_processor):
        processor, _ = setup_processor
        
        df = pd.DataFrame({
            'query_id': ['q1', 'q2', 'q3', 'q4'],
            'eChr': ['chr1', 'chr1', 'chr2', 'chr1'],
            'eStart': [100, 100, 200, 100],
            'eEnd': [150, 150, 250, 150],
            'copyNum': [5, 3, 2, 7]
        })
        
        result = processor.cluster_by_signature(df)
        
        # q1, q2, q4 should be clustered (same location)
        assert len(result) == 2  # One cluster + one unclustered
        clustered = result[result['cluster_size'] == 3]
        assert len(clustered) == 1
        assert clustered['copyNum'].iloc[0] == 15  # 5 + 3 + 7
    
    def test_calculate_maturity_degree(self, setup_processor):
        processor, _ = setup_processor
        
        df = pd.DataFrame({
            'Gap_Percentage': [10.5, 25.0, None, 0.0]
        })
        
        result = processor.calculate_maturity_degree(df)
        
        assert result['MatDegree'].iloc[0] == 89.5
        assert result['MatDegree'].iloc[1] == 75.0
        assert pd.isna(result['MatDegree'].iloc[2])
        assert result['MatDegree'].iloc[3] == 100.0
    
    def test_process_full_pipeline(self, setup_processor):
        processor, tmp_path = setup_processor
        
        # Create test CSV
        csv_file = tmp_path / "mecc.csv"
        csv_file.write_text("""query_id,eChr,eStart,eEnd,q_start,eLength,copyNum,Gap_Percentage
read1,chr1,100,150,5,8,5,10.0
read2,chr1,100,150,3,6,3,15.0
read3,chr2,200,250,1,10,2,5.0
""")
        
        output_dir = tmp_path / "output"
        result = processor.process([csv_file], output_dir, "test", cluster=True)
        
        assert result is not None
        assert len(result) == 2  # Two clusters
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
            'start': [100, 200, 300],
            'end': [150, 250, 350],
            'segment_in_circle': [1, 2, 3]
        })
        
        signature = processor.generate_cecc_signature(df)
        assert signature == "chr1:100-150;chr2:200-250;chr3:300-350"
    
    def test_cluster_by_signature_with_segments(self, setup_processor):
        processor, _ = setup_processor
        
        # Complex segments for two circular DNAs
        df = pd.DataFrame({
            'query_id': ['q1', 'q1', 'q2', 'q2', 'q3'],
            'chr': ['chr1', 'chr2', 'chr1', 'chr2', 'chr3'],
            'start': [100, 200, 100, 200, 300],
            'end': [150, 250, 150, 250, 350],
            'segment_in_circle': [1, 2, 1, 2, 1],
            'copyNum': [5, 5, 3, 3, 2]
        })
        
        result = processor.cluster_by_signature(df)
        
        # q1 and q2 should be clustered (same segments)
        assert len(result) == 3  # 2 segments for cluster + 1 for q3
        clustered = result[result['cluster_size'] > 1]
        assert len(clustered) == 2  # Two segments in the cluster
        assert all(clustered['copyNum'] == 8)  # 5 + 3
    
    def test_process_full_pipeline(self, setup_processor):
        processor, tmp_path = setup_processor
        
        # Create test CSV with multi-segment circular DNA
        csv_file = tmp_path / "cecc.csv"
        csv_file.write_text("""query_id,chr,start,end,segment_in_circle,q_start,consLen
read1,chr1,100,150,1,5,12
read1,chr2,200,250,2,5,12
read2,chr1,300,350,1,3,10
""")
        
        output_dir = tmp_path / "output"
        result = processor.process([csv_file], output_dir, "test", cluster=True)
        
        assert result is not None
        assert (output_dir / "test_CeccDNA_processed.csv").exists()
        assert (output_dir / "test_CeccDNA_pre.fasta").exists()
        
        # Check sequences were generated
        assert len(processor.fasta_records) > 0


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
        
        # Create test CSVs
        uecc_csv = tmp_path / "uecc.csv"
        uecc_csv.write_text("""query_id,eName,q_start,eLength
uecc1,repeat1,5,8
""")
        
        mecc_csv = tmp_path / "mecc.csv"
        mecc_csv.write_text("""query_id,eChr,eStart,eEnd,q_start,eLength
mecc1,chr1,100,150,3,6
""")
        
        cecc_csv = tmp_path / "cecc.csv"
        cecc_csv.write_text("""query_id,chr,start,end,q_start,consLen
cecc1,chr1,100,150,1,10
""")
        
        return fasta_file, uecc_csv, mecc_csv, cecc_csv, tmp_path
    
    def test_process_all_with_xecc(self, setup_pipeline):
        fasta_file, uecc_csv, mecc_csv, cecc_csv, tmp_path = setup_pipeline
        
        config = UMCProcessConfig(
            process_xecc=True,
            xecc_first=True,
            cluster_mecc=True,
            cluster_cecc=True,
            aggregate_uecc=True
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
    
    def test_process_all_xecc_last(self, setup_pipeline):
        fasta_file, uecc_csv, mecc_csv, cecc_csv, tmp_path = setup_pipeline
        
        config = UMCProcessConfig(
            process_xecc=True,
            xecc_first=False
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
        # XeccDNA should still be created
        assert results["xecc"].exists()
    
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
    
    def test_run_pipeline_backward_compatibility(self, setup_pipeline):
        fasta_file, uecc_csv, mecc_csv, cecc_csv, tmp_path = setup_pipeline
        
        processor = UMCProcess()
        output_dir = tmp_path / "output"
        
        # Use legacy method
        results = processor.run_pipeline(
            fasta_file=fasta_file,
            uecc_files=[uecc_csv],
            mecc_files=[mecc_csv],
            cecc_files=[cecc_csv],
            output_dir=output_dir,
            prefix="test",
            process_xecc=False
        )
        
        # Should only return DataFrames, not XeccDNA path
        assert all(isinstance(v, pd.DataFrame) for v in results.values())
        assert len(results) == 3
    
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


class TestBackwardCompatibility:
    def test_menagerie_alias(self):
        # Test that backward compatibility aliases exist
        assert hasattr(umc_module, 'Menagerie')
        assert hasattr(umc_module, 'MenagerieConfig')
        
        # They should be subclasses/aliases
        assert issubclass(umc_module.Menagerie, UMCProcess)
        assert issubclass(umc_module.MenagerieConfig, UMCProcessConfig)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])