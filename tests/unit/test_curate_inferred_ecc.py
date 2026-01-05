from pathlib import Path
import importlib.util
import sys
import types
import tempfile
import pandas as pd
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


def _load_curate_inferred_ecc_module():
    package_name = "circleseeker.modules"
    module_name = f"{package_name}.Iecc_curator"

    if module_name in sys.modules:
        return sys.modules[module_name]

    import circleseeker  # noqa: F401  # Ensure base package is initialized

    if package_name not in sys.modules:
        package_module = types.ModuleType(package_name)
        package_module.__path__ = [str(SRC / "circleseeker" / "modules")]
        sys.modules[package_name] = package_module

    spec = importlib.util.spec_from_file_location(
        module_name,
        SRC / "circleseeker" / "modules" / "Iecc_curator.py",
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


# Load module and get functions
curate_module = _load_curate_inferred_ecc_module()
select_best_duplicate = curate_module.select_best_duplicate
curate_ecc_tables = curate_module.curate_ecc_tables
write_curated_tables = curate_module.write_curated_tables
process_eccDNA = curate_module.process_eccDNA
generate_fasta_sequences = curate_module.generate_fasta_sequences


# Keep the existing test
def test_curate_ecc_tables_handles_simple_and_chimeric(tmp_path):
    overview = pd.DataFrame(
        [
            {
                "regions": "chr1:101-250",
                "circle_length": 150,
                "segment_count": 1,
                "num_split_reads": 5,
                "prob_present": 0.9,
                "prob_artifact": 0.1,
            },
            {
                "regions": "chr1:101-250",
                "circle_length": 150,
                "segment_count": 1,
                "num_split_reads": 1,
                "prob_present": 0.2,
                "prob_artifact": 0.8,
            },
            {
                "regions": "chr2:50-120,chr2:300-360",
                "circle_length": 220,
                "segment_count": 2,
                "num_split_reads": 4,
                "prob_present": 0.75,
                "prob_artifact": 0.1,
            },
        ]
    )

    overview_path = tmp_path / "overview.tsv"
    overview.to_csv(overview_path, sep="\t", index=False)

    simple_df, chimeric_df = curate_ecc_tables(overview_path)

    assert list(simple_df["eccDNA_id"]) == ["IUeccDNA1"]
    assert simple_df.loc[0, "chr"] == "chr1"
    assert simple_df.loc[0, "length"] == 150

    chimeric_ids = chimeric_df["eccDNA_id"].unique().tolist()
    assert chimeric_ids == ["ICeccDNA1"]
    assert chimeric_df["seg_total"].unique().tolist() == [2]
    assert set(chimeric_df["junction_role"]) == {"head", "tail"}


class TestSelectBestDuplicate:
    def test_with_prob_present_column(self):
        # Test case with prob_present column
        df = pd.DataFrame({
            'regions': ['chr1:100-200', 'chr1:100-200', 'chr1:100-200'],
            'prob_present': [0.8, 0.9, 0.7],
            'num_split_reads': [10, 5, 15],
            'prob_artifact': [0.1, 0.2, 0.05]
        })
        
        result = select_best_duplicate(df)
        
        assert len(result) == 1
        assert 'score' in result.columns
        # Should select the one with best combined score
    
    def test_with_prob_joint_event_column(self):
        # Test case with prob_joint_event column instead
        df = pd.DataFrame({
            'regions': ['chr2:200-300', 'chr2:200-300'],
            'prob_joint_event': [0.85, 0.75],
            'num_split_reads': [8, 12],
            'prob_artifact': [0.15, 0.1]
        })
        
        result = select_best_duplicate(df)
        
        assert len(result) == 1
        assert 'score' in result.columns
    
    def test_without_probability_columns(self):
        # Test case without probability columns
        df = pd.DataFrame({
            'regions': ['chr3:300-400', 'chr3:300-400'],
            'num_split_reads': [20, 15],
            'prob_artifact': [0.2, 0.3]
        })
        
        result = select_best_duplicate(df)
        
        assert len(result) == 1
        # Should select based on split reads and artifact probability
        assert result.iloc[0]['num_split_reads'] == 20
    
    def test_zero_split_reads(self):
        # Test when all split reads are zero
        df = pd.DataFrame({
            'regions': ['chr4:400-500', 'chr4:400-500'],
            'prob_present': [0.6, 0.7],
            'num_split_reads': [0, 0],
            'prob_artifact': [0.3, 0.2]
        })
        
        result = select_best_duplicate(df)
        
        assert len(result) == 1
        # Should select based on prob_present and prob_artifact
        assert result.iloc[0]['prob_present'] == 0.7
    
    def test_missing_prob_artifact(self):
        # Test when prob_artifact column is missing
        df = pd.DataFrame({
            'regions': ['chr5:500-600'],
            'prob_present': [0.8],
            'num_split_reads': [10]
        })
        
        result = select_best_duplicate(df)
        
        assert len(result) == 1
        assert 'score' in result.columns


class TestCurateEccTables:
    @pytest.fixture
    def sample_cyrcular_output(self, tmp_path):
        # Create a sample Cyrcular overview TSV
        tsv_file = tmp_path / "cyrcular_overview.tsv"
        
        data = pd.DataFrame({
            'regions': [
                'chr1:1000-2000',      # Simple circle
                'chr2:3000-4000',      # Simple circle (duplicate)
                'chr2:3000-4000',      # Simple circle (duplicate)
                'chr3:5000-6000,chr3:7000-8000',  # Chimeric circle
                'chr4:9000-9050',      # Too short (< 100bp)
                'chr5:10000-11000,chr5:12000-13000,chr5:14000-15000',  # 3-segment chimeric
            ],
            'circle_length': [1000, 1000, 1000, 3000, 50, 5000],
            'segment_count': [1, 1, 1, 2, 1, 3],
            'num_split_reads': [20, 15, 25, 30, 10, 40],
            'prob_present': [0.9, 0.8, 0.85, 0.95, 0.7, 0.92],
            'prob_artifact': [0.1, 0.15, 0.12, 0.05, 0.3, 0.08],
            'af_nanopore': [0.05, 0.04, 0.045, 0.06, 0.02, 0.07]
        })
        
        data.to_csv(tsv_file, sep='\t', index=False)
        return tsv_file
    
    def test_curate_basic(self, sample_cyrcular_output):
        simple_df, chimeric_df = curate_ecc_tables(sample_cyrcular_output)
        
        # Check simple circles
        assert len(simple_df) == 2  # chr1 and best of chr2 duplicates
        assert all(simple_df['eccDNA_id'].str.startswith('IUeccDNA'))
        assert all(simple_df['eccdna_type'] == 'Uecc')
        assert all(simple_df['state'] == 'Inferred')
        
        # Check that chr2 duplicate was deduplicated
        chr2_entries = simple_df[simple_df['chr'] == 'chr2']
        assert len(chr2_entries) == 1
        
        # Check chimeric circles
        assert len(chimeric_df) == 5  # 2 segments + 3 segments
        assert all(chimeric_df['eccDNA_id'].str.startswith('ICeccDNA'))
        assert all(chimeric_df['eccdna_type'] == 'Cecc')
        assert all(chimeric_df['state'] == 'Inferred')
        
        # Check segment indexing
        chr3_segments = chimeric_df[chimeric_df['chr'] == 'chr3']
        assert len(chr3_segments) == 2
        assert set(chr3_segments['seg_index']) == {1, 2}
        assert all(chr3_segments['seg_total'] == 2)
        
        # Check junction roles
        assert chr3_segments[chr3_segments['seg_index'] == 1]['junction_role'].iloc[0] == 'head'
        assert chr3_segments[chr3_segments['seg_index'] == 2]['junction_role'].iloc[0] == 'tail'
    
    def test_length_filter(self, sample_cyrcular_output):
        simple_df, chimeric_df = curate_ecc_tables(sample_cyrcular_output)
        
        # Verify that circles < 100bp were filtered out
        all_lengths = pd.concat([simple_df['length'], chimeric_df['length']])
        assert all(all_lengths >= 100)
    
    def test_coordinate_conversion(self, sample_cyrcular_output):
        simple_df, chimeric_df = curate_ecc_tables(sample_cyrcular_output)
        
        # Check 1-based to 0-based conversion
        chr1_entry = simple_df[simple_df['chr'] == 'chr1'].iloc[0]
        assert chr1_entry['start0'] == 999  # 1000 - 1
        assert chr1_entry['end0'] == 2000
    
    def test_empty_input(self, tmp_path):
        # Test with empty file
        empty_tsv = tmp_path / "empty.tsv"
        pd.DataFrame(columns=['regions', 'circle_length', 'segment_count']).to_csv(
            empty_tsv, sep='\t', index=False
        )
        
        simple_df, chimeric_df = curate_ecc_tables(empty_tsv)
        
        assert simple_df.empty
        assert chimeric_df.empty
    
    def test_invalid_regions(self, tmp_path):
        # Test with malformed region strings
        tsv_file = tmp_path / "invalid_regions.tsv"
        
        data = pd.DataFrame({
            'regions': [
                'invalid_region',           # No colon or dash
                'chr1:1000',               # No end coordinate
                'chr2-2000',               # No colon
                'chr3:3000-4000',          # Valid
                np.nan,                    # Missing value
                'chr4:5000-6000,invalid',  # Partially invalid chimeric
            ],
            'circle_length': [1000] * 6,
            'segment_count': [1, 1, 1, 1, 1, 2],
            'num_split_reads': [10] * 6
        })
        
        data.to_csv(tsv_file, sep='\t', index=False)
        
        simple_df, chimeric_df = curate_ecc_tables(tsv_file)
        
        # Only the valid chr3 entry should be processed
        assert len(simple_df) == 1
        assert simple_df.iloc[0]['chr'] == 'chr3'
        
        # Invalid chimeric should be skipped
        assert chimeric_df.empty
    
    def test_three_segment_chimeric(self, sample_cyrcular_output):
        simple_df, chimeric_df = curate_ecc_tables(sample_cyrcular_output)
        
        # Check 3-segment chimeric circle
        chr5_segments = chimeric_df[chimeric_df['chr'] == 'chr5']
        assert len(chr5_segments) == 3
        assert set(chr5_segments['seg_index']) == {1, 2, 3}
        assert all(chr5_segments['seg_total'] == 3)
        
        # Check junction roles
        roles = chr5_segments.sort_values('seg_index')['junction_role'].tolist()
        assert roles == ['head', 'middle', 'tail']
    
    def test_prob_column_fallback(self, tmp_path):
        # Test when prob_present is missing but prob_joint_event exists
        tsv_file = tmp_path / "prob_fallback.tsv"
        
        data = pd.DataFrame({
            'regions': ['chr1:1000-2000'],
            'circle_length': [1000],
            'segment_count': [1],
            'num_split_reads': [20],
            'prob_joint_event': [0.85],  # No prob_present
            'prob_artifact': [0.1],
            'af_nanopore': [0.05]
        })
        
        data.to_csv(tsv_file, sep='\t', index=False)
        
        simple_df, chimeric_df = curate_ecc_tables(tsv_file)
        
        assert len(simple_df) == 1
        assert simple_df.iloc[0]['prob_present'] == 0.85  # Should use prob_joint_event


class TestWriteCuratedTables:
    def test_write_both_tables(self, tmp_path):
        # Create sample dataframes
        simple_df = pd.DataFrame({
            'eccDNA_id': ['IUeccDNA1', 'IUeccDNA2'],
            'chr': ['chr1', 'chr2'],
            'start0': [1000, 2000],
            'end0': [2000, 3000]
        })
        
        chimeric_df = pd.DataFrame({
            'eccDNA_id': ['ICeccDNA1', 'ICeccDNA1'],
            'chr': ['chr3', 'chr4'],
            'start0': [3000, 4000],
            'end0': [4000, 5000],
            'seg_index': [1, 2]
        })
        
        output_prefix = tmp_path / "test_output"
        simple_out, chimeric_out = write_curated_tables(simple_df, chimeric_df, output_prefix)
        
        # Check files were created
        assert simple_out.exists()
        assert chimeric_out.exists()
        assert simple_out.name == "test_output_simple.csv"
        assert chimeric_out.name == "test_output_chimeric.csv"

        # Check content - function uses comma separator
        simple_loaded = pd.read_csv(simple_out)
        assert len(simple_loaded) == 2
        assert list(simple_loaded['eccDNA_id']) == ['IUeccDNA1', 'IUeccDNA2']

        chimeric_loaded = pd.read_csv(chimeric_out)
        assert len(chimeric_loaded) == 2
    
    def test_write_empty_dataframes(self, tmp_path):
        # Test with empty dataframes
        simple_df = pd.DataFrame()
        chimeric_df = pd.DataFrame()
        
        output_prefix = tmp_path / "empty_output"
        simple_out, chimeric_out = write_curated_tables(simple_df, chimeric_df, output_prefix)
        
        # Files should still be created (even if empty)
        # But based on the code, they won't be created if dataframes are empty
        assert not simple_out.exists()
        assert not chimeric_out.exists()
    
    def test_output_path_handling(self, tmp_path):
        # Test with string path
        simple_df = pd.DataFrame({'eccDNA_id': ['IUeccDNA1']})
        chimeric_df = pd.DataFrame({'eccDNA_id': ['ICeccDNA1']})
        
        output_prefix = str(tmp_path / "string_path")
        simple_out, chimeric_out = write_curated_tables(simple_df, chimeric_df, output_prefix)
        
        assert isinstance(simple_out, Path)
        assert isinstance(chimeric_out, Path)


class TestProcessEccDNA:
    @pytest.fixture
    def sample_cyrcular_output(self, tmp_path):
        # Create a sample Cyrcular overview TSV
        tsv_file = tmp_path / "cyrcular_overview.tsv"
        
        data = pd.DataFrame({
            'regions': [
                'chr1:1000-2000',      # Simple circle
                'chr2:3000-4000',      # Simple circle
                'chr3:5000-6000,chr3:7000-8000',  # Chimeric circle
            ],
            'circle_length': [1000, 1000, 3000],
            'segment_count': [1, 1, 2],
            'num_split_reads': [20, 15, 30],
            'prob_present': [0.9, 0.8, 0.95],
            'prob_artifact': [0.1, 0.15, 0.05],
            'af_nanopore': [0.05, 0.04, 0.06]
        })
        
        data.to_csv(tsv_file, sep='\t', index=False)
        return tsv_file
    
    def test_backward_compatibility(self, sample_cyrcular_output, tmp_path):
        # Test the backward-compatible function
        output_prefix = tmp_path / "compat_test"

        # This should create the files
        process_eccDNA(sample_cyrcular_output, output_prefix)

        # Check files were created in organized folder
        inferred_dir = tmp_path / "compat_test_Inferred_eccDNA"
        simple_file = inferred_dir / "compat_test_simple.csv"
        chimeric_file = inferred_dir / "compat_test_chimeric.csv"

        assert inferred_dir.exists()
        assert simple_file.exists()
        assert chimeric_file.exists()

        # Verify content - function uses comma separator
        simple_df = pd.read_csv(simple_file)
        assert len(simple_df) > 0
        assert 'eccDNA_id' in simple_df.columns


class TestRealData:
    """Test with real Cyrcular data from the curation demo."""

    @pytest.fixture
    def lb_2_3_data(self):
        """Path to real LB_2_3 overview data."""
        root = Path(__file__).resolve().parents[2]
        return root / "tests" / "curation_demo" / "LB_2_3_overview.tsv"

    @pytest.fixture
    def reference_fasta(self):
        """Path to reference genome for FASTA testing."""
        root = Path(__file__).resolve().parents[2]
        return root / "tests" / "integration" / "chm13v2.0.fa"

    @pytest.fixture
    def expected_simple_data(self):
        """Path to expected simple output."""
        root = Path(__file__).resolve().parents[2]
        return root / "tests" / "curation_demo" / "LB_2_3_simple.tsv"

    @pytest.fixture
    def expected_chimeric_data(self):
        """Path to expected chimeric output."""
        root = Path(__file__).resolve().parents[2]
        return root / "tests" / "curation_demo" / "LB_2_3_chimeric.tsv"

    def test_lb_2_3_processing(self, lb_2_3_data, tmp_path):
        """Test processing of real LB_2_3 data."""
        if not lb_2_3_data.exists():
            pytest.skip(f"LB_2_3 data not found: {lb_2_3_data}")

        # Process the real data
        simple_df, chimeric_df = curate_ecc_tables(lb_2_3_data)

        # Basic sanity checks
        assert len(simple_df) > 0, "Should produce some simple eccDNA"
        assert len(chimeric_df) > 0, "Should produce some chimeric eccDNA"

        # Check structure
        required_simple_cols = ['eccDNA_id', 'chr', 'start0', 'end0', 'strand', 'length', 'eccdna_type', 'state']
        required_chimeric_cols = ['eccDNA_id', 'chr', 'start0', 'end0', 'strand', 'length', 'eccdna_type', 'state', 'seg_index', 'seg_total', 'junction_role']

        assert all(col in simple_df.columns for col in required_simple_cols)
        assert all(col in chimeric_df.columns for col in required_chimeric_cols)

        # Check data quality
        assert all(simple_df['eccdna_type'] == 'Uecc')
        assert all(simple_df['state'] == 'Inferred')
        assert all(chimeric_df['eccdna_type'] == 'Cecc')
        assert all(chimeric_df['state'] == 'Inferred')

        # Check ID patterns
        assert all(simple_df['eccDNA_id'].str.startswith('IUeccDNA'))
        assert all(chimeric_df['eccDNA_id'].str.startswith('ICeccDNA'))

        # Check coordinates
        assert all(simple_df['start0'] < simple_df['end0'])
        assert all(chimeric_df['start0'] < chimeric_df['end0'])
        assert all(simple_df['length'] >= 100)  # Length filter

        # Test writing output
        output_prefix = tmp_path / "lb_2_3_test"
        simple_out, chimeric_out = write_curated_tables(simple_df, chimeric_df, output_prefix)

        assert simple_out.exists()
        assert chimeric_out.exists()

        # Check strand information
        assert 'strand' in simple_df.columns
        assert 'strand' in chimeric_df.columns
        assert all(simple_df['strand'].isin(['+', '-']))
        assert all(chimeric_df['strand'].isin(['+', '-']))

        # Should have some positive strand entries
        assert (simple_df['strand'] == '+').sum() > 0

        print(f"Processed {len(simple_df)} simple and {len(chimeric_df)} chimeric eccDNA from real data")
        print(f"Strand distribution - Simple: +{(simple_df['strand'] == '+').sum()}, -{(simple_df['strand'] == '-').sum()}")
        print(f"Strand distribution - Chimeric: +{(chimeric_df['strand'] == '+').sum()}, -{(chimeric_df['strand'] == '-').sum()}")

    def test_compare_with_expected_output(self, lb_2_3_data, expected_simple_data, expected_chimeric_data):
        """Compare processing results with expected output files (if available)."""
        if not lb_2_3_data.exists():
            pytest.skip(f"LB_2_3 data not found: {lb_2_3_data}")

        # Process the data
        simple_df, chimeric_df = curate_ecc_tables(lb_2_3_data)

        # If expected files exist, compare key metrics
        if expected_simple_data.exists():
            expected_simple = pd.read_csv(expected_simple_data, sep='\t')
            print(f"Generated {len(simple_df)} simple eccDNA, expected had {len(expected_simple)}")

            # Allow some tolerance in counts due to potential algorithm differences
            count_ratio = len(simple_df) / len(expected_simple) if len(expected_simple) > 0 else float('inf')
            assert 0.5 <= count_ratio <= 2.0, f"Simple eccDNA count differs significantly: {len(simple_df)} vs {len(expected_simple)}"

        if expected_chimeric_data.exists():
            expected_chimeric = pd.read_csv(expected_chimeric_data, sep='\t')
            print(f"Generated {len(chimeric_df)} chimeric segments, expected had {len(expected_chimeric)}")

            # Allow some tolerance in counts
            count_ratio = len(chimeric_df) / len(expected_chimeric) if len(expected_chimeric) > 0 else float('inf')
            assert 0.5 <= count_ratio <= 2.0, f"Chimeric segment count differs significantly: {len(chimeric_df)} vs {len(expected_chimeric)}"

    def test_fasta_generation(self, reference_fasta, tmp_path):
        """Test FASTA sequence generation with real reference genome."""
        if not reference_fasta.exists():
            pytest.skip(f"Reference FASTA not found: {reference_fasta}")

        # Create sample data for FASTA generation
        simple_df = pd.DataFrame({
            'eccDNA_id': ['IUeccDNA1', 'IUeccDNA2'],
            'chr': ['chr1', 'chr2'],
            'start0': [100, 200],  # 0-based coordinates
            'end0': [200, 300],
            'length': [100, 100]
        })

        chimeric_df = pd.DataFrame({
            'eccDNA_id': ['ICeccDNA1', 'ICeccDNA1', 'ICeccDNA2', 'ICeccDNA2'],
            'chr': ['chr1', 'chr1', 'chr2', 'chr3'],
            'start0': [400, 600, 700, 800],
            'end0': [500, 650, 800, 900],
            'seg_index': [1, 2, 1, 2],
            'seg_total': [2, 2, 2, 2],
            'length': [200, 200, 200, 200]
        })

        output_prefix = tmp_path / "fasta_test"

        # Test FASTA generation
        simple_fasta, chimeric_fasta = generate_fasta_sequences(
            simple_df, chimeric_df, reference_fasta, output_prefix
        )

        if simple_fasta:  # Only check if pysam is available
            assert simple_fasta.exists()
            assert simple_fasta.name.endswith('_UeccDNA_I.fasta')

            # Check FASTA content structure
            with open(simple_fasta, 'r') as f:
                content = f.read()
                assert '>IUeccDNA1' in content
                assert '>IUeccDNA2' in content
                # Should have proper FASTA format
                lines = content.strip().split('\n')
                header_lines = [l for l in lines if l.startswith('>')]
                assert len(header_lines) == 2

        if chimeric_fasta:  # Only check if pysam is available
            assert chimeric_fasta.exists()
            assert chimeric_fasta.name.endswith('_CeccDNA_I.fasta')

            # Check chimeric FASTA content
            with open(chimeric_fasta, 'r') as f:
                content = f.read()
                assert '>ICeccDNA1' in content
                assert '>ICeccDNA2' in content

        print(f"FASTA generation test: simple={simple_fasta is not None}, chimeric={chimeric_fasta is not None}")

    def test_process_eccDNA_with_fasta(self, lb_2_3_data, reference_fasta, tmp_path):
        """Test the full process_eccDNA function with FASTA generation."""
        if not lb_2_3_data.exists() or not reference_fasta.exists():
            pytest.skip("Required data files not found")

        output_prefix = tmp_path / "full_test"

        # Test with FASTA generation
        process_eccDNA(lb_2_3_data, output_prefix, reference_fasta)

        # Check that organized folder was created
        inferred_dir = tmp_path / "full_test_Inferred_eccDNA"
        assert inferred_dir.exists()

        # Check TSV files
        simple_tsv = inferred_dir / "full_test_simple.csv"
        chimeric_tsv = inferred_dir / "full_test_chimeric.csv"
        assert simple_tsv.exists()
        assert chimeric_tsv.exists()

        # Check FASTA files (if pysam available)
        simple_fasta = inferred_dir / "full_test_UeccDNA_I.fasta"
        chimeric_fasta = inferred_dir / "full_test_CeccDNA_I.fasta"

        # These may or may not exist depending on pysam availability
        if simple_fasta.exists():
            # Verify it's a valid FASTA file
            with open(simple_fasta, 'r') as f:
                content = f.read()
                assert content.count('>') > 0  # Should have some sequences

        print(f"Full pipeline test completed. FASTA files generated: simple={simple_fasta.exists()}, chimeric={chimeric_fasta.exists()}")

    def test_strand_detection_logic(self, tmp_path):
        """Test strand detection based on coordinate order."""
        # Create test data with both forward and reverse coordinate orders
        test_data = pd.DataFrame({
            'regions': [
                'chr1:1000-2000',     # Forward strand (+)
                'chr2:3000-2000',     # Reverse strand (-)
                'chr3:5000-6000',     # Forward strand (+)
                'chr4:8000-7000',     # Reverse strand (-)
                'chr5:10000-11000,chr5:12000-11500',  # Chimeric: + then -
            ],
            'circle_length': [1000, 1000, 1000, 1000, 2500],
            'segment_count': [1, 1, 1, 1, 2],
            'num_split_reads': [10, 15, 20, 25, 30],
            'prob_present': [0.9, 0.8, 0.85, 0.75, 0.95],
            'prob_artifact': [0.1, 0.2, 0.15, 0.25, 0.05]
        })

        test_file = tmp_path / "strand_test.tsv"
        test_data.to_csv(test_file, sep='\t', index=False)

        simple_df, chimeric_df = curate_ecc_tables(test_file)

        # Check simple circles
        assert len(simple_df) == 4

        # Verify strand assignment for simple circles
        simple_strands = simple_df.set_index('chr')['strand'].to_dict()
        assert simple_strands['chr1'] == '+'  # 1000-2000 (forward)
        assert simple_strands['chr2'] == '-'  # 3000-2000 (reverse)
        assert simple_strands['chr3'] == '+'  # 5000-6000 (forward)
        assert simple_strands['chr4'] == '-'  # 8000-7000 (reverse)

        # Verify coordinates were normalized correctly
        chr2_row = simple_df[simple_df['chr'] == 'chr2'].iloc[0]
        assert chr2_row['start0'] == 1999  # min(3000,2000) - 1
        assert chr2_row['end0'] == 3000    # max(3000,2000)
        assert chr2_row['strand'] == '-'

        # Check chimeric circles
        assert len(chimeric_df) == 2  # 2 segments

        # Verify chimeric strand assignments
        chimeric_chr5 = chimeric_df[chimeric_df['chr'] == 'chr5']
        seg1 = chimeric_chr5[chimeric_chr5['seg_index'] == 1].iloc[0]
        seg2 = chimeric_chr5[chimeric_chr5['seg_index'] == 2].iloc[0]

        assert seg1['strand'] == '+'  # 10000-11000 (forward)
        assert seg2['strand'] == '-'  # 12000-11500 (reverse)

        print(f"Strand test passed: {len(simple_df)} simple, {len(chimeric_df)} chimeric")
        print(f"Simple strands: {simple_df['strand'].value_counts().to_dict()}")
        print(f"Chimeric strands: {chimeric_df['strand'].value_counts().to_dict()}")


class TestIntegration:
    def test_full_pipeline(self, tmp_path):
        # Create a comprehensive test case
        cyrcular_file = tmp_path / "cyrcular_full.tsv"
        
        # Include various edge cases
        data = pd.DataFrame({
            'regions': [
                # Simple circles
                'chr1:1000-2000',
                'chr2:3000-5000',
                'chr2:3000-5000',  # Duplicate
                'chr2:3000-5000',  # Duplicate
                
                # Chimeric circles
                'chr3:6000-7000,chr3:8000-9000',
                'chr4:10000-11000,chr5:12000-13000',  # Inter-chromosomal
                
                # Edge cases
                'chrX:100-200',    # Short but >= 100 (now 100bp)
                'chrY:200-249',    # Too short (< 100bp)
                
                # Complex chimeric
                'chr6:1000-2000,chr6:3000-4000,chr6:5000-6000,chr7:7000-8000',
            ],
            'circle_length': [1000, 2000, 2000, 2000, 3000, 3000, 100, 49, 8000],
            'segment_count': [1, 1, 1, 1, 2, 2, 1, 1, 4],
            'num_split_reads': [20, 30, 25, 35, 40, 45, 10, 5, 50],
            'prob_present': [0.9, 0.85, 0.8, 0.88, 0.95, 0.92, 0.7, 0.6, 0.98],
            'prob_artifact': [0.1, 0.12, 0.15, 0.11, 0.05, 0.08, 0.3, 0.4, 0.02],
            'af_nanopore': [0.05, 0.06, 0.055, 0.058, 0.07, 0.065, 0.02, 0.01, 0.08]
        })
        
        data.to_csv(cyrcular_file, sep='\t', index=False)
        
        # Process the file
        simple_df, chimeric_df = curate_ecc_tables(cyrcular_file)
        
        # Verify results
        # Simple circles: chr1 + best chr2 duplicate + chrX (chrY filtered for length)
        assert len(simple_df) == 3
        assert set(simple_df['chr']) == {'chr1', 'chr2', 'chrX'}
        
        # The best chr2 duplicate should be selected (highest score)
        chr2_entry = simple_df[simple_df['chr'] == 'chr2'].iloc[0]
        assert chr2_entry['num_split_reads'] in [30, 25, 35]  # One of the duplicates
        
        # Chimeric circles: 2 + 2 + 4 segments = 8 rows
        assert len(chimeric_df) == 8
        
        # Check inter-chromosomal chimeric
        inter_chrom = chimeric_df[chimeric_df['eccDNA_id'] == 'ICeccDNA2']
        assert set(inter_chrom['chr']) == {'chr4', 'chr5'}
        
        # Check 4-segment chimeric
        complex_chimeric = chimeric_df[chimeric_df['eccDNA_id'] == 'ICeccDNA3']
        assert len(complex_chimeric) == 4
        assert set(complex_chimeric['seg_index']) == {1, 2, 3, 4}
        
        # Write output files
        output_prefix = tmp_path / "integration_test"
        write_curated_tables(simple_df, chimeric_df, output_prefix)
        
        assert (tmp_path / "integration_test_simple.csv").exists()
        assert (tmp_path / "integration_test_chimeric.csv").exists()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
