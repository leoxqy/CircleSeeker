"""Integration tests for eccDNA curation with real data.

These tests require external data files from tests/curation_demo/ and optionally
a reference genome. They are skipped if the required files are not available.
"""

from pathlib import Path

import pandas as pd
import pytest

from circleseeker.modules.iecc_curator import (
    curate_ecc_tables,
    generate_fasta_sequences,
    process_eccDNA,
    write_curated_tables,
)


@pytest.mark.integration
@pytest.mark.requires_data
class TestRealDataCuration:
    """Test eccDNA curation with real Cyrcular data from the curation demo."""

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

    @pytest.mark.requires_reference
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

    @pytest.mark.requires_reference
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
