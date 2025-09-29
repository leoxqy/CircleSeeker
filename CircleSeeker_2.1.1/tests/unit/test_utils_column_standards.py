"""Tests for utils column_standards module."""

from pathlib import Path
import sys
import pytest
import pandas as pd
import warnings
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.utils.column_standards import (
    ColumnStandard, ColumnStandardizer, standardize_columns,
    LEGACY_MAPPINGS, COORDINATE_CONVERSIONS, SCHEMAS
)


class TestColumnStandard:
    """Test cases for ColumnStandard dataclass."""

    def test_column_standard_attributes(self):
        """Test that ColumnStandard has expected attributes."""
        # Core identity fields
        assert ColumnStandard.ECCDNA_ID == 'eccDNA_id'
        assert ColumnStandard.ORIG_ECCDNA_ID == 'orig_eccDNA_id'
        assert ColumnStandard.CLUSTER_ID == 'cluster_id'

        # Genomic coordinates
        assert ColumnStandard.CHR == 'chr'
        assert ColumnStandard.START0 == 'start0'
        assert ColumnStandard.END0 == 'end0'
        assert ColumnStandard.STRAND == 'strand'
        assert ColumnStandard.LENGTH == 'length'

        # Sequence and read information
        assert ColumnStandard.READS == 'reads'
        assert ColumnStandard.SEQUENCE == 'sequence'
        assert ColumnStandard.READ_COUNT == 'read_count'

        # Quantitative measures
        assert ColumnStandard.COPY_NUMBER == 'copy_number'
        assert ColumnStandard.MATCH_DEGREE == 'match_degree'
        assert ColumnStandard.REPEAT_NUMBER == 'repeat_number'

        # Type and classification
        assert ColumnStandard.ECCDNA_TYPE == 'eccdna_type'
        assert ColumnStandard.STATE == 'state'

    def test_column_standard_is_string_values(self):
        """Test that all ColumnStandard values are strings."""
        for attr in dir(ColumnStandard):
            if not attr.startswith('_'):
                value = getattr(ColumnStandard, attr)
                assert isinstance(value, str), f"{attr} should be a string, got {type(value)}"


class TestLegacyMappings:
    """Test cases for legacy mapping dictionaries."""

    def test_legacy_mappings_exists(self):
        """Test that LEGACY_MAPPINGS dictionary exists and has content."""
        assert isinstance(LEGACY_MAPPINGS, dict)
        assert len(LEGACY_MAPPINGS) > 0

    def test_legacy_mappings_key_types(self):
        """Test that all legacy mapping keys are strings."""
        for key in LEGACY_MAPPINGS.keys():
            assert isinstance(key, str)

    def test_legacy_mappings_value_types(self):
        """Test that all legacy mapping values are strings."""
        for value in LEGACY_MAPPINGS.values():
            assert isinstance(value, str)

    def test_specific_legacy_mappings(self):
        """Test specific legacy mappings."""
        assert LEGACY_MAPPINGS['eChr'] == ColumnStandard.CHR
        assert LEGACY_MAPPINGS['eStart0'] == ColumnStandard.START0
        assert LEGACY_MAPPINGS['eEnd0'] == ColumnStandard.END0
        assert LEGACY_MAPPINGS['readName'] == ColumnStandard.READS
        assert LEGACY_MAPPINGS['copyNum'] == ColumnStandard.COPY_NUMBER

    def test_coordinate_conversions_exists(self):
        """Test that COORDINATE_CONVERSIONS dictionary exists."""
        assert isinstance(COORDINATE_CONVERSIONS, dict)
        assert 'eStart' in COORDINATE_CONVERSIONS
        assert 's_start' in COORDINATE_CONVERSIONS

    def test_coordinate_conversion_functions(self):
        """Test coordinate conversion functions."""
        # Test eStart conversion (1-based to 0-based)
        estart_func = COORDINATE_CONVERSIONS['eStart']
        assert estart_func(1) == 0
        assert estart_func(100) == 99

        # Test eEnd conversion (should remain the same)
        eend_func = COORDINATE_CONVERSIONS['eEnd']
        assert eend_func(100) == 100
        assert eend_func(200) == 200


class TestColumnStandardizer:
    """Test cases for ColumnStandardizer class."""

    def test_column_standardizer_initialization(self):
        """Test ColumnStandardizer initialization."""
        # Default initialization
        standardizer = ColumnStandardizer()
        assert standardizer.strict_mode is False

        # Strict mode initialization
        strict_standardizer = ColumnStandardizer(strict_mode=True)
        assert strict_standardizer.strict_mode is True

    def test_get_all_standard_columns(self):
        """Test _get_all_standard_columns method."""
        columns = ColumnStandardizer._get_all_standard_columns()
        assert isinstance(columns, list)
        assert len(columns) > 0
        assert ColumnStandard.ECCDNA_ID in columns
        assert ColumnStandard.CHR in columns
        assert ColumnStandard.START0 in columns

    def test_standardize_dataframe_empty(self):
        """Test standardize_dataframe with empty DataFrame."""
        standardizer = ColumnStandardizer()
        empty_df = pd.DataFrame()
        result = standardizer.standardize_dataframe(empty_df)
        assert result.empty
        assert isinstance(result, pd.DataFrame)

    def test_standardize_dataframe_basic(self):
        """Test basic DataFrame standardization."""
        test_df = pd.DataFrame({
            'eChr': ['chr1', 'chr2'],
            'eStart0': [100, 200],
            'readName': ['read1', 'read2'],
            'copyNum': [1, 2]
        })

        standardizer = ColumnStandardizer()
        result = standardizer.standardize_dataframe(test_df)

        # Check that columns were renamed
        assert ColumnStandard.CHR in result.columns
        assert ColumnStandard.START0 in result.columns
        assert ColumnStandard.READS in result.columns
        assert ColumnStandard.COPY_NUMBER in result.columns

        # Check that old columns are gone
        assert 'eChr' not in result.columns
        assert 'readName' not in result.columns

        # Check data integrity
        assert result[ColumnStandard.CHR].tolist() == ['chr1', 'chr2']
        assert result[ColumnStandard.READS].tolist() == ['read1', 'read2']

    def test_standardize_dataframe_coordinate_conversion(self):
        """Test DataFrame standardization with coordinate conversion."""
        test_df = pd.DataFrame({
            'eChr': ['chr1', 'chr2'],
            'eStart': [1, 101],  # 1-based coordinates
            'eEnd': [100, 200]
        })

        standardizer = ColumnStandardizer()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # Suppress coordinate conversion warnings
            result = standardizer.standardize_dataframe(test_df)

        # Check coordinate conversion
        assert result[ColumnStandard.START0].tolist() == [0, 100]  # Converted to 0-based
        assert result[ColumnStandard.END0].tolist() == [100, 200]  # eEnd remains same

    def test_standardize_dataframe_column_conflict(self):
        """Test standardization with column conflicts."""
        test_df = pd.DataFrame({
            'eChr': ['chr1', 'chr2'],
            'chr': ['chrX', 'chrY']  # Conflict: both map to same standard column
        })

        standardizer = ColumnStandardizer()
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = standardizer.standardize_dataframe(test_df)

            # Should warn about conflict
            assert len(w) > 0
            assert "Column conflict" in str(w[0].message)

        # Should prefer existing standard column
        assert result[ColumnStandard.CHR].tolist() == ['chrX', 'chrY']

    def test_standardize_dataframe_unknown_columns_permissive(self):
        """Test standardization with unknown columns in permissive mode."""
        test_df = pd.DataFrame({
            'eChr': ['chr1', 'chr2'],
            'unknown_column': ['x', 'y']
        })

        standardizer = ColumnStandardizer(strict_mode=False)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = standardizer.standardize_dataframe(test_df)

            # Should warn about unknown columns
            assert any("Unknown columns" in str(warning.message) for warning in w)

        # Unknown column should be preserved
        assert 'unknown_column' in result.columns
        assert result['unknown_column'].tolist() == ['x', 'y']

    def test_standardize_dataframe_unknown_columns_strict(self):
        """Test standardization with unknown columns in strict mode."""
        test_df = pd.DataFrame({
            'eChr': ['chr1', 'chr2'],
            'unknown_column': ['x', 'y']
        })

        standardizer = ColumnStandardizer(strict_mode=True)
        with pytest.raises(ValueError) as exc_info:
            standardizer.standardize_dataframe(test_df)

        assert "Unknown columns in strict mode" in str(exc_info.value)

    def test_to_legacy_format_e_prefixed(self):
        """Test conversion to e_prefixed legacy format."""
        test_df = pd.DataFrame({
            ColumnStandard.CHR: ['chr1', 'chr2'],
            ColumnStandard.START0: [100, 200],
            ColumnStandard.READS: ['read1', 'read2']
        })

        standardizer = ColumnStandardizer()
        result = standardizer.to_legacy_format(test_df, 'e_prefixed')

        assert 'eChr' in result.columns
        assert 'eStart0' in result.columns
        assert 'eReads' in result.columns
        assert result['eChr'].tolist() == ['chr1', 'chr2']

    def test_to_legacy_format_camel_case(self):
        """Test conversion to camel_case legacy format."""
        test_df = pd.DataFrame({
            ColumnStandard.READS: ['read1', 'read2'],
            ColumnStandard.COPY_NUMBER: [1, 2],
            ColumnStandard.MATCH_DEGREE: [95.5, 88.2]
        })

        standardizer = ColumnStandardizer()
        result = standardizer.to_legacy_format(test_df, 'camel_case')

        assert 'readName' in result.columns
        assert 'copyNum' in result.columns
        assert 'MatDegree' in result.columns

    def test_to_legacy_format_mixed(self):
        """Test conversion to mixed legacy format."""
        test_df = pd.DataFrame({
            ColumnStandard.CHR: ['chr1', 'chr2'],
            ColumnStandard.READS: ['read1', 'read2'],
            ColumnStandard.COPY_NUMBER: [1, 2]
        })

        standardizer = ColumnStandardizer()
        result = standardizer.to_legacy_format(test_df, 'mixed')

        assert 'eChr' in result.columns
        assert 'readName' in result.columns
        assert 'copyNum' in result.columns

    def test_to_legacy_format_empty(self):
        """Test legacy format conversion with empty DataFrame."""
        standardizer = ColumnStandardizer()
        empty_df = pd.DataFrame()
        result = standardizer.to_legacy_format(empty_df)
        assert result.empty

    def test_validate_schema_success(self):
        """Test successful schema validation."""
        test_df = pd.DataFrame({
            ColumnStandard.ECCDNA_ID: ['U001', 'U002'],
            ColumnStandard.CHR: ['chr1', 'chr2'],
            ColumnStandard.START0: [100, 200]
        })

        standardizer = ColumnStandardizer()
        required_cols = [ColumnStandard.ECCDNA_ID, ColumnStandard.CHR]

        result = standardizer.validate_schema(test_df, required_cols)
        assert result is True

    def test_validate_schema_failure(self):
        """Test schema validation failure."""
        test_df = pd.DataFrame({
            ColumnStandard.ECCDNA_ID: ['U001', 'U002']
        })

        standardizer = ColumnStandardizer()
        required_cols = [ColumnStandard.ECCDNA_ID, ColumnStandard.CHR, ColumnStandard.START0]

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = standardizer.validate_schema(test_df, required_cols)

            assert result is False
            assert len(w) > 0
            assert "Missing required columns" in str(w[0].message)


class TestStandardizeColumnsFunction:
    """Test cases for standardize_columns convenience function."""

    def test_standardize_columns_basic(self):
        """Test basic usage of standardize_columns function."""
        test_df = pd.DataFrame({
            'eChr': ['chr1', 'chr2'],
            'readName': ['read1', 'read2']
        })

        result = standardize_columns(test_df)

        assert ColumnStandard.CHR in result.columns
        assert ColumnStandard.READS in result.columns
        assert 'eChr' not in result.columns

    def test_standardize_columns_with_kwargs(self):
        """Test standardize_columns with keyword arguments."""
        test_df = pd.DataFrame({
            'eChr': ['chr1', 'chr2'],
            'unknown_col': ['x', 'y']
        })

        # Should raise error in strict mode
        with pytest.raises(ValueError):
            standardize_columns(test_df, strict_mode=True)

        # Should work in permissive mode
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = standardize_columns(test_df, strict_mode=False)
            assert 'unknown_col' in result.columns


class TestSchemas:
    """Test cases for SCHEMAS dictionary."""

    def test_schemas_exists(self):
        """Test that SCHEMAS dictionary exists and has expected schemas."""
        assert isinstance(SCHEMAS, dict)
        assert 'uecc' in SCHEMAS
        assert 'mecc' in SCHEMAS
        assert 'cecc' in SCHEMAS
        assert 'confirmed_table' in SCHEMAS

    def test_schema_contents(self):
        """Test that schemas contain expected columns."""
        uecc_schema = SCHEMAS['uecc']
        assert ColumnStandard.ECCDNA_ID in uecc_schema
        assert ColumnStandard.CHR in uecc_schema
        assert ColumnStandard.START0 in uecc_schema
        assert ColumnStandard.COPY_NUMBER in uecc_schema

        mecc_schema = SCHEMAS['mecc']
        assert ColumnStandard.HIT_INDEX in mecc_schema
        assert ColumnStandard.HIT_COUNT in mecc_schema

        cecc_schema = SCHEMAS['cecc']
        assert ColumnStandard.SEG_INDEX in cecc_schema
        assert ColumnStandard.SEG_TOTAL in cecc_schema
        assert ColumnStandard.JUNCTION_ROLE in cecc_schema

    def test_schema_types(self):
        """Test that all schema values are lists of strings."""
        for schema_name, schema in SCHEMAS.items():
            assert isinstance(schema, list), f"Schema {schema_name} should be a list"
            for column in schema:
                assert isinstance(column, str), f"Column {column} in {schema_name} should be a string"


class TestColumnStandardsIntegration:
    """Integration tests for column standards functionality."""

    def test_complete_standardization_workflow(self):
        """Test complete workflow from legacy to standard to legacy."""
        # Create DataFrame with mixed legacy naming
        original_df = pd.DataFrame({
            'eccDNA_id': ['U001', 'U002'],
            'eChr': ['chr1', 'chr2'],
            'eStart': [1, 101],  # 1-based
            'eEnd0': [100, 200],
            'readName': ['read1', 'read2'],
            'copyNum': [1, 2],
            'MatDegree': [95.5, 88.2]
        })

        standardizer = ColumnStandardizer()

        # Step 1: Standardize
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            std_df = standardizer.standardize_dataframe(original_df)

        # Check standardized columns
        assert ColumnStandard.ECCDNA_ID in std_df.columns
        assert ColumnStandard.CHR in std_df.columns
        assert ColumnStandard.START0 in std_df.columns
        assert ColumnStandard.READS in std_df.columns

        # Check coordinate conversion
        assert std_df[ColumnStandard.START0].tolist() == [0, 100]  # Converted to 0-based

        # Step 2: Convert back to legacy
        legacy_df = standardizer.to_legacy_format(std_df, 'mixed')

        # Should have legacy column names
        assert 'eChr' in legacy_df.columns
        assert 'readName' in legacy_df.columns
        assert 'copyNum' in legacy_df.columns

        # Data should be preserved
        assert legacy_df['eChr'].tolist() == ['chr1', 'chr2']
        assert legacy_df['copyNum'].tolist() == [1, 2]

    def test_schema_validation_workflow(self):
        """Test schema validation workflow."""
        # Create DataFrame matching UECC schema
        test_df = pd.DataFrame({
            ColumnStandard.ECCDNA_ID: ['U001', 'U002'],
            ColumnStandard.CHR: ['chr1', 'chr2'],
            ColumnStandard.START0: [100, 200],
            ColumnStandard.END0: [200, 300],
            ColumnStandard.STRAND: ['+', '-'],
            ColumnStandard.LENGTH: [100, 100],
            ColumnStandard.READS: ['read1', 'read2'],
            ColumnStandard.COPY_NUMBER: [1, 2],
            ColumnStandard.MATCH_DEGREE: [95.5, 88.2]
        })

        standardizer = ColumnStandardizer()
        uecc_schema = SCHEMAS['uecc']

        # Should validate successfully
        result = standardizer.validate_schema(test_df, uecc_schema)
        assert result is True

        # Remove required column and test failure
        incomplete_df = test_df.drop(columns=[ColumnStandard.CHR])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = standardizer.validate_schema(incomplete_df, uecc_schema)
            assert result is False