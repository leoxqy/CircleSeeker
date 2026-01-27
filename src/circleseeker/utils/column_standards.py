"""
CircleSeeker Column Naming Standards and Standardization Module

This module defines the unified column naming standard and provides tools for
converting between different naming conventions used throughout CircleSeeker.

Naming Philosophy:
- Simple, clear names without unnecessary prefixes
- Consistent snake_case for multi-word fields
- Preserve important domain-specific conventions (e.g., eccDNA_id)
- 0-based coordinates clearly indicated
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Optional
import pandas as pd
import warnings


@dataclass
class ColumnStandard:
    """Defines the unified column naming standard for CircleSeeker."""

    # === Core Identity Fields ===
    ECCDNA_ID = "eccDNA_id"  # Main identifier (preserve existing convention)
    ORIG_ECCDNA_ID = "orig_eccDNA_id"  # Original ID before renaming
    CLUSTER_ID = "cluster_id"  # Cluster assignment

    # === Genomic Coordinates (0-based) ===
    CHR = "chr"  # Chromosome name
    START0 = "start0"  # 0-based start coordinate
    END0 = "end0"  # 0-based end coordinate
    STRAND = "strand"  # Strand (+/-)
    LENGTH = "length"  # Sequence/region length

    # === Sequence and Read Information ===
    READS = "reads"  # Read names (semicolon separated)
    SEQUENCE = "sequence"  # DNA sequence
    READ_COUNT = "read_count"  # Number of reads

    # === Quantitative Measures ===
    COPY_NUMBER = "copy_number"  # Copy number estimate
    MATCH_DEGREE = "match_degree"  # Match percentage (0-100)
    GAP_PERCENTAGE = "gap_percentage"  # Gap percentage (0-100), complement of match_degree
    REPEAT_NUMBER = "repeat_number"  # Number of repeats

    # === Confidence / Evidence (0-1 unless otherwise noted) ===
    CONFIDENCE_SCORE = "confidence_score"
    MAPQ_BEST = "mapq_best"  # max MAPQ among supporting alignments
    MAPQ_MIN = "mapq_min"  # min MAPQ among supporting alignments
    IDENTITY_BEST = "identity_best"  # max identity (%) among supporting alignments
    IDENTITY_MIN = "identity_min"  # min identity (%) among supporting alignments
    QUERY_COV_BEST = "query_cov_best"  # best-locus (or best-chain) ring coverage
    QUERY_COV_2ND = "query_cov_2nd"  # 2nd-best locus/chain ring coverage
    LOW_MAPQ = "low_mapq"  # boolean flag
    LOW_IDENTITY = "low_identity"  # boolean flag

    # === Type and Classification ===
    ECCDNA_TYPE = "eccdna_type"  # UeccDNA, MeccDNA, CeccDNA, etc.
    STATE = "state"  # Confirmed, Inferred, etc.

    # === Segmentation (for CeccDNA) ===
    SEG_INDEX = "seg_index"  # Segment index (1-based)
    SEG_TOTAL = "seg_total"  # Total segments
    JUNCTION_ROLE = "junction_role"  # head, body, tail

    # === Aggregation Metadata ===
    NUM_MERGED = "num_merged"  # Number of entries merged
    MERGED_FROM_IDS = "merged_from_ids"  # IDs that were merged
    HIT_COUNT = "hit_count"  # Number of alignment hits
    HIT_INDEX = "hit_index"  # Hit index for multi-hit entries


class OutputColumns:
    """Column names for user-facing output files (backward compatible format).

    These names are used in merged_output.csv and other user-facing outputs
    to maintain backward compatibility with existing downstream tools.
    """

    # Core columns (same as internal)
    ECCDNA_ID = "eccDNA_id"
    ORIGINAL_ID = "original_id"

    # PascalCase columns for user output compatibility
    REGIONS = "Regions"  # Multiple regions joined by semicolon
    STRAND = "Strand"
    LENGTH = "Length"
    ECCDNA_TYPE = "eccDNA_type"  # Keep camelCase for backward compat
    STATE = "State"
    SEG_TOTAL = "Seg_total"
    HIT_COUNT = "Hit_count"

    # Evidence columns (same as internal)
    CONFIDENCE_SCORE = ColumnStandard.CONFIDENCE_SCORE
    QUERY_COV_BEST = ColumnStandard.QUERY_COV_BEST
    QUERY_COV_2ND = ColumnStandard.QUERY_COV_2ND
    MAPQ_BEST = ColumnStandard.MAPQ_BEST
    IDENTITY_BEST = ColumnStandard.IDENTITY_BEST
    LOW_MAPQ = ColumnStandard.LOW_MAPQ
    LOW_IDENTITY = ColumnStandard.LOW_IDENTITY


# Legacy mapping dictionaries for backward compatibility
LEGACY_MAPPINGS = {
    # E-prefixed legacy format
    "eChr": ColumnStandard.CHR,
    "eStart": ColumnStandard.START0,  # Note: eStart is 1-based, need conversion
    "eStart0": ColumnStandard.START0,
    "eEnd": ColumnStandard.END0,  # Note: eEnd is 1-based, need conversion
    "eEnd0": ColumnStandard.END0,
    "eStrand": ColumnStandard.STRAND,
    "eLength": ColumnStandard.LENGTH,
    "eReads": ColumnStandard.READS,
    "eRepeatNum": ColumnStandard.REPEAT_NUMBER,
    "eClass": ColumnStandard.ECCDNA_TYPE,
    # CamelCase variants
    "readName": ColumnStandard.READS,
    "copyNum": ColumnStandard.COPY_NUMBER,
    "MatDegree": ColumnStandard.MATCH_DEGREE,
    # Snake_case variants
    "read_name": ColumnStandard.READS,
    "start_0based": ColumnStandard.START0,
    "end_0based": ColumnStandard.END0,
    "copy_number": ColumnStandard.COPY_NUMBER,
    "match_degree": ColumnStandard.MATCH_DEGREE,
    "repeat_number": ColumnStandard.REPEAT_NUMBER,
    "eccdna_id": ColumnStandard.ECCDNA_ID,
    # Alternative naming
    "chr_name": ColumnStandard.CHR,
    "chromosome": ColumnStandard.CHR,
    "reads_count": ColumnStandard.READ_COUNT,
    "seq": ColumnStandard.SEQUENCE,
    "eSeq": ColumnStandard.SEQUENCE,
    "data_type": ColumnStandard.ECCDNA_TYPE,
    "State": ColumnStandard.STATE,
    "Length": ColumnStandard.LENGTH,
    "Hit_count": ColumnStandard.HIT_COUNT,
    "Seg_total": ColumnStandard.SEG_TOTAL,
    # BLAST-specific
    "subject_id": ColumnStandard.CHR,
    "s_start": ColumnStandard.START0,  # Note: BLAST coords need conversion
    "s_end": ColumnStandard.END0,
    "q_start": "q_start",  # Keep query coords separate
    "q_end": "q_end",
    "query_id": ColumnStandard.ECCDNA_ID,
    "alignment_length": "alignment_length",
}

# Coordinate conversion requirements
COORDINATE_CONVERSIONS = {
    "eStart": lambda x: x - 1,  # 1-based to 0-based
    "eEnd": lambda x: x,  # 1-based inclusive end -> 0-based end (half-open)
    "s_start": lambda x: x - 1,  # BLAST-like 1-based -> 0-based
    "s_end": lambda x: x,  # BLAST-like end -> 0-based end (half-open)
}


class ColumnStandardizer:
    """Central column standardization tool for CircleSeeker."""

    def __init__(self, strict_mode: bool = False):
        """
        Initialize standardizer.

        Args:
            strict_mode: If True, raise errors for unknown columns
        """
        self.strict_mode = strict_mode

    def standardize_dataframe(
        self, df: pd.DataFrame, source_format: Optional[str] = None
    ) -> pd.DataFrame:
        """
        Convert DataFrame to standard column format.

        Args:
            df: Input DataFrame
            source_format: Hint about source format ('blast', 'legacy', etc.)

        Returns:
            DataFrame with standardized column names
        """
        if df.empty:
            return df

        df = df.copy()

        # Track conversions for logging
        conversions = {}

        # Apply column name mappings
        for old_col, new_col in LEGACY_MAPPINGS.items():
            if old_col in df.columns:
                if new_col in df.columns and old_col != new_col:
                    # Handle conflicts - prefer existing standard column
                    warnings.warn(
                        f"Column conflict: {old_col} maps to {new_col} but {new_col} already exists"
                    )
                    continue

                df = df.rename(columns={old_col: new_col})
                conversions[old_col] = new_col

                # Apply coordinate conversions if needed
                if old_col in COORDINATE_CONVERSIONS:
                    try:
                        df[new_col] = df[new_col].apply(COORDINATE_CONVERSIONS[old_col])
                        conversions[f"{old_col}_converted"] = "Applied coordinate conversion"
                    except (KeyError, TypeError, ValueError) as e:
                        warnings.warn(f"Failed to convert coordinates for {old_col}: {e}")

        # Handle unknown columns
        standard_columns = self._get_all_standard_columns()
        unknown_columns = [col for col in df.columns if col not in standard_columns]

        if unknown_columns:
            if self.strict_mode:
                raise ValueError(f"Unknown columns in strict mode: {unknown_columns}")
            else:
                warnings.warn(f"Unknown columns (keeping as-is): {unknown_columns}")

        # Log conversions if any were made
        if conversions:
            print(f"Column standardization applied: {len(conversions)} conversions")
            for old, new in conversions.items():
                if not old.endswith("_converted"):
                    print(f"  {old} → {new}")

        return df

    def to_legacy_format(self, df: pd.DataFrame, target_format: str = "mixed") -> pd.DataFrame:
        """
        Convert from standard format to legacy format for backward compatibility.

        Args:
            df: DataFrame with standard column names
            target_format: Target legacy format ('e_prefixed', 'camel_case', 'mixed')

        Returns:
            DataFrame with legacy column names
        """
        if df.empty:
            return df

        df = df.copy()

        # Reverse mapping based on target format
        if target_format == "e_prefixed":
            reverse_mapping = {
                ColumnStandard.CHR: "eChr",
                ColumnStandard.START0: "eStart0",
                ColumnStandard.END0: "eEnd0",
                ColumnStandard.STRAND: "eStrand",
                ColumnStandard.LENGTH: "eLength",
                ColumnStandard.READS: "eReads",
                ColumnStandard.REPEAT_NUMBER: "eRepeatNum",
                ColumnStandard.ECCDNA_TYPE: "eClass",
            }
        elif target_format == "camel_case":
            reverse_mapping = {
                ColumnStandard.READS: "readName",
                ColumnStandard.COPY_NUMBER: "copyNum",
                ColumnStandard.MATCH_DEGREE: "MatDegree",
            }
        else:  # mixed - preserve most common legacy names
            reverse_mapping = {
                ColumnStandard.CHR: "eChr",
                ColumnStandard.START0: "eStart0",
                ColumnStandard.END0: "eEnd0",
                ColumnStandard.STRAND: "eStrand",
                ColumnStandard.LENGTH: "eLength",
                ColumnStandard.READS: "readName",
                ColumnStandard.COPY_NUMBER: "copyNum",
                ColumnStandard.MATCH_DEGREE: "MatDegree",
                ColumnStandard.REPEAT_NUMBER: "eRepeatNum",
                ColumnStandard.ECCDNA_TYPE: "eClass",
            }

        # Apply reverse mapping
        for std_col, legacy_col in reverse_mapping.items():
            if std_col in df.columns:
                df = df.rename(columns={std_col: legacy_col})

        return df

    def validate_schema(self, df: pd.DataFrame, required_columns: list[str]) -> bool:
        """
        Validate that DataFrame has required standard columns.

        Args:
            df: DataFrame to validate
            required_columns: List of required standard column names

        Returns:
            True if all required columns present
        """
        missing = [col for col in required_columns if col not in df.columns]
        if missing:
            warnings.warn(f"Missing required columns: {missing}")
            return False
        return True

    @staticmethod
    def _get_all_standard_columns() -> list[str]:
        """Get list of all defined standard column names."""
        return [
            getattr(ColumnStandard, attr)
            for attr in dir(ColumnStandard)
            if not attr.startswith("_") and isinstance(getattr(ColumnStandard, attr), str)
        ]


# Convenience function for quick standardization
def standardize_columns(df: pd.DataFrame, **kwargs) -> pd.DataFrame:
    """Quick function to standardize DataFrame columns."""
    standardizer = ColumnStandardizer(**kwargs)
    return standardizer.standardize_dataframe(df)


# Schema definitions for different data types
SCHEMAS = {
    "uecc": [
        ColumnStandard.ECCDNA_ID,
        ColumnStandard.CHR,
        ColumnStandard.START0,
        ColumnStandard.END0,
        ColumnStandard.STRAND,
        ColumnStandard.LENGTH,
        ColumnStandard.READS,
        ColumnStandard.COPY_NUMBER,
        ColumnStandard.MATCH_DEGREE,
    ],
    "mecc": [
        ColumnStandard.ECCDNA_ID,
        ColumnStandard.CHR,
        ColumnStandard.START0,
        ColumnStandard.END0,
        ColumnStandard.STRAND,
        ColumnStandard.LENGTH,
        ColumnStandard.READS,
        ColumnStandard.COPY_NUMBER,
        ColumnStandard.HIT_INDEX,
        ColumnStandard.HIT_COUNT,
    ],
    "cecc": [
        ColumnStandard.ECCDNA_ID,
        ColumnStandard.CHR,
        ColumnStandard.START0,
        ColumnStandard.END0,
        ColumnStandard.STRAND,
        ColumnStandard.LENGTH,
        ColumnStandard.SEG_INDEX,
        ColumnStandard.SEG_TOTAL,
        ColumnStandard.JUNCTION_ROLE,
        ColumnStandard.READS,
        ColumnStandard.COPY_NUMBER,
    ],
    "confirmed_table": [
        ColumnStandard.ECCDNA_ID,
        "regions",  # Multiple regions joined by semicolon
        ColumnStandard.LENGTH,
        ColumnStandard.ECCDNA_TYPE,
        ColumnStandard.STATE,
        ColumnStandard.SEG_TOTAL,
        ColumnStandard.HIT_COUNT,
    ],
}


# Output column names for user-facing files (may differ from internal standard)
# These are the column names used in merged_output.csv and other user-facing outputs
OUTPUT_COLUMNS = {
    "eccDNA_id": ColumnStandard.ECCDNA_ID,
    "Regions": "regions",  # PascalCase for user output compatibility
    "Strand": ColumnStandard.STRAND,
    "Length": ColumnStandard.LENGTH,
    "eccDNA_type": ColumnStandard.ECCDNA_TYPE,
    "State": ColumnStandard.STATE,
    "Seg_total": ColumnStandard.SEG_TOTAL,
    "Hit_count": ColumnStandard.HIT_COUNT,
}


if __name__ == "__main__":
    # Example usage and testing
    # Test DataFrame with mixed legacy naming
    test_df = pd.DataFrame(
        {
            "eccDNA_id": ["U001", "U002"],
            "eChr": ["chr1", "chr2"],
            "eStart": [101, 201],  # 1-based
            "eEnd0": [200, 300],
            "readName": ["readA", "readB"],
            "copyNum": [1, 2],
            "unknown_col": ["x", "y"],
        }
    )

    print("Original DataFrame:")
    print(test_df)
    print("\nColumns:", list(test_df.columns))

    # Standardize
    standardizer = ColumnStandardizer()
    std_df = standardizer.standardize_dataframe(test_df)

    print("\nStandardized DataFrame:")
    print(std_df)
    print("\nColumns:", list(std_df.columns))

    # Convert back to legacy
    legacy_df = standardizer.to_legacy_format(std_df, "mixed")
    print("\nLegacy format:")
    print(legacy_df)
    print("\nColumns:", list(legacy_df.columns))
