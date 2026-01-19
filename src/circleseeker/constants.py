"""Unified constants for CircleSeeker.

This module centralizes commonly used constants that are shared across
multiple modules. This avoids duplication and ensures consistency.
"""

# ================== Quality/Confidence Thresholds ==================
# These thresholds are used for flagging low-confidence calls.
# They do NOT gate calls - only add warning flags to output.

# MAPQ threshold below which reads are flagged as low confidence
MAPQ_LOW_THRESHOLD: int = 20

# Alignment identity threshold (percent) below which reads are flagged
IDENTITY_LOW_THRESHOLD: float = 95.0


# ================== Deduplication Thresholds ==================
# Used for coordinate-based deduplication in ecc_dedup

# Sequence identity threshold for considering sequences as duplicates (0-1)
DEDUP_IDENTITY_THRESHOLD: float = 0.99

# Position tolerance in base pairs for coordinate overlap detection
DEDUP_POSITION_TOLERANCE: int = 10


# ================== Output Constants ==================
# Default decimal precision for floating point values in output
OUTPUT_DECIMAL_PRECISION: int = 2
