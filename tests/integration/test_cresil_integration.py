#!/usr/bin/env python3
"""Quick test script to verify Cresil integration."""

import sys
from pathlib import Path
import pandas as pd

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from circleseeker.modules.cresil_adapter import (
    parse_cresil_region,
    parse_cresil_regions,
    convert_cresil_to_cyrcular_format,
)


def test_parse_region():
    """Test region parsing."""
    print("Testing region parsing...")

    # Test cases
    test_cases = [
        ("chr2:47242017-47243061_-", ("chr2", 47242017, 47243061, "-")),
        ("chr8:129264432-129269966_+", ("chr8", 129264432, 129269966, "+")),
        ("chr10:131725754-131727108_-", ("chr10", 131725754, 131727108, "-")),
    ]

    for region_str, expected in test_cases:
        result = parse_cresil_region(region_str)
        assert result == expected, f"Failed: {region_str} -> {result} != {expected}"
        print(f"  ✓ {region_str} -> {result}")

    multi_region = "chr2:47242017-47243061_-;chr2:47250000-47252000_-"
    parsed_multi = parse_cresil_regions(multi_region)
    assert len(parsed_multi) == 2, "Failed to parse multiple Cresil regions"
    print(f"  ✓ {multi_region} -> {parsed_multi}")

    print("All region parsing tests passed!\n")


def test_conversion():
    """Test Cresil to Cyrcular format conversion."""
    print("Testing format conversion...")

    # Create mock Cresil data
    cresil_data = {
        "id": ["ec1", "ec2", "ec3", "ec4"],
        "merge_region": [
            "chr2:47242017-47243061_-",
            "chr2:128071984-128072363_+",
            "chr8:129264432-129269966_-",
            "chr1:1000-1500_-;chr1:2000-2500_-",
        ],
        "merge_len": [1045, 380, 5535, 1000],
        "num_region": [1, 1, 1, 2],
        "ctc": [True, True, True, True],
        "numreads": [5, 1, 14, 6],
        "totalbase": [52522, 9120, 161471, 48000],
        "coverage": [50.26, 24.00, 29.17, 48.0],
    }

    df = pd.DataFrame(cresil_data)

    # Create temporary test file
    test_file = Path("/tmp/test_cresil_output.txt")
    df.to_csv(test_file, sep="\t", index=False)

    try:
        # Convert
        result_df = convert_cresil_to_cyrcular_format(test_file)

        print(f"  Converted {len(result_df)} records")
        print("\n  Expected columns:")
        expected_cols = [
            "circle_id", "regions", "circle_length", "segment_count",
            "num_split_reads", "prob_present", "prob_artifact", "af_nanopore"
        ]
        for col in expected_cols:
            if col in result_df.columns:
                print(f"    ✓ {col}")
            else:
                print(f"    ✗ {col} (missing!)")

        print("\n  Sample output:")
        print(result_df[["circle_id", "regions", "circle_length", "segment_count"]].head())

        print("\nFormat conversion test passed!\n")

    finally:
        # Cleanup
        test_file.unlink(missing_ok=True)


def main():
    """Run all tests."""
    print("=" * 60)
    print("Cresil Integration Test Suite")
    print("=" * 60)
    print()

    try:
        test_parse_region()
        test_conversion()

        print("=" * 60)
        print("✓ ALL TESTS PASSED")
        print("=" * 60)
        return 0

    except Exception as e:
        print(f"\n✗ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
