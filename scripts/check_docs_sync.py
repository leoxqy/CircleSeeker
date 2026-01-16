#!/usr/bin/env python3
"""
CircleSeeker Bilingual Documentation Sync Checker

This script checks that bilingual documentation pairs are in sync.
It verifies:
1. Each Chinese document has a corresponding English version
2. Modification times are reasonably close
3. Reports any out-of-sync documents
"""

import os
import sys
from pathlib import Path
from datetime import datetime, timedelta

# Define bilingual document pairs (Chinese -> English)
BILINGUAL_PAIRS = {
    "CLI_Reference.md": "CLI_Reference_en.md",
    "Pipeline_Modules.md": "Pipeline_Modules_en.md",
    "Configuration_Reference.md": "Configuration_Reference_en.md",
    "Output_Format_Reference.md": "Output_Format_Reference_en.md",
    "Validation_Methodology.md": "Validation_Methodology_en.md",
}

# Documents that are intentionally single-language (no sync check needed)
SINGLE_LANGUAGE_DOCS = {
    "UMC_Classification_Model.md",
    "Validation_Report.md",
    "Code_Review_Report_v0.10.0.md",
    "Implementation_Plan.md",
    "Bilingual_Documentation_Guide.md",
}

# Maximum allowed time difference between paired documents (in days)
MAX_TIME_DIFF_DAYS = 7


def get_docs_dir() -> Path:
    """Get the docs directory path."""
    script_dir = Path(__file__).parent
    project_root = script_dir.parent
    return project_root / "docs"


def check_file_exists(docs_dir: Path, filename: str) -> bool:
    """Check if a file exists in the docs directory."""
    return (docs_dir / filename).exists()


def get_mtime(docs_dir: Path, filename: str) -> datetime:
    """Get the modification time of a file."""
    filepath = docs_dir / filename
    return datetime.fromtimestamp(filepath.stat().st_mtime)


def check_bilingual_sync() -> tuple[list[str], list[str], list[str]]:
    """
    Check bilingual documentation sync status.

    Returns:
        Tuple of (missing_english, out_of_sync, synced) document lists
    """
    docs_dir = get_docs_dir()
    missing_english = []
    out_of_sync = []
    synced = []

    for zh_doc, en_doc in BILINGUAL_PAIRS.items():
        # Check if Chinese doc exists
        if not check_file_exists(docs_dir, zh_doc):
            continue  # Skip if Chinese doc doesn't exist

        # Check if English doc exists
        if not check_file_exists(docs_dir, en_doc):
            missing_english.append(zh_doc)
            continue

        # Check modification time difference
        zh_mtime = get_mtime(docs_dir, zh_doc)
        en_mtime = get_mtime(docs_dir, en_doc)
        time_diff = abs(zh_mtime - en_mtime)

        if time_diff > timedelta(days=MAX_TIME_DIFF_DAYS):
            out_of_sync.append(
                f"{zh_doc} <-> {en_doc} "
                f"(diff: {time_diff.days} days)"
            )
        else:
            synced.append(f"{zh_doc} <-> {en_doc}")

    return missing_english, out_of_sync, synced


def main() -> int:
    """Main entry point."""
    print("=" * 60)
    print("CircleSeeker Bilingual Documentation Sync Check")
    print("=" * 60)
    print()

    docs_dir = get_docs_dir()
    if not docs_dir.exists():
        print(f"ERROR: docs directory not found: {docs_dir}")
        return 1

    missing_english, out_of_sync, synced = check_bilingual_sync()

    # Report synced documents
    if synced:
        print("Synced documents:")
        for doc in synced:
            print(f"  [OK] {doc}")
        print()

    # Report out-of-sync documents (warning)
    if out_of_sync:
        print("WARNING: Out-of-sync documents (modification time differs):")
        for doc in out_of_sync:
            print(f"  [WARN] {doc}")
        print()

    # Report missing English documents (error)
    if missing_english:
        print("ERROR: Missing English versions:")
        for doc in missing_english:
            print(f"  [MISSING] {doc} -> {BILINGUAL_PAIRS[doc]}")
        print()

    # Summary
    print("-" * 60)
    print(f"Summary: {len(synced)} synced, {len(out_of_sync)} warnings, "
          f"{len(missing_english)} errors")

    # Return error code if any English docs are missing
    if missing_english:
        print("\nFAILED: Some English documentation is missing!")
        return 1

    if out_of_sync:
        print("\nPASSED with warnings: Some documents may need sync update.")
        return 0

    print("\nPASSED: All bilingual documentation is in sync!")
    return 0


if __name__ == "__main__":
    sys.exit(main())
