"""
ECC Dedup - eccDNA Cluster Deduplication Module with Confirmed Tables Generation

This module performs independent deduplication for each eccDNA type using their respective
CD-HIT clusters. No cross-type priority judgments - each type is processed independently.

Key features:
- Each type (Uecc/Mecc/Cecc) uses its own CD-HIT cluster mapping (.clstr or
  exported CSV) for deduplication
- Uecc: keeps 1 representative row per cluster, aggregates metadata
- Mecc/Cecc: keeps ALL rows of representative ID, aggregates metadata from cluster members
- Outputs 11 standardized files maintaining compatibility with existing pipeline
- Generates confirmed tables with standardized formats for downstream analysis
- Optional file organization into type-specific folders

Integrated functionality from build_core_tables module.
Migrated from step8_harmonizer.py to the CircleSeeker architecture.
"""

from __future__ import annotations

import csv
import sys
import logging
from circleseeker.utils.logging import get_logger
import re
import shutil
from pathlib import Path
from typing import Optional
import pandas as pd
from circleseeker.utils.column_standards import ColumnStandard

# Increase CSV field size limit to handle large sequence fields
csv.field_size_limit(sys.maxsize)


# ================== Configuration ==================
FIXED_FRONT_COLUMNS = [
    ColumnStandard.ECCDNA_ID,
    ColumnStandard.CHR,
    ColumnStandard.START0,
    ColumnStandard.END0,
    ColumnStandard.STRAND,
    ColumnStandard.LENGTH,
    ColumnStandard.MATCH_DEGREE,
    ColumnStandard.COPY_NUMBER,
]

# Mapping between legacy mixed-case column names and new snake_case names.
LEGACY_TO_SNAKE_CASE = {
    "eccDNA_id": "eccdna_id",
    "eChr": "chr",
    "eStart0": "start_0based",
    "eEnd0": "end_0based",
    "eStrand": "strand",
    "eLength": "length",
    "MatDegree": "match_degree",
    "copyNum": "copy_number",
    "eReads": "reads",
    "readName": "read_name",
    "eRepeatNum": "repeat_number",
    "orig_eccDNA_id": "orig_eccdna_id",
    "data_type": "eccdna_type",
    "eClass": "eccdna_type",
    "State": "state",
}

STRAND_POSITIVE = {"+", "plus", "positive", "forward", "f", "pos"}
STRAND_NEGATIVE = {"-", "minus", "negative", "reverse", "rev", "r", "neg"}

# Confidence/evidence heuristics (used for flags only; does not gate calls).
CONF_MAPQ_LOW_THRESHOLD = 20
CONF_IDENTITY_LOW_THRESHOLD = 95.0


def _normalize_strand(value: object) -> Optional[str]:
    """Normalise strand annotations to '+'/'-' where possible."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    text = str(value).strip()
    if not text:
        return None
    lowered = text.lower()
    if lowered in STRAND_POSITIVE:
        return "+"
    if lowered in STRAND_NEGATIVE:
        return "-"
    if text in {"+", "-"}:
        return text
    return text


# ================== Helper Functions ==================
def to_numeric_safe(series, default=0):
    """Safely convert series to numeric, filling NaN with default."""
    try:
        return pd.to_numeric(series, errors="coerce").fillna(default)
    except Exception:
        return pd.Series([default] * len(series))


def format_two_decimals(value) -> str:
    """Format value to 2 decimal places."""
    try:
        if pd.isna(value):
            return ""
        return f"{float(value):.2f}"
    except Exception:
        return ""


def concat_unique_semicolon(series: pd.Series) -> str:
    """Concatenate unique values with semicolon."""
    unique_vals = series.dropna().astype(str).unique()
    return ";".join(unique_vals) if len(unique_vals) else ""


def natural_sort_eccdna_id(df: pd.DataFrame) -> pd.DataFrame:
    """Sort DataFrame by eccDNA_id using natural sorting (e.g., Uecc1, Uecc2, Uecc10, Uecc11)."""
    if df.empty or "eccDNA_id" not in df.columns:
        return df

    # Extract type prefix and numeric part for natural sorting
    df = df.copy()
    df["_type_prefix"] = df["eccDNA_id"].str.extract(r"([A-Za-z]+)", expand=False)
    df["_id_number"] = to_numeric_safe(
        df["eccDNA_id"].str.extract(r"(\d+)", expand=False),
        default=-1,
    )

    # Sort by type prefix first, then by numeric ID
    df = df.sort_values(["_type_prefix", "_id_number"]).drop(columns=["_type_prefix", "_id_number"])

    return df


def merge_read_lists(series: pd.Series) -> str:
    """Merge semicolon-separated read lists, removing duplicates."""
    all_reads = set()
    for reads_str in series.dropna().astype(str):
        if reads_str:
            tokens = [t.strip() for t in reads_str.split(";") if t.strip()]
            all_reads.update(tokens)
    return ";".join(sorted(all_reads)) if all_reads else ""


def majority_vote(series: pd.Series) -> Optional[str]:
    """Return the most common value in series."""
    clean_series = series.dropna().astype(str)
    if clean_series.empty:
        return None
    return str(clean_series.value_counts().idxmax())


def count_reads_from_string(s: str) -> int:
    """Count unique reads from semicolon-separated string."""
    if pd.isna(s) or not isinstance(s, str) or not s.strip():
        return 0
    return len({x.strip() for x in s.split(";") if x.strip()})


def ensure_int64_column(df: pd.DataFrame, col: str) -> None:
    """Ensure column is Int64 type."""
    if col not in df.columns:
        df[col] = pd.Series([pd.NA] * len(df), dtype="Int64")
    else:
        try:
            df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")
        except Exception:
            df[col] = pd.Series([pd.NA] * len(df), dtype="Int64")


# ================== Normalisation Helpers ==================
def clamp_match_degree(series: pd.Series) -> pd.Series:
    """Convert to numeric and clamp to the [0, 100] interval."""
    numeric = pd.to_numeric(series, errors="coerce")
    numeric = numeric.clip(lower=0, upper=100)
    return numeric


# ================== Confirmed Tables Generation Functions ==================
def read_csv_auto(path: str | Path, sep: str | None = None) -> pd.DataFrame:
    """Read CSV with automatic delimiter detection."""
    path = str(path)
    if sep is None:
        try:
            # Try comma-separated first
            df = pd.read_csv(path)
            # Check if it parsed correctly - if only one column and contains tabs, it's likely TSV
            if len(df.columns) == 1 and "\t" in str(df.columns[0]):
                df = pd.read_csv(path, sep="\t")
            return df
        except Exception:
            return pd.read_csv(path, sep="\t")
    else:
        return pd.read_csv(path, sep=sep)


def region_str(chr_: str, start0: int, end: int) -> str:
    """Format genomic region string."""
    return f"{chr_}:{int(start0)}-{int(end)}"


def normalize_column_names(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize column names to expected format.

    Maps old column names to new standard:
    - eChr -> chr
    - eStart0 -> start0
    - eEnd0 -> end0
    - eLength -> length
    """
    column_mapping = {
        "eChr": "chr",
        "eStart0": "start0",
        "eEnd0": "end0",
        "eLength": "length",
        "eClass": "eccdna_type",
        "State": "state",
        "Length": "length",
        "Regions": "regions",
    }

    # Create a copy to avoid modifying original
    df = df.copy()

    # Apply mapping
    for old_name, new_name in column_mapping.items():
        if old_name in df.columns and new_name not in df.columns:
            df.rename(columns={old_name: new_name}, inplace=True)

    return df


def build_u_confirmed_table(u_df: pd.DataFrame) -> pd.DataFrame:
    """Build UeccDNA confirmed table from segments.

    Args:
        u_df: DataFrame with columns: eccDNA_id, chr, start0, end0, strand, length

    Returns:
        DataFrame with standardized confirmed columns including strand
    """
    # Normalize column names
    u_df = normalize_column_names(u_df)

    req = [
        "eccDNA_id",
        ColumnStandard.CHR,
        ColumnStandard.START0,
        ColumnStandard.END0,
        ColumnStandard.LENGTH,
    ]
    for c in req:
        if c not in u_df.columns:
            raise ValueError(f"U confirmed table missing column: {c}")

    df = u_df.copy()

    if "eccDNA_id" not in df.columns:
        for alt in ("eccdna_id", "ecc_id", "id"):
            if alt in df.columns:
                df["eccDNA_id"] = df[alt]
                if alt != "eccDNA_id":
                    df = df.drop(columns=[alt])
                break

    # Handle empty DataFrame case
    if len(df) == 0:
        return pd.DataFrame(
            columns=[
                "eccDNA_id",
                "Regions",
                "Strand",
                "Length",
                "eccDNA_type",
                "State",
                "Seg_total",
                "Hit_count",
            ]
        )

    df["regions"] = df.apply(
        lambda r: region_str(
            r[ColumnStandard.CHR], r[ColumnStandard.START0], r[ColumnStandard.END0]
        ),
        axis=1,
    )

    # Group by eccDNA_id and aggregate, maintaining order
    agg_dict = {"regions": "first", ColumnStandard.LENGTH: "first"}

    # Include strand if available
    if ColumnStandard.STRAND in df.columns:
        agg_dict[ColumnStandard.STRAND] = "first"

    # Optional confidence/evidence columns (max over merged members).
    for col in (
        ColumnStandard.CONFIDENCE_SCORE,
        ColumnStandard.QUERY_COV_BEST,
        ColumnStandard.QUERY_COV_2ND,
        ColumnStandard.MAPQ_BEST,
        ColumnStandard.IDENTITY_BEST,
        ColumnStandard.LOW_MAPQ,
        ColumnStandard.LOW_IDENTITY,
    ):
        if col in df.columns:
            agg_dict[col] = "max"

    g = df.groupby("eccDNA_id", as_index=False, sort=False).agg(agg_dict)

    # Add standard fields
    g = g.rename(columns={"regions": "Regions", ColumnStandard.LENGTH: "Length"})
    if ColumnStandard.STRAND in g.columns:
        g = g.rename(columns={ColumnStandard.STRAND: "Strand"})
    else:
        g["Strand"] = "."  # Default strand if not available

    g["eccDNA_type"] = "UeccDNA"
    g["State"] = "Confirmed"
    g["Seg_total"] = 1
    g["Hit_count"] = 1

    # Apply natural sorting by eccDNA_id
    g = natural_sort_eccdna_id(g)

    cols = [
        "eccDNA_id",
        "Regions",
        "Strand",
        "Length",
        "eccDNA_type",
        "State",
        "Seg_total",
        "Hit_count",
    ]
    for col in (
        ColumnStandard.CONFIDENCE_SCORE,
        ColumnStandard.QUERY_COV_BEST,
        ColumnStandard.QUERY_COV_2ND,
        ColumnStandard.MAPQ_BEST,
        ColumnStandard.IDENTITY_BEST,
        ColumnStandard.LOW_MAPQ,
        ColumnStandard.LOW_IDENTITY,
    ):
        if col in g.columns:
            cols.append(col)
    return g[cols]


def build_m_confirmed_table(m_df: pd.DataFrame) -> pd.DataFrame:
    """Build MeccDNA confirmed table from segments.

    Args:
        m_df: DataFrame with columns: eccDNA_id, chr, start0, end0, strand, length
        Optional: Hit_count or hit_count

    Returns:
        DataFrame with standardized confirmed columns including strand
    """
    # Normalize column names
    m_df = normalize_column_names(m_df)

    req = ["eccDNA_id", "chr", "start0", "end0", "length"]
    for c in req:
        if c not in m_df.columns:
            raise ValueError(f"Mecc confirmed table missing column: {c}")

    df = m_df.copy()

    if "eccDNA_id" not in df.columns:
        for alt in ("eccdna_id", "ecc_id", "id"):
            if alt in df.columns:
                df["eccDNA_id"] = df[alt]
                if alt != "eccDNA_id":
                    df = df.drop(columns=[alt])
                break

    # Handle empty DataFrame case
    if len(df) == 0:
        return pd.DataFrame(
            columns=[
                "eccDNA_id",
                "Regions",
                "Strand",
                "Length",
                "eccDNA_type",
                "State",
                "Seg_total",
                "Hit_count",
            ]
        )

    # Create region strings
    df["region"] = df.apply(lambda r: region_str(r["chr"], r["start0"], r["end0"]), axis=1)

    # Sort by chromosome and position for consistent ordering
    df["__chr"] = df["chr"].astype(str)
    df["__start"] = df["start0"].astype(int)
    df = df.sort_values(["eccDNA_id", "__chr", "__start"])

    # Calculate hit count for each eccDNA_id (number of sites)
    hit_counts = df.groupby("eccDNA_id").size().rename("Hit_count")
    df = df.merge(hit_counts, on="eccDNA_id", how="left")

    # Aggregate by eccDNA_id, maintaining order
    agg: dict[str, object] = {
        "region": lambda x: "|".join(list(x)),  # Use | separator for MeccDNA
        "length": "first",
        "Hit_count": "first",  # Use the calculated hit count
    }

    # Include strand if available
    if ColumnStandard.STRAND in df.columns:
        agg[ColumnStandard.STRAND] = lambda x: "|".join(list(x))  # Join multiple strands with |

    for col in (
        ColumnStandard.CONFIDENCE_SCORE,
        ColumnStandard.QUERY_COV_BEST,
        ColumnStandard.QUERY_COV_2ND,
        ColumnStandard.MAPQ_BEST,
        ColumnStandard.IDENTITY_BEST,
        ColumnStandard.LOW_MAPQ,
        ColumnStandard.LOW_IDENTITY,
    ):
        if col in df.columns:
            agg[col] = "max"

    out = df.groupby("eccDNA_id", sort=False).agg(agg).reset_index()

    # Add standard fields
    out = out.rename(columns={"region": "Regions", "length": "Length"})
    if ColumnStandard.STRAND in out.columns:
        out = out.rename(columns={ColumnStandard.STRAND: "Strand"})
    else:
        out["Strand"] = "."  # Default strand if not available

    out["eccDNA_type"] = "MeccDNA"
    out["State"] = "Confirmed"
    out["Seg_total"] = 1  # MeccDNA is counted as single unit despite multiple sites

    # Apply natural sorting by eccDNA_id
    out = natural_sort_eccdna_id(out)

    cols = [
        "eccDNA_id",
        "Regions",
        "Strand",
        "Length",
        "eccDNA_type",
        "State",
        "Seg_total",
        "Hit_count",
    ]
    for col in (
        ColumnStandard.CONFIDENCE_SCORE,
        ColumnStandard.QUERY_COV_BEST,
        ColumnStandard.QUERY_COV_2ND,
        ColumnStandard.MAPQ_BEST,
        ColumnStandard.IDENTITY_BEST,
        ColumnStandard.LOW_MAPQ,
        ColumnStandard.LOW_IDENTITY,
    ):
        if col in out.columns:
            cols.append(col)
    return out[cols]


def build_c_confirmed_table(c_df: pd.DataFrame) -> pd.DataFrame:
    """Build CeccDNA confirmed table from segments.

    Args:
        c_df: DataFrame with columns: eccDNA_id, chr, start0, end0, strand,
            length, seg_total, seg_index

    Returns:
        DataFrame with standardized confirmed columns including strand
    """
    # Normalize column names
    c_df = normalize_column_names(c_df)

    req = ["eccDNA_id", "chr", "start0", "end0", "length"]
    for c in req:
        if c not in c_df.columns:
            raise ValueError(f"Cecc segments confirmed table missing column: {c}")

    # Handle seg_total and seg_index - these might be computed from existing data
    if "seg_total" not in c_df.columns:
        c_df["seg_total"] = c_df.groupby("eccDNA_id")["eccDNA_id"].transform("size")
    if "seg_index" not in c_df.columns:
        c_df["seg_index"] = c_df.groupby("eccDNA_id").cumcount() + 1

    df = c_df.copy()

    if "eccDNA_id" not in df.columns:
        for alt in ("eccdna_id",):
            if alt in df.columns:
                df["eccDNA_id"] = df[alt]
                if alt != "eccDNA_id":
                    df = df.drop(columns=[alt])
                break

    # Handle empty DataFrame case
    if len(df) == 0:
        return pd.DataFrame(
            columns=[
                "eccDNA_id",
                "Regions",
                "Strand",
                "Length",
                "eccDNA_type",
                "State",
                "Seg_total",
                "Hit_count",
            ]
        )

    # Create region strings
    df["region"] = df.apply(lambda r: region_str(r["chr"], r["start0"], r["end0"]), axis=1)

    # Sort by segment index to maintain canonical order
    df = df.sort_values(["eccDNA_id", "seg_index"])

    # Aggregate by eccDNA_id
    agg_dict = {
        "region": lambda x: ";".join(list(x)),  # Use ; separator for CeccDNA segments
        "length": "first",
        "seg_total": "max",
    }

    # Include strand if available
    if ColumnStandard.STRAND in df.columns:
        agg_dict[ColumnStandard.STRAND] = lambda x: ";".join(
            list(x)
        )  # Join multiple strands with ;

    for col in (
        ColumnStandard.CONFIDENCE_SCORE,
        ColumnStandard.QUERY_COV_BEST,
        ColumnStandard.QUERY_COV_2ND,
        ColumnStandard.MAPQ_BEST,
        ColumnStandard.IDENTITY_BEST,
        ColumnStandard.LOW_MAPQ,
        ColumnStandard.LOW_IDENTITY,
    ):
        if col in df.columns:
            agg_dict[col] = "max"

    out = df.groupby("eccDNA_id", sort=False).agg(agg_dict).reset_index()

    # Add standard fields
    out = out.rename(columns={"region": "Regions", "length": "Length", "seg_total": "Seg_total"})
    if ColumnStandard.STRAND in out.columns:
        out = out.rename(columns={ColumnStandard.STRAND: "Strand"})
    else:
        out["Strand"] = "."  # Default strand if not available

    out["eccDNA_type"] = "CeccDNA"
    out["State"] = "Confirmed"
    out["Hit_count"] = 1

    # Apply natural sorting by eccDNA_id
    out = natural_sort_eccdna_id(out)

    cols = [
        "eccDNA_id",
        "Regions",
        "Strand",
        "Length",
        "eccDNA_type",
        "State",
        "Seg_total",
        "Hit_count",
    ]
    for col in (
        ColumnStandard.CONFIDENCE_SCORE,
        ColumnStandard.QUERY_COV_BEST,
        ColumnStandard.QUERY_COV_2ND,
        ColumnStandard.MAPQ_BEST,
        ColumnStandard.IDENTITY_BEST,
        ColumnStandard.LOW_MAPQ,
        ColumnStandard.LOW_IDENTITY,
    ):
        if col in out.columns:
            cols.append(col)
    return out[cols]


def detect_related_files(base_dir: Path | str, prefix: str, ecc_type: str) -> list[Path]:
    """Detect related files for a given eccDNA type based on naming patterns.

    Args:
        base_dir: Directory to search in
        prefix: Sample prefix
        ecc_type: eccDNA type (U, M, C)

    Returns:
        List of related file paths
    """
    base_dir = Path(base_dir)
    related_files = []

    # Define patterns for each type
    patterns = {
        "U": [
            f"{prefix}_UeccDNA_C.fasta",
            f"{prefix}_UeccDNA*.fa",
            f"{prefix}_UeccDNA*.bed",
            f"{prefix}_Uecc*.fasta",
            f"{prefix}_Uecc*.fa",
            f"{prefix}_Uecc*.bed",
            f"{prefix}*Uecc*.csv",
            f"{prefix}_UeccDNA.core.csv",
        ],
        "M": [
            f"{prefix}_MeccDNA_C.fasta",
            f"{prefix}_MeccDNA*.fa",
            f"{prefix}_MeccDNA*.bed",
            f"{prefix}_MeccSites*.bed",
            f"{prefix}_Mecc*.fasta",
            f"{prefix}_Mecc*.fa",
            f"{prefix}_Mecc*.bestsite.bed",
            f"{prefix}*Mecc*.csv",
            f"{prefix}_MeccSites.core.csv",
            f"{prefix}_MeccSites.bed",
        ],
        "C": [
            f"{prefix}_CeccDNA_C.fasta",
            f"{prefix}_CeccDNA*.fa",
            f"{prefix}_CeccDNA*.bed",
            f"{prefix}_CeccDNA*.bedpe",
            f"{prefix}_CeccSegments*.bed",
            f"{prefix}_CeccSegments*.bedpe",
            f"{prefix}_Cecc*.junctions.bedpe",
            f"{prefix}_Cecc*.segments.bed",
            f"{prefix}*Cecc*.csv",
            f"{prefix}_CeccSegments.core.csv",
            f"{prefix}_CeccJunctions.bedpe",
        ],
    }

    if ecc_type in patterns:
        for pattern in patterns[ecc_type]:
            for file_path in base_dir.glob(pattern):
                if file_path.is_file():
                    # Don't include the Confirmed.csv files we generate
                    if "Confirmed.csv" not in file_path.name:
                        related_files.append(file_path)

    return list(set(related_files))  # Remove duplicates


def organize_umc_files(
    output_dir: Path | str,
    prefix: str,
    umc_files: Optional[dict[str, list[Path]]] = None,
    auto_detect: bool = True,
) -> dict[str, Path]:
    """Organize UMC output files into type-specific folders.

    Creates folders:
    - {prefix}_Uecc_C/  (C = Confirmed)
    - {prefix}_Mecc_C/  (C = Confirmed)
    - {prefix}_Cecc_C/  (C = Confirmed)

    Args:
        output_dir: Base output directory
        prefix: Sample prefix
        umc_files: Optional dict mapping eccDNA type to list of files to move
        auto_detect: If True, automatically detect related files

    Returns:
        Dictionary mapping eccDNA type to created folder path
    """
    organizer_logger = get_logger("ecc_dedup.organizer")

    output_dir = Path(output_dir)
    folders = {}

    # Create folders (C suffix means Confirmed)
    for ecc_type, folder_suffix in [("U", "Uecc_C"), ("M", "Mecc_C"), ("C", "Cecc_C")]:
        folder_path = output_dir / f"{prefix}_{folder_suffix}"
        folder_path.mkdir(parents=True, exist_ok=True)
        folders[ecc_type] = folder_path

        # Collect files to move
        files_to_move = []

        # Add provided files
        if umc_files and ecc_type in umc_files:
            files_to_move.extend(umc_files[ecc_type])

        # Auto-detect additional files
        if auto_detect:
            detected = detect_related_files(output_dir, prefix, ecc_type)
            files_to_move.extend(detected)

        # Remove duplicates and move files
        moved_files = []
        for file_path in set(files_to_move):
            if file_path.exists():
                # Move all UMC output files (including core.csv)
                dest = folder_path / file_path.name
                try:
                    shutil.move(str(file_path), str(dest))
                    moved_files.append(file_path.name)
                except Exception as exc:
                    organizer_logger.warning(f"Could not move {file_path.name}: {exc}")

        if moved_files:
            organizer_logger.info(f"Moved {len(moved_files)} files to {folder_path.name}/")

    return folders


# ================== CD-HIT Cluster Parser ==================
class CDHitClusters:
    """Parse and manage CD-HIT cluster information from .clstr or CSV exports."""

    def __init__(self) -> None:
        self.members: dict[str, list[str]] = {}  # cluster_id -> member_ids
        self.rep_of: dict[str, str] = {}  # cluster_id -> representative_id

    def parse(self, cluster_file: Path) -> None:
        """Parse cluster assignments from .clstr or CSV mapping."""
        cluster_path = Path(cluster_file)
        if not cluster_path.exists():
            raise FileNotFoundError(f"Cluster file not found: {cluster_file}")

        # Reset state so the same parser instance can be reused safely
        self.members.clear()
        self.rep_of.clear()

        suffix = cluster_path.suffix.lower()
        try:
            if suffix == ".csv":
                self._parse_csv(cluster_path)
            else:
                self._parse_csv(cluster_path)
        except (ValueError, csv.Error):
            # Fallback to legacy .clstr parsing if CSV parsing fails
            self._parse_clstr(cluster_path)

    def _parse_clstr(self, clstr_file: Path) -> None:
        current_cluster_id: Optional[str] = None
        current_members: list[str] = []
        current_rep: Optional[str] = None

        with clstr_file.open("r") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line:
                    continue

                # New cluster header
                if line.startswith(">Cluster"):
                    # Save previous cluster if exists
                    if current_cluster_id is not None and current_members:
                        self.members[current_cluster_id] = current_members.copy()
                        if current_rep:
                            self.rep_of[current_cluster_id] = current_rep

                    # Parse cluster ID
                    tokens = line.split()
                    try:
                        current_cluster_id = str(int(tokens[1]))
                    except (IndexError, ValueError):
                        current_cluster_id = line.split(">Cluster", 1)[1].strip()

                    current_members = []
                    current_rep = None
                    continue

                # Member line - extract ID
                # Pattern: "0       12634nt, >U23805... *"
                match = re.search(r">\s*([^\.>\s][^\.]*?)\.\.\.", line)
                if not match:
                    match = re.search(r">\s*([^\s,]+)", line)

                if match:
                    member_id = match.group(1)
                    current_members.append(member_id)

                    # Check if this is the representative (marked with *)
                    if line.endswith("*"):
                        current_rep = member_id

        # Save last cluster
        if current_cluster_id is not None and current_members:
            self.members[current_cluster_id] = current_members.copy()
            if current_rep:
                self.rep_of[current_cluster_id] = current_rep

        # Ensure every cluster has a representative fallback
        for cid, ids in self.members.items():
            if ids and cid not in self.rep_of:
                self.rep_of[cid] = ids[0]

    def _parse_csv(self, csv_file: Path) -> None:
        with csv_file.open("r", newline="") as handle:
            reader = csv.DictReader(handle)
            if (
                not reader.fieldnames
                or "id" not in reader.fieldnames
                or "cluster" not in reader.fieldnames
            ):
                raise ValueError(f"Cluster CSV must contain 'id' and 'cluster' columns: {csv_file}")

            for row in reader:
                seq_id = (row.get("id") or "").strip()
                if not seq_id:
                    continue

                raw_cluster = row.get("cluster")
                cluster_id = str(raw_cluster).strip() if raw_cluster is not None else ""
                if not cluster_id:
                    cluster_id = f"singleton:{seq_id}"

                members = self.members.setdefault(cluster_id, [])
                members.append(seq_id)

                rep_flag = row.get("is_representative")
                if rep_flag is None:
                    rep_flag = row.get("representative")

                if rep_flag is not None:
                    if str(rep_flag).strip().lower() in {"true", "1", "yes", "y"}:
                        self.rep_of.setdefault(cluster_id, seq_id)

        # Fallback: if no representative flagged, pick first member per cluster
        for cid, ids in self.members.items():
            if ids and cid not in self.rep_of:
                self.rep_of[cid] = ids[0]

    def get_id_to_cluster_map(self) -> dict[str, str]:
        """Create mapping from member ID to cluster ID."""
        mapping = {}
        for cluster_id, member_list in self.members.items():
            for member_id in member_list:
                mapping[member_id] = cluster_id
        return mapping


# ================== Main ECC Dedup Class ==================
class eccDedup:
    """eccDNA cluster deduplication processor with confirmed tables generation."""

    def __init__(self, logger: Optional[logging.Logger] = None):
        """Initialize eccDedup."""
        self.logger = logger or get_logger(self.__class__.__name__)

    def normalize_coordinates(self, df: pd.DataFrame) -> pd.DataFrame:
        """Normalize input columns to standard naming and canonical types."""
        df = df.copy()

        # Rename straightforward legacy columns to standard
        rename_map = {
            legacy: snake
            for legacy, snake in LEGACY_TO_SNAKE_CASE.items()
            if legacy in df.columns and snake not in df.columns
        }
        if rename_map:
            df = df.rename(columns=rename_map)

        # eccDNA identifier - use ColumnStandard
        if ColumnStandard.ECCDNA_ID not in df.columns:
            for cand in ("eccDNA_id", "eccdna_id", "ecc_id", "id"):
                if cand in df.columns:
                    df[ColumnStandard.ECCDNA_ID] = df[cand].astype(str)
                    if cand != ColumnStandard.ECCDNA_ID:
                        df = df.drop(columns=[cand])
                    break

        # Chromosome / contig name - use ColumnStandard
        if ColumnStandard.CHR not in df.columns:
            for cand in ("chr", "eChr", "subject_id", "chrom", "chromosome"):
                if cand in df.columns:
                    df[ColumnStandard.CHR] = df[cand]
                    if cand != ColumnStandard.CHR:
                        df = df.drop(columns=[cand])
                    break

        # Strand standardisation - use ColumnStandard
        if ColumnStandard.STRAND not in df.columns:
            for cand in ("strand", "eStrand", "sstrand", "strand_info"):
                if cand in df.columns:
                    df[ColumnStandard.STRAND] = df[cand]
                    if cand != ColumnStandard.STRAND:
                        df = df.drop(columns=[cand])
                    break
        if ColumnStandard.STRAND in df.columns:
            df[ColumnStandard.STRAND] = df[ColumnStandard.STRAND].apply(_normalize_strand)

        # Reads metadata - use ColumnStandard
        if ColumnStandard.READS not in df.columns:
            for cand in ("reads", "eReads", "read_name"):
                if cand in df.columns:
                    df[ColumnStandard.READS] = df[cand]
                    if cand != ColumnStandard.READS:
                        df = df.drop(columns=[cand])
                    break

        # Clean reads metadata
        if ColumnStandard.READS in df.columns:
            reads_stripped = df[ColumnStandard.READS].astype(str).str.strip()
            if "query_id" in df.columns:
                query_stripped = df["query_id"].astype(str).str.strip()
                mask_same_reads = reads_stripped.eq(query_stripped)
                df.loc[mask_same_reads, ColumnStandard.READS] = ""
            mask_empty_reads = reads_stripped.eq("") | reads_stripped.eq("nan")
            df.loc[mask_empty_reads, ColumnStandard.READS] = ""

        # Copy number - use ColumnStandard
        if ColumnStandard.COPY_NUMBER not in df.columns:
            for cand in ("copy_number", "copyNum"):
                if cand in df.columns:
                    df[ColumnStandard.COPY_NUMBER] = df[cand]
                    if cand != ColumnStandard.COPY_NUMBER:
                        df = df.drop(columns=[cand])
                    break
            if ColumnStandard.COPY_NUMBER not in df.columns:
                df[ColumnStandard.COPY_NUMBER] = pd.NA
        ensure_int64_column(df, ColumnStandard.COPY_NUMBER)

        # Handle repeat_number separately to preserve it for tests
        if "repeat_number" not in df.columns:
            for cand in ("eRepeatNum", "repeatNum", "repeat_num"):
                if cand in df.columns:
                    df["repeat_number"] = df[cand]
                    # Don't drop the original column yet
                    break
            else:
                # If no explicit repeat number, use copy_number as fallback
                df["repeat_number"] = df[ColumnStandard.COPY_NUMBER]

        # Coordinate handling (0-based) - use ColumnStandard
        if ColumnStandard.START0 not in df.columns:
            for cand, offset in [
                ("start_0based", 0),
                ("eStart0", 0),
                ("start0", 0),
                ("eStart", -1),
                ("s_start", -1),
                ("start", 0),
            ]:
                if cand in df.columns:
                    df[ColumnStandard.START0] = to_numeric_safe(df[cand], 0) + offset
                    if cand != ColumnStandard.START0:
                        df = df.drop(columns=[cand])
                    break
        ensure_int64_column(df, ColumnStandard.START0)

        if ColumnStandard.END0 not in df.columns:
            for cand in ("end_0based", "eEnd0", "end0", "eEnd", "s_end", "end"):
                if cand in df.columns:
                    df[ColumnStandard.END0] = to_numeric_safe(df[cand], 0)
                    if cand != ColumnStandard.END0:
                        df = df.drop(columns=[cand])
                    break
        ensure_int64_column(df, ColumnStandard.END0)

        # Length - use ColumnStandard
        if ColumnStandard.LENGTH not in df.columns:
            for cand in ("length", "eLength", "consLen"):
                if cand in df.columns:
                    df[ColumnStandard.LENGTH] = to_numeric_safe(df[cand], 0)
                    if cand != ColumnStandard.LENGTH:
                        df = df.drop(columns=[cand])
                    break
            if ColumnStandard.LENGTH not in df.columns:
                df[ColumnStandard.LENGTH] = pd.to_numeric(
                    df.get(ColumnStandard.END0), errors="coerce"
                ) - pd.to_numeric(df.get(ColumnStandard.START0), errors="coerce")
        ensure_int64_column(df, ColumnStandard.LENGTH)

        # Match degree within 0-100 - use ColumnStandard
        if ColumnStandard.MATCH_DEGREE not in df.columns:
            for cand in ("match_degree", "MatDegree", "identity"):
                if cand in df.columns:
                    df[ColumnStandard.MATCH_DEGREE] = df[cand]
                    if cand != ColumnStandard.MATCH_DEGREE:
                        df = df.drop(columns=[cand])
                    break
            if ColumnStandard.MATCH_DEGREE not in df.columns and "Gap_Percentage" in df.columns:
                df[ColumnStandard.MATCH_DEGREE] = 100 - to_numeric_safe(df["Gap_Percentage"], 0)
            elif ColumnStandard.MATCH_DEGREE not in df.columns:
                df[ColumnStandard.MATCH_DEGREE] = pd.NA
        df[ColumnStandard.MATCH_DEGREE] = clamp_match_degree(df[ColumnStandard.MATCH_DEGREE])

        # eccDNA type bookkeeping column
        if ColumnStandard.ECCDNA_TYPE not in df.columns:
            df[ColumnStandard.ECCDNA_TYPE] = pd.NA

        # Handle original eccDNA ID
        if "orig_eccdna_id" not in df.columns and "orig_eccDNA_id" in df.columns:
            df["orig_eccdna_id"] = df["orig_eccDNA_id"]
            df = df.drop(columns=["orig_eccDNA_id"])

        return df

    def reorder_columns_for_output(self, df: pd.DataFrame) -> pd.DataFrame:
        """Reorder columns: fixed front, others, then reads_count and read_name at end."""
        df = df.copy()

        # Map standard columns to legacy aliases for output compatibility
        alias_map = {
            "eccDNA_id": [ColumnStandard.ECCDNA_ID, "eccdna_id"],
            "chr": [ColumnStandard.CHR, "eChr"],
            "start0": [ColumnStandard.START0, "start_0based", "eStart0"],
            "end0": [ColumnStandard.END0, "end_0based", "eEnd0"],
            "strand": [ColumnStandard.STRAND, "eStrand"],
            "length": [ColumnStandard.LENGTH, "eLength"],
            "match_degree": [ColumnStandard.MATCH_DEGREE, "MatDegree"],
            "copy_number": [ColumnStandard.COPY_NUMBER, "copyNum"],
        }

        for target, aliases in alias_map.items():
            if target in df.columns:
                continue
            for alias in aliases:
                if alias in df.columns:
                    df = df.rename(columns={alias: target})
                    break

        # Ensure reads_count column exists
        if "reads_count" not in df.columns:
            df["reads_count"] = pd.NA

        front_cols = [c for c in FIXED_FRONT_COLUMNS if c in df.columns]

        # Use read_name if available, otherwise fall back to reads column
        # Don't create empty reads column if it doesn't exist
        if "read_name" in df.columns:
            reads_col = "read_name"
        elif ColumnStandard.READS in df.columns:
            reads_col = ColumnStandard.READS
        else:
            # No reads column at all - that's okay, just use reads_count
            reads_col = None

        if reads_col:
            middle_cols = [
                c for c in df.columns if c not in front_cols and c not in ("reads_count", reads_col)
            ]
            ordered_cols = front_cols + middle_cols + ["reads_count", reads_col]
        else:
            middle_cols = [c for c in df.columns if c not in front_cols and c != "reads_count"]
            ordered_cols = front_cols + middle_cols + ["reads_count"]

        return df[ordered_cols]

    def _finalize_dataframe(self, df: pd.DataFrame, dtype: str) -> pd.DataFrame:
        """Ensure canonical column naming, typing and metadata."""
        df = self.normalize_coordinates(df)

        df[ColumnStandard.ECCDNA_TYPE] = dtype

        if "orig_eccdna_id" not in df.columns:
            df["orig_eccdna_id"] = df[ColumnStandard.ECCDNA_ID]

        ensure_int64_column(df, ColumnStandard.START0)
        ensure_int64_column(df, ColumnStandard.END0)
        ensure_int64_column(df, ColumnStandard.LENGTH)

        df[ColumnStandard.COPY_NUMBER] = pd.to_numeric(
            df.get(ColumnStandard.COPY_NUMBER), errors="coerce"
        )

        # Ensure repeat_number column exists and is properly set
        if "repeat_number" not in df.columns:
            # Check for legacy repeat number columns
            for cand in ("eRepeatNum", "repeatNum", "repeat_num"):
                if cand in df.columns:
                    df["repeat_number"] = pd.to_numeric(df[cand], errors="coerce")
                    break
            else:
                # Fallback to copy_number if no repeat number found
                df["repeat_number"] = df[ColumnStandard.COPY_NUMBER]

        if "repeat_number" in df.columns:
            df["repeat_number"] = pd.to_numeric(df["repeat_number"], errors="coerce")
            df[ColumnStandard.COPY_NUMBER] = df[ColumnStandard.COPY_NUMBER].where(
                df[ColumnStandard.COPY_NUMBER].notna(), df["repeat_number"]
            )
            df["repeat_number"] = df["repeat_number"].round().astype("Int64")
            # Add legacy alias for tests
            df["eRepeatNum"] = df["repeat_number"]
        df[ColumnStandard.COPY_NUMBER] = (
            df[ColumnStandard.COPY_NUMBER].round().astype("Int64")
            if len(df)
            else pd.Series(dtype="Int64")
        )

        df[ColumnStandard.MATCH_DEGREE] = clamp_match_degree(df.get(ColumnStandard.MATCH_DEGREE))

        if ColumnStandard.STRAND in df.columns:
            df[ColumnStandard.STRAND] = df[ColumnStandard.STRAND].apply(_normalize_strand)

        if "num_merged" not in df.columns:
            df["num_merged"] = 1
        df["num_merged"] = (
            pd.to_numeric(df["num_merged"], errors="coerce").fillna(1).astype("Int64")
        )

        if "merged_from_ids" in df.columns:
            df["merged_from_ids"] = (
                df["merged_from_ids"].fillna(df[ColumnStandard.ECCDNA_ID]).astype(str)
            )
        else:
            df["merged_from_ids"] = df[ColumnStandard.ECCDNA_ID]

        if "reads_count" in df.columns:
            df["reads_count"] = pd.to_numeric(df["reads_count"], errors="coerce").astype("Int64")

        # Ensure reads column exists
        if ColumnStandard.READS in df.columns:
            df[ColumnStandard.READS] = df[ColumnStandard.READS].fillna("").astype(str)

        # Populate legacy mixed-case aliases for downstream compatibility
        legacy_aliases = {
            "eccDNA_id": ColumnStandard.ECCDNA_ID,
            "eChr": ColumnStandard.CHR,
            "eStart0": ColumnStandard.START0,
            "eEnd0": ColumnStandard.END0,
            "eStrand": ColumnStandard.STRAND,
            "eLength": ColumnStandard.LENGTH,
            "MatDegree": ColumnStandard.MATCH_DEGREE,
            "copyNum": ColumnStandard.COPY_NUMBER,
            "eClass": ColumnStandard.ECCDNA_TYPE,
            "data_type": ColumnStandard.ECCDNA_TYPE,
            # readName removed - only using read_name now
            "eReads": ColumnStandard.READS,
            "orig_eccDNA_id": "orig_eccdna_id",
        }
        for legacy, modern in legacy_aliases.items():
            if modern in df.columns:
                df[legacy] = df[modern]

        # Add 1-based coordinate aliases (eStart/eEnd) for legacy compatibility
        if ColumnStandard.START0 in df.columns:
            df["eStart"] = df[ColumnStandard.START0] + 1  # Convert 0-based to 1-based
        if ColumnStandard.END0 in df.columns:
            df["eEnd"] = df[ColumnStandard.END0]  # End coordinate is the same in both systems

        # Add test-expected legacy snake_case columns
        test_legacy_columns = {
            "eccdna_id": ColumnStandard.ECCDNA_ID,
            "chr": ColumnStandard.CHR,
            "start_0based": ColumnStandard.START0,
            "end_0based": ColumnStandard.END0,
            "strand": ColumnStandard.STRAND,
            "length": ColumnStandard.LENGTH,
            "match_degree": ColumnStandard.MATCH_DEGREE,
            "copy_number": ColumnStandard.COPY_NUMBER,
            "eccdna_type": ColumnStandard.ECCDNA_TYPE,
            "reads": ColumnStandard.READS,
        }
        for test_col, standard_col in test_legacy_columns.items():
            if standard_col in df.columns and test_col not in df.columns:
                df[test_col] = df[standard_col]

        # Add repeat_number column if available
        if "repeat_number" in df.columns and "repeat_number" not in test_legacy_columns:
            pass  # repeat_number is already available

        return df

    def process_uecc(self, df: pd.DataFrame, clusters: CDHitClusters, dtype: str) -> pd.DataFrame:
        """Process Uecc data: one representative row per cluster."""
        if df.empty:
            return df

        # Track if the original input had copyNum column
        had_copy_num_input = "copyNum" in df.columns

        df = self.normalize_coordinates(df)

        id_to_cluster = clusters.get_id_to_cluster_map()
        df["cluster_id"] = df[ColumnStandard.ECCDNA_ID].map(id_to_cluster)

        # Handle singletons
        singleton_mask = df["cluster_id"].isna()
        df.loc[singleton_mask, "cluster_id"] = "singleton:" + df.loc[
            singleton_mask, ColumnStandard.ECCDNA_ID
        ].astype(str)

        # Select representative for each cluster
        rep_map = clusters.rep_of

        # Calculate metrics for selection
        if ColumnStandard.LENGTH in df.columns:
            metric = df.groupby(["cluster_id", ColumnStandard.ECCDNA_ID])[
                ColumnStandard.LENGTH
            ].max()
        else:
            metric = df.groupby(["cluster_id", ColumnStandard.ECCDNA_ID]).size()

        metric_df = metric.reset_index().rename(columns={metric.name: "metric"})

        # Add representative flag
        rep_df = pd.DataFrame(
            {"cluster_id": list(rep_map.keys()), "rep_id": list(rep_map.values())}
        )
        metric_df = metric_df.merge(rep_df, on="cluster_id", how="left")
        metric_df["is_rep"] = (metric_df[ColumnStandard.ECCDNA_ID] == metric_df["rep_id"]).astype(
            int
        )

        # Sort to prioritize: representative > metric
        metric_df = metric_df.sort_values(
            ["cluster_id", "is_rep", "metric"], ascending=[True, False, False]
        )

        # Get chosen representative for each cluster
        chosen = metric_df.groupby("cluster_id").head(1)[["cluster_id", ColumnStandard.ECCDNA_ID]]
        chosen = chosen.rename(columns={ColumnStandard.ECCDNA_ID: "chosen_id"})

        # Aggregate metadata from all cluster members
        agg_dict = {}

        # Handle reads column
        if ColumnStandard.READS in df.columns:
            agg_dict[ColumnStandard.READS] = df.groupby("cluster_id")[ColumnStandard.READS].apply(
                merge_read_lists
            )

        # Confidence/evidence aggregation (numeric; conservative defaults).
        numeric_aggs = {
            ColumnStandard.CONFIDENCE_SCORE: "max",
            ColumnStandard.MAPQ_BEST: "max",
            ColumnStandard.MAPQ_MIN: "min",
            ColumnStandard.IDENTITY_BEST: "max",
            ColumnStandard.IDENTITY_MIN: "min",
            ColumnStandard.QUERY_COV_BEST: "max",
            ColumnStandard.QUERY_COV_2ND: "max",
        }
        for col, how in numeric_aggs.items():
            if col not in df.columns:
                continue
            series = pd.to_numeric(df[col], errors="coerce")
            grouped = series.groupby(df["cluster_id"])
            agg_dict[col] = grouped.max() if how == "max" else grouped.min()

        if "repeat_number" in df.columns:
            unique_by_id = (
                df[["cluster_id", ColumnStandard.ECCDNA_ID, "repeat_number"]]
                .dropna()
                .drop_duplicates()
            )
            agg_dict["repeat_number"] = unique_by_id.groupby("cluster_id")["repeat_number"].apply(
                lambda s: to_numeric_safe(s, 0).sum()
            )

        if ColumnStandard.STRAND in df.columns:
            agg_dict[ColumnStandard.STRAND] = df.groupby("cluster_id")[ColumnStandard.STRAND].apply(
                majority_vote
            )

        if agg_dict:
            agg_df = pd.DataFrame(agg_dict)
        else:
            agg_df = pd.DataFrame(index=df["cluster_id"].unique())

        # Join with chosen representatives
        result = df.merge(chosen, on="cluster_id", how="left")
        result = result[result[ColumnStandard.ECCDNA_ID] == result["chosen_id"]].drop_duplicates(
            subset=["cluster_id", ColumnStandard.ECCDNA_ID]
        )

        # Merge aggregated data
        result = result.merge(agg_df, on="cluster_id", how="left", suffixes=("", "_agg"))

        for col in agg_dict.keys():
            col_agg = f"{col}_agg"
            if col_agg in result.columns:
                result[col] = result[col_agg].where(result[col_agg].notna(), result.get(col))
                result = result.drop(columns=[col_agg])

        # Recompute evidence flags from aggregated maxima when possible.
        if ColumnStandard.MAPQ_BEST in result.columns:
            mapq_best = (
                pd.to_numeric(result[ColumnStandard.MAPQ_BEST], errors="coerce")
                .fillna(0)
                .astype(int)
            )
            result[ColumnStandard.LOW_MAPQ] = mapq_best < int(CONF_MAPQ_LOW_THRESHOLD)

        if ColumnStandard.IDENTITY_BEST in result.columns:
            identity_best = pd.to_numeric(
                result[ColumnStandard.IDENTITY_BEST], errors="coerce"
            ).fillna(0.0)
            result[ColumnStandard.LOW_IDENTITY] = identity_best < float(CONF_IDENTITY_LOW_THRESHOLD)

        # Add merge tracking
        cluster_counts = df.groupby("cluster_id")[ColumnStandard.ECCDNA_ID].nunique()
        result["num_merged"] = result["cluster_id"].map(cluster_counts).fillna(1).astype(int)

        merged_ids = df.groupby("cluster_id")[ColumnStandard.ECCDNA_ID].apply(
            lambda s: ";".join(sorted(s.unique()))
        )
        result["merged_from_ids"] = (
            result["cluster_id"].map(merged_ids).fillna(result[ColumnStandard.ECCDNA_ID])
        )
        result["orig_eccdna_id"] = result[ColumnStandard.ECCDNA_ID]

        # Handle copy_number based on input data:
        # If input had copyNum values, they're preserved from the representative
        # If input didn't have copyNum, use the aggregated repeat_number
        if not had_copy_num_input and "repeat_number" in result.columns:
            result[ColumnStandard.COPY_NUMBER] = result["repeat_number"]

        result = self._finalize_dataframe(result, dtype)
        return result

    def process_mecc_cecc(
        self, df: pd.DataFrame, clusters: CDHitClusters, dtype: str
    ) -> pd.DataFrame:
        """Process Mecc/Cecc data: keep all rows of representative ID."""
        if df.empty:
            return df

        full_df = self.normalize_coordinates(df)
        id_to_cluster = clusters.get_id_to_cluster_map()
        full_df["cluster_id"] = full_df[ColumnStandard.ECCDNA_ID].map(id_to_cluster)

        singleton_mask = full_df["cluster_id"].isna()
        full_df.loc[singleton_mask, "cluster_id"] = "singleton:" + full_df.loc[
            singleton_mask, ColumnStandard.ECCDNA_ID
        ].astype(str)

        if dtype == "Mecc":
            metric_col = (
                ColumnStandard.COPY_NUMBER
                if ColumnStandard.COPY_NUMBER in full_df.columns
                else None
            )
            tiebreak_col = (
                ColumnStandard.LENGTH if ColumnStandard.LENGTH in full_df.columns else None
            )
        else:
            metric_col = ColumnStandard.LENGTH
            tiebreak_col = None

        if metric_col and metric_col in full_df.columns:
            metric = full_df.groupby(["cluster_id", ColumnStandard.ECCDNA_ID])[metric_col].max()
        else:
            metric = full_df.groupby(["cluster_id", ColumnStandard.ECCDNA_ID]).size()

        metric_df = metric.reset_index().rename(columns={metric.name: "metric"})

        if tiebreak_col and tiebreak_col in full_df.columns:
            tie = full_df.groupby(["cluster_id", ColumnStandard.ECCDNA_ID])[tiebreak_col].max()
            tie_df = tie.reset_index().rename(columns={tie.name: "tie"})
            metric_df = metric_df.merge(
                tie_df, on=["cluster_id", ColumnStandard.ECCDNA_ID], how="left"
            )
        else:
            metric_df["tie"] = 0

        rep_map = clusters.rep_of
        rep_df = pd.DataFrame(
            {"cluster_id": list(rep_map.keys()), "rep_id": list(rep_map.values())}
        )
        metric_df = metric_df.merge(rep_df, on="cluster_id", how="left")
        metric_df["is_rep"] = (metric_df[ColumnStandard.ECCDNA_ID] == metric_df["rep_id"]).astype(
            int
        )

        metric_df = metric_df.sort_values(
            ["cluster_id", "is_rep", "metric", "tie"], ascending=[True, False, False, False]
        )

        chosen = metric_df.groupby("cluster_id").head(1)[["cluster_id", ColumnStandard.ECCDNA_ID]]
        chosen = chosen.rename(columns={ColumnStandard.ECCDNA_ID: "chosen_id"})

        result = full_df.merge(chosen, on="cluster_id", how="left")
        result = result[result[ColumnStandard.ECCDNA_ID] == result["chosen_id"]].copy()

        agg_dict: dict[str, pd.Series] = {}

        # Handle reads and other metadata aggregation
        if ColumnStandard.READS in full_df.columns:
            agg_dict[ColumnStandard.READS] = full_df.groupby("cluster_id")[
                ColumnStandard.READS
            ].apply(merge_read_lists)
            # readName column removed - only using read_name now

        numeric_aggs = {
            ColumnStandard.CONFIDENCE_SCORE: "max",
            ColumnStandard.MAPQ_BEST: "max",
            ColumnStandard.MAPQ_MIN: "min",
            ColumnStandard.IDENTITY_BEST: "max",
            ColumnStandard.IDENTITY_MIN: "min",
            ColumnStandard.QUERY_COV_BEST: "max",
            ColumnStandard.QUERY_COV_2ND: "max",
        }
        for col, how in numeric_aggs.items():
            if col not in full_df.columns:
                continue
            series = pd.to_numeric(full_df[col], errors="coerce")
            grouped = series.groupby(full_df["cluster_id"])
            agg_dict[col] = grouped.max() if how == "max" else grouped.min()

        if ColumnStandard.COPY_NUMBER in full_df.columns:
            unique_by_id = (
                full_df[["cluster_id", ColumnStandard.ECCDNA_ID, ColumnStandard.COPY_NUMBER]]
                .dropna()
                .drop_duplicates()
            )
            agg_dict[ColumnStandard.COPY_NUMBER] = unique_by_id.groupby("cluster_id")[
                ColumnStandard.COPY_NUMBER
            ].apply(lambda s: to_numeric_safe(s, 0).sum())
        if ColumnStandard.STRAND in full_df.columns:
            agg_dict[ColumnStandard.STRAND] = full_df.groupby("cluster_id")[
                ColumnStandard.STRAND
            ].apply(majority_vote)

        if agg_dict:
            agg_df = pd.DataFrame(agg_dict)
            result = result.merge(agg_df, on="cluster_id", how="left", suffixes=("", "_agg"))
            for col in agg_dict.keys():
                col_agg = f"{col}_agg"
                if col_agg in result.columns:
                    # Always use aggregated value, don't fallback to original
                    result[col] = result[col_agg]
                    result = result.drop(columns=[col_agg])

        if ColumnStandard.MAPQ_BEST in result.columns:
            mapq_best = (
                pd.to_numeric(result[ColumnStandard.MAPQ_BEST], errors="coerce")
                .fillna(0)
                .astype(int)
            )
            result[ColumnStandard.LOW_MAPQ] = mapq_best < int(CONF_MAPQ_LOW_THRESHOLD)

        if ColumnStandard.IDENTITY_BEST in result.columns:
            identity_best = pd.to_numeric(
                result[ColumnStandard.IDENTITY_BEST], errors="coerce"
            ).fillna(0.0)
            result[ColumnStandard.LOW_IDENTITY] = identity_best < float(CONF_IDENTITY_LOW_THRESHOLD)

        cluster_counts = full_df.groupby("cluster_id")[ColumnStandard.ECCDNA_ID].nunique()
        result["num_merged"] = result["cluster_id"].map(cluster_counts).fillna(1).astype(int)

        merged_ids = full_df.groupby("cluster_id")[ColumnStandard.ECCDNA_ID].apply(
            lambda s: ";".join(sorted(s.unique()))
        )
        result["merged_from_ids"] = (
            result["cluster_id"].map(merged_ids).fillna(result[ColumnStandard.ECCDNA_ID])
        )
        result["orig_eccdna_id"] = result[ColumnStandard.ECCDNA_ID]
        result = self._finalize_dataframe(result, dtype)

        return result

    def renumber_eccdna_ids(self, df: pd.DataFrame, dtype: str) -> pd.DataFrame:
        """Renumber eccDNA IDs sequentially by type."""
        if df.empty:
            return df

        df = df.copy()
        if ColumnStandard.ECCDNA_ID in df.columns:
            previous_ids = df[ColumnStandard.ECCDNA_ID].astype(str)
        elif "eccDNA_id" in df.columns:
            previous_ids = df["eccDNA_id"].astype(str)
        else:
            previous_ids = df["cluster_id"].astype(str)

        unique_clusters = sorted(df["cluster_id"].astype(str).unique())

        id_mapping = {
            cluster_id: f"{dtype}DNA{i+1}" for i, cluster_id in enumerate(unique_clusters)
        }

        df[ColumnStandard.ECCDNA_ID] = df["cluster_id"].astype(str).map(id_mapping)
        df["eccDNA_id"] = df[ColumnStandard.ECCDNA_ID]

        df["orig_eccdna_id"] = previous_ids
        df["orig_eccDNA_id"] = previous_ids

        return df

    def write_uecc_outputs(
        self, df: pd.DataFrame, output_dir: Path, prefix: Optional[str], drop_seq: bool = False
    ) -> None:
        """Generate Uecc output files."""
        if df.empty:
            self.logger.info("Uecc data empty, skipping outputs")
            return

        df = self.normalize_coordinates(df)
        df = df.copy()
        df[ColumnStandard.ECCDNA_TYPE] = "Uecc"

        df["type_prefix"] = df[ColumnStandard.ECCDNA_ID].str.extract(r"([A-Za-z]+)", expand=False)
        df["id_number"] = to_numeric_safe(
            df[ColumnStandard.ECCDNA_ID].str.extract(r"(\d+)", expand=False),
            default=-1,
        )
        df = df.sort_values(["type_prefix", "id_number"])
        df = df.drop(columns=["type_prefix", "id_number"])

        read_series = None
        if ColumnStandard.READS in df.columns and df[ColumnStandard.READS].notna().any():
            read_series = df[ColumnStandard.READS].fillna("")
        else:
            read_series = pd.Series([""] * len(df), index=df.index)

        reads_count_series = read_series.apply(count_reads_from_string)

        copy_number = df.get(
            ColumnStandard.COPY_NUMBER, pd.Series([pd.NA] * len(df), index=df.index)
        )
        repeat_number = df.get("repeat_number", copy_number)

        core = pd.DataFrame(
            {
                "eccDNA_id": df[ColumnStandard.ECCDNA_ID],
                "chr": df.get(ColumnStandard.CHR, pd.NA),
                "start0": df.get(ColumnStandard.START0, pd.NA),
                "end0": df.get(ColumnStandard.END0, pd.NA),
                "strand": df.get(ColumnStandard.STRAND, pd.NA),
                "length": df.get(ColumnStandard.LENGTH, pd.NA),
                "match_degree": df.get(ColumnStandard.MATCH_DEGREE, pd.NA).apply(
                    format_two_decimals
                ),
                ColumnStandard.CONFIDENCE_SCORE: df.get(ColumnStandard.CONFIDENCE_SCORE, pd.NA),
                ColumnStandard.QUERY_COV_BEST: df.get(ColumnStandard.QUERY_COV_BEST, pd.NA),
                ColumnStandard.QUERY_COV_2ND: df.get(ColumnStandard.QUERY_COV_2ND, pd.NA),
                ColumnStandard.MAPQ_BEST: df.get(ColumnStandard.MAPQ_BEST, pd.NA),
                ColumnStandard.IDENTITY_BEST: df.get(ColumnStandard.IDENTITY_BEST, pd.NA),
                ColumnStandard.LOW_MAPQ: df.get(ColumnStandard.LOW_MAPQ, pd.NA),
                ColumnStandard.LOW_IDENTITY: df.get(ColumnStandard.LOW_IDENTITY, pd.NA),
                "copy_number": copy_number,
                "repeat_number": repeat_number,
                "eccdna_type": "Uecc",
                "num_merged": df.get("num_merged", 1),
                "merged_from_ids": df.get("merged_from_ids", df[ColumnStandard.ECCDNA_ID]),
                "reads_count": reads_count_series,
                "read_name": read_series,  # Use read_name instead of readName
            }
        )

        core = self.reorder_columns_for_output(core)

        output_file = (
            output_dir / f"{prefix}_UeccDNA.core.csv" if prefix else output_dir / "UeccDNA.core.csv"
        )
        core_df = core.copy()
        if drop_seq and "eSeq" in core_df.columns:
            core_df = core_df.drop(columns=["eSeq"])
        core_df.to_csv(output_file, index=False)
        self.logger.info(f"Wrote {output_file}")

        # Write BED file
        bed_repeat = repeat_number
        bed = pd.DataFrame(
            {
                "chrom": core["chr"],
                "chromStart": core["start0"],
                "chromEnd": core["end0"],
                "name": core["eccDNA_id"],
                "score": core["reads_count"],
                "strand": core["strand"],
                "length": core["length"],
                "eccdna_type": "Uecc",
                "repeat_number": bed_repeat,
                "match_degree": core["match_degree"],
            }
        )

        bed_file = output_dir / f"{prefix}_UeccDNA.bed" if prefix else output_dir / "UeccDNA.bed"
        bed.to_csv(bed_file, sep="\t", header=False, index=False)
        self.logger.info(f"Wrote {bed_file}")

        # Write FASTA if sequence available
        if "eSeq" in df.columns:
            fa_file = (
                output_dir / f"{prefix}_UeccDNA_C.fasta"
                if prefix
                else output_dir / "UeccDNA_C.fasta"
            )
            seq_map = df.set_index(ColumnStandard.ECCDNA_ID)["eSeq"].to_dict()
            copy_map = pd.Series(copy_number.values, index=df[ColumnStandard.ECCDNA_ID]).to_dict()

            with open(fa_file, "w") as f:
                for _, row in core.iterrows():
                    seq = seq_map.get(row["eccDNA_id"])
                    repeats = copy_map.get(row["eccDNA_id"])
                    repeats_str = str(int(repeats)) if pd.notna(repeats) else "NA"
                    reads_str = (
                        str(int(row["reads_count"])) if pd.notna(row["reads_count"]) else "NA"
                    )

                    coord = f"{row['chr']}:{row['start0']}-{row['end0']}"
                    header = (
                        f">{row['eccDNA_id']}|{coord}({row['strand']})|"
                        f"length={row['length']}|repeats={repeats_str}|reads={reads_str}"
                    )
                    f.write(header + "\n")

                    if isinstance(seq, str) and seq:
                        for i in range(0, len(seq), 60):
                            f.write(seq[i : i + 60] + "\n")
                    else:
                        f.write("N\n")

            self.logger.info(f"Wrote {fa_file}")
        else:
            self.logger.warning("No eSeq column in Uecc data, skipping FASTA output")

    def write_mecc_outputs(
        self, df: pd.DataFrame, output_dir: Path, prefix: Optional[str], drop_seq: bool = False
    ) -> None:
        """Generate Mecc output files."""
        if df.empty:
            self.logger.info("Mecc data empty, skipping outputs")
            return

        df = self.normalize_coordinates(df).copy()

        # Sort by eccDNA_id first (natural sort), then by position for multi-site entries
        df["type_prefix"] = df["eccDNA_id"].str.extract(r"([A-Za-z]+)", expand=False)
        df["id_number"] = to_numeric_safe(
            df["eccDNA_id"].str.extract(r"(\d+)", expand=False),
            default=-1,
        )
        df = df.sort_values(["type_prefix", "id_number", ColumnStandard.START0])
        df = df.drop(columns=["type_prefix", "id_number"])

        # Add hit indexing
        df["hit_index"] = df.groupby("eccDNA_id").cumcount() + 1
        hit_counts = df.groupby("eccDNA_id").size().rename("hit_count")
        df = df.merge(hit_counts, on="eccDNA_id", how="left")

        # Calculate reads_count per eccDNA
        if ColumnStandard.READS in df.columns and df[ColumnStandard.READS].notna().any():
            source_col = ColumnStandard.READS
        elif "read_name" in df.columns and df["read_name"].notna().any():
            source_col = "read_name"
        elif "eReads" in df.columns and df["eReads"].notna().any():
            source_col = "eReads"
        else:
            source_col = None

        if source_col:
            reads_count = (
                df.groupby("eccDNA_id")[source_col]
                .apply(lambda s: count_reads_from_string(";".join(s.dropna().astype(str))))
                .rename("reads_count")
            )
            read_names = df.groupby("eccDNA_id")[source_col].apply(merge_read_lists)
        else:
            reads_count = df.groupby("eccDNA_id").size().rename("reads_count")
            if ColumnStandard.READS in df.columns and df[ColumnStandard.READS].notna().any():
                read_names = df.groupby("eccDNA_id")[ColumnStandard.READS].apply(merge_read_lists)
            else:
                read_names = pd.Series(dtype=str)

        # Build core dataframe
        core = pd.DataFrame(
            {
                "eccDNA_id": df["eccDNA_id"],
                "chr": df.get(ColumnStandard.CHR, pd.NA),
                "start0": df[ColumnStandard.START0],
                "end0": df[ColumnStandard.END0],
                "strand": df.get(ColumnStandard.STRAND, pd.NA),
                "length": df[ColumnStandard.LENGTH],
                "match_degree": df.get(ColumnStandard.MATCH_DEGREE, pd.NA).apply(
                    format_two_decimals
                ),
                ColumnStandard.CONFIDENCE_SCORE: df.get(ColumnStandard.CONFIDENCE_SCORE, pd.NA),
                ColumnStandard.QUERY_COV_BEST: df.get(ColumnStandard.QUERY_COV_BEST, pd.NA),
                ColumnStandard.QUERY_COV_2ND: df.get(ColumnStandard.QUERY_COV_2ND, pd.NA),
                ColumnStandard.MAPQ_BEST: df.get(ColumnStandard.MAPQ_BEST, pd.NA),
                ColumnStandard.IDENTITY_BEST: df.get(ColumnStandard.IDENTITY_BEST, pd.NA),
                ColumnStandard.LOW_MAPQ: df.get(ColumnStandard.LOW_MAPQ, pd.NA),
                ColumnStandard.LOW_IDENTITY: df.get(ColumnStandard.LOW_IDENTITY, pd.NA),
                "copy_number": df.get(ColumnStandard.COPY_NUMBER, pd.NA),
                "eccdna_type": "Mecc",
                "hit_index": df["hit_index"],
                "hit_count": df["hit_count"],
                "num_merged": df.get("num_merged", 1),
                "merged_from_ids": df.get("merged_from_ids", df["eccDNA_id"]),
            }
        )

        core = core.merge(reads_count, on="eccDNA_id", how="left")
        if not read_names.empty:
            core["read_name"] = core["eccDNA_id"].map(read_names).fillna("").astype(str)
        else:
            core["read_name"] = ""

        core = self.reorder_columns_for_output(core)

        # Write core CSV
        core_df = core.copy()
        # readName column removed - only using read_name now
        if drop_seq and "eSeq" in df.columns:
            pass  # eSeq might be in original df but not in core

        output_file = (
            output_dir / f"{prefix}_MeccSites.core.csv"
            if prefix
            else output_dir / "MeccSites.core.csv"
        )
        core_df.to_csv(output_file, index=False)
        self.logger.info(f"Wrote {output_file}")

        # Write Sites BED
        bed = pd.DataFrame(
            {
                "chrom": core["chr"],
                "chromStart": core["start0"],
                "chromEnd": core["end0"],
                "name": core["eccDNA_id"]
                + "|site"
                + core["hit_index"].astype(str)
                + "/"
                + core["hit_count"].astype(str),
                "score": 1,
                "strand": core["strand"],
                "length": core["length"],
                "eccdna_type": core["eccdna_type"],
                "copy_number": core["copy_number"],
            }
        )

        bed_file = (
            output_dir / f"{prefix}_MeccSites.bed" if prefix else output_dir / "MeccSites.bed"
        )
        bed.to_csv(bed_file, sep="\t", header=False, index=False)
        self.logger.info(f"Wrote {bed_file}")

        # Write BestSite BED
        # Select best alignment per eccDNA using available quality metrics
        # Priority: bit_score (if non-zero) > alignment_length > first row
        if "bit_score" in df.columns and (df["bit_score"] != 0).any():
            best_idx = df.groupby("eccDNA_id")["bit_score"].idxmax()
            best_df = df.loc[best_idx].copy()
        elif "alignment_length" in df.columns:
            best_idx = df.groupby("eccDNA_id")["alignment_length"].idxmax()
            best_df = df.loc[best_idx].copy()
        else:
            best_df = df.groupby("eccDNA_id").head(1).copy()

        best_df = best_df.merge(reads_count, on="eccDNA_id", how="left")
        if not read_names.empty:
            best_df["read_name"] = best_df["eccDNA_id"].map(read_names).fillna("").astype(str)

        best_bed = pd.DataFrame(
            {
                "chrom": best_df[ColumnStandard.CHR],
                "chromStart": best_df[ColumnStandard.START0],
                "chromEnd": best_df[ColumnStandard.END0],
                "name": best_df["eccDNA_id"],
                "score": best_df["reads_count"].fillna(1).astype(int),
                "strand": best_df.get(ColumnStandard.STRAND, pd.NA),
                "eLength": best_df.get(ColumnStandard.LENGTH, pd.NA),
                "eClass": "Mecc",
                "copyNum": best_df.get(ColumnStandard.COPY_NUMBER, pd.NA),
            }
        )

        best_bed_file = (
            output_dir / f"{prefix}_MeccBestSite.bed" if prefix else output_dir / "MeccBestSite.bed"
        )
        best_bed.to_csv(best_bed_file, sep="\t", header=False, index=False)
        self.logger.info(f"Wrote {best_bed_file}")

        # Write FASTA
        fa_file = (
            output_dir / f"{prefix}_MeccDNA_C.fasta" if prefix else output_dir / "MeccDNA_C.fasta"
        )

        with open(fa_file, "w") as f:
            meta = pd.DataFrame({"eccDNA_id": df["eccDNA_id"].unique()})
            meta = meta.merge(
                df.groupby("eccDNA_id")[ColumnStandard.LENGTH].max().rename("length"),
                on="eccDNA_id",
                how="left",
            )
            meta = meta.merge(
                df.groupby("eccDNA_id")[ColumnStandard.COPY_NUMBER].max().rename("copies"),
                on="eccDNA_id",
                how="left",
            )
            meta = meta.merge(reads_count, on="eccDNA_id", how="left")
            meta = meta.merge(
                df.groupby("eccDNA_id")["hit_count"].max().rename("sites"),
                on="eccDNA_id",
                how="left",
            )

            if "eSeq" in df.columns:
                seq_pick = (
                    df.dropna(subset=["eSeq"]).groupby("eccDNA_id").head(1)[["eccDNA_id", "eSeq"]]
                )
                meta = meta.merge(seq_pick, on="eccDNA_id", how="left")

            for _, row in meta.iterrows():
                reads_n = int(row["reads_count"]) if pd.notna(row["reads_count"]) else 1
                sites_n = int(row["sites"]) if pd.notna(row["sites"]) else 1
                length_n = int(row["length"]) if pd.notna(row["length"]) else "NA"
                copies_n = row["copies"] if pd.notna(row["copies"]) else "NA"

                header = (
                    f">{row['eccDNA_id']}|multi_loci:{sites_n}_sites|"
                    f"length={length_n}|copies={copies_n}|reads={reads_n}"
                )
                f.write(header + "\n")

                seq = row.get("eSeq")
                if isinstance(seq, str) and seq:
                    for i in range(0, len(seq), 60):
                        f.write(seq[i : i + 60] + "\n")
                else:
                    f.write("N\n")

        self.logger.info(f"Wrote {fa_file}")

    def write_cecc_outputs(
        self, df: pd.DataFrame, output_dir: Path, prefix: Optional[str], drop_seq: bool = False
    ) -> None:
        """Generate Cecc output files."""
        if df.empty:
            self.logger.info("Cecc data empty, skipping outputs")
            return

        df = self.normalize_coordinates(df).copy()

        # Sort by eccDNA_id first using natural sort
        df["type_prefix"] = df["eccDNA_id"].str.extract(r"([A-Za-z]+)", expand=False)
        df["id_number"] = to_numeric_safe(
            df["eccDNA_id"].str.extract(r"(\d+)", expand=False),
            default=-1,
        )
        df = df.sort_values(["type_prefix", "id_number"])
        df = df.drop(columns=["type_prefix", "id_number"])

        # Handle segment indexing
        if "segment_order" in df.columns:
            try:
                df["seg_index"] = pd.to_numeric(df["segment_order"], errors="coerce").astype(
                    "Int64"
                )
            except Exception:
                df["seg_index"] = (
                    df.groupby("eccDNA_id")[ColumnStandard.START0].rank(method="first").astype(int)
                )
        else:
            df["seg_index"] = (
                df.groupby("eccDNA_id")[ColumnStandard.START0].rank(method="first").astype(int)
            )

        seg_total = df.groupby("eccDNA_id").size().rename("seg_total")
        df = df.merge(seg_total, on="eccDNA_id", how="left")

        # Assign junction roles
        if "junction_role" not in df.columns:

            def assign_role(row):
                if row["seg_index"] == 1:
                    return "head"
                elif row["seg_index"] == row["seg_total"]:
                    return "tail"
                else:
                    return "middle"

            df["junction_role"] = df.apply(assign_role, axis=1)

        # Calculate reads_count
        if ColumnStandard.READS in df.columns and df[ColumnStandard.READS].notna().any():
            read_source = ColumnStandard.READS
        elif "read_name" in df.columns and df["read_name"].notna().any():
            read_source = "read_name"
        elif "eReads" in df.columns and df["eReads"].notna().any():
            read_source = "eReads"
        else:
            read_source = None

        if read_source:
            reads_per_id = (
                df.groupby("eccDNA_id")[read_source]
                .apply(lambda s: count_reads_from_string(";".join(s.dropna().astype(str))))
                .rename("reads_count")
            )
            read_names = df.groupby("eccDNA_id")[read_source].apply(merge_read_lists)
        else:
            reads_per_id = df.groupby("eccDNA_id").size().rename("reads_count")
            read_names = pd.Series(dtype=str)

        df = df.merge(reads_per_id, on="eccDNA_id", how="left")
        if not read_names.empty:
            mapped_reads = df["eccDNA_id"].map(read_names).fillna("")
        else:
            mapped_reads = pd.Series([""] * len(df), index=df.index)

        # Build core dataframe
        core = pd.DataFrame(
            {
                "eccDNA_id": df["eccDNA_id"],
                "chr": df.get(ColumnStandard.CHR, pd.NA),
                "start0": df[ColumnStandard.START0],
                "end0": df[ColumnStandard.END0],
                "strand": df.get(ColumnStandard.STRAND, pd.NA),
                "seg_index": df["seg_index"],
                "seg_total": df["seg_total"],
                "junction_role": df.get("junction_role", pd.NA),
                "length": df[ColumnStandard.LENGTH],
                "match_degree": df.get(ColumnStandard.MATCH_DEGREE, pd.NA).apply(
                    format_two_decimals
                ),
                ColumnStandard.CONFIDENCE_SCORE: df.get(ColumnStandard.CONFIDENCE_SCORE, pd.NA),
                ColumnStandard.QUERY_COV_BEST: df.get(ColumnStandard.QUERY_COV_BEST, pd.NA),
                ColumnStandard.QUERY_COV_2ND: df.get(ColumnStandard.QUERY_COV_2ND, pd.NA),
                ColumnStandard.MAPQ_BEST: df.get(ColumnStandard.MAPQ_BEST, pd.NA),
                ColumnStandard.IDENTITY_BEST: df.get(ColumnStandard.IDENTITY_BEST, pd.NA),
                ColumnStandard.LOW_MAPQ: df.get(ColumnStandard.LOW_MAPQ, pd.NA),
                ColumnStandard.LOW_IDENTITY: df.get(ColumnStandard.LOW_IDENTITY, pd.NA),
                "copy_number": df.get(ColumnStandard.COPY_NUMBER, pd.NA),
                "num_merged": df.get("num_merged", 1),
                "merged_from_ids": df.get("merged_from_ids", df["eccDNA_id"]),
                "reads_count": df["reads_count"],
                "read_name": mapped_reads,
            }
        )

        core = self.reorder_columns_for_output(core)

        # Write core CSV
        core_df = core.copy()
        # readName column removed - only using read_name now

        output_file = (
            output_dir / f"{prefix}_CeccSegments.core.csv"
            if prefix
            else output_dir / "CeccSegments.core.csv"
        )
        core_df.to_csv(output_file, index=False)
        self.logger.info(f"Wrote {output_file}")

        # Write Segments BED
        name_col = (
            core["eccDNA_id"]
            + "|seg"
            + core["seg_index"].astype(str)
            + "/"
            + core["seg_total"].astype(str)
            + "|"
            + core["junction_role"].astype(str)
        )

        bed = pd.DataFrame(
            {
                "chrom": core["chr"],
                "chromStart": core["start0"],
                "chromEnd": core["end0"],
                "name": name_col,
                "score": core["reads_count"],
                "strand": core["strand"],
                "length": core["length"],
                "eccdna_type": "Cecc",
                "copy_number": core["copy_number"],
                "match_degree": core["match_degree"],
            }
        )

        bed_file = (
            output_dir / f"{prefix}_CeccSegments.bed" if prefix else output_dir / "CeccSegments.bed"
        )
        bed.to_csv(bed_file, sep="\t", header=False, index=False)
        self.logger.info(f"Wrote {bed_file}")

        # Write Junctions BEDPE
        junction_rows = []
        for eid, sub in core.sort_values(["eccDNA_id", "seg_index"]).groupby("eccDNA_id"):
            rows = list(sub.to_dict("records"))
            for i in range(len(rows) - 1):
                a, b = rows[i], rows[i + 1]
                junction_rows.append(
                    {
                        "chrom1": a["chr"],
                        "start1": int(a["end0"]) - 1 if pd.notna(a["end0"]) else pd.NA,
                        "end1": a["end0"],
                        "chrom2": b["chr"],
                        "start2": b["start0"],
                        "end2": int(b["start0"]) + 1 if pd.notna(b["start0"]) else pd.NA,
                        "name": f"{eid}|seg{a['seg_index']}->seg{b['seg_index']}",
                        "score": a["reads_count"],
                        "strand1": a["strand"],
                        "strand2": b["strand"],
                    }
                )

        if junction_rows:
            bedpe = pd.DataFrame(junction_rows)
            bedpe_file = (
                output_dir / f"{prefix}_CeccJunctions.bedpe"
                if prefix
                else output_dir / "CeccJunctions.bedpe"
            )
            bedpe.to_csv(bedpe_file, sep="\t", header=False, index=False)
            self.logger.info(f"Wrote {bedpe_file}")
        else:
            self.logger.info("No junctions to output for Cecc")

        # Write FASTA
        fa_file = (
            output_dir / f"{prefix}_CeccDNA_C.fasta" if prefix else output_dir / "CeccDNA_C.fasta"
        )

        with open(fa_file, "w") as f:
            seg_n = core.groupby("eccDNA_id")["seg_total"].max().rename("N")
            copies = core.groupby("eccDNA_id")["copy_number"].max().rename("copies")
            reads = core.groupby("eccDNA_id")["reads_count"].max().rename("reads")
            length = core.groupby("eccDNA_id")["length"].max().rename("length")

            junctions = (
                core.sort_values(["eccDNA_id", "seg_index"])
                .groupby("eccDNA_id")["chr"]
                .apply(lambda s: "-".join(s.astype(str).tolist()))
                .rename("junctions")
            )

            meta = pd.concat([seg_n, copies, reads, length, junctions], axis=1).reset_index()

            if "eSeq" in df.columns:
                seq_pick = (
                    df.dropna(subset=["eSeq"]).groupby("eccDNA_id").head(1)[["eccDNA_id", "eSeq"]]
                )
                meta = meta.merge(seq_pick, on="eccDNA_id", how="left")

            for _, row in meta.iterrows():
                segments_n = int(row["N"]) if pd.notna(row["N"]) else 1
                length_n = int(row["length"]) if pd.notna(row["length"]) else "NA"
                copies_n = row["copies"] if pd.notna(row["copies"]) else "NA"
                reads_n = int(row["reads"]) if pd.notna(row["reads"]) else "NA"

                header = (
                    f">{row['eccDNA_id']}|segments:{segments_n}|"
                    f"junctions:{row['junctions']}|length={length_n}|"
                    f"copies={copies_n}|reads={reads_n}"
                )
                f.write(header + "\n")

                seq = row.get("eSeq")
                if isinstance(seq, str) and seq:
                    for i in range(0, len(seq), 60):
                        f.write(seq[i : i + 60] + "\n")
                else:
                    f.write("N\n")

        self.logger.info(f"Wrote {fa_file}")

    def _generate_unified_confirmed_table(
        self, results: dict[str, pd.DataFrame], output_dir: Path, prefix: Optional[str]
    ) -> pd.DataFrame:
        """Generate unified confirmed table from all eccDNA types.

        Args:
            results: Results from deduplication process
            output_dir: Output directory
            prefix: Sample prefix

        Returns:
            Unified confirmed table DataFrame
        """
        confirmed_tables = []

        # Generate confirmed tables for each type
        for ecc_type, df in results.items():
            if df.empty:
                continue

            try:
                if ecc_type == "Uecc":
                    confirmed_df = build_u_confirmed_table(df)
                elif ecc_type == "Mecc":
                    confirmed_df = build_m_confirmed_table(df)
                elif ecc_type == "Cecc":
                    confirmed_df = build_c_confirmed_table(df)
                else:
                    continue

                confirmed_tables.append(confirmed_df)
                self.logger.info(
                    f"Generated confirmed entries for {ecc_type}: {len(confirmed_df)} entries"
                )

            except Exception as e:
                self.logger.warning(f"Failed to generate confirmed table for {ecc_type}: {e}")

        # Combine all confirmed tables
        if confirmed_tables:
            unified_df = pd.concat(confirmed_tables, ignore_index=True)

            # Apply natural sorting by eccDNA_id for the unified table
            unified_df = natural_sort_eccdna_id(unified_df)

            # Write unified confirmed table
            filename = f"{prefix}_eccDNA_Confirmed.csv" if prefix else "eccDNA_Confirmed.csv"
            filepath = output_dir / filename
            unified_df.to_csv(filepath, index=False)

            self.logger.info(
                f"Generated unified confirmed table: {filename} ({len(unified_df)} total entries)"
            )
            return unified_df
        else:
            self.logger.warning("No confirmed tables generated")
            return pd.DataFrame()

    def _organize_output_files(
        self, results: dict[str, pd.DataFrame], output_dir: Path, prefix: Optional[str]
    ) -> None:
        """Organize output files into type-specific folders.

        Args:
            results: Results from deduplication process
            output_dir: Output directory
            prefix: Sample prefix
        """
        umc_files = {}

        # Collect files for each type that was processed
        for ecc_type in results.keys():
            if ecc_type == "Confirmed":
                continue  # Skip the unified confirmed table

            if ecc_type == "Uecc":
                files = (
                    [
                        output_dir / f"{prefix}_UeccDNA.core.csv",
                        output_dir / f"{prefix}_UeccDNA.bed",
                        output_dir / f"{prefix}_UeccDNA_C.fasta",
                    ]
                    if prefix
                    else [
                        output_dir / "UeccDNA.core.csv",
                        output_dir / "UeccDNA.bed",
                        output_dir / "UeccDNA_C.fasta",
                    ]
                )
                umc_files["U"] = [f for f in files if f.exists()]

            elif ecc_type == "Mecc":
                files = (
                    [
                        output_dir / f"{prefix}_MeccSites.core.csv",
                        output_dir / f"{prefix}_MeccSites.bed",
                        output_dir / f"{prefix}_MeccBestSite.bed",
                        output_dir / f"{prefix}_MeccDNA_C.fasta",
                    ]
                    if prefix
                    else [
                        output_dir / "MeccSites.core.csv",
                        output_dir / "MeccSites.bed",
                        output_dir / "MeccBestSite.bed",
                        output_dir / "MeccDNA_C.fasta",
                    ]
                )
                umc_files["M"] = [f for f in files if f.exists()]

            elif ecc_type == "Cecc":
                files = (
                    [
                        output_dir / f"{prefix}_CeccSegments.core.csv",
                        output_dir / f"{prefix}_CeccSegments.bed",
                        output_dir / f"{prefix}_CeccJunctions.bedpe",
                        output_dir / f"{prefix}_CeccDNA_C.fasta",
                    ]
                    if prefix
                    else [
                        output_dir / "CeccSegments.core.csv",
                        output_dir / "CeccSegments.bed",
                        output_dir / "CeccJunctions.bedpe",
                        output_dir / "CeccDNA_C.fasta",
                    ]
                )
                umc_files["C"] = [f for f in files if f.exists()]

        # Organize files into folders
        if umc_files:
            organize_umc_files(output_dir, prefix or "sample", umc_files, auto_detect=True)

    def run_deduplication(
        self,
        output_dir: Path,
        prefix: str,
        uecc_input: Optional[Path] = None,
        uecc_cluster: Optional[Path] = None,
        mecc_input: Optional[Path] = None,
        mecc_cluster: Optional[Path] = None,
        cecc_input: Optional[Path] = None,
        cecc_cluster: Optional[Path] = None,
        generate_confirmed_tables: bool = True,
        organize_output_files: bool = False,
    ) -> dict[str, pd.DataFrame]:
        """Run deduplication pipeline for all eccDNA types."""
        return self.process_all_types(
            uecc_csv=uecc_input,
            uecc_clstr=uecc_cluster,
            mecc_csv=mecc_input,
            mecc_clstr=mecc_cluster,
            cecc_csv=cecc_input,
            cecc_clstr=cecc_cluster,
            output_dir=output_dir,
            prefix=prefix,
            drop_seq=True,  # Default to dropping sequence data
            generate_confirmed_tables=generate_confirmed_tables,
            organize_output_files=organize_output_files,
        )

    def process_all_types(
        self,
        uecc_csv: Optional[Path],
        uecc_clstr: Optional[Path],
        mecc_csv: Optional[Path],
        mecc_clstr: Optional[Path],
        cecc_csv: Optional[Path],
        cecc_clstr: Optional[Path],
        output_dir: Path,
        prefix: Optional[str],
        drop_seq: bool,
        generate_confirmed_tables: bool = True,
        organize_output_files: bool = True,
    ) -> dict[str, pd.DataFrame]:
        """Process all eccDNA types independently and generate outputs."""

        # Create output directory
        output_dir.mkdir(parents=True, exist_ok=True)

        results = {}

        # Process Uecc
        if uecc_csv and uecc_clstr:
            self.logger.info("Processing Uecc data...")
            if not uecc_csv.exists():
                self.logger.warning(f"Uecc CSV not found: {uecc_csv}")
            elif not uecc_clstr.exists():
                self.logger.warning(f"Uecc cluster file not found: {uecc_clstr}")
            else:
                uecc_clusters = CDHitClusters()
                uecc_clusters.parse(uecc_clstr)

                uecc_df = pd.read_csv(uecc_csv)
                if "eccDNA_id" not in uecc_df.columns and "eccdna_id" not in uecc_df.columns:
                    raise ValueError(f"Uecc CSV missing 'eccDNA_id' column: {uecc_csv}")
                if "eccDNA_id" not in uecc_df.columns and "eccdna_id" in uecc_df.columns:
                    uecc_df["eccDNA_id"] = uecc_df["eccdna_id"]

                processed_uecc = self.process_uecc(uecc_df, uecc_clusters, "Uecc")
                processed_uecc = self.renumber_eccdna_ids(processed_uecc, "Uecc")

                results["Uecc"] = processed_uecc
                self.write_uecc_outputs(processed_uecc, output_dir, prefix, drop_seq)

        # Process Mecc
        if mecc_csv and mecc_clstr:
            self.logger.info("Processing Mecc data...")
            if not mecc_csv.exists():
                self.logger.warning(f"Mecc CSV not found: {mecc_csv}")
            elif not mecc_clstr.exists():
                self.logger.warning(f"Mecc cluster file not found: {mecc_clstr}")
            else:
                mecc_clusters = CDHitClusters()
                mecc_clusters.parse(mecc_clstr)

                mecc_df = pd.read_csv(mecc_csv)
                if "eccDNA_id" not in mecc_df.columns and "eccdna_id" not in mecc_df.columns:
                    raise ValueError(f"Mecc CSV missing 'eccDNA_id' column: {mecc_csv}")
                if "eccDNA_id" not in mecc_df.columns and "eccdna_id" in mecc_df.columns:
                    mecc_df["eccDNA_id"] = mecc_df["eccdna_id"]

                processed_mecc = self.process_mecc_cecc(mecc_df, mecc_clusters, "Mecc")
                processed_mecc = self.renumber_eccdna_ids(processed_mecc, "Mecc")

                results["Mecc"] = processed_mecc
                self.write_mecc_outputs(processed_mecc, output_dir, prefix, drop_seq)

        # Process Cecc
        if cecc_csv and cecc_clstr:
            self.logger.info("Processing Cecc data...")
            if not cecc_csv.exists():
                self.logger.warning(f"Cecc CSV not found: {cecc_csv}")
            elif not cecc_clstr.exists():
                self.logger.warning(f"Cecc cluster file not found: {cecc_clstr}")
            else:
                cecc_clusters = CDHitClusters()
                cecc_clusters.parse(cecc_clstr)

                cecc_df = pd.read_csv(cecc_csv)
                if "eccDNA_id" not in cecc_df.columns and "eccdna_id" not in cecc_df.columns:
                    raise ValueError(f"Cecc CSV missing 'eccDNA_id' column: {cecc_csv}")
                if "eccDNA_id" not in cecc_df.columns and "eccdna_id" in cecc_df.columns:
                    cecc_df["eccDNA_id"] = cecc_df["eccdna_id"]

                processed_cecc = self.process_mecc_cecc(cecc_df, cecc_clusters, "Cecc")
                processed_cecc = self.renumber_eccdna_ids(processed_cecc, "Cecc")

                results["Cecc"] = processed_cecc
                self.write_cecc_outputs(processed_cecc, output_dir, prefix, drop_seq)

        # Print summary
        self.logger.info("\n" + "=" * 50)
        self.logger.info("PROCESSING COMPLETE")
        self.logger.info("=" * 50)

        for dtype in ["Uecc", "Mecc", "Cecc"]:
            if dtype in results:
                df = results[dtype]
                unique_ids = (
                    df["eccdna_id"].nunique()
                    if "eccdna_id" in df.columns
                    else df["eccDNA_id"].nunique()
                )
                total_rows = len(df)
                merged_clusters = df[df["num_merged"] > 1]["cluster_id"].nunique()
                self.logger.info(
                    f"{dtype}: {unique_ids} unique eccDNAs, {total_rows} total rows, "
                    f"{merged_clusters} merged clusters"
                )
            else:
                self.logger.info(f"{dtype}: Not processed")

        self.logger.info("=" * 50)

        # Generate unified confirmed table if requested
        if generate_confirmed_tables and results:
            self.logger.info("Generating unified confirmed table...")
            unified_confirmed = self._generate_unified_confirmed_table(results, output_dir, prefix)
            results["Confirmed"] = unified_confirmed

        # Organize files if requested
        if organize_output_files and results:
            self.logger.info("Organizing output files into type-specific folders...")
            self._organize_output_files(results, output_dir, prefix)

        return results


def _parse_args():
    """Parse CLI arguments for direct script execution"""
    import argparse

    parser = argparse.ArgumentParser(description="ECC Dedup - eccDNA clustering and deduplication")
    parser.add_argument(
        "-i",
        "--input-prefix",
        required=True,
        help="Input prefix (no extension), will find CSV and clstr files",
    )
    parser.add_argument("-o", "--output-prefix", required=True, help="Output prefix (no extension)")
    parser.add_argument(
        "--uecc-csv", help="Uecc CSV file (optional, default: input-prefix.uecc.csv)"
    )
    parser.add_argument(
        "--uecc-clstr", help="Uecc cluster file (optional, default: input-prefix.uecc.clstr)"
    )
    parser.add_argument(
        "--mecc-csv", help="Mecc CSV file (optional, default: input-prefix.mecc.csv)"
    )
    parser.add_argument(
        "--mecc-clstr", help="Mecc cluster file (optional, default: input-prefix.mecc.clstr)"
    )
    parser.add_argument(
        "--cecc-csv", help="Cecc CSV file (optional, default: input-prefix.cecc.csv)"
    )
    parser.add_argument(
        "--cecc-clstr", help="Cecc cluster file (optional, default: input-prefix.cecc.clstr)"
    )
    parser.add_argument(
        "--drop-seq", action="store_true", help="Drop sequence column to reduce file size"
    )
    parser.add_argument(
        "--no-confirmed-tables", action="store_true", help="Skip confirmed table generation"
    )
    parser.add_argument("--no-organize", action="store_true", help="Don't organize output by type")
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
        help="Log level (default: INFO)",
    )
    return parser.parse_args()


def main():
    """Main function for CLI execution"""
    from pathlib import Path
    import logging as _logging
    from circleseeker.utils.logging import get_logger

    args = _parse_args()

    input_prefix = Path(args.input_prefix)
    output_prefix = Path(args.output_prefix)

    # Set root log level
    _logging.getLogger().setLevel(getattr(_logging, args.log_level))

    # Create logger
    logger = get_logger("ECCDedup")

    # Determine input file paths
    uecc_csv = Path(args.uecc_csv) if args.uecc_csv else input_prefix.with_suffix(".uecc.csv")
    uecc_clstr = (
        Path(args.uecc_clstr) if args.uecc_clstr else input_prefix.with_suffix(".uecc.clstr")
    )
    mecc_csv = Path(args.mecc_csv) if args.mecc_csv else input_prefix.with_suffix(".mecc.csv")
    mecc_clstr = (
        Path(args.mecc_clstr) if args.mecc_clstr else input_prefix.with_suffix(".mecc.clstr")
    )
    cecc_csv = Path(args.cecc_csv) if args.cecc_csv else input_prefix.with_suffix(".cecc.csv")
    cecc_clstr = (
        Path(args.cecc_clstr) if args.cecc_clstr else input_prefix.with_suffix(".cecc.clstr")
    )

    # Determine output directory and prefix
    output_dir = output_prefix.parent
    prefix = output_prefix.name

    # Create and run ECC dedup
    dedup = eccDedup(logger=logger)

    try:
        # Run deduplication
        results = dedup.process_all_types(
            uecc_csv=uecc_csv,
            uecc_clstr=uecc_clstr,
            mecc_csv=mecc_csv,
            mecc_clstr=mecc_clstr,
            cecc_csv=cecc_csv,
            cecc_clstr=cecc_clstr,
            output_dir=output_dir,
            prefix=prefix,
            drop_seq=args.drop_seq,
            generate_confirmed_tables=not args.no_confirmed_tables,
            organize_output_files=not args.no_organize,
        )

        logger.info(f"ECC dedup complete! Output: {output_dir}")
        for file_type, df in results.items():
            if not df.empty:
                logger.info(f"  {file_type}: {len(df)} records")

    except Exception as e:
        logger.error(f"ECC dedup failed: {e}")
        raise


if __name__ == "__main__":
    main()
