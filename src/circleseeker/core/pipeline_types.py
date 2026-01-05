"""Shared pipeline types.

This module intentionally contains only lightweight dataclasses/constants so it can
be imported by step definitions/contracts without pulling in the full pipeline
implementation (and its heavier dependencies).
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Dict, List, Optional


class ResultKeys:
    """Canonical keys for pipeline step results to avoid typos.

    Note: Some keys are dynamically generated using f-strings:
    - "{ecc_type}_fasta" where ecc_type is uecc/mecc/cecc
    - "{ecc_type}_processed" where ecc_type is uecc/mecc/cecc
    - "{ecc_type}_nr99" where ecc_type is uecc/mecc/cecc
    - "{ecc_type}_clusters" where ecc_type is uecc/mecc/cecc
    - "ecc_dedup_{ecc_type}_dir" where ecc_type is uecc/mecc/cecc
    - "ecc_dedup_{ecc_type}_fasta" where ecc_type is uecc/mecc/cecc
    """

    T2R_CSV = "tandem_to_ring_csv"
    T2R_FASTA = "tandem_to_ring_fasta"
    ALIGNMENT_OUTPUT = "alignment_output"
    UECC_CSV = "uecc_csv"
    MECC_CSV = "mecc_csv"
    UNCLASSIFIED_CSV = "unclassified_csv"
    READ_FILTER_OUTPUT_FASTA = "read_filter_output_fasta"
    MINIMAP2_BAM = "minimap2_bam"
    CYRCULAR_OVERVIEW = "cyrcular_overview_table"
    INFERRED_SIMPLE_TSV = "inferred_simple_tsv"
    INFERRED_CHIMERIC_TSV = "inferred_chimeric_tsv"
    INFERRED_SIMPLE_FASTA = "inferred_simple_fasta"
    INFERRED_CHIMERIC_FASTA = "inferred_chimeric_fasta"
    CONFIRMED_CSV = "confirmed_csv"
    UNIFIED_CSV = "unified_csv"
    SUMMARY_DIR = "summary_dir"
    FINAL_RESULTS = "final_results"
    TIDEHUNTER_OUTPUT = "tidehunter_output"
    REFERENCE_MMI = "reference_mmi"
    UECC_COUNT = "uecc_count"
    MECC_COUNT = "mecc_count"
    UECC_PROCESSED = "uecc_processed"
    MECC_PROCESSED = "mecc_processed"
    CECC_PROCESSED = "cecc_processed"
    UECC_CLUSTERS = "uecc_clusters"
    MECC_CLUSTERS = "mecc_clusters"
    CECC_CLUSTERS = "cecc_clusters"
    CYRCULAR_RESULT_COUNT = "cyrcular_result_count"
    OVERLAP_REPORT = "overlap_report"
    OVERLAP_STATS = "overlap_stats"
    INFERENCE_INPUT_EMPTY = "inference_input_empty"
    INFERENCE_FAILED = "inference_failed"
    INFERENCE_ERROR = "inference_error"
    INFERENCE_ERROR_FILE = "inference_error_file"
    UECC_CORE_CSV = "uecc_core_csv"
    MECC_CORE_CSV = "mecc_core_csv"
    CECC_CORE_CSV = "cecc_core_csv"
    CECC_BUILD_OUTPUT = "cecc_build_output"
    CECC_BUILD_COUNT = "cecc_build_count"
    UECC_FASTA = "uecc_fasta"
    MECC_FASTA = "mecc_fasta"
    CECC_FASTA = "cecc_fasta"
    UECC_NR99 = "uecc_nr99"
    MECC_NR99 = "mecc_nr99"
    CECC_NR99 = "cecc_nr99"
    UECC_HARMONIZED = "uecc_harmonized"
    MECC_HARMONIZED = "mecc_harmonized"
    CECC_HARMONIZED = "cecc_harmonized"
    # Read filter statistics
    READ_FILTER_TOTAL = "read_filter_total_reads"
    READ_FILTER_FILTERED = "read_filter_filtered_reads"
    READ_FILTER_RETAINED = "read_filter_retained_reads"
    READ_FILTER_CTCR = "read_filter_ctcr_reads"
    ALL_FILTERED_FASTA = "all_filtered_fasta"
    # Inferred eccDNA directory
    INFERRED_DIR = "inferred_dir"


@dataclass
class PipelineStep:
    """Represents a pipeline step."""

    name: str
    description: str
    display_name: Optional[str] = None
    required: bool = True
    skip_condition: Optional[str] = None


@dataclass
class StepMetadata:
    """Metadata for a pipeline step execution."""

    step_name: str
    start_time: float
    end_time: Optional[float] = None
    duration: Optional[float] = None
    status: str = "running"  # running, completed, failed
    error_message: Optional[str] = None
    input_files: List[str] = field(default_factory=list)
    output_files: List[str] = field(default_factory=list)
    file_checksums: Dict[str, str] = field(default_factory=dict)


@dataclass
class PipelineState:
    """Enhanced pipeline execution state with recovery support."""

    completed_steps: List[str]
    current_step: Optional[str] = None
    failed_step: Optional[str] = None
    results: Dict[str, Any] = field(default_factory=dict)
    step_metadata: Dict[str, StepMetadata] = field(default_factory=dict)
    pipeline_start_time: Optional[float] = None
    last_checkpoint_time: Optional[float] = None
    config_hash: Optional[str] = None
    version: str = "2.0"

    def add_step_metadata(self, step_name: str, **kwargs) -> None:
        """Add or update step metadata."""
        if step_name not in self.step_metadata:
            self.step_metadata[step_name] = StepMetadata(
                step_name=step_name, start_time=time.time()
            )

        # Update metadata
        metadata = self.step_metadata[step_name]
        for key, value in kwargs.items():
            if hasattr(metadata, key):
                setattr(metadata, key, value)

    def complete_step(self, step_name: str, output_files: List[str] | None = None) -> None:
        """Mark step as completed and calculate duration."""
        if step_name in self.step_metadata:
            metadata = self.step_metadata[step_name]
            metadata.end_time = time.time()
            metadata.duration = metadata.end_time - metadata.start_time
            metadata.status = "completed"
            if output_files:
                metadata.output_files = output_files

        if step_name not in self.completed_steps:
            self.completed_steps.append(step_name)

    def fail_step(self, step_name: str, error_message: str) -> None:
        """Mark step as failed."""
        if step_name in self.step_metadata:
            metadata = self.step_metadata[step_name]
            metadata.end_time = time.time()
            metadata.duration = metadata.end_time - metadata.start_time
            metadata.status = "failed"
            metadata.error_message = error_message

        self.failed_step = step_name

    def get_total_runtime(self) -> Optional[float]:
        """Get total pipeline runtime so far."""
        if not self.pipeline_start_time:
            return None
        return time.time() - self.pipeline_start_time

    def get_step_summary(self) -> Dict[str, Any]:
        """Get summary of step execution times and status."""
        summary: Dict[str, Any] = {}
        for step_name, metadata in self.step_metadata.items():
            summary[step_name] = {
                "status": metadata.status,
                "duration": metadata.duration,
                "start_time": (
                    datetime.fromtimestamp(metadata.start_time).isoformat()
                    if metadata.start_time
                    else None
                ),
                "error": metadata.error_message,
            }
        return summary

