"""Main pipeline orchestrator for CircleSeeker"""

from __future__ import annotations

import logging
import json
import time
import hashlib
import gzip
import shutil
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Optional, Dict, Any, List, Union, Set, Iterable
from dataclasses import dataclass, asdict, field

import pandas as pd

from circleseeker.config import Config
from circleseeker.exceptions import PipelineError
from circleseeker.utils.progress import iter_progress
from circleseeker.utils.logging import get_logger

# Centralized imports for modules and external tools used across steps

# Step 0 - make_blastdb (blast)
from circleseeker.external.blast import MakeBlastDB
# Step 1 - tidehunter
from circleseeker.external.tidehunter import TideHunter
# Step 2 - tandem_to_ring
from circleseeker.modules.tandem_to_ring import TandemToRing
# Step 3 - run_blast (blast)
from circleseeker.external.blast import BlastRunner
# Step 4 - um_classify
from circleseeker.modules.um_classify import UMeccClassifier
# Step 5 - cecc_build
from circleseeker.modules.cecc_build import CeccBuild
# Step 6 - umc_process
from circleseeker.modules.umc_process import UMCProcess
# Step 7 - cd_hit
from circleseeker.external.cd_hit import CDHitEst, CDHitConfig
# Step 8 - ecc_dedup
from circleseeker.modules.ecc_dedup import eccDedup, organize_umc_files
# Step 9 - read_filter
from circleseeker.modules.read_filter import Sieve
# Step 10 - minimap2
from circleseeker.external.minimap2 import Minimap2, MiniMapConfig
# Step 11 - cyrcular_calling
from circleseeker.modules.cyrcular_calling import PipelineConfig as CCConfig
from circleseeker.modules.cyrcular_calling import CyrcularCallingPipeline
# Step 12 - curate_inferred_ecc
from circleseeker.modules.iecc_curator import (
    curate_ecc_tables,
    write_curated_tables,
    write_curated_tables_with_fasta,
)
# Step 13 - ecc_unify
from circleseeker.modules.ecc_unify import merge_eccdna_tables
# Step 14 - ecc_summary
from circleseeker.modules.ecc_summary import EccSummary
# Step 15 - ecc_packager
from circleseeker.modules import ecc_packager as ecc_packager_module


class ResultKeys:
    """Canonical keys for pipeline step results to avoid typos."""
    T2R_CSV = "tandem_to_ring_csv"
    T2R_FASTA = "tandem_to_ring_fasta"
    BLAST_DB = "blast_db"
    BLAST_OUTPUT = "blast_output"
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


@dataclass
class PipelineStep:
    """Represents a pipeline step."""
    name: str
    description: str
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
                step_name=step_name,
                start_time=time.time()
            )
        
        # Update metadata
        metadata = self.step_metadata[step_name]
        for key, value in kwargs.items():
            if hasattr(metadata, key):
                setattr(metadata, key, value)
    
    def complete_step(self, step_name: str, output_files: List[str] = None) -> None:
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
        summary = {}
        for step_name, metadata in self.step_metadata.items():
            summary[step_name] = {
                'status': metadata.status,
                'duration': metadata.duration,
                'start_time': datetime.fromtimestamp(metadata.start_time).isoformat() if metadata.start_time else None,
                'error': metadata.error_message
            }
        return summary


class Pipeline:
    """Main pipeline orchestrator with complete fixes."""
    
    # Define pipeline steps (final order; deprecated/unused steps removed)
    STEPS = [
        PipelineStep("make_blastdb", "Build BLAST database"),
        PipelineStep("tidehunter", "Run TideHunter for tandem repeat detection"),
        PipelineStep("tandem_to_ring", "Process TideHunter output"),
        PipelineStep("run_blast", "Run BLAST alignment"),
        PipelineStep("um_classify", "Classify eccDNA types"),
        PipelineStep("cecc_build", "Process Cecc candidates"),
        PipelineStep("umc_process", "Generate FASTA files"),
        PipelineStep("cd_hit", "Remove redundant sequences"),
        PipelineStep("ecc_dedup", "Coordinate results"),
        PipelineStep("read_filter", "Filter sequences"),
        PipelineStep("minimap2", "Map with Minimap2"),
        PipelineStep("cyrcular_calling", "Detect circular DNA via Cyrcular/Varlociraptor"),
        PipelineStep("curate_inferred_ecc", "Curate inferred eccDNA tables from Cyrcular"),
        PipelineStep("ecc_unify", "Merge eccDNA tables into unified output"),
        PipelineStep("ecc_summary", "Generate summary report"),
        PipelineStep("ecc_packager", "Package output files"),
    ]

    def _set_result(self, key: str, value: Any) -> None:
        """Set a result key (new canonical names only)."""
        self.state.results[key] = value

    def _serialize_path_for_state(self, path: Path | str) -> str:
        """Serialize a filesystem path with relativisation hints for checkpoint storage."""
        candidate = Path(path)
        try:
            resolved = candidate.resolve(strict=False)
        except Exception:
            resolved = candidate

        base_candidates: List[tuple[str, Optional[Path]]] = [
            ("output", getattr(self.config, "output_dir", None)),
            ("final", getattr(self, "final_output_dir", None)),
        ]

        for tag, base in base_candidates:
            if base is None:
                continue
            base_path = Path(base)
            try:
                base_resolved = base_path.resolve(strict=False)
            except Exception:
                base_resolved = base_path
            try:
                rel = resolved.relative_to(base_resolved)
                return f"{tag}::{rel.as_posix()}"
            except ValueError:
                continue

        return resolved.as_posix()

    def _resolve_stored_path(
        self,
        stored_path: Optional[str],
        fallback_names: Iterable[Union[str, Path]] = (),
    ) -> Optional[Path]:
        """Resolve stored or fallback paths to an existing Path if possible."""

        candidates: List[Path] = []
        seen: Set[Path] = set()

        def add_candidate(candidate: Union[str, Path]) -> None:
            path_candidate = Path(candidate)
            if path_candidate in seen:
                return
            seen.add(path_candidate)
            candidates.append(path_candidate)

        if stored_path:
            if "::" in stored_path:
                prefix, rel_fragment = stored_path.split("::", 1)
                rel_path = Path(rel_fragment)
                base_lookup = {
                    "output": getattr(self.config, "output_dir", None),
                    "final": getattr(self, "final_output_dir", None),
                }
                base_dir = base_lookup.get(prefix)
                if base_dir is not None:
                    add_candidate(Path(base_dir) / rel_path)
                add_candidate(rel_path)
            else:
                add_candidate(stored_path)

        for name in fallback_names:
            fallback_path = Path(name)
            if fallback_path.is_absolute():
                add_candidate(fallback_path)
                continue
            add_candidate(fallback_path)
            add_candidate(Path(self.config.output_dir) / fallback_path)
            final_dir = getattr(self, "final_output_dir", None)
            if final_dir is not None:
                add_candidate(Path(final_dir) / fallback_path)

        for candidate in list(candidates):
            if not candidate.is_absolute():
                add_candidate(Path.cwd() / candidate)

        for candidate in candidates:
            try:
                resolved = candidate.resolve(strict=False)
            except Exception:
                resolved = candidate
            if resolved.exists():
                return resolved

        return None

    def _get_result(self, key: str, default: Any = None) -> Any:
        """Get a result by canonical key only."""
        return self.state.results.get(key, default)

    def _create_combined_fasta(self, output_file: Path) -> None:
        """Combine all FASTA outputs into a single file."""
        fasta_files = []

        # Collect all FASTA files
        patterns = [
            f"{self.config.prefix}_UeccDNA*.fasta",
            f"{self.config.prefix}_MeccDNA*.fasta",
            f"{self.config.prefix}_CeccDNA*.fasta",
            "tandem_to_ring.fasta",
            f"{self.config.prefix}_circular.fasta"
        ]

        for pattern in patterns:
            fasta_files.extend(self.config.output_dir.glob(pattern))

        # Also check state results for FASTA files
        for key in self.state.results:
            if key.endswith("_fasta"):
                fasta_path = Path(self.state.results[key])
                if fasta_path.exists() and fasta_path not in fasta_files:
                    fasta_files.append(fasta_path)

        # Combine all FASTA files
        with open(output_file, 'w') as out_f:
            for fasta in fasta_files:
                if fasta.exists() and fasta.stat().st_size > 0:
                    with open(fasta, 'r') as in_f:
                        out_f.write(in_f.read())
                        if not in_f.read().endswith('\n'):
                            out_f.write('\n')

    def _safe_copy_file(self, src: Optional[Path], dest: Path) -> bool:
        """Copy a single file into the final output directory if it exists."""
        if not src:
            return False

        src_path = Path(src)
        if not src_path.exists():
            return False

        try:
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src_path, dest)
            self.logger.debug(f"Copied {src_path} -> {dest}")
            return True
        except Exception as exc:
            self.logger.warning(f"Failed to copy {src_path} to {dest}: {exc}")
            return False

    def _copy_directory_contents(self, src: Optional[Path], dest: Path) -> bool:
        """Copy all contents of a directory into destination if available."""
        if not src:
            return False

        src_path = Path(src)
        if not src_path.exists() or not src_path.is_dir():
            return False

        try:
            dest.mkdir(parents=True, exist_ok=True)
            for item in src_path.iterdir():
                target = dest / item.name
                if item.is_dir():
                    shutil.copytree(item, target, dirs_exist_ok=True)
                elif item.is_file():
                    shutil.copy2(item, target)
            self.logger.debug(f"Copied directory contents {src_path} -> {dest}")
            return True
        except Exception as exc:
            self.logger.warning(f"Failed to copy directory {src_path} to {dest}: {exc}")
            return False

    def _create_basic_output_structure(
        self,
        final_root: Path,
        *,
        uecc_dir: Optional[Path],
        mecc_dir: Optional[Path],
        cecc_dir: Optional[Path],
        inferred_dir: Optional[Path],
        merged_csv: Optional[Path],
        html_report: Optional[Path],
        text_summary: Optional[Path],
    ) -> None:
        """Fallback organizer when ecc_packager cannot run."""
        prefix = self.config.prefix
        final_root.mkdir(parents=True, exist_ok=True)

        directory_specs = [
            (uecc_dir, f"{prefix}_Confirmed_UeccDNA"),
            (mecc_dir, f"{prefix}_Confirmed_MeccDNA"),
            (cecc_dir, f"{prefix}_Confirmed_CeccDNA"),
            (inferred_dir, f"{prefix}_Inferred_eccDNA"),
        ]

        for src_dir, dest_name in directory_specs:
            dest_path = final_root / dest_name
            if self._copy_directory_contents(src_dir, dest_path):
                self.logger.info(f"Copied {dest_name} into final output")

        self._rename_inferred_simple_file(final_root)

        if merged_csv:
            merged_name = f"{prefix}_merged_output.csv"
            if not self._safe_copy_file(merged_csv, final_root / merged_name):
                self.logger.debug("Merged CSV not available for fallback packaging")

        if html_report:
            self._safe_copy_file(html_report, final_root / f"{prefix}_report.html")

        if text_summary:
            self._safe_copy_file(text_summary, final_root / f"{prefix}_summary.txt")

    def _rename_inferred_simple_file(self, final_root: Path) -> None:
        """Ensure inferred simple CSV follows the *_UeccDNA_I.csv naming convention."""
        prefix = self.config.prefix
        inferred_dir = final_root / f"{prefix}_Inferred_eccDNA"
        if not inferred_dir.exists():
            return

        legacy_simple = inferred_dir / f"{prefix}_simple.csv"
        target = inferred_dir / f"{prefix}_UeccDNA_I.csv"
        if not legacy_simple.exists():
            return

        try:
            if target.exists():
                target.unlink()
            legacy_simple.rename(target)
            self.logger.debug("Renamed inferred simple CSV to %s", target.name)
        except Exception as exc:
            self.logger.warning(
                "Could not rename %s to %s: %s",
                legacy_simple.name,
                target.name,
                exc,
            )

    def __init__(self, config: Config):
        self.config = config
        self.logger = get_logger(self.__class__.__name__)
        
        # Setup directory structure
        # Check if output_dir is already a temp directory from previous run
        if config.output_dir.name == ".tmp_work":
            # We're already in temp dir, go back to parent
            self.final_output_dir = config.output_dir.parent
            self.temp_dir = config.output_dir
        else:
            # Normal case: set up temp directory
            self.final_output_dir = config.output_dir
            self.temp_dir = config.output_dir / ".tmp_work"
        
        # Create directories
        self.final_output_dir.mkdir(parents=True, exist_ok=True)
        self.temp_dir.mkdir(parents=True, exist_ok=True)
        
        # Override config to use temp directory for all intermediate files
        self.config.output_dir = self.temp_dir
        
        # Keep checkpoint in final directory for persistence across runs
        self.state_file = self.final_output_dir / f"{config.prefix}.checkpoint"
        self.state = self._load_state()
        self._temp_files = []  # Track temporary files for cleanup
    
    def _compute_config_hash(self) -> str:
        """Compute hash of current configuration for validation."""
        config_str = json.dumps(asdict(self.config), sort_keys=True, default=str)
        return hashlib.sha256(config_str.encode()).hexdigest()[:16]
    
    def _load_state(self) -> PipelineState:
        """Load pipeline state from checkpoint with validation."""
        if self.state_file.exists():
            try:
                with open(self.state_file, 'r') as f:
                    data = json.load(f)
                
                # Check version compatibility
                if data.get('version', '1.0') != '2.0':
                    self.logger.warning("Checkpoint from older version, starting fresh")
                    return self._create_fresh_state()
                
                # Check config compatibility
                current_hash = self._compute_config_hash()
                saved_hash = data.get('config_hash')
                if saved_hash and saved_hash != current_hash:
                    self.logger.warning("Configuration changed since last checkpoint")
                    policy = getattr(self.config.runtime, 'checkpoint_policy', 'continue')
                    if policy == 'reset':
                        self.logger.info("Checkpoint policy = reset; starting fresh state")
                        return self._create_fresh_state()
                    elif policy == 'fail':
                        raise PipelineError("Config changed and checkpoint_policy=fail")
                    else:
                        self.logger.info("Checkpoint policy = continue; resuming with existing checkpoint")
                
                # Reconstruct step metadata
                step_metadata = {}
                for step_name, meta_data in data.get('step_metadata', {}).items():
                    step_metadata[step_name] = StepMetadata(**meta_data)
                
                state = PipelineState(
                    completed_steps=data.get('completed_steps', []),
                    current_step=data.get('current_step'),
                    failed_step=data.get('failed_step'),
                    results=data.get('results', {}),
                    step_metadata=step_metadata,
                    pipeline_start_time=data.get('pipeline_start_time'),
                    last_checkpoint_time=data.get('last_checkpoint_time'),
                    config_hash=current_hash,
                    version=data.get('version', '2.0')
                )
                
                self.logger.info(f"Loaded checkpoint with {len(state.completed_steps)} completed steps")
                if state.failed_step:
                    self.logger.warning(f"Previous run failed at step: {state.failed_step}")
                
                return state
                
            except Exception as e:
                self.logger.error(f"Could not load checkpoint: {e}")
                self.logger.info("Starting with fresh state")
        
        return self._create_fresh_state()
    
    def _create_fresh_state(self) -> PipelineState:
        """Create a fresh pipeline state."""
        return PipelineState(
            completed_steps=[],
            config_hash=self._compute_config_hash(),
            pipeline_start_time=time.time()
        )
    
    def _save_state(self) -> None:
        """Save pipeline state to checkpoint with backup."""
        self.state_file.parent.mkdir(parents=True, exist_ok=True)
        
        # Update checkpoint time
        self.state.last_checkpoint_time = time.time()
        
        # Create backup of existing checkpoint
        backup_file = self.state_file.with_suffix('.checkpoint.bak')
        if self.state_file.exists():
            try:
                backup_file.write_text(self.state_file.read_text())
            except Exception as e:
                self.logger.warning(f"Could not create checkpoint backup: {e}")
        
        # Save new checkpoint
        try:
            checkpoint_data = asdict(self.state)
            # Add timestamp for debugging
            checkpoint_data['_saved_at'] = datetime.now().isoformat()
            
            with open(self.state_file, 'w') as f:
                json.dump(checkpoint_data, f, indent=2, default=str)
            
            self.logger.debug(f"Checkpoint saved to {self.state_file}")
            
        except Exception as e:
            self.logger.error(f"Failed to save checkpoint: {e}")
            # Try to restore from backup
            if backup_file.exists():
                try:
                    self.state_file.write_text(backup_file.read_text())
                    self.logger.info("Restored checkpoint from backup")
                except Exception:
                    self.logger.error("Could not restore from backup")
            raise PipelineError(f"Failed to save checkpoint: {e}")
    
    def _cleanup_temp_files(self) -> None:
        """Clean up any temporary files created during pipeline execution."""
        for temp_file in self._temp_files:
            try:
                if temp_file.exists():
                    temp_file.unlink()
                    self.logger.debug(f"Cleaned up temp file: {temp_file}")
            except Exception as e:
                self.logger.warning(f"Could not clean up temp file {temp_file}: {e}")
        self._temp_files.clear()
    
    def _finalize_outputs(self) -> None:
        """Finalize outputs: only checkpoint and temp cleanup; packaging handled by ecc_packager."""
        # Copy log files if they exist
        for log_file in self.temp_dir.glob("*.log"):
            shutil.copy2(log_file, self.final_output_dir / log_file.name)

        # Remove checkpoint/config artifacts from final output
        for cleanup_path in [
            self.state_file,
            self.state_file.with_suffix(self.state_file.suffix + ".bak"),
            self.final_output_dir / "config.yaml",
        ]:
            try:
                if cleanup_path and cleanup_path.exists():
                    cleanup_path.unlink()
                    self.logger.debug(f"Removed auxiliary file: {cleanup_path}")
            except Exception as exc:
                self.logger.warning(f"Could not remove auxiliary file {cleanup_path}: {exc}")

        # Clean up temp directory based on configuration
        if not self.config.keep_tmp:
            self.logger.info(f"Cleaning up temporary directory: {self.temp_dir}")
            try:
                shutil.rmtree(self.temp_dir)
            except Exception as e:
                self.logger.warning(f"Could not completely remove temp directory: {e}")
        else:
            self.logger.info(f"Temporary files retained in: {self.temp_dir}")

    def show_steps(self, detailed: bool = False) -> None:
        """Show pipeline steps with enhanced status information."""
        import click
        click.echo("CircleSeeker Pipeline Steps:")
        click.echo("=" * 70)
        
        name_width = max((len(s.name) for s in self.STEPS), default=0)
        name_width = max(name_width, 16)
        idx_width = len(str(len(self.STEPS)))

        for i, step in enumerate(self.STEPS, 1):
            if step.name in self.state.completed_steps:
                status = "✓"
                status_color = "completed"
            elif step.name == self.state.current_step:
                status = "◉"
                status_color = "running"
            elif step.name == self.state.failed_step:
                status = "✗"
                status_color = "failed"
            else:
                status = "○"
                status_color = "pending"
            
            # Simplified output: aligned columns, no skippable info
            click.echo(f"{status} Step {i:{idx_width}d}: {step.name:<{name_width}} - {step.description}")
            
            # Detailed info (optional)
            if detailed and step.name in self.state.step_metadata:
                meta = self.state.step_metadata[step.name]
                if meta.duration:
                    click.echo(f"    Duration: {meta.duration:.1f}s")
                if meta.error_message:
                    click.echo(f"    Error: {meta.error_message}")
                if meta.output_files:
                    click.echo(f"    Output files: {len(meta.output_files)} files")
        
        click.echo("=" * 70)
        
        # Show runtime summary (only if detailed)
        if detailed and self.state.pipeline_start_time:
            total_time = self.state.get_total_runtime()
            if total_time:
                click.echo(f"Total runtime: {total_time:.1f}s")
        
        if detailed and self.state.failed_step:
            click.echo(f"⚠️  Previous run failed at: {self.state.failed_step}")
            click.echo("   Use --force to restart from beginning, or resume will continue from failure point")
    
    def run(
        self,
        start_from: Optional[int] = None,
        stop_at: Optional[int] = None,
        force: bool = False
    ) -> Dict[str, Any]:
        """Run the pipeline with automatic cleanup and temp directory management."""
        if force:
            self.state = self._create_fresh_state()
            self._save_state()
        
        try:
            # Determine step range
            start_idx = (start_from - 1) if start_from else 0
            end_idx = stop_at if stop_at else len(self.STEPS)
            
            steps_slice = self.STEPS[start_idx:end_idx]
            iterator = iter_progress(
                steps_slice,
                total=len(steps_slice),
                desc="Process",
                enabled=getattr(self.config.runtime, "enable_progress", True),
            )
            for i, step in enumerate(iterator, start_idx + 1):
                if step.name in self.state.completed_steps and not force:
                    self.logger.info(f"Skipping completed step {i}: {step.name}")
                    continue
                
                # Check skip condition
                if step.skip_condition and getattr(self.config, step.skip_condition, False):
                    self.logger.info(f"Skipping step {i}: {step.name} (user requested)")
                    continue
                
                self.logger.info(f"Running step {i}: {step.name}")
                self.state.current_step = step.name
                self.state.add_step_metadata(step.name)
                self._save_state()
                
                # Execute step
                step_start_time = time.time()
                try:
                    output_files = self._execute_step(step)
                    
                    # Mark as completed
                    self.state.complete_step(step.name, output_files or [])
                    self.state.current_step = None
                    self.state.failed_step = None  # Clear any previous failure
                    
                    step_duration = time.time() - step_start_time
                    self.logger.info(f"Completed step {i}: {step.name} ({step_duration:.1f}s)")
                    
                except Exception as e:
                    # Handle step failure
                    error_msg = str(e)
                    self.state.fail_step(step.name, error_msg)
                    self.state.current_step = None
                    
                    step_duration = time.time() - step_start_time
                    self.logger.error(f"Step {i} failed: {step.name} ({step_duration:.1f}s)")
                    self.logger.error(f"Error: {error_msg}")
                    
                    # Save state with failure info
                    self._save_state()
                    
                    # Re-raise the exception
                    raise PipelineError(f"Pipeline failed at step {step.name}: {error_msg}") from e
                
                # Save successful completion
                self._save_state()
            
            # After successful completion, finalize outputs unless skipped
            if not getattr(self.config, 'skip_organize', False):
                self._finalize_outputs()
                self.logger.info("Outputs organized to final directories")
            else:
                self.logger.info("Skipping output organization (--skip-organize specified)")
            
            self.logger.info("Pipeline completed successfully")
            return self.state.results
            
        except PipelineError as e:
            # Keep temp directory for debugging on failure
            self.logger.info(f"Pipeline failed. Temporary files retained in: {self.temp_dir}")
            raise
        except Exception as e:
            if self.state.current_step:
                self.state.fail_step(self.state.current_step, str(e))
            self._save_state()
            self.logger.info(f"Unexpected failure. Temporary files retained in: {self.temp_dir}")
            raise PipelineError(f"Unexpected pipeline failure: {e}") from e
        finally:
            # Always clean up temporary files
            self._cleanup_temp_files()
    
    def _execute_step(self, step: PipelineStep) -> Optional[List[str]]:
        """Execute a pipeline step and return output files."""
        method_name = f"_step_{step.name}"
        if hasattr(self, method_name):
            method = getattr(self, method_name)
            result = method()
            
            # If method returns a list of output files, use it
            if isinstance(result, list):
                return result
            
            return None
        else:
            raise PipelineError(f"Step implementation not found: {method_name}")
    
    # ===================== STEP IMPLEMENTATIONS =====================

    # New step name wrappers for backward compatibility
    # Removed alias methods - using direct step names now
    
    def _step_make_blastdb(self) -> None:
        """Step 0: Build BLAST database."""
        db_prefix = self.config.output_dir / f"{self.config.prefix}_blast_db"
        
        makeblastdb = MakeBlastDB()
        if makeblastdb.verify_database(db_prefix):
            self.logger.info(f"BLAST database already exists: {db_prefix}")
        else:
            makeblastdb.build_database(
                input_file=self.config.reference,
                output_db=db_prefix,
                dbtype=self.config.tools.blast.get("dbtype", "nucl"),
                title=f"{self.config.prefix}_reference"
            )
        
        self.state.results[ResultKeys.BLAST_DB] = str(db_prefix)
    
    def _step_tidehunter(self) -> None:
        """Step 1: Run TideHunter."""
        
        output_file = self.config.output_dir / f"{self.config.prefix}.TH.ecc_candidates.txt"
        
        tidehunter = TideHunter(threads=self.config.threads)
        tidehunter.run_analysis(
            input_file=self.config.input_file,
            output_file=output_file,
            **self.config.tools.tidehunter
        )
        
        self.state.results[ResultKeys.TIDEHUNTER_OUTPUT] = str(output_file)
    
    def _step_tandem_to_ring(self) -> None:
        """Step 2: Process TideHunter output using tandem_to_ring module."""
        tidehunter_output = self.state.results.get(ResultKeys.TIDEHUNTER_OUTPUT)
        if not tidehunter_output:
            raise PipelineError("TideHunter output not found in state")

        csv_output = self.config.output_dir / "tandem_to_ring.csv"
        fasta_output = self.config.output_dir / "tandem_to_ring.fasta"

        # Use Python import and class
        processor = TandemToRing(
            input_file=tidehunter_output,
            output_file=csv_output,
            circular_fasta=fasta_output,
            logger=self.logger.getChild("tandem_to_ring")
        )

        self.logger.info("Running tandem_to_ring module")
        processor.run()

        self._set_result(ResultKeys.T2R_CSV, str(csv_output))
        self._set_result(ResultKeys.T2R_FASTA, str(fasta_output))
    
    def _step_run_blast(self) -> None:
        """Step 3: Run BLAST alignment."""
        blast_output = self.config.output_dir / f"{self.config.prefix}_blast_results.tsv"
        
        query_file = Path(self._get_result(ResultKeys.T2R_FASTA, default=self.config.input_file))
        db_prefix = Path(self.state.results[ResultKeys.BLAST_DB])
        
        blast_runner = BlastRunner(
            num_threads=self.config.threads,
            **self.config.tools.blast
        )
        blast_runner.run(
            database=db_prefix,
            query_file=query_file,
            output_file=blast_output
        )
        
        self.state.results[ResultKeys.BLAST_OUTPUT] = str(blast_output)
    
    def _step_um_classify(self) -> None:
        """Step 4: Classify eccDNA types using um_classify module."""
        import pandas as pd

        blast_output = Path(self.state.results[ResultKeys.BLAST_OUTPUT])
        output_prefix = self.config.output_dir / "um_classify"

        # Use Python import and class
        classifier = UMeccClassifier(
            logger=self.logger.getChild("um_classify")
        )

        self.logger.info("Running um_classify module")
        # Read blast results
        import pandas as pd
        blast_df = pd.read_csv(blast_output, sep='\t', header=None)
        uecc_df, mecc_df, unclassified_df = classifier.classify_blast_results(
            blast_df, output_prefix
        )

        # Check for output files
        uecc_output = self.config.output_dir / "um_classify.uecc.csv"
        mecc_output = self.config.output_dir / "um_classify.mecc.csv"
        unclass_output = self.config.output_dir / "um_classify.unclassified.csv"

        if uecc_output.exists():
            uecc_df = pd.read_csv(uecc_output)
            self.state.results[ResultKeys.UECC_CSV] = str(uecc_output)
            self.state.results[ResultKeys.UECC_COUNT] = len(uecc_df)

        if mecc_output.exists():
            mecc_df = pd.read_csv(mecc_output)
            self.state.results[ResultKeys.MECC_CSV] = str(mecc_output)
            self.state.results[ResultKeys.MECC_COUNT] = len(mecc_df)

        if unclass_output.exists():
            unclassified_df = pd.read_csv(unclass_output)
            self.state.results[ResultKeys.UNCLASSIFIED_CSV] = str(unclass_output)
            # Don't prepare trapeze data here, let cecc_build handle it
    
    def _prepare_trapeze_data(self, unclassified_df) -> Any:
        """Prepare unclassified BLAST data for trapeze step."""
        import pandas as pd
        
        cecc_df = unclassified_df.copy()
        
        if 'sstrand' in cecc_df.columns:
            cecc_df['strand'] = cecc_df['sstrand']
        
        carousel_file = self.config.output_dir / f"{self.config.prefix}_processed.csv"
        if carousel_file.exists():
            try:
                carousel_df = pd.read_csv(carousel_file)
                cecc_df['readName'] = cecc_df['query_id'].str.replace('_reprep0', '', regex=False)
                
                available_cols = ['readName']
                if 'consLen' in carousel_df.columns:
                    available_cols.append('consLen')
                if 'copyNum' in carousel_df.columns:
                    available_cols.append('copyNum')
                
                carousel_subset = carousel_df[available_cols].drop_duplicates(subset=['readName'])
                cecc_df = cecc_df.merge(carousel_subset, on='readName', how='left')
                
                if 'consLen' not in cecc_df.columns:
                    cecc_df['consLen'] = cecc_df['alignment_length']
                if 'copyNum' not in cecc_df.columns:
                    cecc_df['copyNum'] = 1.0
                    
            except Exception as e:
                self.logger.warning(f"Could not merge carousel data: {e}")
                cecc_df['readName'] = cecc_df['query_id'].str.replace('_reprep0', '', regex=False)
                cecc_df['consLen'] = cecc_df['alignment_length']
                cecc_df['copyNum'] = 1.0
        else:
            cecc_df['readName'] = cecc_df['query_id'].str.replace('_reprep0', '', regex=False)
            cecc_df['consLen'] = cecc_df['alignment_length'] 
            cecc_df['copyNum'] = 1.0
        
        required_cols = ["query_id", "subject_id", "q_start", "q_end", "s_start", "s_end",
                        "strand", "alignment_length", "consLen", "readName", "copyNum"]
        
        for col in required_cols:
            if col not in cecc_df.columns:
                cecc_df[col] = 'plus' if col == 'strand' else 0
        
        return cecc_df[required_cols]
    
    def _step_cecc_build(self) -> None:
        """Step 6: Process Cecc candidates using cecc_build module."""
        import pandas as pd

        unclass_input = self.config.output_dir / "um_classify.unclassified.csv"

        if not unclass_input.exists():
            self.logger.info("No unclassified candidates found for CECC")
            output_file = self.config.output_dir / "cecc_build.csv"
            pd.DataFrame().to_csv(output_file, index=False)
            self._set_result(ResultKeys.CECC_BUILD_OUTPUT, str(output_file))
            self._set_result(ResultKeys.CECC_BUILD_COUNT, 0)
            return

        output_file = self.config.output_dir / "cecc_build.csv"

        # Use Python import and class
        builder = CeccBuild(
            logger=self.logger.getChild("cecc_build")
        )

        self.logger.info("Running cecc_build module")
        try:
            # Call run_pipeline directly with proper parameters
            result_df = builder.run_pipeline(
                input_csv=unclass_input,
                output_csv=output_file,
                overlap_threshold=0.95,
                min_segments=2,
                edge_tolerance=20,
                position_tolerance=50
            )
            if output_file.exists():
                df = pd.read_csv(output_file)
                self._set_result(ResultKeys.CECC_BUILD_OUTPUT, str(output_file))
                self._set_result(ResultKeys.CECC_BUILD_COUNT, len(df))
            else:
                self._set_result(ResultKeys.CECC_BUILD_OUTPUT, str(output_file))
                self._set_result(ResultKeys.CECC_BUILD_COUNT, 0)
        except Exception as e:
            self.logger.warning(f"cecc_build failed: {e}")
            # Create empty file and continue
            pd.DataFrame().to_csv(output_file, index=False)
            self._set_result(ResultKeys.CECC_BUILD_OUTPUT, str(output_file))
            self._set_result(ResultKeys.CECC_BUILD_COUNT, 0)
    
    def _step_umc_process(self) -> None:
        """Step 7: Generate FASTA files using umc_process module."""
        circular_fasta = self.config.output_dir / "tandem_to_ring.fasta"
        if not circular_fasta.exists():
            # Try alternative name
            circular_fasta = self.config.output_dir / f"{self.config.prefix}_circular.fasta"
            if not circular_fasta.exists():
                self.logger.warning("Circular FASTA file not found, skipping umc_process step")
                return

        uecc_csv = self.config.output_dir / "um_classify.uecc.csv"
        mecc_csv = self.config.output_dir / "um_classify.mecc.csv"
        cecc_csv = self.config.output_dir / "cecc_build.csv"
        output_prefix = self.config.output_dir / "umc_process"

        processor = UMCProcess(
            logger=self.logger.getChild("umc_process")
        )

        self.logger.info("Running umc_process module")
        processor.run(
            fasta_file=circular_fasta,
            uecc_csv=uecc_csv if uecc_csv.exists() else None,
            mecc_csv=mecc_csv if mecc_csv.exists() else None,
            cecc_csv=cecc_csv if cecc_csv.exists() else None,
            output_dir=self.config.output_dir,
            prefix=self.config.prefix
        )

        # Check for output files
        type_mapping = {'Uecc': 'uecc', 'Mecc': 'mecc', 'Cecc': 'cecc'}

        for ecc_name, ecc_type in type_mapping.items():
            fasta_file = (self.config.output_dir / f"{self.config.prefix}_{ecc_name}DNA_pre.fasta")
            if not fasta_file.exists():
                continue

            self.state.results[f"{ecc_type}_fasta"] = str(fasta_file)

            # UMCProcess writes *_processed.csv files; fall back to *.csv for legacy runs
            processed_csv = (self.config.output_dir / f"{self.config.prefix}_{ecc_name}DNA_processed.csv")
            if not processed_csv.exists():
                fallback_csv = fasta_file.with_suffix('.csv')
                processed_csv = fallback_csv if fallback_csv.exists() else processed_csv

            if processed_csv.exists():
                self.state.results[f"{ecc_type}_processed"] = str(processed_csv)

            self.logger.debug(
                "Found %s outputs: fasta=%s csv=%s",
                ecc_type,
                fasta_file,
                processed_csv if processed_csv.exists() else "missing"
            )
    
    def _step_cd_hit(self) -> None:
        """Step 8: Remove redundant sequences using CD-HIT-EST wrapper."""

        # Process each eccDNA type independently
        fasta_specs = [
            ("uecc", f"{self.config.prefix}_UeccDNA_pre.fasta", f"{self.config.prefix}_U"),
            ("mecc", f"{self.config.prefix}_MeccDNA_pre.fasta", f"{self.config.prefix}_M"),
            ("cecc", f"{self.config.prefix}_CeccDNA_pre.fasta", f"{self.config.prefix}_C"),
        ]

        # Initialize CD-HIT-EST tool
        cd_hit = CDHitEst(logger=self.logger.getChild("cd_hit"), threads=self.config.threads)

        for ecc_type, input_name, output_prefix in fasta_specs:
            input_fasta = (self.config.output_dir / input_name)
            if not input_fasta.exists():
                self.logger.warning(f"Skipping cd-hit for {input_name}: file not found")
                continue

            output_path = (self.config.output_dir / output_prefix)

            try:
                clstr_path = cd_hit.cluster_sequences(input_fasta, output_path)
            except Exception as e:
                self.logger.warning(f"cd-hit failed for {input_name}: {e}")
                continue

            # Store results
            output_fasta = output_path.with_suffix('.fasta')
            cluster_file = clstr_path if clstr_path and Path(clstr_path).exists() else output_path.with_suffix('.clstr')

            if output_fasta.exists():
                self.state.results[f"{ecc_type}_nr99"] = str(output_fasta)
            if cluster_file and Path(cluster_file).exists():
                self.state.results[f"{ecc_type}_clusters"] = str(cluster_file)
    
    def _step_ecc_dedup(self) -> None:
        """Step 9: Coordinate and deduplicate results using ecc_dedup module."""
        
        harmonizer = eccDedup(self.logger.getChild("ecc_dedup"))
        
        uecc_input = Path(self.state.results[ResultKeys.UECC_PROCESSED]) if ResultKeys.UECC_PROCESSED in self.state.results else None
        uecc_cluster = Path(self.state.results[ResultKeys.UECC_CLUSTERS]) if ResultKeys.UECC_CLUSTERS in self.state.results else None
        mecc_input = Path(self.state.results[ResultKeys.MECC_PROCESSED]) if ResultKeys.MECC_PROCESSED in self.state.results else None
        mecc_cluster = Path(self.state.results[ResultKeys.MECC_CLUSTERS]) if ResultKeys.MECC_CLUSTERS in self.state.results else None
        cecc_input = Path(self.state.results[ResultKeys.CECC_PROCESSED]) if ResultKeys.CECC_PROCESSED in self.state.results else None
        cecc_cluster = Path(self.state.results[ResultKeys.CECC_CLUSTERS]) if ResultKeys.CECC_CLUSTERS in self.state.results else None
        
        results = harmonizer.run_deduplication(
            output_dir=self.config.output_dir,
            prefix=self.config.prefix,
            uecc_input=uecc_input if uecc_input and uecc_input.exists() else None,
            uecc_cluster=uecc_cluster if uecc_cluster and uecc_cluster.exists() else None,
            mecc_input=mecc_input if mecc_input and mecc_input.exists() else None,
            mecc_cluster=mecc_cluster if mecc_cluster and mecc_cluster.exists() else None,
            cecc_input=cecc_input if cecc_input and cecc_input.exists() else None,
            cecc_cluster=cecc_cluster if cecc_cluster and cecc_cluster.exists() else None
        )

        # Fix inconsistent naming by renaming files if needed
        rename_mappings = [
            (f"{self.config.prefix}_Mecc.fa", f"{self.config.prefix}_MeccDNA_C.fasta"),
            (f"{self.config.prefix}_Cecc.fa", f"{self.config.prefix}_CeccDNA_C.fasta"),
            (f"{self.config.prefix}_UeccDNA.fa", f"{self.config.prefix}_UeccDNA_C.fasta"),
        ]

        for old_name, new_name in rename_mappings:
            old_path = self.config.output_dir / old_name
            new_path = self.config.output_dir / new_name
            if old_path.exists() and not new_path.exists():
                old_path.rename(new_path)
                self.logger.debug(f"Renamed {old_name} to {new_name}")

        umc_dirs: dict[str, Path] = {}
        if results:
            prefix_for_dirs = self.config.prefix or "sample"
            try:
                base_lists = {
                    "U": [
                        "UeccDNA_C.fasta",
                        "UeccDNA.bed",
                        "UeccDNA.core.csv",
                    ],
                    "M": [
                        "MeccDNA_C.fasta",
                        "MeccSites.bed",
                        "MeccBestSite.bed",
                        "MeccSites.core.csv",
                    ],
                    "C": [
                        "CeccDNA_C.fasta",
                        "CeccSegments.bed",
                        "CeccSegments.core.csv",
                        "CeccJunctions.bedpe",
                    ],
                }

                def with_prefix(name: str) -> str:
                    return f"{prefix_for_dirs}_{name}" if prefix_for_dirs else name

                files_to_move: dict[str, list[Path]] = {}
                for ecc_code, file_list in base_lists.items():
                    collected: list[Path] = []
                    for base_name in file_list:
                        candidate = self.config.output_dir / with_prefix(base_name)
                        if candidate.exists():
                            collected.append(candidate)
                        elif not prefix_for_dirs:
                            fallback = self.config.output_dir / base_name
                            if fallback.exists():
                                collected.append(fallback)
                    if collected:
                        files_to_move[ecc_code] = collected

                umc_dirs = organize_umc_files(
                    self.config.output_dir,
                    prefix_for_dirs,
                    umc_files=files_to_move if files_to_move else None,
                    auto_detect=False,
                )
                dir_key_map = {'U': 'uecc', 'M': 'mecc', 'C': 'cecc'}
                for key, directory in umc_dirs.items():
                    mapped = dir_key_map.get(key.upper())
                    if mapped and directory.exists():
                        self.state.results[f"ecc_dedup_{mapped}_dir"] = str(directory)
            except Exception as exc:
                self.logger.warning(f"Failed to organize eccDedup outputs: {exc}")
                umc_dirs = {}

        def _locate_output(filename: str, ecc_code: str | None = None) -> Optional[Path]:
            candidate = self.config.output_dir / filename
            if candidate.exists():
                return candidate
            if ecc_code:
                code = ecc_code.upper()
                target_dir = umc_dirs.get(code)
                if target_dir:
                    fallback = target_dir / filename
                    if fallback.exists():
                        return fallback
            for fallback in self.config.output_dir.rglob(filename):
                if fallback.exists():
                    return fallback
            return None

        folder_key_map = {'uecc': 'U', 'mecc': 'M', 'cecc': 'C'}

        # Store harmonizer output paths with consistent naming
        fasta_mappings = {
            'uecc': f"{self.config.prefix}_UeccDNA_C.fasta",
            'mecc': f"{self.config.prefix}_MeccDNA_C.fasta",
            'cecc': f"{self.config.prefix}_CeccDNA_C.fasta"
        }

        for ecc_type, filename in fasta_mappings.items():
            fasta_file = _locate_output(filename, folder_key_map.get(ecc_type))
            if fasta_file:
                self._set_result(f"ecc_dedup_{ecc_type}_fasta", str(fasta_file))
                self.logger.debug(f"Found harmonizer output: {fasta_file}")
            else:
                self.logger.warning(f"eccDedup output not found: {filename}")
        
        # Store harmonized CSV files
        for ecc_type in ['Uecc', 'Mecc', 'Cecc']:
            harmonized_csv = self.config.output_dir / f"{self.config.prefix}_{ecc_type}_harmonized.csv"
            if harmonized_csv.exists():
                key = {
                    'Uecc': ResultKeys.UECC_HARMONIZED,
                    'Mecc': ResultKeys.MECC_HARMONIZED,
                    'Cecc': ResultKeys.CECC_HARMONIZED,
                }[ecc_type]
                self.state.results[key] = str(harmonized_csv)
        
        # Store other output files (BED, BEDPE, etc.)
        additional_files = {
            'uecc': [
                (f"{self.config.prefix}_UeccDNA.bed", "uecc_bed"),
                (f"{self.config.prefix}_UeccDNA.core.csv", "uecc_core_csv"),
            ],
            'mecc': [
                (f"{self.config.prefix}_MeccSites.bed", "mecc_sites_bed"),
                (f"{self.config.prefix}_MeccBestSite.bed", "mecc_bestsite_bed"),
                (f"{self.config.prefix}_MeccSites.core.csv", "mecc_core_csv"),
            ],
            'cecc': [
                (f"{self.config.prefix}_CeccJunctions.bedpe", "cecc_junctions_bedpe"),
                (f"{self.config.prefix}_CeccSegments.bed", "cecc_segments_bed"),
                (f"{self.config.prefix}_CeccSegments.core.csv", "cecc_core_csv"),
            ]
        }
        
        for ecc_type, file_list in additional_files.items():
            for filename, key in file_list:
                file_path = _locate_output(filename, folder_key_map.get(ecc_type))
                if file_path:
                    self.state.results[key] = str(file_path)
                    self.logger.debug(f"Found {key}: {file_path}")

        confirmed_file = self.config.output_dir / f"{self.config.prefix}_eccDNA_Confirmed.csv"
        if confirmed_file.exists():
            self.state.results[ResultKeys.CONFIRMED_CSV] = self._serialize_path_for_state(confirmed_file)

    def _step_read_filter(self) -> None:
        """Step 10: Filter sequences using read_filter module."""
        
        classification_csv = Path(self._get_result(ResultKeys.T2R_CSV, default=self.config.output_dir / f"{self.config.prefix}_processed.csv"))
        
        if not classification_csv.exists() or classification_csv.stat().st_size == 0:
            self.logger.warning("TandemToRing classification CSV not found; skipping sieve step")
            return
        
        original_fasta = self.config.input_file
        if not original_fasta.exists():
            self.logger.error("Original input FASTA file not found; cannot run sieve step")
            return
        
        output_fasta = self.config.output_dir / f"{self.config.prefix}_all_filtered.fasta"
        
        sieve = Sieve(self.logger.getChild("read_filter"))
        stats = sieve.run_sieve(
            carousel_csv=classification_csv,
            input_fastas=[original_fasta],
            output_fasta=output_fasta
        )
        
        self._set_result(ResultKeys.READ_FILTER_OUTPUT_FASTA, str(output_fasta))
        self._set_result("read_filter_total_reads", stats.total_reads)
        self._set_result("read_filter_filtered_reads", stats.filtered_reads)
        self._set_result("read_filter_retained_reads", stats.retained_reads)
        
        self.logger.info(
            f"Sieve complete: {stats.filtered_reads} filtered, "
            f"{stats.retained_reads} retained ({stats.filtered_percentage:.2f}% filtered)"
        )
    
    def _step_minimap2(self) -> None:
        """Step 11: Map filtered sequences to reference."""
        
        # Use consistent configuration
        config = MiniMapConfig(
            preset="map-hifi",
            threads=self.config.threads,
            output_format="bam",
            allow_secondary=True,
            build_index=True,
            sort_bam=True,
            index_bam=True
        )
        
        minimap2 = Minimap2(config=config, logger=self.logger.getChild("minimap2"))
        
        # Determine input FASTA with robust fallback
        sieve_output = self._get_result(ResultKeys.READ_FILTER_OUTPUT_FASTA)
        if sieve_output and Path(sieve_output).exists() and Path(sieve_output).stat().st_size > 0:
            input_fasta = Path(sieve_output)
        else:
            # Create combined FASTA as fallback
            input_fasta = self.config.output_dir / f"{self.config.prefix}_all_filtered.fasta"
            try:
                seen_fastas: set[Path] = set()
                with open(input_fasta, 'w') as out_f:
                    for ecc_type in ['uecc', 'mecc', 'cecc']:
                        candidate_keys = [
                            f"ecc_dedup_{ecc_type}_fasta",
                            f"{ecc_type}_fasta",
                            f"harmonizer_{ecc_type}_fasta",
                        ]
                        for key in candidate_keys:
                            fasta_path_str = self.state.results.get(key)
                            if not fasta_path_str:
                                continue
                            fasta_path = Path(fasta_path_str)
                            if fasta_path in seen_fastas:
                                continue
                            if fasta_path.exists() and fasta_path.stat().st_size > 0:
                                with open(fasta_path, 'r') as in_f:
                                    out_f.write(in_f.read())
                                seen_fastas.add(fasta_path)
            except Exception as e:
                self.logger.error(f"Failed to create combined FASTA: {e}")
                raise

        output_bam = self.config.output_dir / f"{self.config.prefix}_sorted.bam"
        
        # Run alignment
        minimap2.align(
            reference=self.config.reference,
            query=input_fasta,
            output_file=output_bam
        )
        
        self._set_result(ResultKeys.MINIMAP2_BAM, str(output_bam))
        self._set_result("all_filtered_fasta", str(input_fasta))

    def _step_cyrcular_calling(self) -> None:
        """Post-alignment circular DNA detection via Cyrcular and Varlociraptor."""
        bam_file = self.state.results.get(ResultKeys.MINIMAP2_BAM)
        if not bam_file or not Path(bam_file).exists():
            self.logger.warning("No BAM file from minimap2, skipping Cyrcular-Calling step")
            return

        cfg = CCConfig(
            bam_file=Path(bam_file),
            reference=self.config.reference,
            output_dir=self.config.output_dir,
            sample_name=self.config.prefix,
            threads=self.config.threads,
        )

        pipeline = CyrcularCallingPipeline(cfg, logger=self.logger.getChild("cyrcular_calling"))
        results = pipeline.run()

        # Store key outputs
        overview = pipeline.file_paths.get("overview_table")
        if overview and overview.exists():
            serialized = self._serialize_path_for_state(overview)
            self.state.results[ResultKeys.CYRCULAR_OVERVIEW] = serialized
        else:
            self.logger.warning("Cyrcular overview table not found; inferred eccDNA curation may be skipped")
        self.state.results[ResultKeys.CYRCULAR_RESULT_COUNT] = len(results) if results else 0

    def _step_curate_inferred_ecc(self) -> list[str] | None:
        """Step 12: Curate inferred eccDNA simple/chimeric tables from Cyrcular overview TSV."""
        import pandas as pd

        overview_path_str = self.state.results.get(ResultKeys.CYRCULAR_OVERVIEW)

        overview_path_obj = self._resolve_stored_path(
            overview_path_str,
            [f"{self.config.prefix}_overview.tsv"],
        )
        if not overview_path_obj:
            raise PipelineError("Cyrcular overview table missing; cannot curate inferred eccDNA")

        # Update stored path using relative serialization for checkpoints
        self.state.results[ResultKeys.CYRCULAR_OVERVIEW] = self._serialize_path_for_state(overview_path_obj)

        simple_df, chimeric_df = curate_ecc_tables(overview_path_obj)

        inferred_prefix = self.config.output_dir / self.config.prefix
        reference_for_inferred: Optional[Path] = None
        if self.config.reference:
            ref_path = Path(self.config.reference)
            if ref_path.exists():
                reference_for_inferred = ref_path
            else:
                self.logger.warning(f"Reference FASTA not found for inferred eccDNA FASTA generation: {ref_path}")

        simple_csv_path: Optional[Path] = None
        chimeric_csv_path: Optional[Path] = None
        simple_fasta_path: Optional[Path] = None
        chimeric_fasta_path: Optional[Path] = None

        try:
            simple_csv_path, chimeric_csv_path, simple_fasta_path, chimeric_fasta_path = write_curated_tables_with_fasta(
                simple_df,
                chimeric_df,
                inferred_prefix,
                reference_fasta=reference_for_inferred,
                organize_files=True,
            )
        except Exception as exc:
            self.logger.warning(f"Failed to generate inferred eccDNA tables with FASTA: {exc}")
            simple_csv_path, chimeric_csv_path = write_curated_tables(
                simple_df,
                chimeric_df,
                inferred_prefix,
            )

        written: List[str] = []
        if simple_csv_path and simple_csv_path.exists() and simple_csv_path.stat().st_size > 0:
            self.state.results[ResultKeys.INFERRED_SIMPLE_TSV] = self._serialize_path_for_state(simple_csv_path)
            written.append(str(simple_csv_path))
        else:
            self.logger.info("No inferred simple eccDNA records; simple TSV not created")

        if chimeric_csv_path and chimeric_csv_path.exists() and chimeric_csv_path.stat().st_size > 0:
            self.state.results[ResultKeys.INFERRED_CHIMERIC_TSV] = self._serialize_path_for_state(chimeric_csv_path)
            written.append(str(chimeric_csv_path))
        else:
            self.logger.info("No inferred chimeric eccDNA records; chimeric TSV not created")

        if simple_fasta_path and simple_fasta_path.exists() and simple_fasta_path.stat().st_size > 0:
            self.state.results[ResultKeys.INFERRED_SIMPLE_FASTA] = self._serialize_path_for_state(simple_fasta_path)

        if chimeric_fasta_path and chimeric_fasta_path.exists() and chimeric_fasta_path.stat().st_size > 0:
            self.state.results[ResultKeys.INFERRED_CHIMERIC_FASTA] = self._serialize_path_for_state(chimeric_fasta_path)

        inferred_dir_path = None
        if simple_csv_path:
            inferred_dir_path = simple_csv_path.parent
        elif chimeric_csv_path:
            inferred_dir_path = chimeric_csv_path.parent
        if inferred_dir_path and inferred_dir_path.exists():
            self.state.results["inferred_dir"] = self._serialize_path_for_state(inferred_dir_path)

        if not simple_df.empty and not (simple_fasta_path and simple_fasta_path.exists() and simple_fasta_path.stat().st_size > 0):
            raise PipelineError(
                "Inferred simple eccDNA FASTA is required but was not generated. Please ensure pysam is installed and the reference genome is accessible."
            )

        if not chimeric_df.empty and not (chimeric_fasta_path and chimeric_fasta_path.exists() and chimeric_fasta_path.stat().st_size > 0):
            raise PipelineError(
                "Inferred chimeric eccDNA FASTA is required but was not generated. Please ensure pysam is installed and the reference genome is accessible."
            )

        return written or None

    def _step_ecc_unify(self) -> None:
        """Step 13: Merge eccDNA tables into unified output."""
        import pandas as pd

        # Prepare input files
        confirmed_csv_path = self._resolve_stored_path(
            self.state.results.get(ResultKeys.CONFIRMED_CSV),
            [
                f"{self.config.prefix}_confirmed.csv",
                f"{self.config.prefix}_eccDNA_Confirmed.csv",
            ],
        )

        if not confirmed_csv_path:
            self.logger.info("No confirmed eccDNA file found, skipping ecc_unify step")
            return

        # Persist resolved confirmed path
        self.state.results[ResultKeys.CONFIRMED_CSV] = self._serialize_path_for_state(confirmed_csv_path)

        simple_csv_path = self._resolve_stored_path(
            self.state.results.get(ResultKeys.INFERRED_SIMPLE_TSV),
            [
                f"{self.config.prefix}_inferred_simple.csv",
                Path(f"{self.config.prefix}_Inferred_eccDNA") / f"{self.config.prefix}_simple.csv",
            ],
        )

        chimeric_csv_path = self._resolve_stored_path(
            self.state.results.get(ResultKeys.INFERRED_CHIMERIC_TSV),
            [
                f"{self.config.prefix}_inferred_chimeric.csv",
                Path(f"{self.config.prefix}_Inferred_eccDNA") / f"{self.config.prefix}_chimeric.csv",
            ],
        )

        # Output files
        unified_csv = self.config.output_dir / f"{self.config.prefix}_unified.csv"
        overlap_report = self.config.output_dir / f"{self.config.prefix}_overlap_report.txt"
        overlap_stats = self.config.output_dir / f"{self.config.prefix}_overlap_stats.json"

        self.logger.info("Running ecc_unify module")
        try:
            # Call the merge function directly
            merged_df, report_text, stats_dict = merge_eccdna_tables(
                confirmed_file=str(confirmed_csv_path),
                inferred_simple=str(simple_csv_path) if simple_csv_path and simple_csv_path.exists() else None,
                inferred_chimeric=str(chimeric_csv_path) if chimeric_csv_path and chimeric_csv_path.exists() else None,
                overlap_report_file=overlap_report,
                overlap_stats_json=overlap_stats
            )

            # Save unified CSV
            merged_df.to_csv(unified_csv, index=False)
            self.logger.info(f"Unified {len(merged_df)} eccDNA entries")

        except Exception as e:
            self.logger.error(f"ecc_unify failed: {e}")
            return

        # Store results
        if unified_csv.exists():
            self.state.results[ResultKeys.UNIFIED_CSV] = self._serialize_path_for_state(unified_csv)
        if overlap_report.exists():
            self.state.results[ResultKeys.OVERLAP_REPORT] = self._serialize_path_for_state(overlap_report)
        if overlap_stats.exists():
            self.state.results[ResultKeys.OVERLAP_STATS] = self._serialize_path_for_state(overlap_stats)

    def _step_ecc_summary(self) -> None:
        """Step 14: Generate summary report."""

        # Check required files
        processed_csv = self.config.output_dir / "tandem_to_ring.csv"
        unified_csv_path = self._resolve_stored_path(
            self.state.results.get(ResultKeys.UNIFIED_CSV),
            [f"{self.config.prefix}_unified.csv"],
        )
        confirmed_csv_path = self._resolve_stored_path(
            self.state.results.get(ResultKeys.CONFIRMED_CSV),
            [
                f"{self.config.prefix}_unified.csv",
                f"{self.config.prefix}_confirmed.csv",
                f"{self.config.prefix}_eccDNA_Confirmed.csv",
            ],
        )
        stats_json_path = self._resolve_stored_path(
            self.state.results.get(ResultKeys.OVERLAP_STATS),
            [f"{self.config.prefix}_overlap_stats.json"],
        )

        main_csv = unified_csv_path if unified_csv_path else confirmed_csv_path

        if not (processed_csv.exists() and main_csv and main_csv.exists()):
            self.logger.info("Missing required files for summary, skipping")
            return

        # Create all.fasta combining all FASTA outputs
        all_fasta = self.config.output_dir / f"{self.config.prefix}_all.fasta"
        self._create_combined_fasta(all_fasta)

        # Create empty stats file if it doesn't exist (validate path)
        stats_json = stats_json_path if stats_json_path else (self.config.output_dir / f"{self.config.prefix}_overlap_stats.json")
        if not stats_json.exists():
            try:
                stats_json.parent.mkdir(parents=True, exist_ok=True)
                import json as _json
                with open(stats_json, 'w') as f:
                    _json.dump({"overlap_stats": {}}, f)
            except Exception as e:
                self.logger.warning(f"Failed to create stats file {stats_json}: {e}")

        # Generate summary - EccSummary class exists and is correct
        output_dir = self.config.output_dir / "summary_output"
        output_dir.mkdir(exist_ok=True)

        summary = EccSummary(
            sample_name=self.config.prefix,
            output_dir=output_dir,
            logger=self.logger.getChild("ecc_summary")
        )

        # Process the data
        summary.process_fasta(self.config.input_file)
        summary.process_processed_csv(processed_csv)
        summary.process_merged_csv(main_csv)
        if stats_json.exists():
            summary.process_overlap_stats(stats_json)

        # Generate reports
        summary.generate_html_report()
        summary.generate_text_summary()

        self.state.results[ResultKeys.SUMMARY_DIR] = str(output_dir)

    def _step_ecc_packager(self) -> None:
        """Step 15: Package output files using ecc_packager module."""
        base_final_root = self.final_output_dir if self.config.output_dir.name == ".tmp_work" else self.config.output_dir
        final_output_root = base_final_root
        final_output_root.mkdir(parents=True, exist_ok=True)

        def resolve_dir(key: str, default_name: str) -> Optional[Path]:
            stored = self.state.results.get(key)
            resolved = self._resolve_stored_path(stored, [default_name])
            if resolved and resolved.exists():
                return resolved
            fallback_path = self.config.output_dir / default_name
            return fallback_path if fallback_path.exists() else None

        uecc_dir = resolve_dir("ecc_dedup_uecc_dir", f"{self.config.prefix}_Uecc_C")
        mecc_dir = resolve_dir("ecc_dedup_mecc_dir", f"{self.config.prefix}_Mecc_C")
        cecc_dir = resolve_dir("ecc_dedup_cecc_dir", f"{self.config.prefix}_Cecc_C")
        inferred_dir = self._resolve_stored_path(
            self.state.results.get("inferred_dir"),
            [f"{self.config.prefix}_Inferred_eccDNA"],
        )
        if inferred_dir and not inferred_dir.exists():
            inferred_dir = None

        # Prepare files
        merged_csv_path = self._resolve_stored_path(
            self.state.results.get(ResultKeys.UNIFIED_CSV),
            [f"{self.config.prefix}_unified.csv"],
        )
        if not merged_csv_path or not merged_csv_path.exists():
            merged_csv_path = self._resolve_stored_path(
                self.state.results.get(ResultKeys.CONFIRMED_CSV),
                [
                    f"{self.config.prefix}_unified.csv",
                    f"{self.config.prefix}_confirmed.csv",
                    f"{self.config.prefix}_eccDNA_Confirmed.csv",
                ],
            )

        merged_csv = merged_csv_path if merged_csv_path and merged_csv_path.exists() else None

        summary_dir = self.config.output_dir / "summary_output"
        html_report = summary_dir / f"{self.config.prefix}_report.html"
        html_report = html_report if html_report.exists() else None
        text_summary = summary_dir / f"{self.config.prefix}_summary.txt"
        text_summary = text_summary if text_summary.exists() else None

        packager_inputs_available = any(
            candidate is not None
            for candidate in (
                uecc_dir,
                mecc_dir,
                cecc_dir,
                inferred_dir,
                merged_csv,
                html_report,
                text_summary,
            )
        )

        if packager_inputs_available:
            self_config = self.config

            class MockArgs:
                def __init__(self) -> None:
                    self.sample_name = self_config.prefix
                    self.output_dir = str(final_output_root)
                    self.out_dir = str(final_output_root)
                    self.uecc_dir = str(uecc_dir) if uecc_dir else None
                    self.mecc_dir = str(mecc_dir) if mecc_dir else None
                    self.cecc_dir = str(cecc_dir) if cecc_dir else None
                    self.inferred_dir = str(inferred_dir) if inferred_dir else None
                    self.merged_csv = str(merged_csv) if merged_csv else None
                    self.html = str(html_report) if html_report else None
                    self.text = str(text_summary) if text_summary else None
                    self.overwrite = True
                    self.dry_run = False
                    runtime = getattr(self_config, "runtime", None)
                    debug_verbose = bool(getattr(runtime, "debug_verbose", False)) if runtime else False
                    self.verbose = debug_verbose

            mock_args = MockArgs()

            self.logger.info("Running ecc_packager module to organize final results")
            packager_success = False
            try:
                result_code = ecc_packager_module.run(mock_args)
                if result_code == 0:
                    packager_success = True
                    self.logger.info(f"Results successfully packaged in: {final_output_root}")
                else:
                    self.logger.warning("ecc_packager completed with warnings")
            except Exception as exc:
                self.logger.warning(f"ecc_packager encountered issues: {exc}")

            if packager_success:
                packaged_root = Path(mock_args.out_dir) / mock_args.sample_name
                if packaged_root.exists() and packaged_root != final_output_root:
                    for item in packaged_root.iterdir():
                        destination = final_output_root / item.name
                        try:
                            if destination.exists():
                                if destination.is_dir():
                                    shutil.rmtree(destination)
                                else:
                                    destination.unlink()
                            shutil.move(str(item), str(destination))
                        except Exception as exc:
                            self.logger.warning(f"Could not relocate {item} to {destination}: {exc}")
                    try:
                        shutil.rmtree(packaged_root)
                    except Exception as exc:
                        self.logger.warning(f"Could not remove intermediate directory {packaged_root}: {exc}")
                self._rename_inferred_simple_file(final_output_root)
            else:
                self._create_basic_output_structure(
                    final_output_root,
                    uecc_dir=uecc_dir,
                    mecc_dir=mecc_dir,
                    cecc_dir=cecc_dir,
                    inferred_dir=inferred_dir,
                    merged_csv=merged_csv,
                    html_report=html_report,
                    text_summary=text_summary,
                )
        else:
            self.logger.info("No inputs available for ecc_packager; using fallback organizer")
            self._create_basic_output_structure(
                final_output_root,
                uecc_dir=uecc_dir,
                mecc_dir=mecc_dir,
                cecc_dir=cecc_dir,
                inferred_dir=inferred_dir,
                merged_csv=merged_csv,
                html_report=html_report,
                text_summary=text_summary,
            )

        self.state.results[ResultKeys.FINAL_RESULTS] = str(final_output_root)

        # Copy Xecc artifacts (if produced) into final output root
        # Note: XeccDNA file pattern uses capital X and DNA suffix
        xecc_patterns = [f"{self.config.prefix}_XeccDNA*", f"{self.config.prefix}_Xecc*"]
        for pattern in xecc_patterns:
            for xecc_file in self.config.output_dir.glob(pattern):
                try:
                    destination = final_output_root / xecc_file.name
                    if xecc_file.is_file():
                        shutil.copy2(xecc_file, destination)
                        self.logger.info(f"Copied Xecc artifact to final output: {destination}")
                except Exception as exc:
                    self.logger.warning(f"Failed to copy Xecc artifact {xecc_file}: {exc}")

    def _step_merge_eccdna(self) -> list[str] | None:
        """Step 14 (deprecated): Merge confirmed cores with inferred tables into unified final table."""
        # This step is deprecated - functionality moved to ecc_unify
        self.logger.info("Merge eccDNA step skipped - replaced by ecc_unify")

        # Store the individual core files that would have been merged
        u_core = self.config.output_dir / f"{self.config.prefix}_UeccDNA.core.csv"
        m_core = self.config.output_dir / f"{self.config.prefix}_MeccSites.core.csv"
        c_core = self.config.output_dir / f"{self.config.prefix}_CeccSegments.core.csv"

        available_cores = []
        if u_core.exists():
            self.state.results[ResultKeys.UECC_CORE_CSV] = str(u_core)
            available_cores.append(str(u_core))
        if m_core.exists():
            self.state.results[ResultKeys.MECC_CORE_CSV] = str(m_core)
            available_cores.append(str(m_core))
        if c_core.exists():
            self.state.results[ResultKeys.CECC_CORE_CSV] = str(c_core)
            available_cores.append(str(c_core))

        if available_cores:
            self.logger.info(f"Core CSV files available: {len(available_cores)} files")
            return available_cores

        return None
    
    
    # playbill/propmaster steps removed
