"""Main pipeline orchestrator for CircleSeeker - Final fixed version."""

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
from typing import Optional, Dict, Any, List, Union
from dataclasses import dataclass, asdict, field

from circleseeker.config import Config
from circleseeker.exceptions import PipelineError
from circleseeker.external.blast import MakeBlastDB, BlastRunner


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
    
    # Define pipeline steps
    STEPS = [
        PipelineStep("make_blastdb", "Build BLAST database", skip_condition="skip_make_db"),
        PipelineStep("tidehunter", "Run TideHunter for tandem repeat detection", skip_condition="skip_tidehunter"),
        PipelineStep("tandem_to_ring", "Process TideHunter output", skip_condition="skip_carousel"),
        PipelineStep("run_blast", "Run BLAST alignment", skip_condition="skip_blast"),
        PipelineStep("um_classify", "Classify eccDNA types", skip_condition="skip_gatekeeper"),
        PipelineStep("cecc_build", "Process Cecc candidates"),
        PipelineStep("umc_process", "Generate FASTA files"),
        PipelineStep("cd_hit", "Remove redundant sequences"),
        PipelineStep("ecc_dedup", "Coordinate results"),
        PipelineStep("read_filter", "Filter sequences"),
        PipelineStep("minimap2", "Map with Minimap2"),
        PipelineStep("mosdepth", "Calculate depth coverage"),
        PipelineStep("split_refine", "Split regions"),
        PipelineStep("circle_validate", "Final processing"),
        PipelineStep("report_generator", "Generate HTML report", skip_condition="skip_report"),
    ]

    # Mapping of legacy state keys to new canonical keys
    _RESULT_KEY_ALIASES = {
        # tandem_to_ring (formerly carousel)
        "carousel_csv": "tandem_to_ring_csv",
        "carousel_fasta": "tandem_to_ring_fasta",
        # cecc_build (formerly trapeze)
        "trapeze_output": "cecc_build_output",
        "trapeze_count": "cecc_build_count",
        # ecc_dedup (formerly harmonizer)
        "harmonizer_uecc_fasta": "ecc_dedup_uecc_fasta",
        "harmonizer_mecc_fasta": "ecc_dedup_mecc_fasta",
        "harmonizer_cecc_fasta": "ecc_dedup_cecc_fasta",
        # read_filter (formerly sieve)
        "sieve_output_fasta": "read_filter_output_fasta",
        "sieve_total_reads": "read_filter_total_reads",
        "sieve_filtered_reads": "read_filter_filtered_reads",
        "sieve_retained_reads": "read_filter_retained_reads",
        # split_refine (formerly contortionist)
        "contortionist_output": "split_refine_output",
        # circle_validate (formerly juggler)
        "juggler_clean": "circle_validate_clean",
        "juggler_fasta": "circle_validate_fasta",
        "juggler_bed": "circle_validate_bed",
        # report_generator (formerly playbill)
        "playbill_report": "report_generator_report",
        "playbill_results": "report_generator_results",
    }

    def _set_result(self, key: str, value: Any) -> None:
        """Set a result key using new canonical name and mirror to legacy name.

        If `key` is a legacy name, it will be upgraded to the new canonical key
        and mirrored back to the legacy key for compatibility. If `key` is a
        new name, the legacy alias (if any) will also be set.
        """
        # If it's a legacy key, upgrade to new key name
        canonical = self._RESULT_KEY_ALIASES.get(key, None)
        if canonical:
            self.state.results[canonical] = value
            self.state.results[key] = value
            return
        # If it's a new key, also set legacy alias if exists (reverse lookup)
        rev = {v: k for k, v in self._RESULT_KEY_ALIASES.items()}
        legacy = rev.get(key)
        self.state.results[key] = value
        if legacy:
            self.state.results[legacy] = value

    def _get_result(self, *keys: str, default: Any = None) -> Any:
        """Get first available result from a list of keys (new preferred)."""
        for k in keys:
            if k in self.state.results:
                return self.state.results[k]
        return default
    
    def __init__(self, config: Config):
        self.config = config
        self.logger = logging.getLogger(self.__class__.__name__)
        
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
                    self.logger.warning("Configuration changed, checkpoint may be invalid")
                    response = input("Continue with existing checkpoint? [y/N]: ")
                    if response.lower() != 'y':
                        return self._create_fresh_state()
                
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
        """Move final results from temp directory to output directory."""
        import pandas as pd
        
        # Files to copy directly to final directory (organized results)
        final_files = [
            # FASTA and BED files (copy directly)
            (f"{self.config.prefix}_UeccDNA.fa", f"{self.config.prefix}_UeccDNA.fa"),
            (f"{self.config.prefix}_UeccDNA.bed", f"{self.config.prefix}_UeccDNA.bed"),
            (f"{self.config.prefix}_MeccDNA.fa", f"{self.config.prefix}_MeccDNA.fa"),
            (f"{self.config.prefix}_MeccSites.bed", f"{self.config.prefix}_MeccDNA.sites.bed"),
            (f"{self.config.prefix}_MeccBestSite.bed", f"{self.config.prefix}_MeccDNA.bestsite.bed"),
            (f"{self.config.prefix}_CeccDNA.fa", f"{self.config.prefix}_CeccDNA.fa"),
            (f"{self.config.prefix}_CeccJunctions.bedpe", f"{self.config.prefix}_CeccDNA.junctions.bedpe"),
            (f"{self.config.prefix}_CeccSegments.bed", f"{self.config.prefix}_CeccDNA.segments.bed"),
            (f"{self.config.prefix}_final_validation.renamed.fasta", f"{self.config.prefix}_InferredUeccDNA.fa"),
            (f"{self.config.prefix}_final_validation.simplified.bed", f"{self.config.prefix}_InferredUeccDNA.bed"),
            (f"{self.config.prefix}_eccDNA_report.html", f"{self.config.prefix}_eccDNA_report.html"),
            ("config.yaml", "config.yaml"),
        ]
        
        # CSV files that need column cleaning
        csv_files = [
            (f"{self.config.prefix}_UeccDNA.core.csv", f"{self.config.prefix}_UeccDNA.core.csv"),
            (f"{self.config.prefix}_MeccSites.core.csv", f"{self.config.prefix}_MeccDNA.core.csv"),
            (f"{self.config.prefix}_CeccSegments.core.csv", f"{self.config.prefix}_CeccDNA.segments.core.csv"),
            (f"{self.config.prefix}_final_validation.clean.csv", f"{self.config.prefix}_InferredUeccDNA.csv"),
        ]
        
        # Add XeccDNA if it exists
        if self.config.enable_xecc:
            xecc_candidates = [
                f"{self.config.prefix}_XeccDNA.fa",
                f"{self.config.prefix}_XeccDNA.fasta",
            ]
            for filename in xecc_candidates:
                if (self.temp_dir / filename).exists():
                    final_files.append((filename, f"{self.config.prefix}_XeccDNA.fa"))
                    break
        
        # Copy non-CSV files directly
        copied_count = 0
        for src_name, dst_name in final_files:
            src = self.temp_dir / src_name
            dst = self.final_output_dir / dst_name
            if src.exists():
                shutil.copy2(src, dst)
                self.logger.debug(f"Finalized: {src_name} → {dst_name}")
                copied_count += 1
        
        # Process CSV files - remove merged_from_ids and num_merged columns
        columns_to_remove = ['merged_from_ids', 'num_merged']
        for src_name, dst_name in csv_files:
            src = self.temp_dir / src_name
            dst = self.final_output_dir / dst_name
            if src.exists():
                try:
                    df = pd.read_csv(src)
                    # Remove unwanted columns if they exist
                    cols_to_drop = [col for col in columns_to_remove if col in df.columns]
                    if cols_to_drop:
                        df = df.drop(columns=cols_to_drop)
                        self.logger.debug(f"Removed columns {cols_to_drop} from {src_name}")
                    df.to_csv(dst, index=False)
                    self.logger.debug(f"Finalized CSV: {src_name} → {dst_name}")
                    copied_count += 1
                except Exception as e:
                    # If CSV processing fails, just copy the file as-is
                    self.logger.warning(f"Could not process CSV {src_name}: {e}")
                    shutil.copy2(src, dst)
                    copied_count += 1
        
        self.logger.info(f"Finalized {copied_count} files to output directory")
        
        # Copy log files if they exist
        for log_file in self.temp_dir.glob("*.log"):
            shutil.copy2(log_file, self.final_output_dir / log_file.name)
        
        # Clean up temp directory based on configuration
        if not self.config.keep_tmp:
            self.logger.info(f"Cleaning up temporary directory: {self.temp_dir}")
            try:
                shutil.rmtree(self.temp_dir)
            except Exception as e:
                self.logger.warning(f"Could not completely remove temp directory: {e}")
        else:
            self.logger.info(f"Temporary files retained in: {self.temp_dir}")
    
    def _prepare_bed_file(self, bed_path: Path) -> Path:
        """Prepare BED file for processing, handling compressed files."""
        if not bed_path.exists():
            raise FileNotFoundError(f"BED file not found: {bed_path}")
        
        # If it's compressed, decompress to a temporary file
        if bed_path.suffix == ".gz":
            temp_bed = Path(tempfile.mkstemp(prefix="per-base.", suffix=".bed")[1])
            self._temp_files.append(temp_bed)  # Track for cleanup
            
            self.logger.debug(f"Decompressing {bed_path} to {temp_bed}")
            with gzip.open(bed_path, "rb") as fin:
                with open(temp_bed, "wb") as fout:
                    shutil.copyfileobj(fin, fout)
            
            return temp_bed
        
        return bed_path
    
    def show_steps(self, detailed: bool = False) -> None:
        """Show pipeline steps with enhanced status information."""
        print("CircleSeeker Pipeline Steps:")
        print("=" * 70)
        
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
            
            skip_info = f" (skippable: --{step.skip_condition})" if step.skip_condition else ""
            
            # Basic info
            print(f"{status} Step {i:2d}: {step.name:<15} - {step.description}{skip_info}")
            
            # Detailed info
            if detailed and step.name in self.state.step_metadata:
                meta = self.state.step_metadata[step.name]
                if meta.duration:
                    print(f"    Duration: {meta.duration:.1f}s")
                if meta.error_message:
                    print(f"    Error: {meta.error_message}")
                if meta.output_files:
                    print(f"    Output files: {len(meta.output_files)} files")
        
        print("=" * 70)
        
        # Show runtime summary
        if self.state.pipeline_start_time:
            total_time = self.state.get_total_runtime()
            if total_time:
                print(f"Total runtime: {total_time:.1f}s")
        
        if self.state.failed_step:
            print(f"⚠️  Previous run failed at: {self.state.failed_step}")
            print("   Use --force to restart from beginning, or resume will continue from failure point")
    
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
            
            for i, step in enumerate(self.STEPS[start_idx:end_idx], start_idx + 1):
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
            
            # After successful completion, finalize outputs
            self._finalize_outputs()
            
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
    def _step_tandem_to_ring(self) -> Optional[List[str]]:
        return self._step_carousel()

    def _step_um_classify(self) -> Optional[List[str]]:
        return self._step_gatekeeper()

    def _step_cecc_build(self) -> Optional[List[str]]:
        return self._step_trapeze()

    def _step_umc_process(self) -> Optional[List[str]]:
        return self._step_menagerie()

    def _step_ecc_dedup(self) -> Optional[List[str]]:
        return self._step_harmonizer()

    def _step_read_filter(self) -> Optional[List[str]]:
        return self._step_sieve()

    def _step_split_refine(self) -> Optional[List[str]]:
        return self._step_contortionist()

    def _step_circle_validate(self) -> Optional[List[str]]:
        return self._step_juggler()

    def _step_report_generator(self) -> Optional[List[str]]:
        return self._step_playbill()
    
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
        
        self.state.results["blast_db"] = str(db_prefix)
    
    def _step_tidehunter(self) -> None:
        """Step 1: Run TideHunter."""
        from circleseeker.external.tidehunter import TideHunter
        
        output_file = self.config.output_dir / f"{self.config.prefix}.TH.ecc_candidates.txt"
        
        tidehunter = TideHunter(threads=self.config.threads)
        tidehunter.run_analysis(
            input_file=self.config.input_file,
            output_file=output_file,
            **self.config.tools.tidehunter
        )
        
        self.state.results["tidehunter_output"] = str(output_file)
    
    def _step_carousel(self) -> None:
        """Step 2: Process TideHunter output."""
        from circleseeker.modules.carousel import Carousel
        
        tidehunter_output = self.state.results.get("tidehunter_output")
        if not tidehunter_output:
            raise PipelineError("TideHunter output not found in state")
        
        csv_output = self.config.output_dir / f"{self.config.prefix}_processed.csv"
        fasta_output = self.config.output_dir / f"{self.config.prefix}_circular.fasta"
        
        carousel = Carousel(
            input_file=tidehunter_output,
            output_file=str(csv_output),
            circular_fasta=str(fasta_output),
            logger=self.logger.getChild("tandem_to_ring")
        )
        carousel.process()
        
        self._set_result("tandem_to_ring_csv", str(csv_output))
        self._set_result("tandem_to_ring_fasta", str(fasta_output))
    
    def _step_run_blast(self) -> None:
        """Step 3: Run BLAST alignment."""
        blast_output = self.config.output_dir / f"{self.config.prefix}_blast_results.tsv"
        
        query_file = Path(self._get_result("tandem_to_ring_fasta", "carousel_fasta", default=self.config.input_file))
        db_prefix = Path(self.state.results["blast_db"])
        
        blast_runner = BlastRunner(
            num_threads=self.config.threads,
            **self.config.tools.blast
        )
        blast_runner.run(
            database=db_prefix,
            query_file=query_file,
            output_file=blast_output
        )
        
        self.state.results["blast_output"] = str(blast_output)
    
    def _step_gatekeeper(self) -> None:
        """Step 4: Classify eccDNA types."""
        from circleseeker.modules.gatekeeper import GatekeeperClassifier
        import pandas as pd
        
        blast_output = Path(self.state.results["blast_output"])
        
        uecc_output = self.config.output_dir / f"{self.config.prefix}_uecc.csv"
        mecc_output = self.config.output_dir / f"{self.config.prefix}_mecc.csv"
        unclass_output = self.config.output_dir / f"{self.config.prefix}_unclassified.csv"
        
        classifier = GatekeeperClassifier(
            gap_threshold=self.config.quality.min_coverage,
            min_full_length_coverage=95.0
        )
        
        uecc_df, mecc_df, unclassified_df = classifier.run(blast_output)
        
        # Save results
        if uecc_df is not None and not uecc_df.empty:
            uecc_df.to_csv(uecc_output, index=False)
            self.state.results["uecc_csv"] = str(uecc_output)
            self.state.results["uecc_count"] = len(uecc_df)
        
        if mecc_df is not None and not mecc_df.empty:
            mecc_df.to_csv(mecc_output, index=False)
            self.state.results["mecc_csv"] = str(mecc_output)
            self.state.results["mecc_count"] = len(mecc_df)
        
        if unclassified_df is not None and not unclassified_df.empty:
            unclassified_df.to_csv(unclass_output, index=False)
            self.state.results["unclassified_csv"] = str(unclass_output)
            
            # Prepare data for trapeze step
            cecc_output = self.config.output_dir / f"{self.config.prefix}_cecc.csv"
            cecc_df = self._prepare_trapeze_data(unclassified_df)
            cecc_df.to_csv(cecc_output, index=False)
            self.state.results["cecc_csv"] = str(cecc_output)
            self.state.results["cecc_count"] = len(cecc_df)
    
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
    
    def _step_trapeze(self) -> None:
        """Step 6: Process Cecc candidates."""
        from circleseeker.modules.trapeze import Trapeze
        import pandas as pd
        
        cecc_input = self.config.output_dir / f"{self.config.prefix}_cecc.csv"
        
        if not cecc_input.exists():
            self.logger.info("No CECC candidates found")
            output_file = self.config.output_dir / f"{self.config.prefix}_CeccDNA_processed.csv"
            pd.DataFrame().to_csv(output_file, index=False)
            self._set_result("cecc_build_output", str(output_file))
            self._set_result("cecc_build_count", 0)
            return
        
        output_file = self.config.output_dir / f"{self.config.prefix}_CeccDNA_processed.csv"
        
        trapeze = Trapeze(self.logger.getChild("cecc_build"))
        results_df = trapeze.run_pipeline(
            cecc_input, 
            output_file,
            overlap_threshold=0.8,
            min_segments=2,
            edge_tolerance=10,
            position_tolerance=100
        )
        
        self._set_result("cecc_build_output", str(output_file))
        self._set_result("cecc_build_count", len(results_df) if results_df is not None else 0)
    
    def _step_menagerie(self) -> None:
        """Step 7: Generate FASTA files."""
        from circleseeker.modules.menagerie import Menagerie, MenagerieConfig
        
        circular_fasta = self.config.output_dir / f"{self.config.prefix}_circular.fasta"
        if not circular_fasta.exists():
            self.logger.warning("Circular FASTA file not found, skipping menagerie step")
            return
        
        uecc_files = []
        mecc_files = []
        cecc_files = []
        
        if "uecc_csv" in self.state.results:
            uecc_path = Path(self.state.results["uecc_csv"])
            if uecc_path.exists():
                uecc_files.append(uecc_path)
        
        if "mecc_csv" in self.state.results:
            mecc_path = Path(self.state.results["mecc_csv"])
            if mecc_path.exists():
                mecc_files.append(mecc_path)
        
        cecc_key = self._get_result("cecc_build_output", "trapeze_output")
        if cecc_key:
            cecc_path = Path(cecc_key)
            if cecc_path.exists() and cecc_path.stat().st_size > 0:
                cecc_files.append(cecc_path)
        
        config = MenagerieConfig()
        config.process_xecc = self.config.enable_xecc
        config.validate_sequences = True
        
        menagerie = Menagerie(config, self.logger.getChild("umc_process"))
        
        results = menagerie.run_pipeline(
            fasta_file=circular_fasta,
            uecc_files=uecc_files,
            mecc_files=mecc_files,
            cecc_files=cecc_files,
            output_dir=self.config.output_dir,
            prefix=self.config.prefix,
            process_xecc=self.config.enable_xecc
        )
        
        type_mapping = {'uecc': 'UeccDNA', 'mecc': 'MeccDNA', 'cecc': 'CeccDNA'}
        
        for ecc_type, df in results.items():
            if df is not None and ecc_type in type_mapping:
                ecc_name = type_mapping[ecc_type]
                processed_file = self.config.output_dir / f"{self.config.prefix}_{ecc_name}_processed.csv"
                fasta_file = self.config.output_dir / f"{self.config.prefix}_{ecc_name}_pre.fasta"
                
                if processed_file.exists():
                    self.state.results[f"{ecc_type}_processed"] = str(processed_file)
                if fasta_file.exists():
                    self.state.results[f"{ecc_type}_fasta"] = str(fasta_file)
        
        # Store XeccDNA file if it exists
        if self.config.enable_xecc:
            xecc_fasta = self.config.output_dir / f"{self.config.prefix}_XeccDNA_pre.fasta"
            if xecc_fasta.exists():
                self.state.results["xecc_fasta"] = str(xecc_fasta)
                self.logger.debug(f"Found XeccDNA file: {xecc_fasta}")
    
    def _step_cd_hit(self) -> None:
        """Step 8: Remove redundant sequences."""
        from circleseeker.external.cd_hit import CDHitEst, CDHitConfig
        
        config = CDHitConfig(
            similarity_threshold=0.99,
            threads=self.config.threads,
            memory_limit=8000
        )
        
        cd_hit = CDHitEst(config, self.logger.getChild("cd_hit"))
        
        for ecc_type in ['uecc', 'mecc', 'cecc']:
            fasta_key = f"{ecc_type}_fasta"
            if fasta_key in self.state.results:
                input_fasta = Path(self.state.results[fasta_key])
                if input_fasta.exists() and input_fasta.stat().st_size > 0:
                    output_fasta = self.config.output_dir / f"{self.config.prefix}_{ecc_type.capitalize()}DNA_pre_nr99"
                    
                    cluster_file = cd_hit.cluster_sequences(input_fasta, output_fasta)
                    self.state.results[f"{ecc_type}_nr99"] = str(output_fasta)
                    self.state.results[f"{ecc_type}_clusters"] = str(cluster_file)
    
    def _step_harmonizer(self) -> None:
        """Step 9: Coordinate and deduplicate results with consistent naming."""
        from circleseeker.modules.harmonizer import Harmonizer
        
        harmonizer = Harmonizer(self.logger.getChild("ecc_dedup"))
        
        uecc_input = Path(self.state.results["uecc_processed"]) if "uecc_processed" in self.state.results else None
        uecc_cluster = Path(self.state.results["uecc_clusters"]) if "uecc_clusters" in self.state.results else None
        mecc_input = Path(self.state.results["mecc_processed"]) if "mecc_processed" in self.state.results else None
        mecc_cluster = Path(self.state.results["mecc_clusters"]) if "mecc_clusters" in self.state.results else None
        cecc_input = Path(self.state.results["cecc_processed"]) if "cecc_processed" in self.state.results else None
        cecc_cluster = Path(self.state.results["cecc_clusters"]) if "cecc_clusters" in self.state.results else None
        
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
            (f"{self.config.prefix}_Mecc.fa", f"{self.config.prefix}_MeccDNA.fa"),
            (f"{self.config.prefix}_Cecc.fa", f"{self.config.prefix}_CeccDNA.fa"),
            # UeccDNA.fa should already be correct
        ]
        
        for old_name, new_name in rename_mappings:
            old_path = self.config.output_dir / old_name
            new_path = self.config.output_dir / new_name
            if old_path.exists() and not new_path.exists():
                old_path.rename(new_path)
                self.logger.debug(f"Renamed {old_name} to {new_name}")
        
        # Store harmonizer output paths with consistent naming
        fasta_mappings = {
            'uecc': f"{self.config.prefix}_UeccDNA.fa",
            'mecc': f"{self.config.prefix}_MeccDNA.fa",
            'cecc': f"{self.config.prefix}_CeccDNA.fa"
        }
        
        for ecc_type, filename in fasta_mappings.items():
            fasta_file = self.config.output_dir / filename
            if fasta_file.exists():
                self._set_result(f"ecc_dedup_{ecc_type}_fasta", str(fasta_file))
                self.logger.debug(f"Found harmonizer output: {fasta_file}")
            else:
                self.logger.warning(f"Harmonizer output not found: {fasta_file}")
        
        # Store harmonized CSV files
        for ecc_type in ['Uecc', 'Mecc', 'Cecc']:
            harmonized_csv = self.config.output_dir / f"{self.config.prefix}_{ecc_type}_harmonized.csv"
            if harmonized_csv.exists():
                self.state.results[f"{ecc_type.lower()}_harmonized"] = str(harmonized_csv)
        
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
                file_path = self.config.output_dir / filename
                if file_path.exists():
                    self.state.results[key] = str(file_path)
                    self.logger.debug(f"Found {key}: {file_path}")
    
    def _step_sieve(self) -> None:
        """Step 10: Filter sequences using Carousel classification."""
        from circleseeker.modules.sieve import Sieve
        
        classification_csv = Path(
            self._get_result(
                "tandem_to_ring_csv",
                "carousel_csv",
                default=self.config.output_dir / f"{self.config.prefix}_processed.csv"
            )
        )
        
        if not classification_csv.exists() or classification_csv.stat().st_size == 0:
            self.logger.warning("Carousel classification CSV not found; skipping sieve step")
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
        
        self._set_result("read_filter_output_fasta", str(output_fasta))
        self._set_result("read_filter_total_reads", stats.total_reads)
        self._set_result("read_filter_filtered_reads", stats.filtered_reads)
        self._set_result("read_filter_retained_reads", stats.retained_reads)
        
        self.logger.info(
            f"Sieve complete: {stats.filtered_reads} filtered, "
            f"{stats.retained_reads} retained ({stats.filtered_percentage:.2f}% filtered)"
        )
    
    def _step_minimap2(self) -> None:
        """Step 11: Map filtered sequences to reference."""
        from circleseeker.external.minimap2 import Minimap2, MiniMapConfig
        
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
        
        # Determine input FASTA
        sieve_output = self._get_result("read_filter_output_fasta", "sieve_output_fasta")
        if sieve_output and Path(sieve_output).exists() and Path(sieve_output).stat().st_size > 0:
            input_fasta = Path(sieve_output)
        else:
            # Create combined FASTA as fallback
            input_fasta = self.config.output_dir / f"{self.config.prefix}_all_filtered.fasta"
            with open(input_fasta, 'w') as out_f:
                for ecc_type in ['uecc', 'mecc', 'cecc']:
                    fasta_key = f"{ecc_type}_fasta"
                    if fasta_key in self.state.results:
                        fasta_path = Path(self.state.results[fasta_key])
                        if fasta_path.exists() and fasta_path.stat().st_size > 0:
                            with open(fasta_path, 'r') as in_f:
                                out_f.write(in_f.read())
        
        output_bam = self.config.output_dir / f"{self.config.prefix}_sorted.bam"
        
        # Run alignment
        minimap2.align(
            reference=self.config.reference,
            query=input_fasta,
            output_file=output_bam
        )
        
        self._set_result("minimap2_bam", str(output_bam))
        self._set_result("all_filtered_fasta", str(input_fasta))
    
    def _step_mosdepth(self) -> None:
        """Step 12: Calculate depth coverage."""
        from circleseeker.external.mosdepth import Mosdepth, MosdepthConfig
        
        bam_file = self.state.results.get("minimap2_bam")
        if not bam_file or not Path(bam_file).exists():
            self.logger.warning("No BAM file from minimap2, skipping depth analysis")
            return
        
        config = MosdepthConfig(
            window_size=50,
            threads=self.config.threads,
            keep_intermediate=False
        )
        
        mosdepth = Mosdepth(config, self.logger.getChild("mosdepth"))
        
        output_prefix = str(self.config.output_dir / f"{self.config.prefix}_depth")
        per_base_bed = mosdepth.calculate_depth(
            bam_file=bam_file,
            output_prefix=output_prefix
        )
        
        self.state.results["depth_bed"] = str(per_base_bed)
    
    def _step_contortionist(self) -> None:
        """Step 13: Detect split regions."""
        from circleseeker.modules.contortionist import Contortionist
        
        bam_file = self.state.results.get("minimap2_bam")
        depth_bed = self.state.results.get("depth_bed")
        
        if not bam_file or not depth_bed:
            self.logger.warning("Missing BAM or depth files, skipping contortionist analysis")
            return
        
        # Prepare BED file (handle compression if needed)
        per_base_bed = self._prepare_bed_file(Path(depth_bed))
        
        # Initialize with exactly the same parameters as independent test
        contortionist = Contortionist(
            min_depth=1,
            merge_distance=5,
            min_mapq=0,
            min_softclip=20,
            breakpoint_distance=50
        )
        
        output_file = self.config.output_dir / f"{self.config.prefix}_split_regions.tsv"
        
        # Run analysis with simple direct call
        result_file = contortionist.run_analysis(
            per_base_bed=per_base_bed,
            bam_path=Path(bam_file),
            output_tsv=output_file
        )
        
        self._set_result("split_refine_output", str(result_file if result_file else output_file))
        
        # Count circular candidates
        circular_count = 0
        try:
            import pandas as pd
            output_path = Path(result_file) if result_file else output_file
            if output_path.exists() and output_path.stat().st_size > 0:
                df = pd.read_csv(output_path, sep='\t')
                if not df.empty and 'is_circular_pattern' in df.columns:
                    circular_count = int((df['is_circular_pattern'] == 'Yes').sum())
        except Exception as e:
            self.logger.debug(f"Could not count circular candidates: {e}")
        
        self.state.results["circular_candidates"] = circular_count
    
    def _step_juggler(self) -> None:
        """Step 14: Final validation."""
        from circleseeker.modules.juggler import Juggler, ValidationConfig
        
        contortionist_output = self._get_result("split_refine_output", "contortionist_output")
        combined_fasta = self._get_result(
            "read_filter_output_fasta",
            "sieve_output_fasta",
            default=self.config.output_dir / f"{self.config.prefix}_all_filtered.fasta"
        )
        
        if not contortionist_output or not Path(combined_fasta).exists():
            self.logger.warning("Missing input files for juggler validation")
            return
        
        config = ValidationConfig(
            min_cross_reads=1,
            min_overhang=20,
            mapq_threshold=20,
            consensus_method="pileup",
            threads=4,
            max_workers=1,
            keep_temp=False,
            only_circular=False
        )
        
        juggler = Juggler(config, self.logger.getChild("circle_validate"))
        
        output_prefix = self.config.output_dir / f"{self.config.prefix}_final_validation"
        
        results_df = juggler.run_validation(
            candidates_tsv=Path(contortionist_output),
            reads_fasta=Path(combined_fasta),
            reference_fasta=self.config.reference,
            output_prefix=output_prefix,
            max_regions=None
        )
        
        if results_df is not None and not results_df.empty:
            # Save outputs
            output_csv = f"{output_prefix}.csv"
            results_df.to_csv(output_csv, index=False)
            
            validated_count = int((results_df.get('validated', False) == True).sum())
            self.state.results["final_validation"] = str(output_csv)
            self.state.results["validated_circular_count"] = validated_count
            
            # Store juggler output paths for later steps
            self._set_result("circle_validate_clean", str(f"{output_prefix}.clean.csv"))
            self._set_result("circle_validate_fasta", str(f"{output_prefix}.renamed.fasta"))
            self._set_result("circle_validate_bed", str(f"{output_prefix}.simplified.bed"))
        else:
            self.state.results["validated_circular_count"] = 0
    
    def _step_playbill(self) -> None:
        """Step 15: Generate HTML report with correct file paths."""
        from circleseeker.modules.playbill import PlaybillReportGenerator
        
        # Use harmonizer output files (correct naming)
        uecc_fasta = Path(self._get_result("ecc_dedup_uecc_fasta", "harmonizer_uecc_fasta",
                                           default=self.config.output_dir / f"{self.config.prefix}_UeccDNA.fa"))
        mecc_fasta = Path(self._get_result("ecc_dedup_mecc_fasta", "harmonizer_mecc_fasta",
                                           default=self.config.output_dir / f"{self.config.prefix}_MeccDNA.fa"))
        cecc_fasta = Path(self._get_result("ecc_dedup_cecc_fasta", "harmonizer_cecc_fasta",
                                           default=self.config.output_dir / f"{self.config.prefix}_CeccDNA.fa"))
        
        # Use juggler output for inferred eccDNA
        inferred_fasta = None
        jf = self._get_result("circle_validate_fasta", "juggler_fasta")
        if jf:
            inferred_fasta = Path(jf)
        
        if not inferred_fasta or not inferred_fasta.exists():
            # Try alternative path
            inferred_fasta = self.config.output_dir / f"{self.config.prefix}_final_validation.renamed.fasta"
        
        if not inferred_fasta.exists():
            # Fallback to sieve output
            self.logger.warning("Juggler renamed.fasta not found, using filtered fasta as fallback")
            inferred_fasta = Path(self._get_result("read_filter_output_fasta", "sieve_output_fasta",
                                                  default=self.config.output_dir / f"{self.config.prefix}_all_filtered.fasta"))
        
        # Other required files
        processed_csv = self.config.output_dir / f"{self.config.prefix}_processed.csv"
        output_html = self.config.output_dir / f"{self.config.prefix}_eccDNA_report.html"
        
        # Create empty files for missing ones to avoid errors
        for f in [uecc_fasta, mecc_fasta, cecc_fasta]:
            if not f.exists():
                self.logger.warning(f"Creating empty file: {f}")
                f.touch()
        
        if not inferred_fasta.exists():
            self.logger.warning(f"Creating empty inferred file: {inferred_fasta}")
            inferred_fasta.touch()
        
        playbill = PlaybillReportGenerator(self.logger.getChild("report_generator"))
        
        results = playbill.generate_report(
            uecc_fasta=uecc_fasta,
            mecc_fasta=mecc_fasta,
            cecc_fasta=cecc_fasta,
            inferred_fasta=inferred_fasta,
            reads_fasta=Path(self.config.input_file),
            processed_csv=processed_csv,
            output_html=output_html,
            sample_name=self.config.prefix
        )
        
        self._set_result("report_generator_report", str(output_html))
        self._set_result("report_generator_results", results)
    
    def _step_propmaster(self) -> None:
        """Step 16: Organize output files with correct file type mapping."""
        organized_dir = self.config.output_dir / "organized_results"
        organized_dir.mkdir(parents=True, exist_ok=True)
        
        # Direct file copying with correct mappings
        file_mappings = [
            # UeccDNA files
            (f"{self.config.prefix}_UeccDNA.fa", f"{self.config.prefix}_UeccDNA.fa"),
            (f"{self.config.prefix}_UeccDNA.bed", f"{self.config.prefix}_UeccDNA.bed"),
            (f"{self.config.prefix}_UeccDNA.core.csv", f"{self.config.prefix}_UeccDNA.core.csv"),
            
            # MeccDNA files
            (f"{self.config.prefix}_MeccDNA.fa", f"{self.config.prefix}_MeccDNA.fa"),
            (f"{self.config.prefix}_MeccSites.bed", f"{self.config.prefix}_MeccDNA.sites.bed"),
            (f"{self.config.prefix}_MeccBestSite.bed", f"{self.config.prefix}_MeccDNA.bestsite.bed"),
            (f"{self.config.prefix}_MeccSites.core.csv", f"{self.config.prefix}_MeccDNA.core.csv"),
            
            # CeccDNA files
            (f"{self.config.prefix}_CeccDNA.fa", f"{self.config.prefix}_CeccDNA.fa"),
            (f"{self.config.prefix}_CeccJunctions.bedpe", f"{self.config.prefix}_CeccDNA.junctions.bedpe"),
            (f"{self.config.prefix}_CeccSegments.bed", f"{self.config.prefix}_CeccDNA.segments.bed"),
            (f"{self.config.prefix}_CeccSegments.core.csv", f"{self.config.prefix}_CeccDNA.segments.core.csv"),
            
            # Inferred files
            (f"{self.config.prefix}_final_validation.renamed.fasta", f"{self.config.prefix}_InferredUeccDNA.fa"),
            (f"{self.config.prefix}_final_validation.simplified.bed", f"{self.config.prefix}_InferredUeccDNA.bed"),
            (f"{self.config.prefix}_final_validation.clean.csv", f"{self.config.prefix}_InferredUeccDNA.csv"),
            
            # Report
            (f"{self.config.prefix}_eccDNA_report.html", f"{self.config.prefix}_eccDNA_report.html"),
        ]
        
        # Add XeccDNA if enabled
        if self.config.enable_xecc:
            xecc_candidates = [
                f"{self.config.prefix}_XeccDNA.fa",
                f"{self.config.prefix}_XeccDNA_pre.fasta",
            ]
            for filename in xecc_candidates:
                if (self.config.output_dir / filename).exists():
                    file_mappings.append((filename, f"{self.config.prefix}_XeccDNA.fa"))
                    break
        
        # Copy files with correct mappings
        organized_count = 0
        for src_name, dst_name in file_mappings:
            src_path = self.config.output_dir / src_name
            dst_path = organized_dir / dst_name
            if src_path.exists():
                shutil.copy2(src_path, dst_path)
                self.logger.debug(f"Organized: {src_name} → {dst_name}")
                organized_count += 1
            else:
                self.logger.debug(f"File not found for organization: {src_name}")
        
        # Write summary file
        summary_path = organized_dir / f"{self.config.prefix}_file_summary.txt"
        summary_lines = [f"CircleSeeker Organized Files - Sample: {self.config.prefix}\n"]
        summary_lines.append(f"Total files organized: {organized_count}\n")
        summary_lines.append("\nFiles in organized_results/:\n")
        
        for file_path in sorted(organized_dir.glob("*")):
            if file_path.is_file() and file_path.name != f"{self.config.prefix}_file_summary.txt":
                size = file_path.stat().st_size
                size_str = f"{size:,} bytes" if size < 1024 else f"{size/1024:.1f} KB" if size < 1024*1024 else f"{size/1024/1024:.1f} MB"
                summary_lines.append(f"  {file_path.name} ({size_str})\n")
        
        summary_path.write_text("".join(summary_lines))
        
        self.logger.info(f"Organized {organized_count} files into {organized_dir}")
        self.state.results["propmaster_output"] = str(organized_dir)
        self.state.results["organized_file_count"] = organized_count
