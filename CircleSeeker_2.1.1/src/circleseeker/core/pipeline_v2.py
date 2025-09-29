"""
Simplified pipeline orchestrator for CircleSeeker v2.
Reduces technical debt and provides a cleaner execution model.
"""

from __future__ import annotations

import json
import logging
import time
from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Type

import yaml

from circleseeker.utils.logging import get_logger
from circleseeker.modules.base import ModuleBase, ModuleResult


@dataclass
class PipelineConfig:
    """Simplified pipeline configuration."""

    # Basic settings
    input_file: str
    output_dir: str
    sample_name: str
    reference_file: Optional[str] = None

    # Execution control
    threads: int = 4
    debug: bool = False
    skip_steps: List[str] = field(default_factory=list)

    # Tool parameters (flat structure)
    tool_params: Dict[str, Dict[str, Any]] = field(default_factory=dict)

    @classmethod
    def from_yaml(cls, yaml_path: Path) -> PipelineConfig:
        """Load configuration from YAML file."""
        with open(yaml_path, 'r') as f:
            data = yaml.safe_load(f)

        # Extract known fields
        config = cls(
            input_file=data['input_file'],
            output_dir=data['output_dir'],
            sample_name=data['sample_name'],
            reference_file=data.get('reference_file'),
            threads=data.get('threads', 4),
            debug=data.get('debug', False),
            skip_steps=data.get('skip_steps', []),
            tool_params=data.get('tool_params', {})
        )

        return config

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return asdict(self)

    def save(self, path: Path) -> None:
        """Save configuration to file."""
        with open(path, 'w') as f:
            yaml.dump(self.to_dict(), f, default_flow_style=False)


@dataclass
class PipelineStep:
    """Definition of a pipeline step."""

    name: str
    module_class: Type[ModuleBase]
    description: str
    inputs: Dict[str, str]  # key: parameter name, value: source (file or previous step)
    outputs: List[str]  # list of output keys
    required: bool = True

    def should_skip(self, skip_list: List[str]) -> bool:
        """Check if this step should be skipped."""
        return self.name in skip_list


@dataclass
class PipelineState:
    """Simplified pipeline state for recovery."""

    completed_steps: List[str] = field(default_factory=list)
    step_results: Dict[str, ModuleResult] = field(default_factory=dict)
    start_time: float = field(default_factory=time.time)

    def mark_completed(self, step_name: str, result: ModuleResult) -> None:
        """Mark a step as completed."""
        if step_name not in self.completed_steps:
            self.completed_steps.append(step_name)
        self.step_results[step_name] = result

    def is_completed(self, step_name: str) -> bool:
        """Check if a step is completed."""
        return step_name in self.completed_steps

    def get_output_file(self, step_name: str, output_key: str) -> Optional[Path]:
        """Get output file from a previous step."""
        if step_name in self.step_results:
            result = self.step_results[step_name]
            return result.output_files.get(output_key)
        return None

    def save(self, path: Path) -> None:
        """Save state to JSON file."""
        state_dict = {
            'completed_steps': self.completed_steps,
            'start_time': self.start_time,
            'step_outputs': {}
        }

        # Save output file paths (not full results to keep it simple)
        for step_name, result in self.step_results.items():
            state_dict['step_outputs'][step_name] = {
                k: str(v) for k, v in result.output_files.items()
            }

        with open(path, 'w') as f:
            json.dump(state_dict, f, indent=2)

    @classmethod
    def load(cls, path: Path) -> PipelineState:
        """Load state from JSON file."""
        with open(path, 'r') as f:
            data = json.load(f)

        state = cls(
            completed_steps=data['completed_steps'],
            start_time=data['start_time']
        )

        # Reconstruct basic results from saved outputs
        for step_name, outputs in data.get('step_outputs', {}).items():
            result = ModuleResult(
                success=True,
                module_name=step_name,
                output_files={k: Path(v) for k, v in outputs.items()}
            )
            state.step_results[step_name] = result

        return state


class SimplePipeline:
    """Simplified pipeline executor."""

    def __init__(self, config: PipelineConfig, logger: Optional[logging.Logger] = None):
        """Initialize pipeline."""
        self.config = config
        self.logger = logger or get_logger("CircleSeekerPipeline")
        self.state = PipelineState()
        self.steps: List[PipelineStep] = []
        self.output_dir = Path(config.output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Setup paths
        self.state_file = self.output_dir / "pipeline_state.json"
        self.config_file = self.output_dir / "pipeline_config.yaml"

        # Save config for reproducibility
        self.config.save(self.config_file)

    def define_steps(self) -> None:
        """Define pipeline steps based on standard workflow."""
        # Import modules as needed
        from circleseeker.modules.tandem_to_ring import TandemToRing
        from circleseeker.modules.read_filter import ReadFilter
        from circleseeker.modules.um_classify import UMeccClassifier
        from circleseeker.modules.cecc_build import CeccBuilder
        from circleseeker.modules.umc_process import UMCProcessor
        from circleseeker.modules.ecc_dedup import EccDeduplicator
        from circleseeker.modules.cyrcular_calling import CyrcularCalling
        from circleseeker.modules.iecc_curator import IeccCurator
        from circleseeker.modules.ecc_unify import EccUnifier
        from circleseeker.modules.ecc_summary import EccSummary
        from circleseeker.modules.ecc_packager import EccPackager

        # Standard CircleSeeker workflow
        self.steps = [
            # Step 1: Process tandem repeats
            PipelineStep(
                name="tandem_to_ring",
                module_class=TandemToRing,
                description="Process tandem repeats to circular sequences",
                inputs={"input_file": "config.input_file"},
                outputs=["csv", "fasta", "ctcr_reads"]
            ),

            # Step 2: Filter reads
            PipelineStep(
                name="read_filter",
                module_class=ReadFilter,
                description="Filter out CtcR variant reads",
                inputs={
                    "fasta_file": "config.input_file",
                    "csv_file": "tandem_to_ring.csv"
                },
                outputs=["filtered_fasta"]
            ),

            # Step 3: BLAST and classify
            PipelineStep(
                name="um_classify",
                module_class=UMeccClassifier,
                description="Classify eccDNA into U/M/unclassified",
                inputs={"blast_file": "external.blast_results"},
                outputs=["uecc_csv", "mecc_csv", "unclassified_csv"]
            ),

            # Step 4: Build complex eccDNA
            PipelineStep(
                name="cecc_build",
                module_class=CeccBuilder,
                description="Build complex eccDNA from unclassified",
                inputs={"input_file": "um_classify.unclassified_csv"},
                outputs=["cecc_csv"]
            ),

            # Step 5: Process UMC types
            PipelineStep(
                name="umc_process",
                module_class=UMCProcessor,
                description="Process and generate type-specific FASTA files",
                inputs={
                    "fasta_file": "tandem_to_ring.fasta",
                    "uecc_csv": "um_classify.uecc_csv",
                    "mecc_csv": "um_classify.mecc_csv",
                    "cecc_csv": "cecc_build.cecc_csv"
                },
                outputs=["uecc_fasta", "mecc_fasta", "cecc_fasta",
                        "uecc_processed", "mecc_processed", "cecc_processed"]
            ),

            # Step 6: Deduplicate (after CD-HIT clustering)
            PipelineStep(
                name="ecc_dedup",
                module_class=EccDeduplicator,
                description="Deduplicate eccDNA using CD-HIT clusters",
                inputs={
                    "uecc_csv": "umc_process.uecc_processed",
                    "mecc_csv": "umc_process.mecc_processed",
                    "cecc_csv": "umc_process.cecc_processed",
                    "uecc_clstr": "external.uecc_clstr",
                    "mecc_clstr": "external.mecc_clstr",
                    "cecc_clstr": "external.cecc_clstr"
                },
                outputs=["confirmed_tables", "output_dirs"]
            ),

            # Step 7: Cyrcular calling (if BAM available)
            PipelineStep(
                name="cyrcular_calling",
                module_class=CyrcularCalling,
                description="Detect circular DNA from alignment",
                inputs={
                    "bam_file": "external.bam_file",
                    "reference": "config.reference_file"
                },
                outputs=["overview_tsv"],
                required=False
            ),

            # Step 8: Curate inferred eccDNA
            PipelineStep(
                name="iecc_curator",
                module_class=IeccCurator,
                description="Curate inferred eccDNA",
                inputs={
                    "overview_tsv": "cyrcular_calling.overview_tsv",
                    "reference": "config.reference_file"
                },
                outputs=["simple_csv", "chimeric_csv", "inferred_dir"],
                required=False
            ),

            # Step 9: Unify results
            PipelineStep(
                name="ecc_unify",
                module_class=EccUnifier,
                description="Merge confirmed and inferred eccDNA",
                inputs={
                    "confirmed_csv": "ecc_dedup.confirmed_tables",
                    "simple_csv": "iecc_curator.simple_csv",
                    "chimeric_csv": "iecc_curator.chimeric_csv"
                },
                outputs=["merged_csv"]
            ),

            # Step 10: Generate summary
            PipelineStep(
                name="ecc_summary",
                module_class=EccSummary,
                description="Generate analysis summary reports",
                inputs={
                    "merged_csv": "ecc_unify.merged_csv",
                    "input_fasta": "config.input_file"
                },
                outputs=["html_report", "text_summary"]
            ),

            # Step 11: Package results
            PipelineStep(
                name="ecc_packager",
                module_class=EccPackager,
                description="Organize and package final results",
                inputs={
                    "uecc_dir": "ecc_dedup.output_dirs",
                    "mecc_dir": "ecc_dedup.output_dirs",
                    "cecc_dir": "ecc_dedup.output_dirs",
                    "inferred_dir": "iecc_curator.inferred_dir",
                    "merged_csv": "ecc_unify.merged_csv",
                    "html_report": "ecc_summary.html_report",
                    "text_summary": "ecc_summary.text_summary"
                },
                outputs=["final_output"]
            )
        ]

    def resolve_input(self, input_spec: str) -> Any:
        """
        Resolve input specification to actual value.

        Examples:
            "config.input_file" -> self.config.input_file
            "tandem_to_ring.csv" -> output from tandem_to_ring step
        """
        if input_spec.startswith("config."):
            # Get from config
            attr_name = input_spec.split(".", 1)[1]
            return getattr(self.config, attr_name)

        elif input_spec.startswith("external."):
            # External file that should exist
            file_key = input_spec.split(".", 1)[1]
            # This would be provided via config or command line
            return self.config.tool_params.get("external", {}).get(file_key)

        elif "." in input_spec:
            # Output from previous step
            step_name, output_key = input_spec.rsplit(".", 1)
            return self.state.get_output_file(step_name, output_key)

        else:
            # Direct value
            return input_spec

    def run_step(self, step: PipelineStep) -> ModuleResult:
        """Run a single pipeline step."""
        self.logger.info(f"Running step: {step.name} - {step.description}")

        # Instantiate module
        module = step.module_class(
            name=step.name,
            logger=self.logger,
            debug=self.config.debug
        )

        # Resolve inputs
        kwargs = {}
        for param_name, input_spec in step.inputs.items():
            value = self.resolve_input(input_spec)
            if value is not None:
                kwargs[param_name] = value

        # Add step-specific parameters from config
        if step.name in self.config.tool_params:
            kwargs.update(self.config.tool_params[step.name])

        # Add common parameters
        kwargs['output_dir'] = self.output_dir / step.name
        kwargs['sample_name'] = self.config.sample_name
        kwargs['threads'] = self.config.threads

        # Run module
        result = module.run(**kwargs)

        # Save state
        self.state.mark_completed(step.name, result)
        self.state.save(self.state_file)

        return result

    def run(self, resume: bool = False) -> bool:
        """
        Run the complete pipeline.

        Args:
            resume: Resume from saved state if available

        Returns:
            True if pipeline completed successfully
        """
        # Load state if resuming
        if resume and self.state_file.exists():
            self.logger.info("Resuming from saved state")
            self.state = PipelineState.load(self.state_file)

        # Define steps
        self.define_steps()

        # Run steps
        success = True
        for step in self.steps:
            # Skip if requested
            if step.should_skip(self.config.skip_steps):
                self.logger.info(f"Skipping step: {step.name}")
                continue

            # Skip if already completed
            if self.state.is_completed(step.name):
                self.logger.info(f"Step already completed: {step.name}")
                continue

            # Check if step is required
            if not step.required:
                # Check if inputs are available
                inputs_available = all(
                    self.resolve_input(input_spec) is not None
                    for input_spec in step.inputs.values()
                )
                if not inputs_available:
                    self.logger.info(f"Skipping optional step: {step.name} (inputs not available)")
                    continue

            # Run step
            try:
                result = self.run_step(step)
                if not result.success:
                    self.logger.error(f"Step {step.name} failed: {result.error_message}")
                    success = False
                    break
            except Exception as e:
                self.logger.error(f"Step {step.name} failed with exception: {e}")
                success = False
                break

        # Final summary
        if success:
            self.logger.info("Pipeline completed successfully")
            self._print_summary()
        else:
            self.logger.error("Pipeline failed")

        return success

    def _print_summary(self) -> None:
        """Print pipeline execution summary."""
        total_time = time.time() - self.state.start_time

        self.logger.info("="*60)
        self.logger.info("Pipeline Execution Summary")
        self.logger.info("="*60)
        self.logger.info(f"Sample: {self.config.sample_name}")
        self.logger.info(f"Output: {self.output_dir}")
        self.logger.info(f"Total time: {total_time:.2f} seconds")
        self.logger.info(f"Completed steps: {len(self.state.completed_steps)}")

        # Print metrics from each step
        for step_name, result in self.state.step_results.items():
            if result.metrics:
                self.logger.info(f"\n{step_name} metrics:")
                for key, value in result.metrics.items():
                    self.logger.info(f"  {key}: {value}")


def main():
    """Main entry point for the pipeline."""
    import argparse

    parser = argparse.ArgumentParser(description="CircleSeeker Pipeline v2")
    parser.add_argument("-c", "--config", required=True, help="Pipeline configuration file")
    parser.add_argument("-r", "--resume", action="store_true", help="Resume from saved state")
    parser.add_argument("-d", "--debug", action="store_true", help="Enable debug logging")

    args = parser.parse_args()

    # Load config
    config = PipelineConfig.from_yaml(Path(args.config))
    if args.debug:
        config.debug = True

    # Create and run pipeline
    pipeline = SimplePipeline(config)
    success = pipeline.run(resume=args.resume)

    return 0 if success else 1


if __name__ == "__main__":
    import sys
    sys.exit(main())