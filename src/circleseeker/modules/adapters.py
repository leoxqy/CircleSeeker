"""
Adapters for existing CircleSeeker modules to work with the new pipeline architecture.
This allows gradual migration while keeping existing modules functional.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from circleseeker.modules.base import ModuleBase, ModuleResult


class CLIModuleAdapter(ModuleBase):
    """
    Adapter to run existing CLI modules through the new pipeline interface.
    This allows us to use existing modules without modification.
    """

    def __init__(self, module_path: str, **kwargs):
        """
        Initialize CLI module adapter.

        Args:
            module_path: Path to the Python module file
        """
        super().__init__(**kwargs)
        self.module_path = Path(module_path)

        if not self.module_path.exists():
            raise FileNotFoundError(f"Module not found: {module_path}")

    def validate_inputs(self, **kwargs) -> bool:
        """Basic validation - check if required files exist."""
        # Check common input files
        for key in ["input_file", "input", "fasta_file", "blast_file"]:
            if key in kwargs and kwargs[key]:
                file_path = Path(kwargs[key])
                if not file_path.exists():
                    raise FileNotFoundError(f"Input file not found: {file_path}")
        return True

    def build_command(self, **kwargs) -> list:
        """Build command line arguments for the module."""
        cmd = [sys.executable, str(self.module_path)]

        # Map common parameters
        param_mapping = {
            "input_file": "-i",
            "output_file": "-o",
            "output_dir": "-o",
            "sample_name": "-p",
            "threads": "--threads",
            "debug": "--debug",
        }

        for key, value in kwargs.items():
            if value is not None:
                if key in param_mapping:
                    cmd.extend([param_mapping[key], str(value)])
                elif key.startswith("--"):
                    cmd.extend([key, str(value)])

        return cmd

    def execute(self, **kwargs) -> ModuleResult:
        """Execute the CLI module and capture results."""
        result = ModuleResult(success=False, module_name=self.name)

        try:
            cmd = self.build_command(**kwargs)
            self.logger.debug(f"Running command: {' '.join(cmd)}")

            # Run the module
            process = subprocess.run(cmd, capture_output=True, text=True, check=False)

            if process.returncode == 0:
                result.success = True
                self.logger.info(f"{self.name} completed successfully")

                # Try to detect output files
                output_dir = kwargs.get("output_dir", Path("."))
                if output_dir:
                    output_dir = Path(output_dir)
                    if output_dir.exists():
                        # List files created in output directory
                        for file in output_dir.iterdir():
                            if file.is_file():
                                result.add_output(file.name, file)

            else:
                result.error_message = process.stderr
                self.logger.error(f"{self.name} failed: {process.stderr}")

        except Exception as e:
            result.error_message = str(e)
            self.logger.error(f"{self.name} exception: {e}")

        return result


class TandemToRingAdapter(CLIModuleAdapter):
    """Specific adapter for tandem_to_ring module."""

    def __init__(self, **kwargs):
        module_path = "src/circleseeker/modules/tandem_to_ring.py"
        super().__init__(module_path=module_path, name="tandem_to_ring", **kwargs)

    def execute(self, **kwargs) -> ModuleResult:
        """Execute tandem_to_ring with specific handling."""
        # Import the module directly for better integration
        try:
            from circleseeker.modules.tandem_to_ring import TandemToRing

            # Create module instance
            module = TandemToRing(
                input_file=kwargs["input_file"],
                output_file=kwargs.get("output_file", "tandem_to_ring.csv"),
                circular_fasta=kwargs.get("circular_fasta", "circular.fasta"),
                logger=self.logger,
            )

            # Run module
            df_result = module.run()

            # Create result
            result = ModuleResult(
                success=True if df_result is not None else False, module_name=self.name
            )

            # Add output files
            output_file = Path(kwargs.get("output_file", "tandem_to_ring.csv"))
            circular_fasta = Path(kwargs.get("circular_fasta", "circular.fasta"))

            if output_file.exists():
                result.add_output("csv", output_file)
            if circular_fasta.exists():
                result.add_output("fasta", circular_fasta)

            # Add metrics
            if df_result is not None:
                result.add_metric("total_reads", len(df_result))
                if "readClass" in df_result.columns:
                    class_counts = df_result["readClass"].value_counts().to_dict()
                    result.add_metric("class_distribution", class_counts)

            return result

        except ImportError:
            # Fallback to CLI mode
            return super().execute(**kwargs)


class UMClassifyAdapter(CLIModuleAdapter):
    """Specific adapter for um_classify module."""

    def __init__(self, **kwargs):
        module_path = "src/circleseeker/modules/um_classify.py"
        super().__init__(module_path=module_path, name="um_classify", **kwargs)

    def execute(self, **kwargs) -> ModuleResult:
        """Execute um_classify with specific handling."""
        try:
            from circleseeker.modules.um_classify import UMeccClassifier

            # Create module instance
            classifier = UMeccClassifier(logger=self.logger)

            # Read BLAST results
            blast_file = Path(kwargs["blast_file"])
            df = classifier.read_blast_results(blast_file)

            # Classify
            uecc_df, mecc_df, unclassified_df = classifier.classify(df)

            # Save outputs
            output_prefix = kwargs.get("output_prefix", "um_classify")

            uecc_file = Path(f"{output_prefix}.uecc.csv")
            mecc_file = Path(f"{output_prefix}.mecc.csv")
            unclass_file = Path(f"{output_prefix}.unclassified.csv")

            uecc_df.to_csv(uecc_file, index=False)
            mecc_df.to_csv(mecc_file, index=False)
            unclassified_df.to_csv(unclass_file, index=False)

            # Create result
            result = ModuleResult(success=True, module_name=self.name)

            result.add_output("uecc_csv", uecc_file)
            result.add_output("mecc_csv", mecc_file)
            result.add_output("unclassified_csv", unclass_file)

            # Add metrics
            result.add_metric("uecc_count", len(uecc_df))
            result.add_metric("mecc_count", len(mecc_df))
            result.add_metric("unclassified_count", len(unclassified_df))

            return result

        except ImportError:
            # Fallback to CLI mode
            return super().execute(**kwargs)


class ReadFilterAdapter(CLIModuleAdapter):
    """Specific adapter for read_filter module."""

    def __init__(self, **kwargs):
        module_path = "src/circleseeker/modules/read_filter.py"
        super().__init__(module_path=module_path, name="read_filter", **kwargs)

    def execute(self, **kwargs) -> ModuleResult:
        """Execute read_filter with specific handling."""
        try:
            from circleseeker.modules.read_filter import ReadFilter

            # Create module instance
            filter_module = ReadFilter(
                fasta_file=kwargs["fasta_file"],
                csv_file=kwargs["csv_file"],
                output_file=kwargs.get("output_file", "filtered.fasta"),
                logger=self.logger,
            )

            # Run filter
            stats = filter_module.filter_reads()

            # Create result
            result = ModuleResult(success=True, module_name=self.name)

            output_file = Path(kwargs.get("output_file", "filtered.fasta"))
            if output_file.exists():
                result.add_output("filtered_fasta", output_file)

            # Add metrics
            if stats:
                result.add_metric("total_reads", stats.total_reads)
                result.add_metric("filtered_reads", stats.filtered_reads)
                result.add_metric("retained_reads", stats.retained_reads)
                result.add_metric("filtered_percentage", stats.filtered_percentage)

            return result

        except ImportError:
            # Fallback to CLI mode
            return super().execute(**kwargs)
