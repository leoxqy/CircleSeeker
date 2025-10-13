"""
Cyrcular-Calling Pipeline - HiFi Circular DNA Detection

A post-alignment pipeline for detecting circular DNA from HiFi sequencing data
using Cyrcular and Varlociraptor tools.

This module integrates with CircleSeeker by consuming the BAM produced by
the minimap2 step and generating annotated tables of circular candidates.

Refactored to call external wrappers for cyrcular/varlociraptor/bcftools/samtools.
"""

from __future__ import annotations

import subprocess
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List, Dict, Any
import pandas as pd

from circleseeker.exceptions import PipelineError
from circleseeker.utils.logging import get_logger
from circleseeker.external.samtools import Samtools
from circleseeker.external.cyrcular import Cyrcular
from circleseeker.external.varlociraptor import Varlociraptor
from circleseeker.external.bcftools import Bcftools
import gzip


@dataclass
class PipelineConfig:
    """Configuration for Cyrcular-Calling pipeline."""

    # Input/Output
    bam_file: Path
    reference: Path
    output_dir: Path
    sample_name: str = "sample"

    # Processing parameters
    threads: int = 4
    min_read_depth: int = 2
    min_split_reads: int = 2
    max_paths_per_component: int = 30
    max_deletion_length: int = 1000

    # Varlociraptor parameters
    fdr: float = 0.2
    vaf_resolution: float = 0.01
    max_depth: int = 200

    # Pipeline options
    keep_intermediate: bool = False
    memory_limit: str = "4G"

    def __post_init__(self) -> None:
        self.bam_file = Path(self.bam_file)
        self.reference = Path(self.reference)
        self.output_dir = Path(self.output_dir)
        if not self.bam_file.exists():
            raise FileNotFoundError(f"BAM file not found: {self.bam_file}")
        if not self.reference.exists():
            raise FileNotFoundError(f"Reference file not found: {self.reference}")


@dataclass
class CircularDNAResult:
    circle_id: str
    circle_length: int
    num_segments: int = 0
    vaf: Optional[float] = None
    read_support: Optional[int] = None

    @classmethod
    def from_dataframe_row(cls, row: pd.Series) -> "CircularDNAResult":
        return cls(
            circle_id=str(row.get("circle_id", "unknown")),
            circle_length=int(row.get("circle_length", 0)),
            num_segments=int(row.get("num_segments", 0)),
            vaf=float(row.get("vaf")) if row.get("vaf") is not None else None,
            read_support=int(row.get("read_support")) if row.get("read_support") is not None else None,
        )


class DependencyChecker:
    REQUIRED_TOOLS = ["samtools", "bcftools", "cyrcular", "varlociraptor"]

    def __init__(self, logger=None) -> None:
        self.logger = logger or get_logger(self.__class__.__name__)

    def check_all(self) -> None:
        self.logger.info("Checking dependencies...")
        missing = [t for t in self.REQUIRED_TOOLS if not shutil.which(t)]
        for t in missing:
            self.logger.error(f"{t} not found in PATH")
        if missing:
            raise PipelineError(f"Missing required tools: {', '.join(missing)}")
        self.logger.info("All dependencies satisfied")


class CyrcularCallingPipeline:
    def __init__(self, config: PipelineConfig, logger=None) -> None:
        self.config = config
        self.logger = logger or get_logger(self.__class__.__name__)
        self.checker = DependencyChecker(self.logger)
        self.file_paths: Dict[str, Path] = {}
        self._initialize_file_paths()
        # External tool wrappers
        self._samtools = Samtools(logger=self.logger.getChild("samtools"))
        self._cyrcular = Cyrcular(logger=self.logger.getChild("cyrcular"))
        self._bcftools = Bcftools(logger=self.logger.getChild("bcftools"))
        self._varloc = Varlociraptor(logger=self.logger.getChild("varlociraptor"))

    def _initialize_file_paths(self) -> None:
        prefix = self.config.sample_name
        output = self.config.output_dir
        self.file_paths = {
            # Cyrcular outputs
            "candidates": output / f"{prefix}.candidates.bcf",
            "candidates_sorted": output / f"{prefix}.candidates.sorted.bcf",
            "graph": output / f"{prefix}.graph",
            "dot_dir": output / f"{prefix}.dot",
            # Varlociraptor outputs
            "alignprops": output / f"{prefix}.alignprops.json",
            "obs": output / f"{prefix}.obs.bcf",
            "obs_sorted": output / f"{prefix}.obs.sorted.bcf",
            "calls": output / f"{prefix}.calls.bcf",
            "calls_fdr": output / f"{prefix}.calls.fdr.bcf",
            # Annotation outputs
            "annotated_graph": output / f"{prefix}.annotated.graph",
            "overview_table": output / f"{prefix}_overview.tsv",
            "details_dir": output / f"{prefix}_tables" / f"{prefix}_details",
            # Temporary files
            "scenario": output / "scenario.yaml",
            "empty_gene": output / "empty_gene.gff3.gz",
            "empty_regulatory": output / "empty_regulatory.gff3.gz",
        }

    def run(self) -> List[CircularDNAResult]:
        self.logger.info("=" * 60)
        self.logger.info(f"Cyrcular-Calling Pipeline - {self.config.sample_name}")
        self.logger.info("=" * 60)
        self.logger.info(f"BAM file: {self.config.bam_file}")
        self.logger.info(f"Reference: {self.config.reference}")
        self.logger.info(f"Output directory: {self.config.output_dir}")

        self.checker.check_all()
        self.config.output_dir.mkdir(parents=True, exist_ok=True)

        self._prepare_indices()
        self._run_cyrcular()
        self._run_varlociraptor()
        self._annotate_results()
        results = self._parse_results()

        if not self.config.keep_intermediate:
            self._cleanup_intermediate_files()

        self._print_summary(results)
        return results

    def _prepare_indices(self) -> None:
        ref_fai = Path(f"{self.config.reference}.fai")
        if not ref_fai.exists():
            self.logger.info("Building reference index")
            self._samtools.faidx(self.config.reference)
        bam_bai = Path(f"{self.config.bam_file}.bai")
        if not bam_bai.exists():
            self.logger.info("Building BAM index")
            self._samtools.index_bam(self.config.bam_file)

    def _run_cyrcular(self) -> None:
        # Run Cyrcular
        self._cyrcular.graph_breakends(
            bam=self.config.bam_file,
            reference=self.config.reference,
            output_candidates=self.file_paths["candidates"],
            output_graph=self.file_paths["graph"],
            dot_dir=self.file_paths["dot_dir"],
            min_read_depth=self.config.min_read_depth,
            min_split_reads=self.config.min_split_reads,
            max_paths_per_component=self.config.max_paths_per_component,
            max_deletion_length=self.config.max_deletion_length,
            threads=self.config.threads,
        )

        # Sort candidates and index
        self.logger.info("Sorting candidate variants")
        self._bcftools.sort(
            input_vcf_bcf=self.file_paths["candidates"],
            output_bcf=self.file_paths["candidates_sorted"],
            memory_limit=self.config.memory_limit,
        )
        self._bcftools.index(self.file_paths["candidates_sorted"])

    def _run_varlociraptor(self) -> None:
        # Estimate alignment properties - NOTE: Uses --bams (plural)
        self.logger.info("Estimating alignment properties")
        self._varloc.estimate_alignment_properties(
            reference=self.config.reference,
            bam=self.config.bam_file,
            output_json=self.file_paths["alignprops"],
        )

        # Preprocess variants - NOTE: Uses --bam (singular)
        self.logger.info("Preprocessing variants with Varlociraptor")
        self._varloc.preprocess_variants(
            reference=self.config.reference,
            candidates_bcf_sorted=self.file_paths["candidates_sorted"],
            alignprops_json=self.file_paths["alignprops"],
            bam=self.config.bam_file,
            output_obs_bcf=self.file_paths["obs"],
            max_depth=self.config.max_depth,
        )

        # Sort observations
        self.logger.info("Sorting observation file")
        self._bcftools.sort(
            input_vcf_bcf=self.file_paths["obs"],
            output_bcf=self.file_paths["obs_sorted"],
            memory_limit=self.config.memory_limit,
        )
        self._bcftools.index(self.file_paths["obs_sorted"], force=True)

        # Create scenario file
        scenario = (
            f"samples:\n  {self.config.sample_name}:\n    resolution: {self.config.vaf_resolution}\n    universe: \"[0.0,1.0]\"\n"
            f"events:\n  PRESENT: \"{self.config.sample_name}:]0.0,1.0[\"\n"
        )
        self.file_paths["scenario"].write_text(scenario)
        self.logger.debug(f"Created scenario file: {self.file_paths['scenario']}")

        # Call variants
        self.logger.info("Calling variants")
        self._varloc.call_variants_generic(
            obs_bcf_sorted=self.file_paths["obs_sorted"],
            sample_name=self.config.sample_name,
            scenario_yaml=self.file_paths["scenario"],
            output_calls_bcf=self.file_paths["calls"],
        )

        # FDR filtering
        self.logger.info(f"FDR filtering (threshold={self.config.fdr})")
        self._varloc.filter_calls_fdr_local_smart(
            input_calls_bcf=self.file_paths["calls"],
            output_calls_fdr_bcf=self.file_paths["calls_fdr"],
            fdr=self.config.fdr,
            memory_limit=self.config.memory_limit,
        )
        self._bcftools.index(self.file_paths["calls_fdr"])

    def _annotate_results(self) -> None:
        # Create empty annotation files (gzip)
        for key in ("empty_gene", "empty_regulatory"):
            p = self.file_paths[key]
            p.parent.mkdir(parents=True, exist_ok=True)
            with gzip.open(p, "wt", encoding="utf-8") as f:
                f.write("##gff-version 3\n")

        # Annotate graph
        self._cyrcular.graph_annotate(
            reference=self.config.reference,
            gene_annotation_gff_gz=self.file_paths["empty_gene"],
            regulatory_annotation_gff_gz=self.file_paths["empty_regulatory"],
            graph_input=self.file_paths["graph"],
            output_graph=self.file_paths["annotated_graph"],
        )

        # Create output dirs and tables
        self.file_paths["details_dir"].mkdir(parents=True, exist_ok=True)
        self._cyrcular.graph_table(
            annotated_graph=self.file_paths["annotated_graph"],
            calls_bcf=self.file_paths["calls_fdr"],
            reference=self.config.reference,
            circle_table=self.file_paths["overview_table"],
            segment_tables_dir=self.file_paths["details_dir"],
        )

    def _parse_results(self) -> List[CircularDNAResult]:
        results: List[CircularDNAResult] = []
        overview = self.file_paths["overview_table"]
        if not overview.exists():
            self.logger.warning(f"Overview file not found: {overview}")
            return results
        try:
            df = pd.read_csv(overview, sep="\t")
            if df.empty:
                self.logger.info("No circular DNA detected")
                return results
            for _, row in df.iterrows():
                results.append(CircularDNAResult.from_dataframe_row(row))
            self.logger.info(f"Parsed {len(results)} circular DNA results")
        except Exception as e:
            self.logger.error(f"Failed to parse results: {e}")
            raise PipelineError(f"Result parsing failed: {e}") from e
        return results

    def _cleanup_intermediate_files(self) -> None:
        """Clean up intermediate files, keeping only essential results."""
        self.logger.info("Cleaning up intermediate files...")
        
        # Remove intermediate BCF files and temporary files
        intermediate_files = [
            "obs", "obs_sorted", "candidates", "candidates_sorted", 
            "calls", "calls_fdr", "scenario", "empty_gene", "empty_regulatory"
        ]
        self._remove_files(intermediate_files)
        
        # Remove graph files (keep only annotated graph)
        graph_files = ["graph", "dot_dir"]
        self._remove_files(graph_files)
        
        # Remove alignment properties file
        alignprops = self.file_paths.get("alignprops")
        if alignprops and alignprops.exists():
            try:
                alignprops.unlink()
                self.logger.debug(f"Removed: {alignprops}")
            except Exception as e:
                self.logger.warning(f"Failed to remove {alignprops}: {e}")

    def _remove_files(self, file_keys: List[str]) -> None:
        """Remove files specified by their keys in file_paths."""
        for key in file_keys:
            p = self.file_paths.get(key)
            if p and p.exists():
                try:
                    # Try to remove potential index files first (e.g., .bcf.csi, .bcf.tbi)
                    try:
                        idx_csi = p.with_suffix(p.suffix + ".csi")
                        if idx_csi.exists():
                            idx_csi.unlink()
                            self.logger.debug(f"Removed index file: {idx_csi}")
                    except Exception:
                        pass
                    try:
                        idx_tbi = p.with_suffix(p.suffix + ".tbi")
                        if idx_tbi.exists():
                            idx_tbi.unlink()
                            self.logger.debug(f"Removed index file: {idx_tbi}")
                    except Exception:
                        pass

                    if p.is_file():
                        p.unlink()
                        self.logger.debug(f"Removed file: {p}")
                    elif p.is_dir():
                        import shutil
                        shutil.rmtree(p)
                        self.logger.debug(f"Removed directory: {p}")
                except Exception as e:
                    self.logger.warning(f"Failed to remove {p}: {e}")

    def _print_summary(self, results: List[CircularDNAResult]) -> None:
        self.logger.info("=" * 60)
        self.logger.info("ANALYSIS COMPLETE")
        self.logger.info("=" * 60)
        self.logger.info(f"Output file: {self.file_paths['overview_table']}")
        self.logger.info(f"Total circular DNAs detected: {len(results)}")
        self.logger.info("=" * 60)


def _parse_args():
    """Parse CLI arguments for direct script execution"""
    import argparse
    parser = argparse.ArgumentParser(
        description="Cyrcular-Calling - HiFi circular DNA detection"
    )
    parser.add_argument(
        "-b", "--bam",
        required=True,
        help="Input BAM file"
    )
    parser.add_argument(
        "-r", "--reference",
        required=True,
        help="Reference FASTA file"
    )
    parser.add_argument(
        "-o", "--output-dir",
        required=True,
        help="Output directory"
    )
    parser.add_argument(
        "-s", "--sample-name",
        default="sample",
        help="Sample name (default: sample)"
    )
    parser.add_argument(
        "--min-read-depth",
        type=int,
        default=2,
        help="Min read depth (default: 2)"
    )
    parser.add_argument(
        "--min-split-reads",
        type=int,
        default=2,
        help="Min split reads (default: 2)"
    )
    parser.add_argument(
        "--max-paths-per-component",
        type=int,
        default=30,
        help="Max paths per component (default: 30)"
    )
    parser.add_argument(
        "--max-deletion-length",
        type=int,
        default=1000,
        help="Max deletion length (default: 1000)"
    )
    parser.add_argument(
        "--fdr",
        type=float,
        default=0.2,
        help="FDR threshold (default: 0.2)"
    )
    parser.add_argument(
        "--max-depth",
        type=int,
        default=200,
        help="Max sequencing depth (default: 200)"
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Threads (default: 4)"
    )
    parser.add_argument(
        "--memory-limit",
        default="6G",
        help="Memory limit (default: 6G)"
    )
    parser.add_argument(
        "--keep-intermediate",
        action="store_true",
        help="Keep intermediate files"
    )
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
        help="Log level (default: INFO)"
    )
    return parser.parse_args()


def main():
    """Main function for CLI execution"""
    from pathlib import Path
    import logging as _logging
    
    args = _parse_args()
    
    bam_file = Path(args.bam)
    reference = Path(args.reference)
    output_dir = Path(args.output_dir)
    
    # Set root log level
    _logging.getLogger().setLevel(getattr(_logging, args.log_level))
    
    # Create config
    config = PipelineConfig(
        bam_file=bam_file,
        reference=reference,
        output_dir=output_dir,
        sample_name=args.sample_name,
        threads=args.threads,
        min_read_depth=args.min_read_depth,
        min_split_reads=args.min_split_reads,
        max_paths_per_component=args.max_paths_per_component,
        max_deletion_length=args.max_deletion_length,
        fdr=args.fdr,
        vaf_resolution=0.01,
        max_depth=args.max_depth,
        keep_intermediate=args.keep_intermediate,
        memory_limit=args.memory_limit
    )
    
    # Create and run pipeline
    pipeline = CyrcularCallingPipeline(config)
    
    try:
        results = pipeline.run()
        print(f"\nDetection complete! Found {len(results)} circular DNAs")
        print(f"Results: {output_dir / f'{args.sample_name}_overview.tsv'}")
        
    except Exception as e:
        print(f"Detection failed: {e}")
        raise


if __name__ == "__main__":
    main()
