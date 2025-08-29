#!/usr/bin/env python3
"""
CircleSeeker2.py - Continuous CircleSeeker2 pipeline orchestrator.

Runs multiple steps in a single Python process to minimize IO and avoid
breaking between modules. Initially wires:
- make-blastdb (build reference database)
- run-blast (align query fasta to reference DB)

Designed to be the single entry point for packaging (e.g., bioconda).
"""

from __future__ import annotations

import argparse
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


def setup_logging(verbosity: int) -> None:
    """Configure root logging level once for the whole pipeline."""
    level = logging.WARNING
    if verbosity >= 2:
        level = logging.DEBUG
    elif verbosity == 1:
        level = logging.INFO
    logging.basicConfig(level=level, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')


@dataclass
class PipelineConfig:
    # Inputs
    reference_fasta: Path
    query_fasta: Path

    # Outputs/prefixes
    db_prefix: Path
    blast_output: Path
    tidehunter_output: Path
    carousel_csv: Path
    carousel_fasta: Path
    gate_uecc_csv: Path
    gate_mecc_csv: Path
    gate_unclass_csv: Path

    # Parameters
    threads: int = 8
    word_size: int = 100
    evalue: str = "1e-50"
    perc_identity: float = 99.0
    dbtype: str = "nucl"  # 'nucl' or 'prot'
    parse_seqids: bool = False
    taxid: Optional[int] = None
    max_file_size: str = "1GB"
    gate_gap_threshold: float = 10.0
    gate_min_coverage: float = 95.0

    # Flow control
    skip_make_db: bool = False
    skip_blast: bool = False
    skip_tidehunter: bool = False
    skip_carousel: bool = False
    skip_gatekeeper: bool = False


class CircleSeeker2:
    def __init__(self, logger: Optional[logging.Logger] = None) -> None:
        self.logger = logger or logging.getLogger(self.__class__.__name__)

    def run_pipeline(self, cfg: PipelineConfig) -> int:
        """Run the continuous pipeline. Returns 0 on success, non-zero on failure."""
        try:
            self.logger.info("Starting CircleSeeker2 pipeline")
            # Step 0: Build BLAST DB (if not skipped)
            if not cfg.skip_make_db:
                self._step_make_blastdb(cfg)
            else:
                self.logger.info("Skipping make-blastdb as requested")

            # Step 1: TideHunter (if not skipped)
            if not cfg.skip_tidehunter:
                self._step_tidehunter(cfg)
            else:
                self.logger.info("Skipping TideHunter as requested")

            # Step 2: Carousel (if not skipped)
            if not cfg.skip_carousel:
                self._step_carousel(cfg)
            else:
                self.logger.info("Skipping Carousel as requested")

            # Step 3: Run BLAST (if not skipped)
            if not cfg.skip_blast:
                self._step_run_blast(cfg)
            else:
                self.logger.info("Skipping run-blast as requested")

            # Step 4: Gatekeeper classification
            if not cfg.skip_gatekeeper:
                self._step_gatekeeper(cfg)
            else:
                self.logger.info("Skipping Gatekeeper as requested")

            self.logger.info("Pipeline completed successfully")
            return 0
        except Exception as e:
            self.logger.exception(f"Pipeline failed: {e}")
            return 1

    def _step_make_blastdb(self, cfg: PipelineConfig) -> None:
        """Build BLAST database from reference FASTA, minimizing redundant work."""
        logger = self.logger.getChild("make_blastdb")
        db_base = cfg.db_prefix

        # If core DB files exist, assume already built to avoid IO
        required_exts = ['.nhr', '.nin', '.nsq'] if cfg.dbtype == 'nucl' else ['.phr', '.pin', '.psq']
        if all((db_base.with_suffix(ext)).exists() for ext in required_exts):
            logger.info(f"BLAST DB already present at prefix: {db_base}")
            return

        # Lazy import to support both -m package and script execution
        try:
            from .step0_make_blastdb import BlastDatabaseBuilder  # type: ignore
        except ImportError:
            from step0_make_blastdb import BlastDatabaseBuilder  # type: ignore

        builder = BlastDatabaseBuilder(dbtype=cfg.dbtype, parse_seqids=cfg.parse_seqids, taxid=cfg.taxid)

        # Quick input validation (in-memory checks only)
        ok, size_mb, seq_count = builder._check_input_file(cfg.reference_fasta)
        if not ok:
            raise RuntimeError(f"Reference FASTA not valid: {cfg.reference_fasta}")
        if size_mb > 1000:
            logger.warning(f"Large reference FASTA ({size_mb:.2f} MB), operation may take a while")

        # Build
        logger.info("Building BLAST database...")
        builder.build_database(
            input_file=cfg.reference_fasta,
            output_db=cfg.db_prefix,
            title=str(cfg.db_prefix.name),
            max_file_size=cfg.max_file_size,
        )

        # Verify
        if not builder.verify_database(cfg.db_prefix):
            raise RuntimeError("BLAST database verification failed")

    def _step_run_blast(self, cfg: PipelineConfig) -> None:
        """Run BLAST against the built database."""
        logger = self.logger.getChild("run_blast")
        # Lazy import
        try:
            from .step3_run_blast import BlastRunner  # type: ignore
        except ImportError:
            from step3_run_blast import BlastRunner  # type: ignore

        runner = BlastRunner(
            num_threads=cfg.threads,
            word_size=cfg.word_size,
            evalue=cfg.evalue,
            perc_identity=cfg.perc_identity,
            outfmt="6 std sstrand",
        )

        # Ensure output directory exists (one IO op)
        cfg.blast_output.parent.mkdir(parents=True, exist_ok=True)

        # Choose query: prefer Carousel circular FASTA if it exists
        query_path = cfg.carousel_fasta if cfg.carousel_fasta.exists() else cfg.query_fasta

        # We avoid re-reading files; BLAST needs file paths only
        logger.info("Running BLAST alignment...")
        runner.run(
            database=cfg.db_prefix,
            query_file=query_path,
            output_file=cfg.blast_output,
            use_time_cmd=False,
        )

    def _step_tidehunter(self, cfg: PipelineConfig) -> None:
        """Run TideHunter on input reads to produce tandem repeat candidates."""
        logger = self.logger.getChild("tidehunter")
        # Ensure output directory exists
        cfg.tidehunter_output.parent.mkdir(parents=True, exist_ok=True)
        # Lazy import
        try:
            from .step1_tidehunter import TideHunterRunner  # type: ignore
        except ImportError:
            from step1_tidehunter import TideHunterRunner  # type: ignore

        runner = TideHunterRunner(num_threads=cfg.threads)
        logger.info("Running TideHunter...")
        runner.run(input_fasta=cfg.query_fasta, output_path=cfg.tidehunter_output)

    def _step_carousel(self, cfg: PipelineConfig) -> None:
        """Process TideHunter output and generate circular sequences and classification."""
        logger = self.logger.getChild("carousel")
        # Ensure output directories exist
        cfg.carousel_csv.parent.mkdir(parents=True, exist_ok=True)
        cfg.carousel_fasta.parent.mkdir(parents=True, exist_ok=True)
        # Validate TideHunter output if TideHunter was skipped
        if cfg.skip_tidehunter and not cfg.tidehunter_output.exists():
            raise FileNotFoundError(
                f"TideHunter output not found: {cfg.tidehunter_output}. Provide --tidehunter-out or run TideHunter step."
            )
        logger.info("Running Carousel processing...")
        # Lazy import; provide helpful message if heavy deps missing
        try:
            from .step2_carousel import Carousel  # type: ignore
        except ImportError:
            try:
                from step2_carousel import Carousel  # type: ignore
            except ImportError as e:
                raise ImportError(
                    "Carousel dependencies not available. Please ensure pandas, numpy, biopython, and networkx are installed.\n"
                    "Install via: pip install pandas numpy biopython networkx (intervaltree optional)."
                ) from e

        car = Carousel(
            input_file=str(cfg.tidehunter_output),
            output_file=str(cfg.carousel_csv),
            circular_fasta=str(cfg.carousel_fasta),
            logger=logger,
        )
        car.process()

    def _step_gatekeeper(self, cfg: PipelineConfig) -> None:
        """Classify eccDNA from BLAST results into Uecc/Mecc and extract unclassified."""
        logger = self.logger.getChild("gatekeeper")
        # Lazy import
        try:
            from .step4_gatekeeper import GatekeeperClassifier  # type: ignore
        except ImportError:
            from step4_gatekeeper import GatekeeperClassifier  # type: ignore

        # Ensure outputs dir exists
        for p in [cfg.gate_uecc_csv, cfg.gate_mecc_csv, cfg.gate_unclass_csv]:
            p.parent.mkdir(parents=True, exist_ok=True)

        logger.info("Running Gatekeeper classification...")
        classifier = GatekeeperClassifier(
            gap_threshold=cfg.gate_gap_threshold,
            min_full_length_coverage=cfg.gate_min_coverage,
        )
        uecc_df, mecc_df, unclassified_df = classifier.run(cfg.blast_output)

        # Write outputs if non-empty
        if uecc_df is not None and not uecc_df.empty:
            uecc_df.to_csv(cfg.gate_uecc_csv, index=False)
            logger.info(f"Saved Uecc to {cfg.gate_uecc_csv}")
        else:
            logger.info("No Uecc results to save")

        if mecc_df is not None and not mecc_df.empty:
            mecc_df.to_csv(cfg.gate_mecc_csv, index=False)
            logger.info(f"Saved Mecc to {cfg.gate_mecc_csv}")
        else:
            logger.info("No Mecc results to save")

        if unclassified_df is not None and not unclassified_df.empty:
            unclassified_df.to_csv(cfg.gate_unclass_csv, index=False)
            logger.info(f"Saved Unclassified to {cfg.gate_unclass_csv}")
        else:
            logger.info("No Unclassified results to save")


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="CircleSeeker2 continuous pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # IO
    p.add_argument("--reference", "-r", type=Path, required=True, help="Reference FASTA for DB build")
    p.add_argument("--db-prefix", "-d", type=Path, required=False, help="BLAST DB prefix (path without extension)")
    p.add_argument("--query", "-q", type=Path, required=False, help="Reads FASTA (used by TideHunter; also used as BLAST query if Carousel is skipped)")
    # Note: free '-o' for output dir in simple mode; keep --blast-out long-only
    p.add_argument("--blast-out", type=Path, required=False, help="BLAST output TSV path")
    p.add_argument("--tidehunter-out", type=Path, required=False, help="TideHunter output TSV path (step1)")
    p.add_argument("--carousel-csv", type=Path, required=False, help="Carousel classification CSV (step2)")
    p.add_argument("--carousel-fasta", type=Path, required=False, help="Carousel circular FASTA (step2)")
    p.add_argument("--gate-uecc", type=Path, required=False, help="Gatekeeper Uecc CSV (step4)")
    p.add_argument("--gate-mecc", type=Path, required=False, help="Gatekeeper Mecc CSV (step4)")
    p.add_argument("--gate-unclass", type=Path, required=False, help="Gatekeeper Unclassified CSV (step4)")

    # Simple-mode (v1-like) convenience interface
    p.add_argument("-i", "--input", type=Path, required=False, help="Input reads FASTA (simple mode)")
    p.add_argument("-p", "--prefix", type=str, required=False, help="Output prefix base name (simple mode)")
    p.add_argument("-o", "--output", type=Path, required=False, help="Output directory (simple mode)")

    # Parameters
    p.add_argument("--threads", "-t", type=int, default=8, help="Number of threads")
    p.add_argument("--word-size", "-w", type=int, default=100, help="BLAST word size")
    p.add_argument("--evalue", "-e", type=str, default="1e-50", help="BLAST e-value")
    p.add_argument("--perc-identity", "-P", type=float, default=99.0, help="BLAST minimum percent identity")
    p.add_argument("--dbtype", choices=["nucl", "prot"], default="nucl", help="Database type")
    p.add_argument("--parse-seqids", action="store_true", help="Parse sequence IDs in DB build")
    p.add_argument("--taxid", type=int, default=None, help="Taxonomy ID applied to DB")
    p.add_argument("--max-file-size", type=str, default="1GB", help="DB volume max file size")
    p.add_argument("--gate-gap-threshold", type=float, default=10.0, help="Gatekeeper max gap percentage")
    p.add_argument("--gate-min-coverage", type=float, default=95.0, help="Gatekeeper min coverage for full-length")

    # Flow control
    p.add_argument("--skip-make-db", action="store_true", help="Skip DB build if already present")
    p.add_argument("--skip-blast", action="store_true", help="Skip BLAST step")
    p.add_argument("--skip-tidehunter", action="store_true", help="Skip TideHunter step")
    p.add_argument("--skip-carousel", action="store_true", help="Skip Carousel step")
    p.add_argument("--skip-gatekeeper", action="store_true", help="Skip Gatekeeper step")

    # Logging
    p.add_argument("-v", action="count", default=0, help="Increase verbosity (-v INFO, -vv DEBUG)")
    return p


def main(argv: Optional[list[str]] = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    setup_logging(args.v)
    logger = logging.getLogger("CircleSeeker2")

    # Simple mode detection: if input and prefix are provided, compute all IO paths
    simple_mode = args.input is not None and args.prefix is not None
    if simple_mode:
        out_dir = args.output or Path.cwd()
        out_dir.mkdir(parents=True, exist_ok=True)
        tmp_dir = out_dir / ".tmp"
        tmp_dir.mkdir(parents=True, exist_ok=True)

        base = args.prefix
        # Derive IO paths following v1-style naming
        db_prefix = tmp_dir / f"{base}_blast_db"
        tide_out = tmp_dir / f"{base}.TH.ecc_candidates.txt"
        car_csv = tmp_dir / f"{base}.carousel_processed_results.csv"
        car_fa = tmp_dir / f"{base}.carousel_circular_sequences.fasta"
        blast_out = tmp_dir / f"{base}.blast1_out_alignment_results.txt"
        gate_uecc = out_dir / f"{base}_Uecc.csv"
        gate_mecc = out_dir / f"{base}_Mecc.csv"
        gate_unclass = out_dir / f"{base}_Unclassified.csv"

        # Map to advanced args if they weren't explicitly provided
        args.reference = args.reference  # already required
        args.db_prefix = db_prefix
        args.blast_out = blast_out
        args.tidehunter_out = tide_out
        args.carousel_csv = car_csv
        args.carousel_fasta = car_fa
        args.gate_uecc = gate_uecc
        args.gate_mecc = gate_mecc
        args.gate_unclass = gate_unclass
        args.query = args.input
    else:
        # Advanced mode: Derive default step1/2/4 outputs if not provided
        # Place alongside BLAST output for coherence
        if args.blast_out is None:
            parser.error("--blast-out is required unless using simple mode (-i and -p)")
        out_dir = args.blast_out.parent
        tide_out = args.tidehunter_out or out_dir / "step1_tidehunter.txt"
        car_csv = args.carousel_csv or out_dir / "step2_processed.csv"
        car_fa = args.carousel_fasta or out_dir / "step2_circular.fasta"
        gate_uecc = args.gate_uecc or out_dir / "step4_uecc.csv"
        gate_mecc = args.gate_mecc or out_dir / "step4_mecc.csv"
        gate_unclass = args.gate_unclass or out_dir / "step4_unclassified.csv"

        # Validate required inputs depending on skipped steps
        if args.query is None and not args.skip_tidehunter:
            parser.error("--query is required unless --skip-tidehunter is set")

    cfg = PipelineConfig(
        reference_fasta=args.reference,
        query_fasta=args.query if args.query is not None else car_fa,
        db_prefix=args.db_prefix,
        blast_output=args.blast_out,
        tidehunter_output=tide_out,
        carousel_csv=car_csv,
        carousel_fasta=car_fa,
        gate_uecc_csv=gate_uecc,
        gate_mecc_csv=gate_mecc,
        gate_unclass_csv=gate_unclass,
        threads=args.threads,
        word_size=args.word_size,
        evalue=args.evalue,
        perc_identity=args.perc_identity,
        dbtype=args.dbtype,
        parse_seqids=args.parse_seqids,
        taxid=args.taxid,
        max_file_size=args.max_file_size,
        gate_gap_threshold=args.gate_gap_threshold,
        gate_min_coverage=args.gate_min_coverage,
        skip_make_db=args.skip_make_db,
        skip_blast=args.skip_blast,
        skip_tidehunter=args.skip_tidehunter,
        skip_carousel=args.skip_carousel,
        skip_gatekeeper=args.skip_gatekeeper,
    )

    engine = CircleSeeker2(logger=logger)
    return engine.run_pipeline(cfg)


if __name__ == "__main__":
    import sys
    sys.exit(main())