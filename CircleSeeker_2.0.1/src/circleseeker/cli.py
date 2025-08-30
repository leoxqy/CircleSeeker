"""Command-line interface for CircleSeeker."""

from __future__ import annotations

import sys
import logging
from pathlib import Path
from typing import Optional
import click

from circleseeker import __version__
from circleseeker.config import Config, load_config, save_config
from circleseeker.exceptions import CircleSeekerError, ConfigurationError


def setup_logging(level: int, log_file: Optional[Path] = None) -> None:
    """Setup logging configuration."""
    handlers = []
    
    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(level)
    handlers.append(console_handler)
    
    # File handler if specified
    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        handlers.append(file_handler)
    
    # Configure logging
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=handlers
    )


@click.group(context_settings=dict(help_option_names=["-h", "--help"]), invoke_without_command=True)
@click.version_option(version=__version__, prog_name="CircleSeeker")
# Direct-run options on the top-level command
@click.option("-i", "--input", "input_file", type=click.Path(exists=True, path_type=Path), required=False, help="Input FASTA file (HiFi reads)")
@click.option("-r", "--reference", type=click.Path(exists=True, path_type=Path), required=False, help="Reference genome FASTA file")
@click.option("-o", "--output", type=click.Path(path_type=Path), default=Path("circleseeker_output"), help="Output directory")
@click.option("-p", "--prefix", default="sample", help="Output file prefix")
@click.option("-c", "--config", type=click.Path(exists=True, path_type=Path), help="Configuration file (YAML)")
@click.option("-t", "--threads", type=int, default=8, help="Number of threads [default: 8]")
@click.option("--keep-tmp", is_flag=True, help="Keep temporary files")
@click.option("--enable-xecc", "--enable-xeccdna", "enable_xecc", is_flag=True, help="Enable XeccDNA detection (alias: --enable-xeccdna)")
@click.option("--enable-uecc", is_flag=True, default=True, help="Enable UeccDNA detection [default: True]")
@click.option("--enable-mecc", is_flag=True, default=True, help="Enable MeccDNA detection [default: True]")
@click.option("--enable-cecc", is_flag=True, default=True, help="Enable CeccDNA detection [default: True]")
@click.option("--start-from", type=int, help="Resume from specific step")
@click.option("--stop-at", type=int, help="Stop at specific step (1-based index)")
@click.option("--stop-after", type=str, help="Stop after the named step (e.g., 'carousel')")
@click.option("--resume", is_flag=True, help="Resume from last checkpoint (compatibility flag)")
@click.option("--generate-config", is_flag=True, help="Print default config YAML to stdout and exit")
@click.option("--show-steps", is_flag=True, help="Show pipeline steps and exit")
@click.option("--force", is_flag=True, help="Force re-run all steps")
@click.option("-v", "--verbose", count=True, help="Increase verbosity (-v/-vv)")
@click.option("--log-file", type=click.Path(path_type=Path), help="Log file path")
@click.option("--dry-run", is_flag=True, help="Show what would be done")
@click.option("--skip-report", is_flag=True, help="Skip report generation (playbill)")
@click.option("--skip-organize", is_flag=True, help="Skip organizing outputs (propmaster)")
@click.option("--debug", is_flag=True, help="Show advanced commands and enable debug logging")
@click.pass_context
def cli(
    ctx,
    input_file: Path | None,
    reference: Path | None,
    output: Path,
    prefix: str,
    config: Optional[Path],
    threads: int,
    keep_tmp: bool,
    enable_xecc: bool,
    enable_uecc: bool,
    enable_mecc: bool,
    enable_cecc: bool,
    start_from: Optional[int],
    stop_at: Optional[int],
    stop_after: Optional[str],
    resume: bool,
    generate_config: bool,
    show_steps: bool,
    force: bool,
    verbose: int,
    log_file: Optional[Path],
    dry_run: bool,
    skip_report: bool,
    skip_organize: bool,
    debug: bool,
):
    """CircleSeeker: Comprehensive eccDNA detection from HiFi sequencing data.

    Run directly as: CircleSeeker -i <reads.fa> -r <ref.fa> [options]
    """
    ctx.ensure_object(dict)
    ctx.obj["debug"] = debug

    # Hide or reveal advanced commands in help based on --debug
    advanced_cmds = {"init-config", "show-checkpoint", "validate", "benchmark", "run"}
    for name, cmd in list(ctx.command.commands.items()):
        if name in advanced_cmds:
            cmd.hidden = not debug

    # If a subcommand was invoked, do not run the pipeline here
    if ctx.invoked_subcommand:
        return

    # Setup logging
    log_level = logging.DEBUG if debug else logging.WARNING
    if verbose == 1 and not debug:
        log_level = logging.INFO
    elif verbose >= 2:
        log_level = logging.DEBUG

    setup_logging(level=log_level, log_file=log_file)
    logger = logging.getLogger(__name__)

    try:
        # Generate default config and exit
        if generate_config:
            try:
                from circleseeker.resources import get_default_config
                click.echo(get_default_config())
            except Exception as e:
                logger.error(f"Failed to generate default config: {e}")
                sys.exit(1)
            return
        # Show steps early if requested (no config needed)
        if show_steps:
            from circleseeker.core.pipeline import Pipeline
            from circleseeker.config import Config as ConfigClass
            minimal_cfg = ConfigClass()
            minimal_cfg.input_file = Path("/tmp/dummy.fasta")
            minimal_cfg.reference = Path("/tmp/dummy.fa")
            pipeline = Pipeline(minimal_cfg)
            pipeline.show_steps(detailed=True)
            return

        # Check required parameters for actual run
        if not input_file or not reference:
            click.echo("Error: Both --input and --reference are required to run the pipeline")
            click.echo("Use --show-steps to see pipeline steps without running")
            ctx = click.get_current_context()
            click.echo(ctx.get_help())
            sys.exit(1)

        # Load configuration
        cfg = load_config(config) if config else Config()

        # Override with CLI arguments
        cfg.input_file = input_file
        cfg.reference = reference
        cfg.output_dir = output
        cfg.prefix = prefix
        cfg.threads = threads
        cfg.keep_tmp = keep_tmp
        cfg.enable_xecc = enable_xecc
        if skip_report:
            cfg.skip_report = True
        if skip_organize:
            cfg.skip_organize = True
        if not enable_uecc:
            cfg.skip_uecc = True
        if not enable_mecc:
            cfg.skip_mecc = True
        if not enable_cecc:
            cfg.skip_cecc = True

        # Import pipeline here to avoid circular imports
        from circleseeker.core.pipeline import Pipeline
        pipeline = Pipeline(cfg)

        # Save effective configuration
        output.mkdir(parents=True, exist_ok=True)
        try:
            save_config(cfg, output / "config.yaml")
        except Exception as e:
            logger.warning(f"Could not save config: {e}")

        if dry_run:
            logger.info("Dry run mode - showing what would be executed:")
            pipeline.show_steps()
            click.echo(f"\nWould process: {input_file}")
            click.echo(f"With reference: {reference}")
            click.echo(f"Output to: {output}")
            click.echo(f"Using {threads} threads")
            return

        # Resolve stop-after by name (support legacy aliases)
        if stop_after:
            try:
                from circleseeker.core.pipeline import Pipeline as _P
                names = [s.name for s in _P.STEPS]
                # Legacy -> New mapping
                alias = {
                    'carousel': 'tandem_to_ring',
                    'gatekeeper': 'um_classify',
                    'trapeze': 'cecc_build',
                    'menagerie': 'umc_process',
                    'harmonizer': 'ecc_dedup',
                    'sieve': 'read_filter',
                    'contortionist': 'split_refine',
                    'juggler': 'circle_validate',
                    'playbill': 'report_generator',
                }
                step_key = alias.get(stop_after, stop_after)
                if step_key not in names:
                    raise ValueError(f"Unknown step name: {stop_after}")
                stop_at = names.index(step_key) + 1  # 1-based
            except Exception as e:
                logger.error(f"Failed to resolve --stop-after: {e}")
                sys.exit(2)

        # Handle --resume flag (compatibility): ensure not forcing rerun
        if resume and force:
            logger.info("--resume specified; ignoring --force to allow resumption")
            force = False

        # Run pipeline
        logger.info(f"Starting CircleSeeker v{__version__}")
        logger.info(f"Input: {input_file}")
        logger.info(f"Reference: {reference}")
        logger.info(f"Output: {output}")

        _ = pipeline.run(start_from=start_from, stop_at=stop_at, force=force)

        click.echo("\n" + "=" * 60)
        click.echo("Pipeline completed successfully!")
        click.echo(f"Finalized outputs saved in: {output}")
        if cfg.keep_tmp:
            click.echo(f"Temporary workdir retained: {output}/.tmp_work")
        else:
            click.echo("Temporary workdir cleaned up")
        click.echo("=" * 60)

    except CircleSeekerError as e:
        logger.error(f"Pipeline error: {e}")
        sys.exit(1)
    except Exception as e:
        logger.exception(f"Unexpected error: {e}")
        sys.exit(2)


@cli.command(hidden=True)
@click.option("-i", "--input", "input_file", type=click.Path(exists=True, path_type=Path), 
              required=False, help="Input FASTA file (HiFi reads)")
@click.option("-r", "--reference", type=click.Path(exists=True, path_type=Path),
              required=False, help="Reference genome FASTA file")
@click.option("-o", "--output", type=click.Path(path_type=Path), 
              default=Path("circleseeker_output"), help="Output directory")
@click.option("-p", "--prefix", default="sample", help="Output file prefix")
@click.option("-c", "--config", type=click.Path(exists=True, path_type=Path),
              help="Configuration file (YAML)")
@click.option("-t", "--threads", type=int, default=8, 
              help="Number of threads [default: 8]")
@click.option("--keep-tmp", is_flag=True, help="Keep temporary files")
@click.option("--enable-xecc", "--enable-xeccdna", "enable_xecc", is_flag=True, help="Enable XeccDNA detection (alias: --enable-xeccdna)")
@click.option("--enable-uecc", is_flag=True, default=True, help="Enable UeccDNA detection [default: True]")
@click.option("--enable-mecc", is_flag=True, default=True, help="Enable MeccDNA detection [default: True]") 
@click.option("--enable-cecc", is_flag=True, default=True, help="Enable CeccDNA detection [default: True]")
@click.option("--start-from", type=int, help="Resume from specific step")
@click.option("--stop-at", type=int, help="Stop at specific step (1-based index)")
@click.option("--stop-after", type=str, help="Stop after the named step (e.g., 'carousel')")
@click.option("--resume", is_flag=True, help="Resume from last checkpoint (compatibility flag)")
@click.option("--show-steps", is_flag=True, help="Show pipeline steps and exit")
@click.option("--force", is_flag=True, help="Force re-run all steps")
@click.option("-v", "--verbose", count=True, help="Increase verbosity")
@click.option("--log-file", type=click.Path(path_type=Path), help="Log file path")
@click.option("--dry-run", is_flag=True, help="Show what would be done")
@click.option("--skip-report", is_flag=True, help="Skip report generation (playbill)")
@click.option("--skip-organize", is_flag=True, help="Skip organizing outputs (propmaster)")
def run(
    input_file: Path,
    reference: Path,
    output: Path,
    prefix: str,
    config: Optional[Path],
    threads: int,
    keep_tmp: bool,
    enable_xecc: bool,
    enable_uecc: bool,
    enable_mecc: bool,
    enable_cecc: bool,
    start_from: Optional[int],
    stop_at: Optional[int],
    stop_after: Optional[str],
    resume: bool,
    show_steps: bool,
    force: bool,
    verbose: int,
    log_file: Optional[Path],
    dry_run: bool,
    skip_report: bool,
    skip_organize: bool,
):
    """Run the CircleSeeker pipeline for eccDNA detection."""
    
    # Setup logging
    log_level = logging.WARNING
    if verbose == 1:
        log_level = logging.INFO
    elif verbose >= 2:
        log_level = logging.DEBUG
    
    setup_logging(level=log_level, log_file=log_file)
    logger = logging.getLogger(__name__)
    
    try:
        # Show steps early if requested (no config needed)
        if show_steps:
            from circleseeker.core.pipeline import Pipeline
            from circleseeker.config import Config as ConfigClass
            # Create a minimal config just for showing steps
            minimal_cfg = ConfigClass()
            minimal_cfg.input_file = Path("/tmp/dummy.fasta")  # Dummy path
            minimal_cfg.reference = Path("/tmp/dummy.fa")      # Dummy path
            pipeline = Pipeline(minimal_cfg)
            pipeline.show_steps(detailed=True)
            return
            
        # Check required parameters for actual run
        if not input_file or not reference:
            click.echo("Error: Both --input and --reference are required to run the pipeline")
            click.echo("Use --show-steps to see pipeline steps without running")
            ctx = click.get_current_context()
            click.echo(ctx.get_help())
            sys.exit(1)
        
        # Load configuration
        from circleseeker.config import Config
        cfg = load_config(config) if config else Config()
        
        # Override with CLI arguments
        cfg.input_file = input_file
        cfg.reference = reference
        cfg.output_dir = output
        cfg.prefix = prefix
        cfg.threads = threads
        cfg.keep_tmp = keep_tmp
        cfg.enable_xecc = enable_xecc
        # Optional skips for report/organize steps
        if skip_report:
            cfg.skip_report = True
        if skip_organize:
            cfg.skip_organize = True
        
        # Set eccDNA detection flags
        if not enable_uecc:
            cfg.skip_uecc = True
        if not enable_mecc:
            cfg.skip_mecc = True
        if not enable_cecc:
            cfg.skip_cecc = True
        
        # Import pipeline here to avoid circular imports
        from circleseeker.core.pipeline import Pipeline
        pipeline = Pipeline(cfg)
        
        # Save effective configuration
        output.mkdir(parents=True, exist_ok=True)
        try:
            save_config(cfg, output / "config.yaml")
        except Exception as e:
            logger.warning(f"Could not save config: {e}")
        
        if dry_run:
            logger.info("Dry run mode - showing what would be executed:")
            pipeline.show_steps()
            click.echo(f"\nWould process: {input_file}")
            click.echo(f"With reference: {reference}")
            click.echo(f"Output to: {output}")
            click.echo(f"Using {threads} threads")
            return
        
        # Resolve stop-after by name to stop_at index
        if stop_after:
            try:
                from circleseeker.core.pipeline import Pipeline as _P
                names = [s.name for s in _P.STEPS]
                if stop_after not in names:
                    raise ValueError(f"Unknown step name: {stop_after}")
                stop_at = names.index(stop_after) + 1  # 1-based
            except Exception as e:
                logger.error(f"Failed to resolve --stop-after: {e}")
                sys.exit(2)

        # Handle --resume flag (compatibility): ensure not forcing rerun
        if resume and force:
            logger.info("--resume specified; ignoring --force to allow resumption")
            force = False

        # Run pipeline
        logger.info(f"Starting CircleSeeker v{__version__}")
        logger.info(f"Input: {input_file}")
        logger.info(f"Reference: {reference}")
        logger.info(f"Output: {output}")
        
        from circleseeker.core.pipeline import Pipeline
        pipeline = Pipeline(cfg)
        results = pipeline.run(
            start_from=start_from,
            stop_at=stop_at,
            force=force
        )
        
        # Print simple completion message
        click.echo("\n" + "="*60)
        click.echo("Pipeline completed successfully!")
        click.echo(f"Finalized outputs saved in: {output}")
        if cfg.keep_tmp:
            click.echo(f"Temporary workdir retained: {output}/.tmp_work")
        else:
            click.echo("Temporary workdir cleaned up")
        click.echo("="*60)
        
    except CircleSeekerError as e:
        logger.error(f"Pipeline error: {e}")
        sys.exit(1)
    except Exception as e:
        logger.exception(f"Unexpected error: {e}")
        sys.exit(2)


@cli.command(name="show-checkpoint", hidden=True)
@click.option("-o", "--output", type=click.Path(path_type=Path), default=Path("circleseeker_output"), help="Output directory containing checkpoint")
@click.option("-p", "--prefix", default="sample", help="Run prefix (used to locate checkpoint)")
def show_checkpoint(output: Path, prefix: str):
    """Show detailed checkpoint information for a run."""
    from circleseeker.config import Config
    from circleseeker.core.pipeline import Pipeline
    
    # Auto-detect prefix if possible when not provided explicitly
    detected_prefix = prefix
    try:
        if prefix == "sample":
            checkpoint_files = list(output.glob("*.checkpoint"))
            if len(checkpoint_files) == 1:
                detected_prefix = checkpoint_files[0].name.replace(".checkpoint", "")
    except Exception:
        pass
    
    cfg = Config()
    cfg.output_dir = output
    cfg.prefix = detected_prefix
    pipeline = Pipeline(cfg)
    pipeline.show_checkpoint_info()


@cli.command(hidden=True)
@click.option("-o", "--output", type=click.Path(path_type=Path),
              default=Path("config.yaml"), help="Output configuration file")
def init_config(output: Path):
    """Generate a template configuration file."""
    from circleseeker.resources import get_default_config
    
    try:
        config_text = get_default_config()
        output.write_text(config_text)
        click.echo(f"Configuration template saved to: {output}")
        click.echo("Edit this file to customize your analysis parameters.")
    except ImportError:
        # Fallback if resources module not available yet
        default_config = """# CircleSeeker Configuration File

# Input files (can be overridden by CLI arguments)
input_file: ~
reference: ~
output_dir: "circleseeker_output"
prefix: "sample"

# Feature flags
enable_xecc: false

# Runtime settings
runtime:
  log_level: "INFO"
  log_file: ~
  tmp_dir: ".tmp"
  checkpoint_interval: 5
  keep_tmp: false

# Performance settings
performance:
  threads: 8
  max_memory: "16G"
  chunk_size: 10000
  parallel_jobs: 4
  stream_buffer_size: 65536
  enable_profiling: false

# Quality control parameters
quality:
  min_quality_score: 0.99
  min_coverage: 10
  min_eccdna_size: 100
  max_eccdna_size: 1000000
  min_alignment_length: 100
  min_identity: 99.0

# External tool parameters
tools:
  tidehunter:
    k: 16
    w: 1
    p: 100
    P: 2000000
    e: 0.1
    f: 2
  blast:
    word_size: 100
    evalue: "1e-50"
    perc_identity: 99.0
  minimap2:
    preset: "map-hifi"
    additional_args: ""
  samtools: {}
  mosdepth:
    min_mapq: 1
"""
        output.write_text(default_config)
        click.echo(f"Configuration template saved to: {output}")


@cli.command(hidden=True)
@click.option("--full", is_flag=True, help="Run full validation including test data")
def validate(full: bool):
    """Validate CircleSeeker installation and dependencies."""
    click.echo("Validating CircleSeeker installation...")
    
    try:
        from circleseeker.utils.validators import validate_installation
        issues = validate_installation(full_check=full)
        
        if not issues:
            click.echo("✓ All checks passed!")
            click.echo(f"  CircleSeeker version: {__version__}")
        else:
            click.echo("✗ Issues found:")
            for issue in issues:
                click.echo(f"  - {issue}")
            sys.exit(1)
    except ImportError:
        # Basic validation without full validator
        click.echo("✓ Basic installation check passed!")
        click.echo(f"  CircleSeeker version: {__version__}")


@cli.command(hidden=True)
@click.argument("input_files", nargs=-1, required=True, 
                type=click.Path(exists=True, path_type=Path))
@click.option("-o", "--output", type=click.Path(path_type=Path),
              default=Path("benchmark_results.csv"))
@click.option("-t", "--threads", type=int, default=8)
def benchmark(input_files: tuple[Path, ...], output: Path, threads: int):
    """Run performance benchmarks on test data."""
    try:
        from circleseeker.utils.profiling import run_benchmarks
        
        click.echo(f"Running benchmarks on {len(input_files)} files...")
        results = run_benchmarks(list(input_files), threads=threads)
        results.to_csv(output, index=False)
        click.echo(f"Benchmark results saved to: {output}")
    except ImportError:
        click.echo("Benchmarking functionality not available yet.")
        sys.exit(1)


def main(argv: list[str] | None = None) -> int:
    """Main entry point."""
    try:
        cli(argv)
        return 0
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())
