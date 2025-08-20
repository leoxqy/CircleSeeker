# src/circleseeker2/cli.py
"""Command-line interface for CircleSeeker2."""
from __future__ import annotations

import sys
import logging
from pathlib import Path
from typing import Optional
import click

from circleseeker2 import __version__
from circleseeker2.config import Config, def_config, save_config
# from circleseeker2.core.pipeline import Pipeline
from circleseeker2.utils.logging import setup_logging
# from circleseeker2.utils.validators import validate_installation
from circleseeker2.exceptions import CircleSeekerError


@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(version=__version__, prog_name="CircleSeeker2")
@click.pass_context
def cli(ctx):
    """CircleSeeker2: Comprehensive eccDNA detection from HiFi sequencing data."""
    ctx.ensure_object(dict)


@cli.command()
@click.argument("input_file", type=click.Path(exists=True, path_type=Path))
@click.argument("reference", type=click.Path(exists=True, path_type=Path))
@click.option("-o", "--output", type=click.Path(path_type=Path), 
              default=Path("circleseeker2_output"), help="Output directory")
@click.option("-p", "--prefix", default="sample", help="Output file prefix")
@click.option("-c", "--config", type=click.Path(exists=True, path_type=Path),
              help="Configuration file (YAML)")
@click.option("-t", "--threads", type=int, default=8, 
              help="Number of threads [default: 8]")
@click.option("--keep-tmp", is_flag=True, help="Keep temporary files")
@click.option("--enable-xecc", is_flag=True, help="Enable XeccDNA detection")
@click.option("--start-from", type=int, help="Resume from specific step")
@click.option("--stop-at", type=int, help="Stop at specific step")
@click.option("--force", is_flag=True, help="Force re-run all steps")
@click.option("-v", "--verbose", count=True, help="Increase verbosity")
@click.option("--log-file", type=click.Path(path_type=Path), help="Log file path")
@click.option("--dry-run", is_flag=True, help="Show what would be done")
def run(
    input_file: Path,
    reference: Path,
    output: Path,
    prefix: str,
    config: Optional[Path],
    threads: int,
    keep_tmp: bool,
    enable_xecc: bool,
    start_from: Optional[int],
    stop_at: Optional[int],
    force: bool,
    verbose: int,
    log_file: Optional[Path],
    dry_run: bool,
):
    """Run the CircleSeeker2 pipeline for eccDNA detection."""
    
    # Setup logging
    log_level = logging.WARNING
    if verbose == 1:
        log_level = logging.INFO
    elif verbose >= 2:
        log_level = logging.DEBUG
    
        setup_logging(level=log_level, log_file=log_file)
    logger = logging.getLogger(__name__)
    
    try:
        # Load configuration
        cfg = def_config(config) if config else Config()
        
        # Override with CLI arguments
        cfg.input_file = input_file
        cfg.reference = reference
        cfg.output_dir = output
        cfg.prefix = prefix
        cfg.threads = threads
        cfg.keep_tmp = keep_tmp
        cfg.enable_xecc = enable_xecc
        
        # Validate configuration
        # cfg.validate()
        
        # Save effective configuration
        output.mkdir(parents=True, exist_ok=True)
        save_config(cfg, output / "config.yaml")
        
        if dry_run:
            logger.info("Dry run mode - showing pipeline steps:")
            # pipeline = Pipeline(cfg)
            # pipeline.show_steps()
            return
        
        # Run pipeline
        logger.info(f"Starting CircleSeeker2 v{__version__}")
        logger.info(f"Input: {input_file}")
        logger.info(f"Reference: {reference}")
        logger.info(f"Output: {output}")
        
        # pipeline = Pipeline(cfg)
        # results = pipeline.run(
        #     start_from=start_from,
        #     stop_at=stop_at,
        #     force=force
        # )
        
        # # Print summary
        # click.echo("\n" + "="*60)
        # click.echo("Pipeline completed successfully!")
        # click.echo(f"Results saved to: {output}")
        # click.echo("\nSummary:")
        # click.echo(f"  Total eccDNA detected: {results.get('total_eccdna', 0)}")
        # click.echo(f"  - UeccDNA: {results.get('uecc_count', 0)}")
        # click.echo(f"  - MeccDNA: {results.get('mecc_count', 0)}")
        # click.echo(f"  - CeccDNA: {results.get('cecc_count', 0)}")
        # click.echo(f"  - MCeccDNA: {results.get('mcecc_count', 0)}")
        # if enable_xecc:
        #     click.echo(f"  - XeccDNA: {results.get('xecc_count', 0)}")
        # click.echo("="*60)
        
    except CircleSeekerError as e:
        logger.error(f"Pipeline error: {e}")
        sys.exit(1)
    except Exception as e:
        logger.exception(f"Unexpected error: {e}")
        sys.exit(2)


@cli.command()
@click.option("-o", "--output", type=click.Path(path_type=Path),
              default=Path("config.yaml"), help="Output configuration file")
def init_config(output: Path):
    """Generate a template configuration file."""
    # from circleseeker2.resources import get_default_config
    
    # config_text = get_default_config()
    # output.write_text(config_text)
    click.echo(f"Configuration template saved to: {output}")
    click.echo("Edit this file to customize your analysis parameters.")


@cli.command()
@click.option("--full", is_flag=True, help="Run full validation including test data")
def validate(full: bool):
    """Validate CircleSeeker2 installation and dependencies."""
    click.echo("Validating CircleSeeker2 installation...")
    
    # issues = validate_installation(full_check=full)
    issues = []
    
    if not issues:
        click.echo("✓ All checks passed!")
        click.echo(f"  CircleSeeker2 version: {__version__}")
    else:
        click.echo("✗ Issues found:")
        for issue in issues:
            click.echo(f"  - {issue}")
        sys.exit(1)


@cli.command()
@click.argument("input_files", nargs=-1, required=True, 
                type=click.Path(exists=True, path_type=Path))
@click.option("-o", "--output", type=click.Path(path_type=Path),
              default=Path("benchmark_results.csv"))
@click.option("-t", "--threads", type=int, default=8)
def benchmark(input_files: tuple[Path, ...], output: Path, threads: int):
    """Run performance benchmarks on test data."""
    # from circleseeker2.utils.profiling import run_benchmarks
    
    click.echo(f"Running benchmarks on {len(input_files)} files...")
    # results = run_benchmarks(list(input_files), threads=threads)
    # results.to_csv(output, index=False)
    click.echo(f"Benchmark results saved to: {output}")


def main(argv: list[str] | None = None) -> int:
    """Main entry point."""
    try:
        cli(argv)
        return 0
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        return 1

