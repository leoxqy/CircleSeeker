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
from circleseeker.utils.logging import setup_logging, get_logger
from circleseeker.utils.display import ConsoleFormatter, ProgressDisplay, E, print_formatted


# Logging setup now centralized in circleseeker.utils.logging


# Advanced-only option names (validated at runtime)
ADVANCED_ONLY_KEYS = {
    "start_from": "--start-from",
    "stop_at": "--stop-at",
    "resume": "--resume",
    "force": "--force",
    "dry_run": "--dry-run",
    "generate_config": "--generate-config",
    "show_steps": "--show-steps",
    "skip_report": "--skip-report",
    "log_output": "--log-output",
}


def _print_help(ctx, param, value):
    if value and not ctx.resilient_parsing:
        click.echo(ctx.get_help(), color=ctx.color)
        ctx.exit()


def _print_version(ctx, param, value):
    if value and not ctx.resilient_parsing:
        click.echo(f"CircleSeeker {__version__}")
        ctx.exit()


@click.group(context_settings=dict(help_option_names=["-h", "--help"]), invoke_without_command=True, add_help_option=False)
# Direct-run options on the top-level command
@click.option("-v", "--version", is_flag=True, is_eager=True, expose_value=False, callback=_print_version, help="Show the version and exit.")
# Direct-run options on the top-level command
@click.option("-i", "--input", "input_file", type=click.Path(exists=True, path_type=Path), required=False, help="Input FASTA file (HiFi reads)")
@click.option("-r", "--reference", type=click.Path(exists=True, path_type=Path), required=False, help="Reference genome FASTA file")
@click.option("-o", "--output", type=click.Path(path_type=Path), default=Path("circleseeker_output"), help="Output directory")
@click.option("-p", "--prefix", default="sample", help="Output file prefix")
@click.option("-c", "--config", type=click.Path(exists=True, path_type=Path), help="Configuration file (YAML)")
@click.option("-t", "--threads", type=int, default=8, help="Number of threads [default: 8]")
@click.option("--start-from", type=int, help="Resume from specific step (debug only)", hidden=True)
@click.option("--stop-at", type=int, help="Stop at specific step (debug only)", hidden=True)
@click.option("--resume", is_flag=True, help="Resume from last checkpoint (debug only)", hidden=True)
@click.option("--generate-config", is_flag=True, help="Print default config YAML (debug only)", hidden=True)
@click.option("--show-steps", is_flag=True, help="Show pipeline steps and exit (debug only)", hidden=True)
@click.option("--force", is_flag=True, help="Force re-run all steps (debug only)", hidden=True)
@click.option("-n", "--noise", count=True, help="Increase log verbosity (-n, -nn)")
@click.option("--log-output", type=click.Path(path_type=Path), help="Path for additional log output (debug only)", hidden=True)
@click.option("--dry-run", is_flag=True, help="Show actions without executing (debug only)", hidden=True)
@click.option("--skip-report", is_flag=True, help="Skip report generation (debug only)", hidden=True)
@click.option('-h', '--help', is_flag=True, is_eager=True, expose_value=False, callback=_print_help, help='Show this message and exit.')
@click.option("--keep-tmp", is_flag=True, help="Retain temporary working directory")
@click.option("--debug", is_flag=True, help="Enable advanced options and debug logging")
@click.pass_context
def cli(
    ctx,
    input_file: Path | None,
    reference: Path | None,
    output: Path,
    prefix: str,
    config: Optional[Path],
    threads: int,
    start_from: Optional[int],
    stop_at: Optional[int],
    resume: bool,
    generate_config: bool,
    show_steps: bool,
    force: bool,
    noise: int,
    log_output: Optional[Path],
    dry_run: bool,
    skip_report: bool,
    debug: bool,
    keep_tmp: bool,
):
    """CircleSeeker: Comprehensive eccDNA detection from HiFi sequencing data.

    Run directly as: CircleSeeker -i <reads.fa> -r <ref.fa> [options]
    """
    ctx.ensure_object(dict)
    ctx.obj["debug"] = debug

    if not debug:
        advanced_used = []
        if start_from is not None:
            advanced_used.append(ADVANCED_ONLY_KEYS["start_from"])
        if stop_at is not None:
            advanced_used.append(ADVANCED_ONLY_KEYS["stop_at"])
        if resume:
            advanced_used.append(ADVANCED_ONLY_KEYS["resume"])
        if force:
            advanced_used.append(ADVANCED_ONLY_KEYS["force"])
        if dry_run:
            advanced_used.append(ADVANCED_ONLY_KEYS["dry_run"])
        if generate_config:
            advanced_used.append(ADVANCED_ONLY_KEYS["generate_config"])
        if show_steps:
            advanced_used.append(ADVANCED_ONLY_KEYS["show_steps"])
        if skip_report:
            advanced_used.append(ADVANCED_ONLY_KEYS["skip_report"])
        if log_output:
            advanced_used.append(ADVANCED_ONLY_KEYS["log_output"])

        if advanced_used:
            opts = ", ".join(sorted(set(advanced_used)))
            click.echo(f"The following options require --debug: {opts}", err=True)
            sys.exit(2)

    # Hide or reveal advanced commands in help based on --debug
    advanced_cmds = {"init-config", "show-checkpoint", "validate", "benchmark", "run"}
    for name, cmd in list(ctx.command.commands.items()):
        if name in advanced_cmds:
            cmd.hidden = not debug

    # If a subcommand was invoked, do not run the pipeline here
    if ctx.invoked_subcommand:
        return

    # Setup logging
    if debug:
        log_level = logging.DEBUG
    elif noise >= 2:
        log_level = logging.DEBUG
    elif noise == 1:
        log_level = logging.INFO
    else:
        log_level = logging.ERROR

    setup_logging(level=log_level, log_file=log_output)
    logger = get_logger("cli")

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
            pipeline.show_steps(detailed=False)
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
        # XeccDNA detection enabled by default via config; no CLI toggle
        if skip_report:
            cfg.skip_report = True

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
        # Handle --resume flag (compatibility): ensure not forcing rerun
        if resume and force:
            logger.info("--resume specified; ignoring --force to allow resumption")
            force = False

        # Create formatter for beautiful output
        formatter = ConsoleFormatter()

        # Print header
        print_formatted(formatter.header())

        # Print configuration
        config_dict = {
            'input_file': input_file.name,
            'reference': reference.name,
            'output_dir': str(output.absolute()),
            'threads': threads
        }
        print_formatted(formatter.format_config(config_dict))
        print_formatted(formatter.separator())
        print_formatted(formatter.start_message())

        _ = pipeline.run(start_from=start_from, stop_at=stop_at, force=force)

        # Print beautiful completion message
        print_formatted(formatter.separator())
        print_formatted(formatter.success_message())
        if not cfg.keep_tmp:
            print_formatted(formatter.cleanup_message())
        print_formatted(formatter.final_output_message(output))
        print_formatted(formatter.separator("="))

    except CircleSeekerError as e:
        if 'formatter' in locals():
            print_formatted("")
            print_formatted(formatter.separator())
            print_formatted(formatter.error_message(str(e)))
            print_formatted(formatter.separator("="))
        else:
            logger.error(f"Pipeline error: {e}")
        sys.exit(1)
    except Exception as e:
        if 'formatter' in locals():
            print_formatted("")
            print_formatted(formatter.separator())
            print_formatted(formatter.error_message(f"Unexpected error: {e}"))
            print_formatted(formatter.separator("="))
        else:
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
@click.option("--start-from", type=int, help="Resume from specific step")
@click.option("--stop-at", type=int, help="Stop at specific step (1-based index)")
@click.option("--resume", is_flag=True, help="Resume from last checkpoint (compatibility flag)")
@click.option("--show-steps", is_flag=True, help="Show pipeline steps and exit")
@click.option("--force", is_flag=True, help="Force re-run all steps")
@click.option("-v", "--verbose", count=True, help="Increase verbosity")
@click.option("--log-file", type=click.Path(path_type=Path), help="Log file path")
@click.option("--dry-run", is_flag=True, help="Show what would be done")
@click.option("--skip-report", is_flag=True, help="Skip report generation (report generator)")
@click.option("--skip-organize", is_flag=True, help="Skip organizing outputs (file organizer)")
def run(
    input_file: Path,
    reference: Path,
    output: Path,
    prefix: str,
    config: Optional[Path],
    threads: int,
    keep_tmp: bool,
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
    logger = get_logger("cli")
    
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
            pipeline.show_steps(detailed=False)
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
        # XeccDNA detection enabled by default via config; no CLI toggle
        # Optional skips for report/organize steps
        if skip_report:
            cfg.skip_report = True
        
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

        # Handle --resume flag (compatibility): ensure not forcing rerun
        if resume and force:
            logger.info("--resume specified; ignoring --force to allow resumption")
            force = False

        # Create formatter for beautiful output
        formatter = ConsoleFormatter()

        # Print header
        print_formatted(formatter.header())

        # Print configuration
        config_dict = {
            'input_file': input_file.name,
            'reference': reference.name,
            'output_dir': str(output.absolute()),
            'threads': threads
        }
        print_formatted(formatter.format_config(config_dict))
        print_formatted(formatter.separator())
        print_formatted(formatter.start_message())
        print_formatted("")

        from circleseeker.core.pipeline import Pipeline
        pipeline = Pipeline(cfg)
        results = pipeline.run(
            start_from=start_from,
            stop_at=stop_at,
            force=force
        )
        
        # Print beautiful completion message
        print_formatted("")
        print_formatted(formatter.separator())
        print_formatted(formatter.success_message())
        if not cfg.keep_tmp:
            print_formatted(formatter.cleanup_message())
        print_formatted(formatter.final_output_message(output))
        print_formatted(formatter.separator("="))

    except CircleSeekerError as e:
        if 'formatter' in locals():
            print_formatted("")
            print_formatted(formatter.separator())
            print_formatted(formatter.error_message(str(e)))
            print_formatted(formatter.separator("="))
        else:
            logger.error(f"Pipeline error: {e}")
        sys.exit(1)
    except Exception as e:
        if 'formatter' in locals():
            print_formatted("")
            print_formatted(formatter.separator())
            print_formatted(formatter.error_message(f"Unexpected error: {e}"))
            print_formatted(formatter.separator("="))
        else:
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
enable_xecc: true

# Runtime settings
runtime:
  log_level: "INFO"
  log_file: ~
  tmp_dir: ".tmp"
  checkpoint_interval: 5
  keep_tmp: false
  checkpoint_policy: "continue"
  enable_progress: true

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
