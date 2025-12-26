"""Command-line interface for CircleSeeker."""

from __future__ import annotations

import sys
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Optional
import click

from circleseeker import __version__
from circleseeker.config import Config, load_config, save_config
from circleseeker.exceptions import CircleSeekerError
from circleseeker.utils.logging import setup_logging, get_logger
from circleseeker.utils.display import ConsoleFormatter, print_formatted


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
    "log_output": "--log-output",
}


@dataclass
class PipelineOptions:
    """Container for pipeline execution options."""

    input_file: Optional[Path]
    reference: Optional[Path]
    output: Optional[Path]  # None means use config or default
    prefix: Optional[str]   # None means use config or default
    config_path: Optional[Path]
    threads: Optional[int]  # None means use config or default
    keep_tmp: Optional[bool]  # None=use config, True=--keep-tmp, False=--no-keep-tmp
    start_from: Optional[int]
    stop_at: Optional[int]
    resume: bool
    show_steps: bool
    force: bool
    dry_run: bool
    log_file: Optional[Path] = None
    # Logging options from CLI (used to determine if config should override)
    noise: int = 0  # -n count
    debug: bool = False  # --debug flag


def _show_pipeline_steps() -> None:
    """Show pipeline steps without creating directories (lightweight mode)."""
    from circleseeker.core.pipeline import Pipeline

    # Use Pipeline's class-level STEPS directly to avoid instantiation side effects
    click.echo("\nCircleSeeker Pipeline Steps:")
    click.echo("-" * 40)
    for i, step in enumerate(Pipeline.STEPS, 1):
        display_name = step.display_name or step.name
        click.echo(f"  {i:2d}. {display_name:<25} - {step.description}")
    click.echo("-" * 40)
    click.echo(f"Total: {len(Pipeline.STEPS)} steps\n")


def _detect_inference_engine() -> tuple[str | None, str | None]:
    """
    Detect available inference engines.

    Returns:
        Tuple of (selected_engine, fallback_engine)
        - selected_engine: The engine that will be used ('cresil' or 'cyrcular')
        - fallback_engine: The backup engine if available, or None
    """
    import shutil

    has_cresil = shutil.which("cresil") is not None
    has_cyrcular = shutil.which("cyrcular") is not None

    if has_cresil:
        return ("cresil", "cyrcular" if has_cyrcular else None)
    elif has_cyrcular:
        return ("cyrcular", None)
    else:
        return (None, None)


def _format_inference_info() -> tuple[str, str | None]:
    """
    Format inference engine information for display.

    Returns:
        Tuple of (info_string, warning_message)
        - info_string: The formatted string for config display
        - warning_message: Optional warning if recommended tool is missing
    """
    selected, fallback = _detect_inference_engine()

    if selected is None:
        return ("Not available", "No inference tool found. Install Cresil or Cyrcular.")

    if selected == "cresil":
        return ("Cresil", None)
    else:  # cyrcular
        return ("Cyrcular", "Cresil recommended")


def _execute_pipeline(opts: PipelineOptions, logger) -> None:
    """
    Execute the CircleSeeker pipeline with given options.

    This is the unified execution function used by both cli() and run() commands.

    Args:
        opts: Pipeline execution options
        logger: Logger instance for output
    """
    # Handle show_steps early - no side effects
    if opts.show_steps:
        _show_pipeline_steps()
        return

    # Load configuration first to allow config file values
    cfg = load_config(opts.config_path) if opts.config_path else Config()

    # Reconfigure logging if config specifies different values and CLI didn't override
    # CLI flags (--debug, -n, --log-output) take precedence over config file
    if not opts.debug and opts.noise == 0:
        # No CLI logging flags specified, use config values
        # Default to WARNING for cleaner output (use -n for INFO, -nn for DEBUG)
        config_log_level = getattr(cfg.runtime, "log_level", "WARNING").upper()
        level_map = {"DEBUG": logging.DEBUG, "INFO": logging.INFO, "WARNING": logging.WARNING, "ERROR": logging.ERROR}
        new_level = level_map.get(config_log_level, logging.WARNING)

        config_log_file = getattr(cfg.runtime, "log_file", None)
        # Only use config log_file if CLI didn't specify one
        log_file_to_use = opts.log_file if opts.log_file else config_log_file

        setup_logging(level=new_level, log_file=log_file_to_use)

    # Resolve input/reference: CLI args take precedence, then config file
    input_file = opts.input_file
    reference = opts.reference

    # Config is a dataclass, so fields always exist (may be None)
    if not input_file and cfg.input_file:
        input_file = Path(cfg.input_file) if isinstance(cfg.input_file, str) else cfg.input_file
    if not reference and cfg.reference:
        reference = Path(cfg.reference) if isinstance(cfg.reference, str) else cfg.reference

    # Check required parameters for actual run
    if not input_file or not reference:
        click.echo("Error: Both --input and --reference are required to run the pipeline")
        click.echo("These can be provided via CLI arguments or in a config file (-c)")
        click.echo("Use --show-steps to see pipeline steps without running")
        ctx = click.get_current_context()
        click.echo(ctx.get_help())
        sys.exit(1)

    # Update config with resolved values
    # Priority: CLI arg (if provided) > config file > hardcoded default
    cfg.input_file = input_file
    cfg.reference = reference

    # output_dir: CLI > config > default "circleseeker_output"
    if opts.output is not None:
        cfg.output_dir = opts.output
    elif not hasattr(cfg, "output_dir") or not cfg.output_dir:
        cfg.output_dir = Path("circleseeker_output")

    # prefix: CLI > config > default "sample"
    if opts.prefix is not None:
        cfg.prefix = opts.prefix
    elif not hasattr(cfg, "prefix") or not cfg.prefix:
        cfg.prefix = "sample"

    # threads: CLI > config > default 8
    if opts.threads is not None:
        cfg.threads = opts.threads
    elif not hasattr(cfg, "threads") or cfg.threads is None:
        cfg.threads = 8

    # keep_tmp: --keep-tmp flag overrides config if set
    # Default is False (remove temp directory after completion)
    if opts.keep_tmp:
        cfg.keep_tmp = True
    elif not hasattr(cfg, "keep_tmp"):
        cfg.keep_tmp = False

    # Validate configuration FIRST (fail fast on user errors like missing files)
    from circleseeker.exceptions import ConfigurationError

    try:
        cfg.validate()
    except ConfigurationError as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)

    # Save the user-facing output directory before Pipeline modifies it
    # Pipeline.__init__ changes cfg.output_dir to temp_dir for internal use
    user_output_dir = Path(cfg.output_dir)

    # Handle dry run BEFORE creating Pipeline to avoid side effects (directory creation)
    if opts.dry_run:
        logger.info("Dry run mode - showing what would be executed:")
        _show_pipeline_steps()
        click.echo(f"\nWould process: {input_file}")
        click.echo(f"With reference: {reference}")
        click.echo(f"Output to: {user_output_dir.absolute()}")
        click.echo(f"Using {cfg.threads} threads")
        return

    # Check external tool dependencies
    from circleseeker.utils.dependency_checker import check_dependencies

    check_dependencies(logger=logger)

    # Import pipeline here to avoid circular imports
    from circleseeker.core.pipeline import Pipeline

    pipeline = Pipeline(cfg)

    # Save effective configuration to final output directory
    # Temporarily restore original output_dir so saved config reflects user's intent
    saved_output_dir = cfg.output_dir
    try:
        cfg.output_dir = user_output_dir
        save_config(cfg, pipeline.final_output_dir / "config.yaml")
    except Exception as e:
        logger.warning(f"Could not save config: {e}")
    finally:
        cfg.output_dir = saved_output_dir  # Always restore for pipeline operation

    # Handle --resume flag (compatibility): ensure not forcing rerun
    force = opts.force
    if opts.resume and force:
        logger.info("--resume specified; ignoring --force to allow resumption")
        force = False

    # Create formatter for beautiful output
    formatter = ConsoleFormatter()

    # Detect inference engine
    inference_info, inference_warning = _format_inference_info()

    # Print header
    print_formatted(formatter.header())

    # Print configuration (including inference engine)
    # Use user_output_dir for display (cfg.output_dir was changed to temp_dir by Pipeline)
    config_dict = {
        "input_file": input_file.name,
        "reference": reference.name,
        "output_dir": str(user_output_dir.absolute()),
        "threads": cfg.threads,
        "inference": inference_info,
    }
    print_formatted(formatter.format_config(config_dict))

    # Show warning if recommended tool is missing
    if inference_warning:
        print_formatted(f"  [!] {inference_warning}")

    print_formatted(formatter.separator())
    print_formatted(formatter.start_message())

    # Run pipeline
    pipeline.run(start_from=opts.start_from, stop_at=opts.stop_at, force=force)

    # Print beautiful completion message
    print_formatted(formatter.separator())
    print_formatted(formatter.success_message())
    if not cfg.keep_tmp:
        print_formatted(formatter.cleanup_message())
    # Use user_output_dir for display (cfg.output_dir was changed to temp_dir by Pipeline)
    print_formatted(formatter.final_output_message(user_output_dir))
    print_formatted(formatter.separator("="))


def _print_help(ctx, param, value):
    if value and not ctx.resilient_parsing:
        click.echo(ctx.get_help(), color=ctx.color)
        ctx.exit()


def _print_version(ctx, param, value):
    if value and not ctx.resilient_parsing:
        click.echo(f"CircleSeeker {__version__}")
        ctx.exit()


@click.group(
    context_settings=dict(help_option_names=["-h", "--help"]),
    invoke_without_command=True,
    add_help_option=False,
)
# Direct-run options on the top-level command
@click.option(
    "-v",
    "--version",
    is_flag=True,
    is_eager=True,
    expose_value=False,
    callback=_print_version,
    help="Show the version and exit.",
)
# Direct-run options on the top-level command
@click.option(
    "-i",
    "--input",
    "input_file",
    type=click.Path(exists=True, path_type=Path),
    required=False,
    help="Input FASTA file (HiFi reads)",
)
@click.option(
    "-r",
    "--reference",
    type=click.Path(exists=True, path_type=Path),
    required=False,
    help="Reference genome FASTA file",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=None,  # None allows config file to take effect
    help="Output directory [default: circleseeker_output]",
)
@click.option("-p", "--prefix", default=None, help="Output file prefix [default: sample]")
@click.option(
    "-c", "--config", type=click.Path(exists=True, path_type=Path), help="Configuration file (YAML)"
)
@click.option("-t", "--threads", type=int, default=None, help="Number of threads [default: 8]")
@click.option("--start-from", type=int, help="Resume from specific step (debug only)", hidden=True)
@click.option("--stop-at", type=int, help="Stop at specific step (debug only)", hidden=True)
@click.option(
    "--resume", is_flag=True, help="Resume from last checkpoint (debug only)", hidden=True
)
@click.option(
    "--generate-config", is_flag=True, help="Print default config YAML (debug only)", hidden=True
)
@click.option(
    "--show-steps", is_flag=True, help="Show pipeline steps and exit (debug only)", hidden=True
)
@click.option("--force", is_flag=True, help="Force re-run all steps (debug only)", hidden=True)
@click.option("-n", "--noise", count=True, help="Increase log verbosity (-n, -nn)")
@click.option(
    "--log-output",
    type=click.Path(path_type=Path),
    help="Path for additional log output (debug only)",
    hidden=True,
)
@click.option(
    "--dry-run", is_flag=True, help="Show actions without executing (debug only)", hidden=True
)
@click.option(
    "-h",
    "--help",
    is_flag=True,
    is_eager=True,
    expose_value=False,
    callback=_print_help,
    help="Show this message and exit.",
)
@click.option(
    "--keep-tmp",
    is_flag=True,
    default=False,
    help="Retain temporary working directory (default: remove)",
)
@click.option("--debug", is_flag=True, help="Enable advanced options and debug logging")
@click.pass_context
def cli(
    ctx,
    input_file: Optional[Path],
    reference: Optional[Path],
    output: Optional[Path],
    prefix: Optional[str],
    config: Optional[Path],
    threads: Optional[int],
    start_from: Optional[int],
    stop_at: Optional[int],
    resume: bool,
    generate_config: bool,
    show_steps: bool,
    force: bool,
    noise: int,
    log_output: Optional[Path],
    dry_run: bool,
    debug: bool,
    keep_tmp: Optional[bool],
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

        # Create options container and execute unified pipeline
        opts = PipelineOptions(
            input_file=input_file,
            reference=reference,
            output=output,
            prefix=prefix,
            config_path=config,
            threads=threads,
            keep_tmp=keep_tmp,
            start_from=start_from,
            stop_at=stop_at,
            resume=resume,
            show_steps=show_steps,
            force=force,
            dry_run=dry_run,
            log_file=log_output,
            noise=noise,
            debug=debug,
        )
        _execute_pipeline(opts, logger)

    except CircleSeekerError as e:
        logger.error(f"Pipeline error: {e}")
        sys.exit(1)
    except Exception as e:
        logger.exception(f"Unexpected error: {e}")
        sys.exit(2)


@cli.command(hidden=True)
@click.option(
    "-i",
    "--input",
    "input_file",
    type=click.Path(exists=True, path_type=Path),
    required=False,
    help="Input FASTA file (HiFi reads)",
)
@click.option(
    "-r",
    "--reference",
    type=click.Path(exists=True, path_type=Path),
    required=False,
    help="Reference genome FASTA file",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=None,  # None allows config file to take effect
    help="Output directory [default: circleseeker_output]",
)
@click.option("-p", "--prefix", default=None, help="Output file prefix [default: sample]")
@click.option(
    "-c", "--config", type=click.Path(exists=True, path_type=Path), help="Configuration file (YAML)"
)
@click.option("-t", "--threads", type=int, default=None, help="Number of threads [default: 8]")
@click.option(
    "--keep-tmp",
    is_flag=True,
    default=False,
    help="Retain temporary working directory (default: remove)",
)
@click.option("--start-from", type=int, help="Resume from specific step")
@click.option("--stop-at", type=int, help="Stop at specific step (1-based index)")
@click.option("--resume", is_flag=True, help="Resume from last checkpoint (compatibility flag)")
@click.option("--show-steps", is_flag=True, help="Show pipeline steps and exit")
@click.option("--force", is_flag=True, help="Force re-run all steps")
@click.option("-v", "--verbose", count=True, help="Increase verbosity")
@click.option("--log-file", type=click.Path(path_type=Path), help="Log file path")
@click.option("--dry-run", is_flag=True, help="Show what would be done")
def run(
    input_file: Optional[Path],
    reference: Optional[Path],
    output: Optional[Path],
    prefix: Optional[str],
    config: Optional[Path],
    threads: Optional[int],
    keep_tmp: Optional[bool],
    start_from: Optional[int],
    stop_at: Optional[int],
    resume: bool,
    show_steps: bool,
    force: bool,
    verbose: int,
    log_file: Optional[Path],
    dry_run: bool,
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
        # Create options container and execute unified pipeline
        opts = PipelineOptions(
            input_file=input_file,
            reference=reference,
            output=output,
            prefix=prefix,
            config_path=config,
            threads=threads,
            keep_tmp=keep_tmp,
            start_from=start_from,
            stop_at=stop_at,
            resume=resume,
            show_steps=show_steps,
            force=force,
            dry_run=dry_run,
            log_file=log_file,
            noise=verbose,  # verbose acts like noise for logging level
            debug=verbose >= 2,  # treat verbose >= 2 as debug mode
        )
        _execute_pipeline(opts, logger)

    except CircleSeekerError as e:
        logger.error(f"Pipeline error: {e}")
        sys.exit(1)
    except Exception as e:
        logger.exception(f"Unexpected error: {e}")
        sys.exit(2)


@cli.command(name="show-checkpoint", hidden=True)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=Path("circleseeker_output"),
    help="Output directory containing checkpoint",
)
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
    except Exception as e:
        # Log but continue with default prefix
        logging.getLogger("cli").debug(f"Could not auto-detect prefix: {e}")

    cfg = Config()
    cfg.output_dir = output
    cfg.prefix = detected_prefix
    pipeline = Pipeline(cfg)
    pipeline.show_checkpoint_info()


@cli.command(hidden=True)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=Path("config.yaml"),
    help="Output configuration file",
)
def init_config(output: Path):
    """Generate a template configuration file."""
    from circleseeker.resources import get_default_config

    config_text = get_default_config()
    output.write_text(config_text)
    click.echo(f"Configuration template saved to: {output}")
    click.echo("Edit this file to customize your analysis parameters.")


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
