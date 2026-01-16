"""Shared pipeline execution helpers for the CLI."""

from __future__ import annotations

import logging
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import click

from circleseeker.config import Config, load_config, save_config
from circleseeker.utils.display import ConsoleFormatter, print_formatted
from circleseeker.utils.logging import setup_logging


# Advanced-only option names (validated at runtime)
ADVANCED_ONLY_KEYS = {
    "start_from": "--start-from",
    "stop_at": "--stop-at",
    "resume": "--resume",
    "force": "--force",
    "dry_run": "--dry-run",
    "generate_config": "--generate-config",
    "show_steps": "--show-steps",
    "log_file": "--log-file",
}


@dataclass
class PipelineOptions:
    """Container for pipeline execution options."""

    input_file: Optional[Path]
    reference: Optional[Path]
    output: Optional[Path]  # None means use config or default
    prefix: Optional[str]  # None means use config or default
    config_path: Optional[Path]
    threads: Optional[int]  # None means use config or default
    keep_tmp: bool  # CLI flag: True if --keep-tmp provided
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


def show_pipeline_steps() -> None:
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
    if has_cyrcular:
        return ("cyrcular", None)
    return (None, None)


def _format_inference_info() -> tuple[str, str | None]:
    """
    Format inference engine information for display.

    Returns:
        Tuple of (info_string, warning_message)
        - info_string: The formatted string for config display
        - warning_message: Optional warning if recommended tool is missing
    """
    selected, _fallback = _detect_inference_engine()

    if selected is None:
        return ("Not available", "No inference tool found. Install Cresil or Cyrcular.")

    if selected == "cresil":
        return ("Cresil", None)
    return ("Cyrcular", "Cresil recommended")


def execute_pipeline(
    opts: PipelineOptions,
    logger: logging.Logger,
    ctx: Optional[click.Context] = None,
) -> None:
    """
    Execute the CircleSeeker pipeline with given options.

    This is the unified execution function used by both the top-level command and `run`.

    Args:
        opts: Pipeline execution options
        logger: Logger instance for output
        ctx: Optional Click context for displaying help on error
    """
    # Handle show_steps early - no side effects
    if opts.show_steps:
        show_pipeline_steps()
        return

    # Load configuration first to allow config file values
    cfg = load_config(opts.config_path) if opts.config_path else Config()

    # Reconfigure logging if config specifies different values and CLI didn't override
    # CLI flags (--debug, -n, --log-output) take precedence over config file
    if not opts.debug and opts.noise == 0:
        # No CLI logging flags specified, use config values
        # Default to WARNING for cleaner output (use -n for INFO, -nn for DEBUG)
        config_log_level = getattr(cfg.runtime, "log_level", "WARNING").upper()
        level_map = {
            "DEBUG": logging.DEBUG,
            "INFO": logging.INFO,
            "WARNING": logging.WARNING,
            "ERROR": logging.ERROR,
        }
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
        input_file = cfg.input_file
    if not reference and cfg.reference:
        reference = cfg.reference

    # Check required parameters for actual run
    if not input_file or not reference:
        click.echo("Error: Both --input and --reference are required to run the pipeline", err=True)
        click.echo("These can be provided via CLI arguments or in a config file (-c)", err=True)
        click.echo("Use --show-steps to see pipeline steps without running", err=True)
        # Try to show help if context is available
        if ctx is None:
            try:
                ctx = click.get_current_context(silent=True)
            except RuntimeError:
                ctx = None
        if ctx is not None:
            click.echo(ctx.get_help())
        sys.exit(1)

    # Update config with resolved values
    # Priority: CLI arg (if provided) > config file > hardcoded default
    cfg.input_file = input_file
    cfg.reference = reference

    # output_dir: CLI > config > default "circleseeker_output"
    if opts.output is not None:
        cfg.output_dir = opts.output
    elif not cfg.output_dir:
        cfg.output_dir = Path("circleseeker_output")

    # prefix: CLI > config > default "sample"
    if opts.prefix is not None:
        cfg.prefix = opts.prefix
    elif not cfg.prefix:
        cfg.prefix = "sample"

    # threads: CLI > config > default 8
    if opts.threads is not None:
        cfg.threads = opts.threads
    elif cfg.threads is None:
        cfg.threads = 8

    # keep_tmp: CLI --keep-tmp flag > config file > default False
    # When CLI flag is set, it always takes precedence
    if opts.keep_tmp:
        cfg.keep_tmp = True
    # Otherwise, keep config file value (cfg.keep_tmp already has default False)

    # Validate configuration FIRST (fail fast on user errors like missing files)
    from circleseeker.exceptions import ConfigurationError

    try:
        cfg.validate()
    except ConfigurationError as exc:
        click.echo(f"Error: {exc}", err=True)
        sys.exit(1)

    # Save the user-facing output directory before Pipeline modifies it
    # Pipeline.__init__ changes cfg.output_dir to temp_dir for internal use
    user_output_dir = Path(cfg.output_dir)

    # Handle dry run BEFORE creating Pipeline to avoid side effects (directory creation)
    if opts.dry_run:
        logger.info("Dry run mode - showing what would be executed:")
        show_pipeline_steps()
        click.echo(f"\nWould process: {input_file}")
        click.echo(f"With reference: {reference}")
        click.echo(f"Output to: {user_output_dir.absolute()}")
        click.echo(f"Using {cfg.threads} threads")
        return

    # Note: Dependency checking is now Step 1 of the pipeline

    # Import pipeline here to avoid circular imports
    from circleseeker.core.pipeline import Pipeline

    pipeline = Pipeline(cfg)

    # Save effective configuration to final output directory
    # Temporarily restore original output_dir so saved config reflects user's intent
    saved_output_dir = cfg.output_dir
    try:
        cfg.output_dir = user_output_dir
        save_config(cfg, pipeline.final_output_dir / "config.yaml")
    except Exception as exc:
        logger.warning(f"Could not save config: {exc}")
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
