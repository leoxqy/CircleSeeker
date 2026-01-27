"""Click application entrypoint for CircleSeeker."""

from __future__ import annotations

import logging
import signal
import sys
from pathlib import Path
from types import FrameType
from typing import Optional

import click

from circleseeker.cli.exit_codes import (
    EXIT_SUCCESS,
    EXIT_ERROR,
    EXIT_USAGE,
    EXIT_SIGINT,
    EXIT_SIGTERM,
)

def _handle_signal(signum: int, frame: Optional[FrameType]) -> None:
    """Handle interrupt signals for graceful shutdown."""
    sig_name = "SIGINT" if signum == signal.SIGINT else "SIGTERM"
    click.echo(f"\n{sig_name} received, initiating graceful shutdown...", err=True)
    # Raise KeyboardInterrupt to propagate through the call stack
    raise KeyboardInterrupt(f"{sig_name} received")

from circleseeker import __version__
from circleseeker.exceptions import CircleSeekerError
from circleseeker.utils.logging import get_logger, setup_logging

from .commands.checkpoint import show_checkpoint
from .commands.config import init_config
from .commands.validate import validate
from .common_options import (
    input_option,
    reference_option,
    output_option,
    prefix_option,
    config_option,
    threads_option,
    keep_tmp_option,
    turbo_option,
)
from .pipeline import ADVANCED_ONLY_KEYS, PipelineOptions, execute_pipeline


def _print_help(ctx: click.Context, _param: click.Parameter, value: bool) -> None:
    if value and not ctx.resilient_parsing:
        click.echo(ctx.get_help(), color=ctx.color)
        ctx.exit()


def _print_version(ctx: click.Context, _param: click.Parameter, value: bool) -> None:
    if value and not ctx.resilient_parsing:
        click.echo(f"CircleSeeker {__version__}")
        ctx.exit()


def _print_advanced_help(ctx: click.Context, _param: click.Parameter, value: bool) -> None:
    """Print help for advanced/debug options."""
    if not value or ctx.resilient_parsing:
        return
    click.echo("CircleSeeker Advanced Options")
    click.echo("=" * 50)
    click.echo()
    click.echo("These options require --debug:")
    click.echo()
    group = ctx.command
    for p in group.params:
        if isinstance(p, click.Option) and p.hidden:
            names = " / ".join(p.opts)
            click.echo(f"  {names:<25} {p.help or ''}")
    click.echo()
    # Dynamically list hidden subcommands
    if isinstance(group, click.Group):
        hidden = [(n, c) for n, c in group.commands.items() if getattr(c, 'hidden', False)]
        if hidden:
            click.echo("Hidden Subcommands (require --debug):")
            for name, cmd in hidden:
                click.echo(f"  {name:<25} {cmd.get_short_help_str()}")
            click.echo()
    click.echo("Example:")
    click.echo("  circleseeker --debug --start-from 5 -i reads.fa -r ref.fa")
    click.echo("  circleseeker --debug -vv --resume -i reads.fa -r ref.fa")
    click.echo()
    ctx.exit()


@click.group(
    context_settings=dict(help_option_names=["-h", "--help"]),
    invoke_without_command=True,
    add_help_option=False,
)
# Version option (use -V to avoid conflict with -v/--verbose)
@click.option(
    "-V",
    "--version",
    is_flag=True,
    is_eager=True,
    expose_value=False,
    callback=_print_version,
    help="Show the version and exit.",
)
# Advanced help option
@click.option(
    "--help-advanced",
    is_flag=True,
    is_eager=True,
    expose_value=False,
    callback=_print_advanced_help,
    help="Show help for advanced/debug options and exit.",
)
# Common pipeline options (using shared definitions)
@input_option
@reference_option
@output_option
@prefix_option
@config_option
@threads_option
@keep_tmp_option
@turbo_option
# Sensitivity preset option
@click.option(
    "--preset",
    type=click.Choice(["relaxed", "balanced", "strict"], case_sensitive=False),
    default=None,
    help="Sensitivity preset: relaxed (high recall), balanced (default), strict (high precision)",
)
# Verbosity option (unified: -v/--verbose)
@click.option(
    "-v",
    "--verbose",
    count=True,
    help="Increase log verbosity (-v for INFO, -vv for DEBUG)",
)
# Advanced options (hidden by default, require --debug)
@click.option("--start-from", type=int, help="Resume from specific step (debug only)", hidden=True)
@click.option("--stop-at", type=int, help="Stop at specific step (debug only)", hidden=True)
@click.option("--resume", is_flag=True, help="Resume from last checkpoint (debug only)", hidden=True)
@click.option(
    "--show-steps", is_flag=True, help="Show pipeline steps and exit",
)
@click.option("--force", is_flag=True, help="Force re-run all steps (debug only)", hidden=True)
@click.option(
    "--log-file",
    type=click.Path(path_type=Path),
    help="Path for log file output (debug only)",
    hidden=True,
)
@click.option("--dry-run", is_flag=True, help="Show actions without executing")
@click.option(
    "-h",
    "--help",
    is_flag=True,
    is_eager=True,
    expose_value=False,
    callback=_print_help,
    help="Show this message and exit.",
)
@click.option("--debug", is_flag=True, help="Unlock advanced options and hidden subcommands")
@click.pass_context
def cli(
    ctx: click.Context,
    input_file: Optional[Path],
    reference: Optional[Path],
    output: Optional[Path],
    prefix: Optional[str],
    config: Optional[Path],
    threads: Optional[int],
    keep_tmp: bool,
    turbo: bool,
    preset: Optional[str],
    verbose: int,
    start_from: Optional[int],
    stop_at: Optional[int],
    resume: bool,
    show_steps: bool,
    force: bool,
    log_file: Optional[Path],
    dry_run: bool,
    debug: bool,
) -> None:
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
        if log_file:
            advanced_used.append(ADVANCED_ONLY_KEYS["log_file"])

        if advanced_used:
            opt_names = ", ".join(sorted(set(advanced_used)))
            click.echo(f"Error: The following options require --debug: {opt_names}", err=True)
            click.echo("Tip: Use --help-advanced to see all advanced options.", err=True)
            sys.exit(EXIT_USAGE)

    # If a subcommand was invoked, do not run the pipeline here
    if ctx.invoked_subcommand:
        return

    # Setup logging based on verbosity
    if verbose >= 2:
        log_level = logging.DEBUG
    elif verbose == 1:
        log_level = logging.INFO
    else:
        log_level = logging.WARNING

    setup_logging(level=log_level, log_file=log_file)
    logger = get_logger("cli")

    try:
        opts = PipelineOptions(
            input_file=input_file,
            reference=reference,
            output=output,
            prefix=prefix,
            config_path=config,
            threads=threads,
            keep_tmp=keep_tmp,
            turbo=turbo,
            start_from=start_from,
            stop_at=stop_at,
            resume=resume,
            show_steps=show_steps,
            force=force,
            dry_run=dry_run,
            log_file=log_file,
            noise=verbose,
            debug=debug,
            preset=preset,
        )
        execute_pipeline(opts, logger)

    except KeyboardInterrupt:
        logger.info("Pipeline interrupted by user")
        sys.exit(EXIT_SIGINT)
    except CircleSeekerError as exc:
        logger.error(f"Pipeline error: {exc}")
        sys.exit(EXIT_ERROR)
    except Exception as exc:
        logger.exception(f"Unexpected error: {exc}")
        sys.exit(EXIT_ERROR)


cli.add_command(show_checkpoint)
cli.add_command(init_config)
cli.add_command(validate)


def main(argv: list[str] | None = None) -> int:
    """Main entry point with signal handling."""
    # Set up signal handlers for graceful shutdown
    signal.signal(signal.SIGINT, _handle_signal)
    signal.signal(signal.SIGTERM, _handle_signal)

    try:
        cli(argv)
        return EXIT_SUCCESS
    except KeyboardInterrupt:
        return EXIT_SIGINT
    except SystemExit as exc:
        # Preserve explicit exit codes from cli()
        return exc.code if isinstance(exc.code, int) else EXIT_ERROR
    except Exception as exc:
        click.echo(f"Error: {exc}", err=True)
        return EXIT_ERROR


if __name__ == "__main__":
    sys.exit(main())
