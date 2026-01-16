"""Click application entrypoint for CircleSeeker."""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Optional, cast

import click

from circleseeker import __version__
from circleseeker.exceptions import CircleSeekerError
from circleseeker.utils.logging import get_logger, setup_logging

from .commands.checkpoint import show_checkpoint
from .commands.config import init_config
from .commands.run import run
from .commands.validate import validate
from .common_options import (
    input_option,
    reference_option,
    output_option,
    prefix_option,
    config_option,
    threads_option,
    keep_tmp_option,
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
    if value and not ctx.resilient_parsing:
        click.echo("CircleSeeker Advanced Options")
        click.echo("=" * 50)
        click.echo()
        click.echo("These options require --debug flag to be enabled:")
        click.echo()
        click.echo("Pipeline Control:")
        click.echo("  --start-from INT     Resume from specific step number (1-16)")
        click.echo("  --stop-at INT        Stop at specific step number (1-16)")
        click.echo("  --resume             Resume from last checkpoint")
        click.echo("  --force              Force re-run all steps (ignore checkpoint)")
        click.echo("  --dry-run            Show actions without executing")
        click.echo()
        click.echo("Configuration:")
        click.echo("  --generate-config    Print default config YAML to stdout")
        click.echo("  --show-steps         Show pipeline steps and exit")
        click.echo("  --log-file PATH      Write logs to file")
        click.echo()
        click.echo("Subcommands (require --debug):")
        click.echo("  run                  Run pipeline (explicit subcommand)")
        click.echo("  init-config          Generate config file interactively")
        click.echo("  show-checkpoint      Display checkpoint information")
        click.echo("  validate             Validate installation and dependencies")
        click.echo()
        click.echo("Example:")
        click.echo("  circleseeker --debug --start-from 5 -i reads.fa -r ref.fa -o out/")
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
    "--generate-config", is_flag=True, help="Print default config YAML (debug only)", hidden=True
)
@click.option(
    "--show-steps", is_flag=True, help="Show pipeline steps and exit (debug only)", hidden=True
)
@click.option("--force", is_flag=True, help="Force re-run all steps (debug only)", hidden=True)
@click.option(
    "--log-file",
    type=click.Path(path_type=Path),
    help="Path for log file output (debug only)",
    hidden=True,
)
@click.option("--dry-run", is_flag=True, help="Show actions without executing (debug only)", hidden=True)
@click.option(
    "-h",
    "--help",
    is_flag=True,
    is_eager=True,
    expose_value=False,
    callback=_print_help,
    help="Show this message and exit.",
)
@click.option("--debug", is_flag=True, help="Enable advanced options and debug logging")
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
    verbose: int,
    start_from: Optional[int],
    stop_at: Optional[int],
    resume: bool,
    generate_config: bool,
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
        if dry_run:
            advanced_used.append(ADVANCED_ONLY_KEYS["dry_run"])
        if generate_config:
            advanced_used.append(ADVANCED_ONLY_KEYS["generate_config"])
        if show_steps:
            advanced_used.append(ADVANCED_ONLY_KEYS["show_steps"])
        if log_file:
            advanced_used.append(ADVANCED_ONLY_KEYS["log_file"])

        if advanced_used:
            opt_names = ", ".join(sorted(set(advanced_used)))
            click.echo(f"Error: The following options require --debug: {opt_names}", err=True)
            click.echo("Tip: Use --help-advanced to see all advanced options.", err=True)
            sys.exit(2)

    # Hide or reveal advanced commands in help based on --debug
    advanced_cmds = {"init-config", "show-checkpoint", "validate", "benchmark", "run"}
    group = cast(click.Group, ctx.command)
    for name, cmd in list(group.commands.items()):
        if name in advanced_cmds:
            cmd.hidden = not debug

    # If a subcommand was invoked, do not run the pipeline here
    if ctx.invoked_subcommand:
        return

    # Setup logging based on verbosity
    if debug:
        log_level = logging.DEBUG
    elif verbose >= 2:
        log_level = logging.DEBUG
    elif verbose == 1:
        log_level = logging.INFO
    else:
        log_level = logging.ERROR

    setup_logging(level=log_level, log_file=log_file)
    logger = get_logger("cli")

    try:
        # Generate default config and exit
        if generate_config:
            try:
                from circleseeker.resources import get_default_config

                click.echo(get_default_config())
            except Exception as exc:
                logger.error(f"Failed to generate default config: {exc}")
                sys.exit(1)
            return

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
            noise=verbose,
            debug=debug,
        )
        execute_pipeline(opts, logger)

    except CircleSeekerError as exc:
        logger.error(f"Pipeline error: {exc}")
        sys.exit(1)
    except Exception as exc:
        logger.exception(f"Unexpected error: {exc}")
        sys.exit(2)


cli.add_command(run)
cli.add_command(show_checkpoint)
cli.add_command(init_config)
cli.add_command(validate)


def main(argv: list[str] | None = None) -> int:
    """Main entry point."""
    try:
        cli(argv)
        return 0
    except Exception as exc:
        click.echo(f"Error: {exc}", err=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())

