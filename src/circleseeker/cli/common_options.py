"""Shared Click options for CircleSeeker CLI commands.

This module defines reusable Click option decorators to ensure consistency
between the main CLI command and subcommands like `run`.
"""

from __future__ import annotations

from functools import wraps
from pathlib import Path
from typing import Callable, TypeVar

import click

F = TypeVar("F", bound=Callable[..., None])


def input_option(func: F) -> F:
    """Input FASTA file option."""
    return click.option(
        "-i",
        "--input",
        "input_file",
        type=click.Path(exists=True, path_type=Path),
        required=False,
        help="Input FASTA file (HiFi reads)",
    )(func)


def reference_option(func: F) -> F:
    """Reference genome option."""
    return click.option(
        "-r",
        "--reference",
        type=click.Path(exists=True, path_type=Path),
        required=False,
        help="Reference genome FASTA file",
    )(func)


def output_option(func: F) -> F:
    """Output directory option."""
    return click.option(
        "-o",
        "--output",
        type=click.Path(path_type=Path),
        default=None,
        help="Output directory [default: circleseeker_output]",
    )(func)


def prefix_option(func: F) -> F:
    """Output prefix option."""
    return click.option(
        "-p",
        "--prefix",
        default=None,
        help="Output file prefix [default: sample]",
    )(func)


def config_option(func: F) -> F:
    """Configuration file option."""
    return click.option(
        "-c",
        "--config",
        type=click.Path(exists=True, path_type=Path),
        help="Configuration file (YAML)",
    )(func)


def threads_option(func: F) -> F:
    """Threads option."""
    return click.option(
        "-t",
        "--threads",
        type=int,
        default=None,
        help="Number of threads [default: 8]",
    )(func)


def keep_tmp_option(func: F) -> F:
    """Keep temporary files option."""
    return click.option(
        "--keep-tmp",
        is_flag=True,
        help="Retain temporary working directory (default: remove)",
    )(func)


def verbose_option(func: F) -> F:
    """Verbosity option (unified: use -v/--verbose everywhere)."""
    return click.option(
        "-v",
        "--verbose",
        count=True,
        help="Increase log verbosity (-v for INFO, -vv for DEBUG)",
    )(func)


def log_file_option(hidden: bool = False) -> Callable[[F], F]:
    """Log file option (unified name: --log-file)."""
    def decorator(func: F) -> F:
        return click.option(
            "--log-file",
            type=click.Path(path_type=Path),
            help="Path for log file output",
            hidden=hidden,
        )(func)
    return decorator


def start_from_option(hidden: bool = False) -> Callable[[F], F]:
    """Start from step option."""
    def decorator(func: F) -> F:
        return click.option(
            "--start-from",
            type=int,
            help="Resume from specific step (1-based index)",
            hidden=hidden,
        )(func)
    return decorator


def stop_at_option(hidden: bool = False) -> Callable[[F], F]:
    """Stop at step option."""
    def decorator(func: F) -> F:
        return click.option(
            "--stop-at",
            type=int,
            help="Stop at specific step (1-based index)",
            hidden=hidden,
        )(func)
    return decorator


def resume_option(hidden: bool = False) -> Callable[[F], F]:
    """Resume from checkpoint option."""
    def decorator(func: F) -> F:
        return click.option(
            "--resume",
            is_flag=True,
            help="Resume from last checkpoint",
            hidden=hidden,
        )(func)
    return decorator


def force_option(hidden: bool = False) -> Callable[[F], F]:
    """Force re-run option."""
    def decorator(func: F) -> F:
        return click.option(
            "--force",
            is_flag=True,
            help="Force re-run all steps",
            hidden=hidden,
        )(func)
    return decorator


def dry_run_option(hidden: bool = False) -> Callable[[F], F]:
    """Dry run option."""
    def decorator(func: F) -> F:
        return click.option(
            "--dry-run",
            is_flag=True,
            help="Show actions without executing",
            hidden=hidden,
        )(func)
    return decorator


def show_steps_option(hidden: bool = False) -> Callable[[F], F]:
    """Show steps option."""
    def decorator(func: F) -> F:
        return click.option(
            "--show-steps",
            is_flag=True,
            help="Show pipeline steps and exit",
            hidden=hidden,
        )(func)
    return decorator


def common_pipeline_options(func: F) -> F:
    """Apply all common pipeline options to a command.

    This decorator applies the standard set of options used by both
    the main CLI command and the `run` subcommand.

    Usage:
        @click.command()
        @common_pipeline_options
        def my_command(input_file, reference, output, ...):
            pass
    """
    # Apply options in reverse order (Click applies them bottom-up)
    decorators = [
        input_option,
        reference_option,
        output_option,
        prefix_option,
        config_option,
        threads_option,
        keep_tmp_option,
        verbose_option,
    ]
    for decorator in reversed(decorators):
        func = decorator(func)
    return func


def advanced_pipeline_options(hidden: bool = True) -> Callable[[F], F]:
    """Apply advanced pipeline options to a command.

    Args:
        hidden: Whether to hide these options in help (default: True)

    Usage:
        @click.command()
        @common_pipeline_options
        @advanced_pipeline_options(hidden=not debug)
        def my_command(...):
            pass
    """
    def decorator(func: F) -> F:
        decorators = [
            start_from_option(hidden=hidden),
            stop_at_option(hidden=hidden),
            resume_option(hidden=hidden),
            force_option(hidden=hidden),
            dry_run_option(hidden=hidden),
            show_steps_option(hidden=hidden),
            log_file_option(hidden=hidden),
        ]
        for dec in reversed(decorators):
            func = dec(func)
        return func
    return decorator
