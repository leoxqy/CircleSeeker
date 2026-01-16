"""`run` subcommand implementation."""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Optional

import click

from circleseeker.exceptions import CircleSeekerError
from circleseeker.utils.logging import get_logger, setup_logging

from ..common_options import (
    common_pipeline_options,
    start_from_option,
    stop_at_option,
    resume_option,
    force_option,
    dry_run_option,
    show_steps_option,
    log_file_option,
)
from ..pipeline import PipelineOptions, execute_pipeline


@click.command(hidden=True)
@common_pipeline_options
@start_from_option(hidden=False)
@stop_at_option(hidden=False)
@resume_option(hidden=False)
@force_option(hidden=False)
@dry_run_option(hidden=False)
@show_steps_option(hidden=False)
@log_file_option(hidden=False)
def run(
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
    force: bool,
    dry_run: bool,
    show_steps: bool,
    log_file: Optional[Path],
) -> None:
    """Run the CircleSeeker pipeline for eccDNA detection."""

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
            start_from=start_from,
            stop_at=stop_at,
            resume=resume,
            show_steps=show_steps,
            force=force,
            dry_run=dry_run,
            log_file=log_file,
            noise=verbose,
            debug=verbose >= 2,
        )
        execute_pipeline(opts, logger)

    except CircleSeekerError as exc:
        logger.error(f"Pipeline error: {exc}")
        sys.exit(1)
    except Exception as exc:
        logger.exception(f"Unexpected error: {exc}")
        sys.exit(2)
