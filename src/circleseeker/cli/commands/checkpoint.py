"""`show-checkpoint` subcommand implementation."""

from __future__ import annotations

import logging
from pathlib import Path

import click


@click.command(name="show-checkpoint", hidden=True)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=Path("circleseeker_output"),
    help="Output directory containing checkpoint",
)
@click.option("-p", "--prefix", default="sample", help="Run prefix (used to locate checkpoint)")
def show_checkpoint(output: Path, prefix: str) -> None:
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
    except OSError as exc:
        # Log but continue with default prefix
        logging.getLogger("cli").debug(f"Could not auto-detect prefix: {exc}")

    cfg = Config()
    cfg.output_dir = output
    cfg.prefix = detected_prefix
    pipeline = Pipeline(cfg)
    pipeline.show_checkpoint_info()

