"""Configuration-related CLI commands."""

from __future__ import annotations

from pathlib import Path

import click


@click.command(name="init-config", hidden=True)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=Path("config.yaml"),
    help="Output configuration file",
)
def init_config(output: Path) -> None:
    """Generate a template configuration file."""
    from circleseeker.resources import get_default_config

    config_text = get_default_config()
    output.write_text(config_text)
    click.echo(f"Configuration template saved to: {output}")
    click.echo("Edit this file to customize your analysis parameters.")

