"""Configuration-related CLI commands."""

from __future__ import annotations

from pathlib import Path

import click


@click.command(name="init-config")
@click.option(
    "--output-file",
    type=click.Path(path_type=Path),
    default=Path("config.yaml"),
    help="Output configuration file path",
)
@click.option(
    "--stdout",
    is_flag=True,
    help="Print default config YAML to stdout instead of writing a file",
)
def init_config(output_file: Path, stdout: bool) -> None:
    """Generate a template configuration file."""
    from circleseeker.resources import get_default_config

    config_text = get_default_config()
    if stdout:
        click.echo(config_text)
    else:
        output_file.write_text(config_text)
        click.echo(f"Configuration template saved to: {output_file}")
        click.echo("Edit this file to customize your analysis parameters.")

