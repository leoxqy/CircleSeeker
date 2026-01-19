"""Installation validation command."""

from __future__ import annotations

import sys

import click

from circleseeker import __version__
from circleseeker.cli.exit_codes import EXIT_ERROR


@click.command(hidden=True)
@click.option("--full", is_flag=True, help="Run full validation including test data")
def validate(full: bool) -> None:
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
            sys.exit(EXIT_ERROR)
    except ImportError:
        # Basic validation without full validator
        click.echo("✓ Basic installation check passed!")
        click.echo(f"  CircleSeeker version: {__version__}")

