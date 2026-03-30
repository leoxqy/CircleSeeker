"""Shared pipeline execution helpers for the CLI."""

from __future__ import annotations

import logging
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import click

from circleseeker.config import Config, load_config, save_config, apply_preset, apply_platform_preset, PRESETS
from circleseeker.cli.exit_codes import EXIT_ERROR, EXIT_USAGE
from circleseeker.utils.display import ConsoleFormatter, print_formatted
from circleseeker.utils.logging import setup_logging


# Advanced-only option names (validated at runtime)
ADVANCED_ONLY_KEYS = {
    "start_from": "--start-from",
    "stop_at": "--stop-at",
    "resume": "--resume",
    "force": "--force",
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
    turbo: bool  # CLI flag: True if --turbo provided (use /dev/shm)
    fast_align: bool  # CLI flag: True if --fast-align provided
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
    # Sequencing platform: hifi (default), ont, ngs
    platform: Optional[str] = None
    # Sensitivity preset: relaxed, balanced (default), strict
    preset: Optional[str] = None
    # R2 FASTQ file for NGS paired-end mode (R1 is input_file)
    r2_file: Optional[Path] = None


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


def execute_pipeline(
    opts: PipelineOptions,
    logger: logging.Logger,
    ctx: Optional[click.Context] = None,
) -> None:
    """
    Execute the CircleSeeker pipeline with given options.

    This is the unified execution function used by the top-level command.

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

    # Apply platform preset if specified
    # Priority: CLI --platform > config file > default (hifi)
    if opts.platform:
        cfg.platform = opts.platform
    if cfg.platform not in ("hifi", "ngs"):
        try:
            apply_platform_preset(cfg)
        except (ValueError, KeyError) as exc:
            click.echo(f"Error: {exc}", err=True)
            sys.exit(EXIT_ERROR)

    # Apply sensitivity preset if specified
    # Preset is applied AFTER platform preset but BEFORE CLI overrides
    # This allows: defaults -> config file -> platform -> preset -> CLI args
    if opts.preset:
        try:
            apply_preset(cfg, opts.preset)  # type: ignore[arg-type]
        except (ValueError, KeyError) as exc:
            click.echo(f"Error: {exc}", err=True)
            sys.exit(EXIT_ERROR)

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
        sys.exit(EXIT_ERROR)

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

    # threads: CLI > config (config already has default 8)
    if opts.threads is not None:
        cfg.threads = opts.threads

    # keep_tmp: CLI --keep-tmp flag > config file > default False
    # When CLI flag is set, it always takes precedence
    if opts.keep_tmp:
        cfg.keep_tmp = True
    # Otherwise, keep config file value (cfg.keep_tmp already has default False)

    # turbo: CLI --turbo flag > config file > default False
    # When CLI flag is set, it always takes precedence
    if opts.turbo:
        cfg.runtime.turbo_mode = True
    # Otherwise, keep config file value (cfg.runtime.turbo_mode already has default False)

    # fast_align: CLI --fast-align flag > config file > default False
    if opts.fast_align:
        cecc_cfg = cfg.tools.cecc_build
        if isinstance(cecc_cfg, dict):
            cecc_cfg["fast_last"] = True
        else:
            cecc_cfg.fast_last = True

    # ── NGS platform: independent execution path ──────────────────────
    # NGS uses paired-end FASTQ (R1 via -i, R2 via --r2) and bypasses
    # the HiFi/ONT long-read pipeline entirely.
    if cfg.platform == "ngs":
        _execute_ngs_pipeline(opts, cfg, input_file, reference, logger)
        return

    # Validate configuration FIRST (fail fast on user errors like missing files)
    from circleseeker.exceptions import ConfigurationError

    try:
        cfg.validate()
    except ConfigurationError as exc:
        click.echo(f"Error: {exc}", err=True)
        sys.exit(EXIT_ERROR)

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
        click.echo(f"Platform: {cfg.platform}")
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
    except OSError as exc:
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

    # Print header
    print_formatted(formatter.header())

    # Print configuration
    # Use user_output_dir for display (cfg.output_dir was changed to temp_dir by Pipeline)
    preset_display = opts.preset.capitalize() if opts.preset else "Balanced"

    # Determine turbo mode status for display
    turbo_requested = opts.turbo or cfg.runtime.turbo_mode
    turbo_active = getattr(pipeline, "_turbo_mode_active", False)
    if turbo_requested:
        if turbo_active:
            turbo_status = "Enabled (RAM-backed)"
        else:
            turbo_status = "Fallback (unavailable)"
    else:
        turbo_status = None  # Don't show if not requested

    config_dict = {
        "input_file": input_file.name,
        "reference": reference.name,
        "output_dir": str(user_output_dir.absolute()),
        "platform": cfg.platform.upper(),
        "threads": cfg.threads,
    }
    if turbo_status:
        config_dict["turbo"] = turbo_status
    print_formatted(formatter.format_config(config_dict))

    # Show warning if turbo mode was requested but not available
    if turbo_requested and not turbo_active:
        print_formatted("  [!] Turbo mode unavailable: /dev/shm not accessible or insufficient space")

    print_formatted(formatter.separator())
    print_formatted(formatter.start_message())

    # Run pipeline
    pipeline.run(start_from=opts.start_from, stop_at=opts.stop_at, force=force)

    # Print beautiful completion message
    print_formatted(formatter.separator())
    print_formatted(formatter.success_message())

    # Surface soft-fail warnings prominently in terminal output
    from circleseeker.core.pipeline_types import ResultKeys

    soft_fail_msgs: list[str] = []
    if pipeline.state.results.get(ResultKeys.CECC_BUILD_FAILED):
        soft_fail_msgs.append(
            "[!] CeccDNA detection (cecc_build) FAILED — 0 CeccDNA in results. Check logs for details."
        )
    if pipeline.state.results.get(ResultKeys.INFERENCE_FAILED):
        soft_fail_msgs.append(
            "[!] SplitReads inference FAILED — 0 Inferred eccDNA in results. Check logs for details."
        )
    for msg in soft_fail_msgs:
        print_formatted(f"  {msg}")

    if not cfg.keep_tmp:
        print_formatted(formatter.cleanup_message())
    # Use user_output_dir for display (cfg.output_dir was changed to temp_dir by Pipeline)
    print_formatted(formatter.final_output_message(user_output_dir))
    print_formatted(formatter.separator("="))


def _execute_ngs_pipeline(
    opts: PipelineOptions,
    cfg: Config,
    input_file: Optional[Path],
    reference: Optional[Path],
    logger: logging.Logger,
) -> None:
    """Execute the NGS (short-read paired-end) eccDNA detection pipeline.

    This is a separate execution path from the HiFi/ONT long-read pipeline.
    It validates NGS-specific inputs and delegates to NGSCircleDetect.
    """
    from circleseeker.utils.display import ConsoleFormatter, print_formatted

    # ── Validate NGS-specific inputs ──────────────────────────────────
    if not input_file:
        click.echo(
            "Error: --input (-i) is required for NGS mode (R1 FASTQ file)",
            err=True,
        )
        sys.exit(EXIT_ERROR)

    if not reference:
        click.echo("Error: --reference (-r) is required for NGS mode", err=True)
        sys.exit(EXIT_ERROR)

    r2_file = opts.r2_file

    # Determine input mode: paired-end FASTQ (R1+R2) or sorted BAM
    is_bam_input = input_file.suffix.lower() in (".bam",)

    if not is_bam_input and r2_file is None:
        click.echo(
            "Error: NGS mode requires paired-end FASTQ input.\n"
            "  Provide R1 via -i and R2 via --r2:\n"
            "    circleseeker --platform ngs -i R1.fq -r ref.fa --r2 R2.fq\n"
            "  Or provide a sorted BAM file directly:\n"
            "    circleseeker --platform ngs -i sorted.bam -r ref.fa",
            err=True,
        )
        sys.exit(EXIT_ERROR)

    if not reference.exists():
        click.echo(f"Error: Reference file not found: {reference}", err=True)
        sys.exit(EXIT_ERROR)

    if not input_file.exists():
        click.echo(f"Error: Input file not found: {input_file}", err=True)
        sys.exit(EXIT_ERROR)

    if r2_file is not None and not r2_file.exists():
        click.echo(f"Error: R2 FASTQ file not found: {r2_file}", err=True)
        sys.exit(EXIT_ERROR)

    # ── Resolve output directory and prefix ───────────────────────────
    output_dir = cfg.output_dir
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    prefix = cfg.prefix

    # ── Handle dry run ────────────────────────────────────────────────
    if opts.dry_run:
        logger.info("Dry run mode - showing what would be executed:")
        click.echo(f"\nNGS eccDNA detection pipeline")
        click.echo(f"  Input (R1): {input_file}")
        if r2_file:
            click.echo(f"  Input (R2): {r2_file}")
        else:
            click.echo(f"  Input (BAM): {input_file}")
        click.echo(f"  Reference:  {reference}")
        click.echo(f"  Output to:  {output_dir.absolute()}")
        click.echo(f"  Platform:   NGS")
        click.echo(f"  Threads:    {cfg.threads}")
        return

    # ── Run NGS detection ─────────────────────────────────────────────
    from circleseeker.modules.ngs_circle_detect import NGSCircleDetect, NGSDetectConfig

    ngs_cfg = NGSDetectConfig(threads=cfg.threads)
    detector = NGSCircleDetect(reference=reference, config=ngs_cfg, logger=logger)

    # Print header
    formatter = ConsoleFormatter()
    print_formatted(formatter.header())

    config_dict = {
        "input_file": input_file.name,
        "reference": reference.name,
        "output_dir": str(output_dir.absolute()),
        "platform": "NGS",
        "threads": cfg.threads,
    }
    if r2_file:
        config_dict["r2_file"] = r2_file.name
    print_formatted(formatter.format_config(config_dict))
    print_formatted(formatter.separator())
    print_formatted(formatter.start_message())

    if is_bam_input:
        logger.info("NGS mode: detecting eccDNA from BAM: %s", input_file)
        df = detector.run_from_bam(input_file, output_dir, prefix=prefix)
    else:
        logger.info(
            "NGS mode: detecting eccDNA from FASTQ pair: %s, %s",
            input_file,
            r2_file,
        )
        df = detector.run_from_fastq(input_file, r2_file, output_dir, prefix=prefix)

    # ── Report results ────────────────────────────────────────────────
    print_formatted(formatter.separator())
    n_total = len(df) if df is not None else 0
    if n_total > 0:
        type_counts = df["eccdna_type"].value_counts().to_dict() if "eccdna_type" in df.columns else {}
        type_summary = ", ".join(f"{k}: {v}" for k, v in sorted(type_counts.items()))
        print_formatted(f"  NGS detection complete: {n_total} eccDNA detected ({type_summary})")
        out_csv = output_dir / f"{prefix}_ngs_circles.csv"
        print_formatted(f"  Results saved to: {out_csv}")
    else:
        print_formatted("  NGS detection complete: 0 eccDNA detected")

    print_formatted(formatter.separator("="))
