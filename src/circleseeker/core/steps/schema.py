"""Schema validation utilities for step IO contracts."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import TYPE_CHECKING, Iterable, Optional

from circleseeker.exceptions import PipelineError

from .contracts import ArtifactSpec, StepContract

if TYPE_CHECKING:
    from circleseeker.core.pipeline import Pipeline


def _format_template(template: str, *, prefix: str) -> str:
    try:
        return template.format(prefix=prefix)
    except Exception:
        return template


def resolve_artifact_path(pipeline: Pipeline, spec: ArtifactSpec) -> Optional[Path]:
    """Resolve an ArtifactSpec to a concrete path for a specific pipeline run."""
    if spec.base == "output":
        base = Path(pipeline.config.output_dir)
        return base / _format_template(spec.template, prefix=pipeline.config.prefix)
    if spec.base == "final":
        base = Path(pipeline.final_output_dir)
        return base / _format_template(spec.template, prefix=pipeline.config.prefix)
    if spec.base == "config":
        value = getattr(pipeline.config, spec.template, None)
        if value is None:
            return None
        return Path(value)
    return Path(_format_template(spec.template, prefix=pipeline.config.prefix))


def _read_header(path: Path, *, delimiter: str) -> list[str]:
    with open(path, newline="", encoding="utf-8-sig") as f:
        reader = csv.reader(f, delimiter=delimiter)
        return next(reader, [])


def _validate_required_columns(path: Path, spec: ArtifactSpec) -> None:
    if not spec.required_columns:
        return

    try:
        header = _read_header(path, delimiter=spec.sep)
    except Exception as exc:
        raise PipelineError(f"Failed to read header for {path}: {exc}") from exc

    if not header:
        return

    missing = set(spec.required_columns) - set(header)
    if missing:
        coord_note = f" ({spec.coordinate_system})" if spec.coordinate_system else ""
        raise PipelineError(
            f"Schema validation failed for {path}: missing columns {sorted(missing)}{coord_note}"
        )


def _validate_tsv_no_header(
    path: Path, *, min_fields: int = 12, sample_lines: int = 1000
) -> None:
    """Validate a TSV file without header by checking multiple lines.

    Args:
        path: Path to the TSV file
        min_fields: Minimum number of fields expected per line
        sample_lines: Maximum number of lines to check (default: 1000)
    """
    errors: list[str] = []
    with open(path, encoding="utf-8", errors="replace") as f:
        for line_num, line in enumerate(f, 1):
            if line_num > sample_lines:
                break
            line = line.strip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < min_fields:
                errors.append(
                    f"Line {line_num}: {len(fields)} fields < {min_fields} required"
                )

    if errors:
        error_summary = "\n".join(errors[:5])
        if len(errors) > 5:
            error_summary += f"\n... and {len(errors) - 5} more errors"
        raise PipelineError(
            f"Schema validation failed for {path}:\n{error_summary}"
        )


def _validate_artifacts(
    pipeline: Pipeline, artifacts: Iterable[ArtifactSpec], *, phase: str, step_name: str
) -> None:
    for spec in artifacts:
        path = resolve_artifact_path(pipeline, spec)
        if path is None:
            if spec.required:
                raise PipelineError(
                    f"Schema validation failed for step '{step_name}' ({phase}): "
                    f"required artifact '{spec.name}' is not configured ({spec.base}:{spec.template})"
                )
            continue

        if spec.kind == "dir":
            if spec.required and not path.exists():
                raise PipelineError(
                    f"Schema validation failed for step '{step_name}' ({phase}): "
                    f"missing directory {path}"
                )
            continue

        if spec.required and not path.exists():
            raise PipelineError(
                f"Schema validation failed for step '{step_name}' ({phase}): missing file {path}"
            )
        if not path.exists():
            continue

        try:
            if path.is_file() and path.stat().st_size == 0:
                continue
        except OSError:
            continue

        if spec.kind == "csv":
            _validate_required_columns(path, spec)
        elif spec.kind == "tsv":
            _validate_required_columns(path, spec)
        elif spec.kind == "tsv_no_header":
            _validate_tsv_no_header(path)


def validate_step_contract(pipeline: Pipeline, contract: StepContract, *, when: str) -> None:
    """Validate inputs/outputs for a step contract.

    `when` is either "pre" (validate inputs) or "post" (validate outputs).
    """
    if when == "pre":
        _validate_artifacts(pipeline, contract.inputs, phase="pre", step_name=contract.step_name)
        return
    if when == "post":
        _validate_artifacts(pipeline, contract.outputs, phase="post", step_name=contract.step_name)
        return
    raise ValueError(f"Invalid when={when!r}; expected 'pre' or 'post'")
