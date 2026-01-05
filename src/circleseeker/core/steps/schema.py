"""Schema validation utilities for step IO contracts."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable, Optional

from circleseeker.exceptions import PipelineError

from .contracts import ArtifactSpec, StepContract


def _format_template(template: str, *, prefix: str) -> str:
    try:
        return template.format(prefix=prefix)
    except Exception:
        return template


def resolve_artifact_path(pipeline, spec: ArtifactSpec) -> Optional[Path]:
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


def _validate_tsv_no_header(path: Path, *, min_fields: int = 12) -> None:
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < min_fields:
                raise PipelineError(
                    f"Schema validation failed for {path}: expected >= {min_fields} fields, got {len(fields)}"
                )
            return


def _validate_artifacts(
    pipeline, artifacts: Iterable[ArtifactSpec], *, phase: str, step_name: str
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


def validate_step_contract(pipeline, contract: StepContract, *, when: str) -> None:
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
