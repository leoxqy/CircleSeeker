"""Main pipeline orchestrator for CircleSeeker"""

from __future__ import annotations

import json
import os
import time
import hashlib
import shutil
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable, Optional, Union
from dataclasses import asdict

from circleseeker.config import Config
from circleseeker.exceptions import ConfigurationError, PipelineError
from circleseeker.utils.progress import iter_progress
from circleseeker.utils.logging import get_logger
from circleseeker.core.pipeline_types import (
    PipelineState,
    PipelineStep,
    ResultKeys,
    StepMetadata,
)
from circleseeker.core.steps.definitions import PIPELINE_STEPS

# Step execution lives in `circleseeker.core.steps.*` and is imported lazily by
# wrapper methods on `Pipeline` to keep imports light.

class Pipeline:
    """Main pipeline orchestrator with complete fixes."""

    # Canonical step ordering and user-facing metadata.
    STEPS = PIPELINE_STEPS

    # Instance attributes with type annotations
    state: PipelineState
    temp_dir: Path
    final_output_dir: Path
    _inference_tool: Optional[str]
    _temp_dir_safe_to_delete: bool

    def _set_result(self, key: str, value: Any) -> None:
        """Set a result key (new canonical names only)."""
        self.state.results[key] = value

    def _serialize_path_for_state(self, path: Path | str) -> str:
        """Serialize a filesystem path with relativisation hints for checkpoint storage."""
        candidate = Path(path)
        try:
            resolved = candidate.resolve(strict=False)
        except OSError:
            resolved = candidate

        base_candidates: list[tuple[str, Optional[Path]]] = [
            ("output", getattr(self.config, "output_dir", None)),
            ("final", getattr(self, "final_output_dir", None)),
        ]

        for tag, base in base_candidates:
            if base is None:
                continue
            base_path = Path(base)
            try:
                base_resolved = base_path.resolve(strict=False)
            except OSError:
                base_resolved = base_path
            try:
                rel = resolved.relative_to(base_resolved)
                return f"{tag}::{rel.as_posix()}"
            except ValueError:
                continue

        return resolved.as_posix()

    def _resolve_stored_path(
        self,
        stored_path: Optional[str],
        fallback_names: Iterable[Union[str, Path]] = (),
    ) -> Optional[Path]:
        """Resolve stored or fallback paths to an existing Path if possible."""

        candidates: list[Path] = []
        seen: set[Path] = set()

        def add_candidate(candidate: Union[str, Path]) -> None:
            path_candidate = Path(candidate)
            if path_candidate in seen:
                return
            seen.add(path_candidate)
            candidates.append(path_candidate)

        if stored_path:
            if "::" in stored_path:
                prefix, rel_fragment = stored_path.split("::", 1)
                rel_path = Path(rel_fragment)
                base_lookup = {
                    "output": getattr(self.config, "output_dir", None),
                    "final": getattr(self, "final_output_dir", None),
                }
                base_dir = base_lookup.get(prefix)
                if base_dir is not None:
                    add_candidate(Path(base_dir) / rel_path)
                add_candidate(rel_path)
            else:
                add_candidate(stored_path)

        for name in fallback_names:
            fallback_path = Path(name)
            if fallback_path.is_absolute():
                add_candidate(fallback_path)
                continue
            add_candidate(fallback_path)
            add_candidate(Path(self.config.output_dir) / fallback_path)
            final_dir = getattr(self, "final_output_dir", None)
            if final_dir is not None:
                add_candidate(Path(final_dir) / fallback_path)

        for candidate in list(candidates):
            if not candidate.is_absolute():
                add_candidate(Path.cwd() / candidate)

        for candidate in candidates:
            try:
                resolved = candidate.resolve(strict=False)
            except OSError:
                resolved = candidate
            if resolved.exists():
                return resolved

        return None

    def _get_result(self, key: str, default: Any = None) -> Any:
        """Get a result by canonical key only."""
        return self.state.results.get(key, default)

    def _del_result(self, key: str) -> None:
        """Delete a result by key if it exists."""
        self.state.results.pop(key, None)

    def _select_inference_tool(self) -> str:
        """Determine which inference tool to use.

        SplitReads-Core is the built-in tool, optimized for HiFi data.
        It is always available as all required dependencies are bundled.
        """
        if self._inference_tool is None:
            # SplitReads-Core is built-in and always available
            self._inference_tool = "splitreads_core"
            self.logger.debug("Selected SplitReads-Core for eccDNA inference (built-in, HiFi optimized)")
        return self._inference_tool

    def _ensure_reference_mmi(self) -> Path:
        """Ensure minimap2 reference index (.mmi) exists and return its path."""
        reference = self.config.reference
        if reference is None:
            raise PipelineError("Reference genome is required")
        reference = Path(reference)

        def _canonical_mmi_path(ref: Path) -> Path:
            """Return the canonical minimap2 index path produced by `minimap2 -d` conventions."""
            ref_str = str(ref)
            if ref_str.endswith(".fasta"):
                return Path(ref_str[:-6] + ".mmi")
            if ref_str.endswith(".fa"):
                return Path(ref_str[:-3] + ".mmi")
            if ref_str.endswith(".fna"):
                return Path(ref_str[:-4] + ".mmi")
            return Path(ref_str + ".mmi")

        # Prefer existing index next to the reference. If missing (or not writable), create an index
        # inside the pipeline temp directory.
        canonical_next_to_ref = _canonical_mmi_path(reference)
        legacy_next_to_ref = Path(str(reference) + ".mmi")
        mmi_in_temp = self.temp_dir / canonical_next_to_ref.name

        for candidate in (canonical_next_to_ref, legacy_next_to_ref, mmi_in_temp):
            if candidate.exists():
                return candidate

        reference_mmi = mmi_in_temp
        self.logger.info(f"Reference MMI index not found; creating: {reference_mmi}")
        minimap2_cmd = ["minimap2", "-d", str(reference_mmi), str(reference)]

        log_file = self.temp_dir / "minimap2_index.log"

        def _read_tail(path: Path, max_bytes: int = 32_000) -> str:
            try:
                with open(path, "rb") as handle:
                    try:
                        handle.seek(0, 2)
                        size = handle.tell()
                        handle.seek(max(0, size - max_bytes))
                    except OSError:
                        pass
                    data = handle.read()
                return data.decode(errors="replace")
            except OSError:
                return ""

        try:
            with open(log_file, "w") as log_handle:
                subprocess.run(
                    minimap2_cmd,
                    check=True,
                    stdout=log_handle,
                    stderr=log_handle,
                    text=True,
                )
        except subprocess.CalledProcessError as exc:
            stderr_tail = _read_tail(log_file)
            message = "Failed to create reference MMI index"
            if stderr_tail:
                message = f"{message}. Last minimap2 output:\n{stderr_tail[-2000:]}"
            raise PipelineError(message) from exc
        self.logger.info(f"Created reference index: {reference_mmi}")
        return reference_mmi

    def _create_combined_fasta(self, output_file: Path) -> None:
        """Combine all FASTA outputs into a single file."""
        fasta_files: list[Path] = []

        # Collect all FASTA files
        patterns = [
            f"{self.config.prefix}_UeccDNA*.fasta",
            f"{self.config.prefix}_MeccDNA*.fasta",
            f"{self.config.prefix}_CeccDNA*.fasta",
            "tandem_to_ring.fasta",
            f"{self.config.prefix}_circular.fasta",
        ]

        for pattern in patterns:
            fasta_files.extend(self.config.output_dir.glob(pattern))

        # Also check state results for FASTA files
        # Use _resolve_stored_path to handle serialized paths like "output::..." or "final::..."
        for key in self.state.results:
            if key.endswith("_fasta"):
                stored_path = self.state.results[key]
                fasta_path = self._resolve_stored_path(stored_path)
                if fasta_path and fasta_path.exists() and fasta_path not in fasta_files:
                    fasta_files.append(fasta_path)

        # Combine all FASTA files
        with open(output_file, "w") as out_f:
            for fasta in fasta_files:
                if fasta.exists() and fasta.stat().st_size > 0:
                    with open(fasta, "r") as in_f:
                        contents = in_f.read()
                        if not contents:
                            continue
                        out_f.write(contents)
                        if not contents.endswith("\n"):
                            out_f.write("\n")

    def _safe_copy_file(self, src: Optional[Path], dest: Path) -> bool:
        """Copy a single file into the final output directory if it exists."""
        if not src:
            return False

        src_path = Path(src)
        if not src_path.exists():
            return False

        try:
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src_path, dest)
            self.logger.debug(f"Copied {src_path} -> {dest}")
            return True
        except Exception as exc:
            self.logger.warning(f"Failed to copy {src_path} to {dest}: {exc}")
            return False

    def _copy_directory_contents(self, src: Optional[Path], dest: Path) -> bool:
        """Copy all contents of a directory into destination if available."""
        if not src:
            return False

        src_path = Path(src)
        if not src_path.exists() or not src_path.is_dir():
            return False

        try:
            dest.mkdir(parents=True, exist_ok=True)
            for item in src_path.iterdir():
                target = dest / item.name
                if item.is_dir():
                    shutil.copytree(item, target, dirs_exist_ok=True)
                elif item.is_file():
                    shutil.copy2(item, target)
            self.logger.debug(f"Copied directory contents {src_path} -> {dest}")
            return True
        except Exception as exc:
            self.logger.warning(f"Failed to copy directory {src_path} to {dest}: {exc}")
            return False

    def _create_basic_output_structure(
        self,
        final_root: Path,
        *,
        uecc_dir: Optional[Path],
        mecc_dir: Optional[Path],
        cecc_dir: Optional[Path],
        inferred_dir: Optional[Path],
        merged_csv: Optional[Path],
        html_report: Optional[Path],
        text_summary: Optional[Path],
    ) -> None:
        """Fallback organizer when ecc_packager cannot run."""
        prefix = self.config.prefix
        final_root.mkdir(parents=True, exist_ok=True)

        directory_specs = [
            (uecc_dir, f"{prefix}_Confirmed_UeccDNA"),
            (mecc_dir, f"{prefix}_Confirmed_MeccDNA"),
            (cecc_dir, f"{prefix}_Confirmed_CeccDNA"),
            (inferred_dir, f"{prefix}_Inferred_eccDNA"),
        ]

        for src_dir, dest_name in directory_specs:
            dest_path = final_root / dest_name
            if self._copy_directory_contents(src_dir, dest_path):
                self.logger.info(f"Copied {dest_name} into final output")

        self._rename_inferred_simple_file(final_root)

        if merged_csv:
            merged_name = f"{prefix}_merged_output.csv"
            if not self._safe_copy_file(merged_csv, final_root / merged_name):
                self.logger.debug("Merged CSV not available for fallback packaging")

        if html_report:
            self._safe_copy_file(html_report, final_root / f"{prefix}_report.html")

        if text_summary:
            self._safe_copy_file(text_summary, final_root / f"{prefix}_summary.txt")

    def _rename_inferred_simple_file(self, final_root: Path) -> None:
        """Ensure inferred simple CSV follows the *_UeccDNA_I.csv naming convention."""
        prefix = self.config.prefix
        inferred_dir = final_root / f"{prefix}_Inferred_eccDNA"
        if not inferred_dir.exists():
            return

        legacy_simple = inferred_dir / f"{prefix}_simple.csv"
        target = inferred_dir / f"{prefix}_UeccDNA_I.csv"
        if not legacy_simple.exists():
            return

        try:
            if target.exists():
                target.unlink()
            legacy_simple.rename(target)
            self.logger.debug("Renamed inferred simple CSV to %s", target.name)
        except Exception as exc:
            self.logger.warning(
                "Could not rename %s to %s: %s",
                legacy_simple.name,
                target.name,
                exc,
            )

    def __init__(self, config: Config):
        self.config = config
        self.logger = get_logger(self.__class__.__name__)

        # Instance-level inference tool cache (thread-safe)
        self._inference_tool = None

        # Turbo mode state
        self._turbo_mode_active = False
        self._turbo_shm_path: Optional[Path] = None
        self._symlink_path: Optional[Path] = None

        # Setup directory structure
        self.final_output_dir = config.output_dir
        self.final_output_dir.mkdir(parents=True, exist_ok=True)

        # Track whether temp_dir is safe to auto-delete
        self._temp_dir_safe_to_delete = False

        # Check if turbo mode should be enabled
        if config.runtime.turbo_mode:
            turbo_result = self._setup_turbo_mode(config)
            if turbo_result:
                self.temp_dir, self._turbo_shm_path, self._symlink_path = turbo_result
                self._turbo_mode_active = True
                self._temp_dir_safe_to_delete = True
                self.logger.info(
                    f"Turbo mode enabled: temp files in {self._turbo_shm_path}, "
                    f"symlink at {self._symlink_path}"
                )
            else:
                # Fallback to normal mode
                self.logger.warning("Turbo mode unavailable, falling back to normal mode")
                self._setup_normal_mode(config)
        else:
            self._setup_normal_mode(config)

        # Create temp directory
        self.temp_dir.mkdir(parents=True, exist_ok=True)

        # Override config to use temp directory for all intermediate files
        self.config.output_dir = self.temp_dir

        # Keep checkpoint in final directory for persistence across runs
        self.state_file = self.final_output_dir / f"{config.prefix}.checkpoint"
        self.state = self._load_state()

    def _check_shm_space(self, min_space_gb: float) -> tuple[bool, float]:
        """Check if /dev/shm has enough space.

        Args:
            min_space_gb: Minimum required space in GB

        Returns:
            Tuple of (has_enough_space, available_space_gb)
        """
        shm_path = Path("/dev/shm")
        if not shm_path.exists():
            return False, 0.0

        try:
            import os
            stat = os.statvfs(shm_path)
            available_bytes = stat.f_bavail * stat.f_frsize
            available_gb = available_bytes / (1024 ** 3)
            return available_gb >= min_space_gb, available_gb
        except (OSError, AttributeError):
            return False, 0.0

    def _setup_turbo_mode(self, config: Config) -> Optional[tuple[Path, Path, Path]]:
        """Setup turbo mode with /dev/shm.

        Returns:
            Tuple of (temp_dir, shm_path, symlink_path) if successful, None otherwise
        """
        # Check /dev/shm availability and space
        min_space = config.runtime.turbo_min_space_gb
        has_space, available_gb = self._check_shm_space(min_space)

        if not has_space:
            self.logger.warning(
                f"/dev/shm has insufficient space: {available_gb:.1f}GB available, "
                f"{min_space:.1f}GB required"
            )
            return None

        # Generate unique directory name
        import os
        job_id = os.environ.get("SLURM_JOB_ID", str(os.getpid()))
        shm_dir_name = f"circleseeker_{config.prefix}_{job_id}"
        shm_path = Path("/dev/shm") / shm_dir_name

        # Create /dev/shm directory
        try:
            shm_path.mkdir(parents=True, exist_ok=True)
        except (OSError, PermissionError) as e:
            self.logger.warning(f"Cannot create directory in /dev/shm: {e}")
            return None

        # Create symlink in output directory
        symlink_path = self.final_output_dir / ".tmp"

        # Remove existing symlink or directory if exists
        if symlink_path.is_symlink():
            symlink_path.unlink()
        elif symlink_path.exists():
            # If it's a real directory, we can't use turbo mode safely
            self.logger.warning(
                f"Cannot enable turbo mode: {symlink_path} exists and is not a symlink"
            )
            try:
                shm_path.rmdir()
            except OSError:
                pass
            return None

        # Create symlink
        try:
            symlink_path.symlink_to(shm_path)
        except OSError as e:
            self.logger.warning(f"Cannot create symlink for turbo mode: {e}")
            try:
                shm_path.rmdir()
            except OSError:
                pass
            return None

        self.logger.info(f"Turbo mode: /dev/shm available space {available_gb:.1f}GB")
        return symlink_path, shm_path, symlink_path

    def _setup_normal_mode(self, config: Config) -> None:
        """Setup normal mode with temp directory under output_dir."""
        configured_tmp = config.runtime.tmp_dir if config.runtime.tmp_dir else Path(".tmp_work")

        if configured_tmp.is_absolute():
            # Absolute path: use as-is, but mark as unsafe to auto-delete
            self.temp_dir = configured_tmp
            self._temp_dir_safe_to_delete = False
        else:
            # Basic path safety (no '.', no '..') is enforced by Config.validate().
            # Here we only check output_dir-specific constraints.

            # Check if output_dir is already a temp directory from previous run
            if config.output_dir.name == configured_tmp.name:
                self.final_output_dir = config.output_dir.parent
                self.temp_dir = config.output_dir
            else:
                self.temp_dir = self.final_output_dir / configured_tmp

            self._temp_dir_safe_to_delete = True

            # Validate temp_dir is under final_output_dir
            try:
                final_resolved = Path(self.final_output_dir).resolve(strict=False)
                temp_resolved = Path(self.temp_dir).resolve(strict=False)
                if temp_resolved == final_resolved:
                    raise ConfigurationError(
                        "Invalid runtime.tmp_dir: resolves to output_dir. "
                        "Choose a subdirectory name such as '.tmp_work'."
                    )
                temp_resolved.relative_to(final_resolved)
            except ValueError as exc:
                raise ConfigurationError(
                    "Invalid runtime.tmp_dir: resolves outside output_dir. "
                    "Use an absolute path if you want an external temp directory."
                ) from exc

    def _compute_config_hash(self) -> str:
        """Compute hash of current configuration for validation.

        Normalizes paths to absolute form and uses 20 characters for reduced collision risk.
        """
        config_copy = asdict(self.config)

        def normalize_value(v: Any) -> Any:
            """Normalize a single value, resolving Path objects to absolute paths."""
            if isinstance(v, Path):
                try:
                    return str(v.resolve())
                except (OSError, ValueError):
                    return str(v)
            if isinstance(v, str) and v and (v.startswith("/") or v.startswith(".")):
                try:
                    return str(Path(v).resolve())
                except (OSError, ValueError):
                    return v
            return v

        def normalize_dict(d: Any) -> Any:
            """Recursively normalize dict values."""
            if isinstance(d, dict):
                return {k: normalize_dict(v) for k, v in d.items()}
            if isinstance(d, list):
                return [normalize_dict(i) for i in d]
            return normalize_value(d)

        normalized = normalize_dict(config_copy)
        config_str = json.dumps(normalized, sort_keys=True, default=str)
        return hashlib.sha256(config_str.encode()).hexdigest()[:20]

    def _load_state(self) -> PipelineState:
        """Load pipeline state from checkpoint with validation."""
        if self.state_file.exists():
            try:
                with open(self.state_file, "r") as f:
                    data = json.load(f)

                # Check version compatibility
                if data.get("version", "1.0") != "2.0":
                    self.logger.warning("Checkpoint from older version, starting fresh")
                    return self._create_fresh_state()

                # Check config compatibility
                current_hash = self._compute_config_hash()
                saved_hash = data.get("config_hash")
                if saved_hash and saved_hash != current_hash:
                    self.logger.warning("Configuration changed since last checkpoint")
                    policy = getattr(self.config.runtime, "checkpoint_policy", "continue")
                    if policy == "reset":
                        self.logger.info("Checkpoint policy = reset; starting fresh state")
                        return self._create_fresh_state()
                    elif policy == "fail":
                        raise PipelineError("Config changed and checkpoint_policy=fail")
                    else:
                        self.logger.info(
                            "Checkpoint policy = continue; resuming with existing checkpoint"
                        )

                # Reconstruct step metadata
                step_metadata = {}
                for step_name, meta_data in data.get("step_metadata", {}).items():
                    metadata = StepMetadata(**meta_data)
                    metadata.step_name = step_name
                    step_metadata[step_name] = metadata

                state = PipelineState(
                    completed_steps=data.get("completed_steps", []),
                    current_step=data.get("current_step"),
                    failed_step=data.get("failed_step"),
                    results=data.get("results", {}),
                    step_metadata=step_metadata,
                    pipeline_start_time=data.get("pipeline_start_time"),
                    last_checkpoint_time=data.get("last_checkpoint_time"),
                    config_hash=current_hash,
                    version=data.get("version", "2.0"),
                )

                self.logger.info(
                    f"Loaded checkpoint with {len(state.completed_steps)} completed steps"
                )
                if state.failed_step:
                    self.logger.warning(f"Previous run failed at step: {state.failed_step}")

                return state

            except json.JSONDecodeError as e:
                self.logger.error(f"Checkpoint file corrupted (invalid JSON): {e}")
                self.logger.info("Starting with fresh state")
            except (KeyError, TypeError, ValueError) as e:
                self.logger.error(f"Checkpoint file has incompatible format: {e}")
                self.logger.info("Starting with fresh state")
            except OSError as e:
                self.logger.error(f"Could not read checkpoint file: {e}")
                self.logger.info("Starting with fresh state")

        return self._create_fresh_state()

    def _create_fresh_state(self) -> PipelineState:
        """Create a fresh pipeline state."""
        return PipelineState(
            completed_steps=[],
            config_hash=self._compute_config_hash(),
            pipeline_start_time=time.time(),
        )

    def _save_state(self) -> None:
        """Save pipeline state to checkpoint with atomic write."""
        self.state_file.parent.mkdir(parents=True, exist_ok=True)

        # Update checkpoint time
        self.state.last_checkpoint_time = time.time()

        # Create backup of existing checkpoint
        backup_file = self.state_file.with_suffix(".checkpoint.bak")
        if self.state_file.exists():
            try:
                shutil.copy2(self.state_file, backup_file)
            except OSError as e:
                self.logger.warning(f"Could not create checkpoint backup: {e}")

        # Save new checkpoint atomically using temp file + rename
        temp_file = self.state_file.with_suffix(".checkpoint.tmp")
        try:
            checkpoint_data = asdict(self.state)
            # Add timestamp for debugging
            checkpoint_data["_saved_at"] = datetime.now().isoformat()

            with open(temp_file, "w") as f:
                json.dump(checkpoint_data, f, indent=2, default=str)
                f.flush()
                os.fsync(f.fileno())  # Ensure data is written to disk

            # Atomic rename (on POSIX systems)
            temp_file.replace(self.state_file)
            self.logger.debug(f"Checkpoint saved to {self.state_file}")

        except Exception as e:
            self.logger.error(f"Failed to save checkpoint: {e}")
            # Clean up temp file if it exists
            if temp_file.exists():
                try:
                    temp_file.unlink()
                except OSError:
                    pass
            # Try to restore from backup
            if backup_file.exists():
                try:
                    shutil.copy2(backup_file, self.state_file)
                    self.logger.info("Restored checkpoint from backup")
                except OSError:
                    self.logger.error("Could not restore from backup")
            raise PipelineError(f"Failed to save checkpoint: {e}")

    def _finalize_outputs(self) -> None:
        """Finalize outputs: only checkpoint and temp cleanup; packaging handled by ecc_packager."""
        # Clean up checkpoint/config artifacts unless keep_tmp is set
        if not self.config.keep_tmp:
            for cleanup_path in [
                self.state_file,
                self.state_file.with_suffix(self.state_file.suffix + ".bak"),
                self.final_output_dir / "config.yaml",
            ]:
                try:
                    if cleanup_path and cleanup_path.exists():
                        cleanup_path.unlink()
                        self.logger.debug(f"Removed auxiliary file: {cleanup_path}")
                except OSError as exc:
                    self.logger.warning(f"Could not remove auxiliary file {cleanup_path}: {exc}")

            # Handle turbo mode cleanup
            if self._turbo_mode_active and self._turbo_shm_path and self._symlink_path:
                self._cleanup_turbo_mode(keep_files=False)
            elif self._temp_dir_safe_to_delete:
                self.logger.info(f"Cleaning up temporary directory: {self.temp_dir}")
                try:
                    shutil.rmtree(self.temp_dir)
                except OSError as e:
                    self.logger.warning(f"Could not completely remove temp directory: {e}")
            else:
                # Absolute path - don't auto-delete for safety
                self.logger.info(
                    f"Temporary files in absolute path retained for safety: {self.temp_dir}\n"
                    f"Please manually remove if no longer needed."
                )
        else:
            # keep_tmp is True
            if self._turbo_mode_active and self._turbo_shm_path and self._symlink_path:
                self._cleanup_turbo_mode(keep_files=True)
            else:
                self.logger.info(
                    f"Temporary files retained in: {self.temp_dir}, "
                    f"checkpoint retained at: {self.state_file}"
                )

    def _cleanup_turbo_mode(self, keep_files: bool) -> None:
        """Cleanup turbo mode: move or delete files from /dev/shm.

        Args:
            keep_files: If True, move files to output_dir/tmp/; if False, delete everything
        """
        if not self._turbo_shm_path or not self._symlink_path:
            return

        if keep_files:
            # Move files from /dev/shm to output_dir/tmp/
            dest_dir = self.final_output_dir / "tmp"
            self.logger.info(
                f"Turbo mode: moving temp files from {self._turbo_shm_path} to {dest_dir}"
            )

            try:
                # Remove symlink first
                if self._symlink_path.is_symlink():
                    self._symlink_path.unlink()

                # Move directory contents
                if self._turbo_shm_path.exists():
                    shutil.move(str(self._turbo_shm_path), str(dest_dir))
                    self.logger.info(f"Temporary files retained in: {dest_dir}")
            except OSError as e:
                self.logger.warning(f"Could not move turbo temp files: {e}")
                self.logger.info(
                    f"Temporary files may remain in: {self._turbo_shm_path}"
                )
        else:
            # Delete everything
            self.logger.info(f"Turbo mode: cleaning up {self._turbo_shm_path}")

            try:
                # Remove symlink
                if self._symlink_path.is_symlink():
                    self._symlink_path.unlink()

                # Remove /dev/shm directory
                if self._turbo_shm_path.exists():
                    shutil.rmtree(self._turbo_shm_path)
                    self.logger.info("Turbo mode temp files cleaned up successfully")
            except OSError as e:
                self.logger.warning(
                    f"Could not completely clean up turbo temp files: {e}\n"
                    f"Please manually remove: {self._turbo_shm_path}"
                )

    def show_steps(self, detailed: bool = False) -> None:
        """Show pipeline steps with enhanced status information."""
        import click

        click.echo("CircleSeeker Pipeline Steps:")
        click.echo("=" * 70)

        step_groups = {
            1: "CtcReads-Caller",
            11: "SplitReads-Caller",
            14: "Integration",
        }

        name_width = max(
            (len((s.display_name or s.name)) for s in self.STEPS),
            default=0,
        )
        name_width = max(name_width, 16)
        idx_width = len(str(len(self.STEPS)))

        for i, step in enumerate(self.STEPS, 1):
            group = step_groups.get(i)
            if group:
                click.echo(f"{group}:")
            step_label = step.display_name or step.name
            if step.name in self.state.completed_steps:
                status = "✓"
            elif step.name == self.state.current_step:
                status = "◉"
            elif step.name == self.state.failed_step:
                status = "✗"
            else:
                status = "○"

            # Simplified output: aligned columns, no skippable info
            click.echo(
                f"{status} Step {i:{idx_width}d}: {step_label:<{name_width}} - {step.description}"
            )

            # Detailed info (optional)
            if detailed and step.name in self.state.step_metadata:
                meta = self.state.step_metadata[step.name]
                if meta.duration:
                    click.echo(f"    Duration: {meta.duration:.1f}s")
                if meta.error_message:
                    click.echo(f"    Error: {meta.error_message}")
                if meta.output_files:
                    click.echo(f"    Output files: {len(meta.output_files)} files")

        click.echo("=" * 70)

        # Show runtime summary (only if detailed)
        if detailed and self.state.pipeline_start_time:
            total_time = self.state.get_total_runtime()
            if total_time:
                click.echo(f"Total runtime: {total_time:.1f}s")

        if detailed and self.state.failed_step:
            click.echo(f"⚠️  Previous run failed at: {self.state.failed_step}")
            click.echo(
                "   Use --force to restart from beginning, "
                "or resume will continue from failure point"
            )

    def show_checkpoint_info(self) -> None:
        """Show detailed checkpoint information for debugging."""
        import click
        from datetime import datetime

        click.echo("\nCheckpoint Information")
        click.echo("=" * 50)

        if not self.state_file.exists():
            click.echo(f"No checkpoint found at: {self.state_file}")
            click.echo("Run the pipeline first to create a checkpoint.")
            return

        click.echo(f"Checkpoint file: {self.state_file}")
        click.echo(f"Output directory: {self.final_output_dir}")
        click.echo(f"Prefix: {self.config.prefix}")
        click.echo("-" * 50)

        # Show completion status
        completed = len(self.state.completed_steps)
        total = len(self.STEPS)
        click.echo(f"Progress: {completed}/{total} steps completed")

        if self.state.completed_steps:
            click.echo(f"Completed: {', '.join(self.state.completed_steps)}")

        if self.state.current_step:
            click.echo(f"Current step: {self.state.current_step}")

        if self.state.failed_step:
            click.echo(f"Failed at: {self.state.failed_step}")

        if self.state.last_checkpoint_time:
            ts = datetime.fromtimestamp(self.state.last_checkpoint_time)
            click.echo(f"Last saved: {ts.isoformat()}")

        if self.state.pipeline_start_time:
            runtime = self.state.get_total_runtime()
            if runtime:
                click.echo(f"Total runtime: {runtime:.1f}s")

        click.echo("=" * 50)

    def run(
        self, start_from: Optional[int] = None, stop_at: Optional[int] = None, force: bool = False
    ) -> dict[str, Any]:
        """Run the pipeline with automatic cleanup and temp directory management."""
        if force:
            self.state = self._create_fresh_state()
            self._save_state()

        try:
            total_steps = len(self.STEPS)
            if start_from is not None and not (1 <= start_from <= total_steps):
                raise PipelineError(
                    f"start_from must be within 1..{total_steps}, got {start_from}"
                )
            if stop_at is not None and not (1 <= stop_at <= total_steps):
                raise PipelineError(f"stop_at must be within 1..{total_steps}, got {stop_at}")
            if start_from is not None and stop_at is not None and start_from > stop_at:
                raise PipelineError(
                    f"start_from ({start_from}) cannot be greater than stop_at ({stop_at})"
                )

            # Determine step range
            start_idx = (start_from - 1) if start_from else 0
            end_idx = stop_at if stop_at else len(self.STEPS)

            steps_slice = self.STEPS[start_idx:end_idx]
            iterator = iter_progress(
                steps_slice,
                total=len(steps_slice),
                desc="Process",
                enabled=getattr(self.config.runtime, "enable_progress", True),
            )
            for step_index, step in enumerate(iterator, start_idx):
                step_number = step_index + 1  # display is 1-based (matches --start-from/--stop-at)
                step_label = step.display_name or step.name
                if step.name in self.state.completed_steps and not force:
                    self.logger.info(f"Skipping completed step {step_number}: {step_label}")
                    continue

                # Check skip condition
                if step.skip_condition and getattr(self.config, step.skip_condition, False):
                    self.logger.info(f"Skipping step {step_number}: {step_label} (user requested)")
                    continue

                self.logger.info(f"Running step {step_number}: {step_label}")
                self.state.current_step = step.name
                self.state.add_step_metadata(step.name)
                self._save_state()

                # Execute step
                step_start_time = time.time()
                try:
                    output_files = self._execute_step(step)

                    # Mark as completed
                    self.state.complete_step(step.name, output_files or [])
                    self.state.current_step = None
                    self.state.failed_step = None  # Clear any previous failure

                    step_duration = time.time() - step_start_time
                    self.logger.info(
                        f"Completed step {step_number}: {step_label} ({step_duration:.1f}s)"
                    )

                except Exception as e:
                    # Handle step failure
                    error_msg = str(e)
                    self.state.fail_step(step.name, error_msg)
                    self.state.current_step = None

                    step_duration = time.time() - step_start_time
                    self.logger.error(f"Step {step_number} failed: {step_label} ({step_duration:.1f}s)")
                    self.logger.error(f"Error: {error_msg}", exc_info=True)

                    # Save state with failure info
                    self._save_state()

                    # Re-raise the exception
                    raise PipelineError(f"Pipeline failed at step {step.name}: {error_msg}") from e

                # Save successful completion
                self._save_state()

            # After successful completion, finalize outputs unless skipped
            if not getattr(self.config, "skip_organize", False):
                self._finalize_outputs()
                self.logger.info("Outputs organized to final directories")
            else:
                self.logger.info("Skipping output organization (--skip-organize specified)")

            self.logger.info("Pipeline completed successfully")
            return self.state.results

        except PipelineError:
            # Keep temp directory for debugging on failure
            self.logger.info(f"Pipeline failed. Temporary files retained in: {self.temp_dir}")
            raise
        except Exception as e:
            if self.state.current_step:
                self.state.fail_step(self.state.current_step, str(e))
            self._save_state()
            self.logger.info(f"Unexpected failure. Temporary files retained in: {self.temp_dir}")
            raise PipelineError(f"Unexpected pipeline failure: {e}") from e

    def _execute_step(self, step: PipelineStep) -> Optional[list[str]]:
        """Execute a pipeline step and return output files."""
        method_name = f"_step_{step.name}"
        if hasattr(self, method_name):
            method = getattr(self, method_name)
            from circleseeker.core.steps.contracts import get_step_contract
            from circleseeker.core.steps.schema import validate_step_contract

            contract = get_step_contract(step.name)
            if contract is not None:
                validate_step_contract(self, contract, when="pre")
            result = method()
            if contract is not None:
                validate_step_contract(self, contract, when="post")

            # If method returns a list of output files, use it
            if isinstance(result, list):
                return result

            return None
        else:
            raise PipelineError(f"Step implementation not found: {method_name}")

    # ===================== STEP IMPLEMENTATIONS =====================

    def _step_check_dependencies(self) -> None:
        """Step 0: Check all required tools and dependencies."""
        from circleseeker.core.steps.preprocess import check_dependencies

        check_dependencies(self)

    def _step_tidehunter(self) -> None:
        """Step 1: Run TideHunter."""
        from circleseeker.core.steps.preprocess import tidehunter

        tidehunter(self)

    def _step_tandem_to_ring(self) -> None:
        """Step 2: Process TideHunter output using tandem_to_ring module."""
        from circleseeker.core.steps.preprocess import tandem_to_ring

        tandem_to_ring(self)

    def _step_run_alignment(self) -> None:
        """Step 3: Run minimap2 alignment."""
        from circleseeker.core.steps.preprocess import run_alignment

        run_alignment(self)

    def _step_um_classify(self) -> None:
        """Step 4: Classify eccDNA types using um_classify module."""
        from circleseeker.core.steps.umc import um_classify

        um_classify(self)

    def _step_cecc_build(self) -> None:
        """Step 6: Process Cecc candidates using cecc_build module."""
        from circleseeker.core.steps.umc import cecc_build

        cecc_build(self)

    def _step_umc_process(self) -> None:
        """Step 7: Generate FASTA files using umc_process module."""
        from circleseeker.core.steps.umc import umc_process

        umc_process(self)

    def _step_cd_hit(self) -> None:
        """Step 8: Remove redundant sequences using CD-HIT-EST wrapper."""
        from circleseeker.core.steps.umc import cd_hit

        cd_hit(self)

    def _step_ecc_dedup(self) -> None:
        """Step 9: Coordinate and deduplicate results using ecc_dedup module."""
        from circleseeker.core.steps.umc import ecc_dedup

        ecc_dedup(self)

    def _step_read_filter(self) -> None:
        """Step 10: Filter sequences using read_filter module.

        Based on TandemToRing classification.
        """
        from circleseeker.core.steps.inference import read_filter

        read_filter(self)

    def _step_minimap2(self) -> None:
        """Step 11: Prepare mapping artifacts (alignment optional for SplitReads-Core)."""
        from circleseeker.core.steps.inference import minimap2

        minimap2(self)

    def _step_ecc_inference(self) -> None:
        """Step 12: Circular DNA detection using SplitReads-Core (built-in, HiFi optimized)."""
        from circleseeker.core.steps.inference import ecc_inference

        ecc_inference(self)

    def _run_splitreads_core_inference(self) -> None:
        """Run SplitReads-Core inference pipeline (built-in, HiFi optimized)."""
        from circleseeker.core.steps.inference import run_splitreads_core_inference

        run_splitreads_core_inference(self)

    def _step_curate_inferred_ecc(self) -> list[str] | None:
        """Step 12: Curate inferred eccDNA simple/chimeric tables from SplitReads-Core output."""
        from circleseeker.core.steps.inference import curate_inferred_ecc

        return curate_inferred_ecc(self)

    def _step_ecc_unify(self) -> None:
        """Step 13: Merge eccDNA tables into unified output."""
        from circleseeker.core.steps.postprocess import ecc_unify

        ecc_unify(self)

    def _step_ecc_summary(self) -> None:
        """Step 14: Generate summary report."""
        from circleseeker.core.steps.postprocess import ecc_summary

        ecc_summary(self)

    def _step_ecc_packager(self) -> None:
        """Step 15: Package output files using ecc_packager module."""
        from circleseeker.core.steps.postprocess import ecc_packager

        ecc_packager(self)
