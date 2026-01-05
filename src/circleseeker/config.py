"""Configuration management for CircleSeeker."""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional, Dict, Any
import yaml
from circleseeker.exceptions import ConfigurationError


@dataclass
class RuntimeConfig:
    """Runtime configuration."""

    log_level: str = "WARNING"
    log_file: Optional[Path] = None
    tmp_dir: Path = Path(".tmp_work")
    keep_tmp: bool = False
    # Policy when config changes vs. checkpoint: 'continue' | 'reset' | 'fail'
    checkpoint_policy: str = "continue"
    # Enable tqdm progress where available
    enable_progress: bool = True


@dataclass
class PerformanceConfig:
    """Performance-related configuration."""

    threads: int = 8


@dataclass
class ToolConfig:
    """External tool configuration."""

    tidehunter: Dict[str, Any] = field(
        default_factory=lambda: {"k": 16, "w": 1, "p": 100, "P": 2000000, "e": 0.1, "f": 2}
    )
    um_classify: Dict[str, Any] = field(
        default_factory=lambda: {
            "gap_threshold": 10.0,
            # Coverage model parameters (preferred, fractions 0-1)
            "theta_full": 0.95,
            # Split thresholds (preferred): U and M can use different cutoffs.
            # If unset, code falls back to theta_full for both.
            "theta_u": 0.95,
            "theta_m": 0.95,
            # U requires the 2nd-best locus coverage to be low (fractions 0-1).
            # Set to 1.0 to disable.
            "theta_u2_max": 0.05,
            # Optional: require minimap2 MAPQ >= this threshold for Uecc (0 disables).
            "mapq_u_min": 0,
            # When attempting to call U, veto if there is significant secondary evidence
            # (other-chr or far-away mappings) beyond these thresholds.
            "u_secondary_min_frac": 0.01,
            "u_secondary_min_bp": 50,
            "u_contig_gap_bp": 1000,
            "theta_locus": 0.95,
            "pos_tol_bp": 50,
            # Ambiguity interception (fractions 0-1)
            "delta_uc": 0.05,
            "epsilon_mc": 0.05,
            # Legacy keys (accepted for backward compatibility)
            "min_full_length_coverage": 95.0,
            "max_identity_gap_for_mecc": 5.0,
        }
    )
    cecc_build: Dict[str, Any] = field(
        default_factory=lambda: {
            "overlap_threshold": 0.95,
            "min_segments": 2,
            # Gap tolerance on query (bp)
            "edge_tolerance": 20,
            "tau_gap": 20,
            # Position tolerance used by legacy closure checks (bp)
            "position_tolerance": 50,
            # Reciprocal-overlap threshold to treat two genomic intervals as the same locus
            # in the final overlap filter (fractions 0-1; accepts percents for convenience)
            "locus_overlap_threshold": 0.95,
            # Chain coverage threshold (preferred, fraction 0-1)
            "theta_chain": 0.95,
            # Legacy key (percent 0-100)
            "min_match_degree": 95.0,
            "max_rotations": 20,
        }
    )
    minimap2_align: Dict[str, Any] = field(
        default_factory=lambda: {"preset": "sr", "max_target_seqs": 200, "additional_args": ""}
    )
    minimap2: Dict[str, Any] = field(
        default_factory=lambda: {"preset": "map-hifi", "additional_args": ""}
    )
    samtools: Dict[str, Any] = field(default_factory=dict)


@dataclass
class Config:
    """Main configuration class."""

    # Required parameters (set via CLI or config file)
    input_file: Optional[Path] = None
    reference: Optional[Path] = None
    output_dir: Path = Path("circleseeker_output")
    prefix: str = "sample"

    # Feature flags
    enable_xecc: bool = True

    # Step skip flags (canonical new names provided as properties below)
    skip_tidehunter: bool = False
    skip_carousel: bool = False  # alias of skip_tandem_to_ring
    skip_organize: bool = False

    # Sub-configurations
    runtime: RuntimeConfig = field(default_factory=RuntimeConfig)
    performance: PerformanceConfig = field(default_factory=PerformanceConfig)
    tools: ToolConfig = field(default_factory=ToolConfig)

    # Convenience properties
    @property
    def threads(self) -> int:
        return self.performance.threads

    @threads.setter
    def threads(self, value: int):
        self.performance.threads = value

    @property
    def keep_tmp(self) -> bool:
        return self.runtime.keep_tmp

    @keep_tmp.setter
    def keep_tmp(self, value: bool):
        self.runtime.keep_tmp = value

    def validate(self) -> None:
        """Validate configuration."""
        if not self.input_file:
            raise ConfigurationError("Input file is required")
        if not self.reference:
            raise ConfigurationError("Reference genome is required")
        if not self.input_file.exists():
            raise ConfigurationError(f"Input file not found: {self.input_file}")
        if not self.reference.exists():
            raise ConfigurationError(f"Reference file not found: {self.reference}")

        # Validate numeric ranges
        if self.performance.threads < 1:
            raise ConfigurationError("Threads must be >= 1")

        # Validate runtime tmp_dir safety.
        # - Relative tmp_dir must be a subdirectory name/path (not '.', not escaping via '..')
        # - If you want temp files outside the output dir, use an absolute path.
        tmp_dir = self.runtime.tmp_dir
        if not isinstance(tmp_dir, Path):
            tmp_dir = Path(tmp_dir)
        if not tmp_dir.is_absolute():
            if tmp_dir == Path(".") or str(tmp_dir).strip() in {"", "."}:
                raise ConfigurationError(
                    "Invalid runtime.tmp_dir: must be a subdirectory (not '.'). "
                    "Use an absolute path if you want an external temp directory."
                )
            if ".." in tmp_dir.parts:
                raise ConfigurationError(
                    "Invalid runtime.tmp_dir: must not contain '..'. "
                    "Use an absolute path if you want an external temp directory."
                )

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""

        def path_to_str(obj):
            if isinstance(obj, Path):
                return str(obj)
            elif isinstance(obj, dict):
                return {k: path_to_str(v) for k, v in obj.items()}
            elif isinstance(obj, (list, tuple)):
                return [path_to_str(item) for item in obj]
            return obj

        return path_to_str(asdict(self))

    # ---- New canonical skip flags via properties (backward-compatible) ----
    @property
    def skip_tandem_to_ring(self) -> bool:
        return self.skip_carousel

    @skip_tandem_to_ring.setter
    def skip_tandem_to_ring(self, value: bool) -> None:
        self.skip_carousel = value

def load_config(path: Path) -> Config:
    """Load configuration from YAML file."""
    with open(path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}

    def build_config(data: Dict[str, Any]) -> Config:
        cfg = Config()

        # Direct attributes
        if "input_file" in data and data["input_file"] is not None:
            cfg.input_file = Path(data["input_file"])
        if "reference" in data and data["reference"] is not None:
            cfg.reference = Path(data["reference"])
        if "output_dir" in data and data["output_dir"] is not None:
            cfg.output_dir = Path(data["output_dir"])
        if "prefix" in data:
            cfg.prefix = data["prefix"]
        if "enable_xecc" in data:
            cfg.enable_xecc = data["enable_xecc"]
        if "threads" in data and data["threads"] is not None:
            cfg.performance.threads = data["threads"]

        unsupported_skip_flags = [
            "skip_run_alignment",
            "skip_um_classify",
            "skip_alignment",
            "skip_gatekeeper",
        ]
        present_unsupported = [flag for flag in unsupported_skip_flags if flag in data]
        if present_unsupported:
            raise ConfigurationError(
                "Unsupported config option(s): " + ", ".join(present_unsupported)
            )

        # Skip flags
        skip_flags = ["skip_tidehunter", "skip_carousel", "skip_organize"]
        for flag in skip_flags:
            if flag in data:
                setattr(cfg, flag, data[flag])

        # Runtime config
        if "runtime" in data:
            for key, value in data["runtime"].items():
                if hasattr(cfg.runtime, key):
                    if key in ["log_file", "tmp_dir"] and value:
                        value = Path(value)
                    setattr(cfg.runtime, key, value)

        # Performance config
        if "performance" in data:
            for key, value in data["performance"].items():
                if hasattr(cfg.performance, key):
                    setattr(cfg.performance, key, value)

        # Tool config
        if "tools" in data:
            for tool, params in data["tools"].items():
                if hasattr(cfg.tools, tool):
                    if params is None:
                        continue
                    setattr(cfg.tools, tool, params)

        return cfg

    return build_config(data)


def save_config(cfg: Config, path: Path) -> None:
    """Save configuration to YAML file."""
    data = cfg.to_dict()
    with open(path, "w", encoding="utf-8") as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False)
