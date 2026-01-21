"""Configuration management for CircleSeeker."""

from __future__ import annotations

import os
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional, Any, cast, Literal, Union, Callable, TYPE_CHECKING
import yaml
from circleseeker.exceptions import ConfigurationError


def _default_threads() -> int:
    """Calculate default threads: min(8, cpu_count * 2) to avoid exceeding resources."""
    cpu_count = os.cpu_count() or 4
    return min(8, cpu_count * 2)


# ================== Preset Definitions ==================
# Three sensitivity presets: relaxed (high recall), balanced (default), strict (high precision)

PresetName = Literal["relaxed", "balanced", "strict"]

PRESETS: dict[PresetName, dict[str, dict[str, Any]]] = {
    # ------------------------------------------------
    # RELAXED: High recall, may have more false positives
    # Use when: exploratory analysis, don't want to miss anything
    # ------------------------------------------------
    "relaxed": {
        "tidehunter": {
            # k=16 unchanged (lower k changes output format)
            "e": 0.15,      # Allow more errors (default: 0.1)
            "c": 2,         # Min copy number (TideHunter requires >=2)
            "p": 30,        # Detect shorter periods (default: 100, TideHunter min: 30)
        },
        "tandem_to_ring": {
            "min_ave_match": 98.0,
        },
        "minimap2_align": {
            "min_identity": 97.0,           # Lower identity threshold
            "identity_decay_per_10kb": 1.0, # More lenient for long reads
            "min_identity_floor": 95.0,     # Lower floor
        },
        "um_classify": {
            "theta_full": 0.90,         # Lower coverage requirement
            "theta_u": 0.90,
            "theta_m": 0.90,
            "theta_u2_max": 0.10,       # Allow more secondary coverage
            "theta_locus": 0.90,
            "u_secondary_max_ratio": 0.10,
        },
        "cecc_build": {
            "overlap_threshold": 0.90,
            "theta_chain": 0.90,
            "min_match_degree": 90.0,
            "edge_tolerance": 50,       # More tolerant of edge mismatches
        },
    },

    # ------------------------------------------------
    # BALANCED: Default settings, good precision-recall trade-off
    # Use when: standard analysis
    # ------------------------------------------------
    "balanced": {
        # Empty dict = use default values from ToolConfig
    },

    # ------------------------------------------------
    # STRICT: High precision, may miss some true positives
    # Use when: validation studies, publication-quality results
    # ------------------------------------------------
    "strict": {
        "tidehunter": {
            # k=16 unchanged (default)
            "e": 0.05,      # Require fewer errors
            "c": 3,         # Require more copies (min copy number)
            "p": 150,       # Require longer periods
        },
        "tandem_to_ring": {
            "min_ave_match": 99.5,
        },
        "minimap2_align": {
            "min_identity": 99.5,           # Higher identity threshold
            "identity_decay_per_10kb": 0.25,# Less decay allowed
            "min_identity_floor": 98.0,     # Higher floor
        },
        "um_classify": {
            "theta_full": 0.98,         # Higher coverage requirement
            "theta_u": 0.98,
            "theta_m": 0.98,
            "theta_u2_max": 0.02,       # Stricter secondary filtering
            "mapq_u_min": 20,           # Require decent MAPQ
            "theta_locus": 0.98,
            "u_secondary_max_ratio": 0.02,
        },
        "cecc_build": {
            "overlap_threshold": 0.98,
            "min_segments": 3,          # Require more supporting segments
            "theta_chain": 0.98,
            "min_match_degree": 98.0,
            "edge_tolerance": 10,       # Less tolerant
        },
    },
}


def apply_preset(cfg: "Config", preset: PresetName) -> None:
    """Apply a preset to the configuration.

    Args:
        cfg: Configuration object to modify in-place
        preset: One of 'relaxed', 'balanced', 'strict'

    Raises:
        ConfigurationError: If preset name is invalid
    """
    if preset not in PRESETS:
        valid = ", ".join(PRESETS.keys())
        raise ConfigurationError(f"Invalid preset '{preset}'. Valid options: {valid}")

    preset_values = PRESETS[preset]

    # Apply each tool's preset values
    for tool_name, params in preset_values.items():
        if hasattr(cfg.tools, tool_name):
            tool_config = getattr(cfg.tools, tool_name)
            # Support both typed config classes and dicts
            if isinstance(tool_config, dict):
                tool_config.update(params)
            elif hasattr(tool_config, '__dataclass_fields__'):
                # Update dataclass fields directly
                for key, value in params.items():
                    if hasattr(tool_config, key):
                        setattr(tool_config, key, value)


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
    # Turbo mode: use /dev/shm for temporary files (RAM-backed, faster I/O)
    turbo_mode: bool = False
    # Minimum required space in /dev/shm for turbo mode (in GB)
    turbo_min_space_gb: float = 10.0


@dataclass
class PerformanceConfig:
    """Performance-related configuration."""

    threads: int = field(default_factory=_default_threads)


# ================== Typed Tool Configuration Classes ==================


class ToolConfigMixin:
    """Mixin class providing common methods for tool configuration dataclasses.

    All tool config classes inherit from this mixin to provide:
    - to_dict(): Convert to dictionary
    - from_dict(): Create from dictionary, ignoring unknown keys
    - get(): Dict-like get method
    - __getitem__(): Dict-like access
    """

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary."""
        return asdict(self)  # type: ignore[call-overload, no-any-return]

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "ToolConfigMixin":
        """Create from dictionary, ignoring unknown keys."""
        valid_keys = {f.name for f in cls.__dataclass_fields__.values()}  # type: ignore[attr-defined]
        filtered = {k: v for k, v in data.items() if k in valid_keys}
        return cls(**filtered)

    def get(self, key: str, default: Any = None) -> Any:
        """Dict-like get method."""
        return getattr(self, key, default)

    def __getitem__(self, key: str) -> Any:
        """Dict-like access."""
        return getattr(self, key)


@dataclass
class TideHunterConfig(ToolConfigMixin):
    """TideHunter tool configuration with type safety."""

    k: int = 16  # k-mer size
    w: int = 1  # window size
    p: int = 100  # minimum period size
    P: int = 2000000  # maximum period size
    e: float = 0.1  # maximum allowed error rate
    f: int = 2  # output format
    c: int = 2  # minimum copy number


@dataclass
class TandemToRingConfig(ToolConfigMixin):
    """TandemToRing tool configuration with type safety."""

    min_ave_match: float = 99.0  # Minimum average match percentage


@dataclass
class UMClassifyConfig(ToolConfigMixin):
    """UMClassify tool configuration with type safety."""

    # Gap threshold for alignment analysis
    gap_threshold: float = 10.0
    # Coverage model parameters (fractions 0-1)
    theta_full: float = 0.95
    theta_u: float = 0.95
    theta_m: float = 0.95
    theta_u2_max: float = 0.05
    theta_locus: float = 0.95
    # MAPQ thresholds
    mapq_u_min: int = 0
    mapq_m_ambiguous_threshold: int = 0
    mecc_identity_gap_threshold: float = 0.0
    # Secondary evidence thresholds
    u_secondary_min_frac: float = 0.01
    u_secondary_min_bp: int = 50
    u_contig_gap_bp: int = 1000
    u_secondary_max_ratio: float = 0.05
    u_high_coverage_threshold: float = 0.98
    u_high_mapq_threshold: int = 50
    pos_tol_bp: int = 50
    # Ambiguity interception
    delta_uc: float = 0.05
    epsilon_mc: float = 0.05
    # Legacy keys (for backward compatibility)
    min_full_length_coverage: float = 95.0
    max_identity_gap_for_mecc: float = 5.0


@dataclass
class CeccBuildConfig(ToolConfigMixin):
    """CeccBuild tool configuration with type safety."""

    overlap_threshold: float = 0.95
    min_segments: int = 2
    edge_tolerance: int = 20  # Gap tolerance on query (bp)
    tau_gap: int = 20
    position_tolerance: int = 50  # Position tolerance for closure checks (bp)
    half_query_buffer: int = 50  # Buffer for doubled-sequence detection (bp)
    locus_overlap_threshold: float = 0.95
    theta_chain: float = 0.95  # Chain coverage threshold (fraction 0-1)
    min_match_degree: float = 95.0  # Legacy key (percent 0-100)
    max_rotations: int = 20


@dataclass
class AlignmentConfig(ToolConfigMixin):
    """General alignment configuration (supports multiple aligners)."""

    aligner: str = "minimap2"  # "minimap2" or "last"
    min_identity: float = 99.0
    min_alignment_length: int = 50
    db_prefix: Optional[str] = None  # LAST-specific: pre-built database


@dataclass
class Minimap2AlignConfig(ToolConfigMixin):
    """Minimap2-specific alignment configuration for run_alignment step."""

    preset: str = "map-hifi"
    max_target_seqs: int = 5
    additional_args: str = ""
    min_identity: float = 99.0  # Base identity threshold (%)
    identity_decay_per_10kb: float = 0.5  # Length-based compensation (%)
    min_identity_floor: float = 97.0  # Identity threshold floor (%)
    split_by_length: bool = False  # Use different presets for short/long
    split_length: int = 5000  # Length threshold in bp
    preset_short: str = "map-hifi"  # Preset for sequences < split_length
    preset_long: str = "map-hifi"  # Preset for sequences >= split_length


@dataclass
class Minimap2Config(ToolConfigMixin):
    """Minimap2 tool configuration for inference step."""

    preset: str = "map-hifi"
    additional_args: str = ""


@dataclass
class SamtoolsConfig(ToolConfigMixin):
    """Samtools tool configuration."""

    # Currently no specific parameters, but structured for future use
    pass


# Type alias for tool config that can be either typed class or dict
ToolConfigValue = Union[
    TideHunterConfig, TandemToRingConfig, UMClassifyConfig, CeccBuildConfig,
    AlignmentConfig, Minimap2AlignConfig, Minimap2Config, SamtoolsConfig,
    dict[str, Any]
]


def _ensure_tool_config(
    value: Optional[ToolConfigValue],
    config_class: type,
    default_factory: Callable[[], Any],
) -> Any:
    """Convert dict to typed config class, or return default if None."""
    if value is None:
        return default_factory()
    if isinstance(value, config_class):
        return value
    if isinstance(value, dict):
        return config_class.from_dict(value)  # type: ignore[attr-defined]
    return value


@dataclass
class ToolConfig:
    """External tool configuration with typed sub-configs."""

    tidehunter: Union[TideHunterConfig, dict[str, Any]] = field(
        default_factory=TideHunterConfig
    )
    tandem_to_ring: Union[TandemToRingConfig, dict[str, Any]] = field(
        default_factory=TandemToRingConfig
    )
    um_classify: Union[UMClassifyConfig, dict[str, Any]] = field(
        default_factory=UMClassifyConfig
    )
    cecc_build: Union[CeccBuildConfig, dict[str, Any]] = field(
        default_factory=CeccBuildConfig
    )
    alignment: Union[AlignmentConfig, dict[str, Any]] = field(
        default_factory=AlignmentConfig
    )
    minimap2_align: Union[Minimap2AlignConfig, dict[str, Any]] = field(
        default_factory=Minimap2AlignConfig
    )
    minimap2: Union[Minimap2Config, dict[str, Any]] = field(
        default_factory=Minimap2Config
    )
    samtools: Union[SamtoolsConfig, dict[str, Any]] = field(
        default_factory=SamtoolsConfig
    )

    def __post_init__(self) -> None:
        """Convert dict values to typed configs after initialization."""
        self.tidehunter = _ensure_tool_config(
            self.tidehunter, TideHunterConfig, TideHunterConfig
        )
        self.tandem_to_ring = _ensure_tool_config(
            self.tandem_to_ring, TandemToRingConfig, TandemToRingConfig
        )
        self.um_classify = _ensure_tool_config(
            self.um_classify, UMClassifyConfig, UMClassifyConfig
        )
        self.cecc_build = _ensure_tool_config(
            self.cecc_build, CeccBuildConfig, CeccBuildConfig
        )
        self.alignment = _ensure_tool_config(
            self.alignment, AlignmentConfig, AlignmentConfig
        )
        self.minimap2_align = _ensure_tool_config(
            self.minimap2_align, Minimap2AlignConfig, Minimap2AlignConfig
        )
        self.minimap2 = _ensure_tool_config(
            self.minimap2, Minimap2Config, Minimap2Config
        )
        self.samtools = _ensure_tool_config(
            self.samtools, SamtoolsConfig, SamtoolsConfig
        )


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
    def threads(self, value: int) -> None:
        self.performance.threads = value

    @property
    def keep_tmp(self) -> bool:
        return self.runtime.keep_tmp

    @keep_tmp.setter
    def keep_tmp(self, value: bool) -> None:
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

        # Validate threads upper bound to prevent resource exhaustion
        import os
        max_threads = (os.cpu_count() or 8) * 2
        if self.performance.threads > max_threads:
            raise ConfigurationError(
                f"Threads must be <= {max_threads} (2x CPU cores), got {self.performance.threads}"
            )

        # Validate tool configuration parameters
        self._validate_tool_params()

        # Validate runtime tmp_dir safety.
        # - Relative tmp_dir must be a subdirectory name/path (not '.', not escaping via '..')
        # - If you want temp files outside the output dir, use an absolute path.
        tmp_dir = self.runtime.tmp_dir
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

    def _validate_tool_params(self) -> None:
        """Validate tool configuration parameter ranges."""
        # Helper for range validation
        def check_range(
            value: float, name: str, min_val: float = 0.0, max_val: float = 1.0
        ) -> None:
            if not min_val <= value <= max_val:
                raise ConfigurationError(
                    f"{name} must be in [{min_val}, {max_val}], got {value}"
                )

        def check_positive(value: float, name: str) -> None:
            if value < 0:
                raise ConfigurationError(f"{name} must be non-negative, got {value}")

        # Validate UMClassify thresholds (fractions in [0, 1])
        # Note: __post_init__ ensures these are typed configs, not dicts
        um = cast(UMClassifyConfig, self.tools.um_classify)
        check_range(um.theta_full, "um_classify.theta_full")
        check_range(um.theta_u, "um_classify.theta_u")
        check_range(um.theta_m, "um_classify.theta_m")
        check_range(um.theta_u2_max, "um_classify.theta_u2_max")
        check_range(um.theta_locus, "um_classify.theta_locus")
        check_range(um.u_secondary_max_ratio, "um_classify.u_secondary_max_ratio")
        check_range(um.u_high_coverage_threshold, "um_classify.u_high_coverage_threshold")
        check_positive(um.gap_threshold, "um_classify.gap_threshold")

        # Validate CeccBuild thresholds
        cecc = cast(CeccBuildConfig, self.tools.cecc_build)
        check_range(cecc.overlap_threshold, "cecc_build.overlap_threshold")
        check_range(cecc.locus_overlap_threshold, "cecc_build.locus_overlap_threshold")
        check_range(cecc.theta_chain, "cecc_build.theta_chain")
        check_range(cecc.min_match_degree, "cecc_build.min_match_degree", 0.0, 100.0)
        if cecc.min_segments < 1:
            raise ConfigurationError(
                f"cecc_build.min_segments must be >= 1, got {cecc.min_segments}"
            )

        # Validate TandemToRing
        ttr = cast(TandemToRingConfig, self.tools.tandem_to_ring)
        check_range(ttr.min_ave_match, "tandem_to_ring.min_ave_match", 0.0, 100.0)

        # Validate Minimap2Align identity thresholds
        mm2 = cast(Minimap2AlignConfig, self.tools.minimap2_align)
        check_range(mm2.min_identity, "minimap2_align.min_identity", 0.0, 100.0)
        check_range(mm2.min_identity_floor, "minimap2_align.min_identity_floor", 0.0, 100.0)
        check_positive(mm2.identity_decay_per_10kb, "minimap2_align.identity_decay_per_10kb")

        # Validate TideHunter
        th = cast(TideHunterConfig, self.tools.tidehunter)
        if th.c < 2:
            raise ConfigurationError(
                f"tidehunter.c (min copy number) must be >= 2, got {th.c}"
            )
        check_range(th.e, "tidehunter.e (error rate)", 0.0, 1.0)

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary."""

        def path_to_str(obj: Any) -> Any:
            if isinstance(obj, Path):
                return str(obj)
            elif isinstance(obj, dict):
                return {k: path_to_str(v) for k, v in obj.items()}
            elif isinstance(obj, (list, tuple)):
                return [path_to_str(item) for item in obj]
            return obj

        converted = path_to_str(asdict(self))
        if not isinstance(converted, dict):
            raise TypeError("Config.to_dict() expected dict")
        return cast(dict[str, Any], converted)

    # ---- New canonical skip flags via properties (backward-compatible) ----
    @property
    def skip_tandem_to_ring(self) -> bool:
        return self.skip_carousel

    @skip_tandem_to_ring.setter
    def skip_tandem_to_ring(self, value: bool) -> None:
        self.skip_carousel = value

def load_config(path: Path) -> Config:
    """Load configuration from YAML file.

    Args:
        path: Path to the YAML configuration file

    Returns:
        Config object populated with values from the file

    Raises:
        ConfigurationError: If the file cannot be read or contains invalid YAML
    """
    try:
        with open(path, "r", encoding="utf-8") as f:
            data = yaml.safe_load(f) or {}
    except yaml.YAMLError as e:
        raise ConfigurationError(f"Invalid YAML syntax in config file {path}: {e}") from e
    except FileNotFoundError:
        raise ConfigurationError(f"Config file not found: {path}") from None
    except PermissionError:
        raise ConfigurationError(f"Permission denied reading config file: {path}") from None
    except OSError as e:
        raise ConfigurationError(f"Error reading config file {path}: {e}") from e

    def build_config(data: dict[str, Any]) -> Config:
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
            "skip_tandem_to_ring",
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

        # Tool config with deep merge support
        if "tools" in data:
            for tool, params in data["tools"].items():
                if hasattr(cfg.tools, tool) and params is not None:
                    existing = getattr(cfg.tools, tool)
                    if isinstance(existing, dict) and isinstance(params, dict):
                        # Deep merge: update existing dict with new values
                        existing.update(params)
                    else:
                        # Non-dict values: replace entirely
                        setattr(cfg.tools, tool, params)

        return cfg

    return build_config(data)


def save_config(cfg: Config, path: Path) -> None:
    """Save configuration to YAML file."""
    data = cfg.to_dict()
    with open(path, "w", encoding="utf-8") as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False)
