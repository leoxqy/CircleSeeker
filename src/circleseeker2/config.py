"""Configuration management for CircleSeeker2."""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional, Dict, Any
import yaml
from circleseeker2.exceptions import ConfigurationError


@dataclass
class RuntimeConfig:
    """Runtime configuration."""
    log_level: str = "INFO"
    log_file: Optional[Path] = None
    tmp_dir: Path = Path(".tmp")
    checkpoint_interval: int = 5
    keep_tmp: bool = False


@dataclass
class PerformanceConfig:
    """Performance-related configuration."""
    threads: int = 8
    max_memory: str = "16G"
    chunk_size: int = 10000
    parallel_jobs: int = 4
    stream_buffer_size: int = 65536
    enable_profiling: bool = False


@dataclass
class QualityConfig:
    """Quality control parameters."""
    min_quality_score: float = 0.99
    min_coverage: int = 10
    min_eccdna_size: int = 100
    max_eccdna_size: int = 1000000
    min_alignment_length: int = 100
    min_identity: float = 99.0


@dataclass
class ToolConfig:
    """External tool configuration."""
    tidehunter: Dict[str, Any] = field(default_factory=lambda: {
        "k": 16, "w": 1, "p": 100, "P": 2000000, "e": 0.1, "f": 2
    })
    blast: Dict[str, Any] = field(default_factory=lambda: {
        "word_size": 100, "evalue": "1e-50", "perc_identity": 99.0
    })
    minimap2: Dict[str, Any] = field(default_factory=lambda: {
        "preset": "map-hifi", "additional_args": ""
    })
    samtools: Dict[str, Any] = field(default_factory=dict)
    mosdepth: Dict[str, Any] = field(default_factory=lambda: {
        "min_mapq": 1
    })


@dataclass
class Config:
    """Main configuration class."""
    # Required parameters (set via CLI or config file)
    input_file: Optional[Path] = None
    reference: Optional[Path] = None
    output_dir: Path = Path("circleseeker2_output")
    prefix: str = "sample"
    
    # Feature flags
    enable_xecc: bool = False
    
    # Step skip flags
    skip_make_db: bool = False
    skip_tidehunter: bool = False
    skip_carousel: bool = False
    skip_blast: bool = False
    skip_gatekeeper: bool = False
    skip_report: bool = False
    skip_organize: bool = False
    
    # Sub-configurations
    runtime: RuntimeConfig = field(default_factory=RuntimeConfig)
    performance: PerformanceConfig = field(default_factory=PerformanceConfig)
    quality: QualityConfig = field(default_factory=QualityConfig)
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
        if self.quality.min_eccdna_size < 1:
            raise ConfigurationError("Minimum eccDNA size must be >= 1")
        if self.quality.min_identity < 0 or self.quality.min_identity > 100:
            raise ConfigurationError("Identity must be between 0 and 100")
    
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


def load_config(path: Path) -> Config:
    """Load configuration from YAML file."""
    with open(path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    
    def build_config(data: Dict[str, Any]) -> Config:
        cfg = Config()
        
        # Direct attributes
        if "input_file" in data:
            cfg.input_file = Path(data["input_file"])
        if "reference" in data:
            cfg.reference = Path(data["reference"])
        if "output_dir" in data:
            cfg.output_dir = Path(data["output_dir"])
        if "prefix" in data:
            cfg.prefix = data["prefix"]
        if "enable_xecc" in data:
            cfg.enable_xecc = data["enable_xecc"]
        
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
        
        # Quality config
        if "quality" in data:
            for key, value in data["quality"].items():
                if hasattr(cfg.quality, key):
                    setattr(cfg.quality, key, value)
        
        # Tool config
        if "tools" in data:
            for tool, params in data["tools"].items():
                if hasattr(cfg.tools, tool):
                    setattr(cfg.tools, tool, params)
        
        return cfg
    
    return build_config(data)


def save_config(cfg: Config, path: Path) -> None:
    """Save configuration to YAML file."""
    data = cfg.to_dict()
    with open(path, "w", encoding="utf-8") as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False)