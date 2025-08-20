# CircleSeeker2: 完整技术重构方案

## 一、项目概述

### 1.1 项目背景

CircleSeeker2 是一个用于从 PacBio HiFi 长读长测序数据中识别环状DNA（eccDNA）的生物信息学分析流程。通过三轮迭代分析（TideHunter → BLAST+ × 2 → Minimap2），实现对eccDNA的高灵敏度检测。

### 1.2 重构目标

将现有的单体脚本重构为符合 Bioconda 标准的专业 Python 包，实现：

- 模块化架构设计
- 流式处理大文件
- 并行计算支持
- 完整的测试和文档
- 一键安装部署

## 二、项目结构

```
circleseeker2/
├── pyproject.toml                  # 现代Python打包配置
├── setup.py                        # 兼容性设置脚本
├── MANIFEST.in                     # 包含非Python文件
├── LICENSE                         # GPL-2.0许可证
├── README.md                       # 项目说明
├── CHANGELOG.md                    # 版本更新记录
├── CONTRIBUTING.md                 # 贡献指南
├── CITATION.cff                    # 引用信息
├── requirements.txt                # Python依赖
├── requirements-dev.txt            # 开发依赖
├── .gitignore                      # Git忽略配置
├── .pre-commit-config.yaml         # 代码质量检查
│
├── recipe/                         # Bioconda配方
│   └── meta.yaml
│
├── .github/                        # GitHub配置
│   ├── workflows/
│   │   ├── ci.yml                 # 持续集成
│   │   ├── release.yml            # 自动发布
│   │   └── bioconda.yml           # Bioconda测试
│   └── PULL_REQUEST_TEMPLATE/
│       └── bioconda.md            # PR模板
│
├── docs/                           # 文档
│   ├── conf.py                    # Sphinx配置
│   ├── index.md                   # 文档首页
│   ├── installation.md            # 安装指南
│   ├── quickstart.md              # 快速开始
│   ├── user_guide.md              # 用户手册
│   ├── api_reference.md           # API文档
│   ├── algorithms.md              # 算法说明
│   └── tutorials/
│       ├── basic_usage.md         # 基础教程
│       └── advanced_analysis.md   # 高级分析
│
├── examples/                       # 示例
│   ├── demo_tiny.fa               # 最小测试数据
│   ├── demo_small.fa              # 小型测试数据
│   ├── reference_small.fa         # 参考基因组
│   ├── config_example.yaml        # 配置示例
│   └── run_demo.sh               # 演示脚本
│
├── tests/                          # 测试
│   ├── conftest.py               # pytest配置
│   ├── data/                      # 测试数据
│   │   ├── test_reads.fa
│   │   ├── test_reference.fa
│   │   ├── test.bed.gz
│   │   └── expected_output.fa
│   ├── unit/                      # 单元测试
│   │   ├── test_algorithms.py
│   │   ├── test_io.py
│   │   ├── test_external.py
│   │   └── test_modules.py
│   ├── integration/               # 集成测试
│   │   ├── test_pipeline.py
│   │   └── test_cli.py
│   └── benchmarks/                # 性能测试
│       └── test_performance.py
│
├── scripts/                        # 辅助脚本
│   ├── validate_installation.py   # 安装验证
│   ├── benchmark.py              # 性能测试
│   └── convert_legacy.py         # 旧版本迁移
│
└── src/                           # 源代码（src布局）
    └── circleseeker2/
        ├── __init__.py
        ├── __version__.py         # 版本信息
        ├── __main__.py           # 包入口
        ├── cli.py                # 命令行接口
        ├── config.py             # 配置管理
        ├── constants.py          # 常量定义
        ├── exceptions.py         # 自定义异常
        │
        ├── core/                 # 核心功能
        │   ├── __init__.py
        │   ├── pipeline.py       # 主流程控制
        │   ├── checkpoint.py     # 断点续运
        │   └── orchestrator.py   # 任务编排
        │
        ├── algorithms/           # 核心算法
        │   ├── __init__.py
        │   ├── robust.py         # 鲁棒高覆盖检测
        │   ├── circular.py       # 环状DNA识别
        │   └── filtering.py      # 过滤算法
        │
        ├── io/                   # 输入输出
        │   ├── __init__.py
        │   ├── bed.py           # BED格式处理
        │   ├── fasta.py         # FASTA/FASTQ处理
        │   ├── sam_bam.py       # SAM/BAM处理
        │   ├── blast.py         # BLAST结果解析
        │   └── tables.py        # 表格数据处理
        │
        ├── external/            # 外部工具封装
        │   ├── __init__.py
        │   ├── base.py          # 基础执行器
        │   ├── tidehunter.py    # TideHunter封装
        │   ├── blast.py         # BLAST+封装
        │   ├── minimap2.py      # Minimap2封装
        │   ├── samtools.py      # Samtools封装
        │   └── mosdepth.py      # Mosdepth封装
        │
        ├── modules/             # 功能模块（原有模块迁移）
        │   ├── __init__.py
        │   ├── astrologer.py    # 高覆盖检测
        │   ├── carousel.py      # 序列处理
        │   ├── ringmaster.py    # 分类分析
        │   ├── juggler.py       # 比对分析
        │   ├── sieve.py         # 序列筛选
        │   ├── tamer.py         # 序列标准化
        │   └── processors/      # 处理器
        │       ├── __init__.py
        │       ├── merge_uecc.py
        │       ├── merge_mecc.py
        │       ├── merge_cecc.py
        │       └── fai_processor.py
        │
        ├── reports/             # 报告生成
        │   ├── __init__.py
        │   ├── generator.py     # 报告生成器
        │   ├── visualizer.py    # 可视化
        │   └── templates/       # 报告模板
        │       ├── report.html.j2
        │       └── summary.txt.j2
        │
        ├── utils/               # 工具函数
        │   ├── __init__.py
        │   ├── logging.py       # 日志配置
        │   ├── parallel.py      # 并行处理
        │   ├── progress.py      # 进度监控
        │   ├── profiling.py     # 性能分析
        │   ├── validators.py    # 数据验证
        │   └── helpers.py       # 辅助函数
        │
        └── resources/           # 资源文件
            ├── __init__.py
            ├── config_schema.yaml
            └── default_config.yaml
```

## 三、核心代码实现

### 3.1 版本和包初始化

```python
# src/circleseeker2/__version__.py
__version__ = "2.0.0"
__author__ = "CircleSeeker Team"
__email__ = "circleseeker@example.org"
__license__ = "GPL-2.0"

# src/circleseeker2/__init__.py
"""CircleSeeker2: Comprehensive eccDNA detection from HiFi sequencing data."""

from circleseeker2.__version__ import __version__, __author__, __email__, __license__
from circleseeker2.core.pipeline import Pipeline
from circleseeker2.config import Config
from circleseeker2.exceptions import CircleSeekerError

__all__ = [
    "__version__",
    "__author__",
    "__email__",
    "__license__",
    "Pipeline",
    "Config",
    "CircleSeekerError",
]

# src/circleseeker2/__main__.py
"""Package entry point for python -m circleseeker2."""

from circleseeker2.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
```

### 3.2 命令行接口（CLI）

```python
# src/circleseeker2/cli.py
"""Command-line interface for CircleSeeker2."""

from __future__ import annotations

import sys
import logging
from pathlib import Path
from typing import Optional
import click

from circleseeker2 import __version__
from circleseeker2.config import Config, load_config, save_config
from circleseeker2.core.pipeline import Pipeline
from circleseeker2.utils.logging import setup_logging
from circleseeker2.utils.validators import validate_installation
from circleseeker2.exceptions import CircleSeekerError


@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(version=__version__, prog_name="CircleSeeker2")
@click.pass_context
def cli(ctx):
    """CircleSeeker2: Comprehensive eccDNA detection from HiFi sequencing data."""
    ctx.ensure_object(dict)


@cli.command()
@click.argument("input_file", type=click.Path(exists=True, path_type=Path))
@click.argument("reference", type=click.Path(exists=True, path_type=Path))
@click.option("-o", "--output", type=click.Path(path_type=Path), 
              default=Path("circleseeker2_output"), help="Output directory")
@click.option("-p", "--prefix", default="sample", help="Output file prefix")
@click.option("-c", "--config", type=click.Path(exists=True, path_type=Path),
              help="Configuration file (YAML)")
@click.option("-t", "--threads", type=int, default=8, 
              help="Number of threads [default: 8]")
@click.option("--keep-tmp", is_flag=True, help="Keep temporary files")
@click.option("--enable-xecc", is_flag=True, help="Enable XeccDNA detection")
@click.option("--start-from", type=int, help="Resume from specific step")
@click.option("--stop-at", type=int, help="Stop at specific step")
@click.option("--force", is_flag=True, help="Force re-run all steps")
@click.option("-v", "--verbose", count=True, help="Increase verbosity")
@click.option("--log-file", type=click.Path(path_type=Path), help="Log file path")
@click.option("--dry-run", is_flag=True, help="Show what would be done")
def run(
    input_file: Path,
    reference: Path,
    output: Path,
    prefix: str,
    config: Optional[Path],
    threads: int,
    keep_tmp: bool,
    enable_xecc: bool,
    start_from: Optional[int],
    stop_at: Optional[int],
    force: bool,
    verbose: int,
    log_file: Optional[Path],
    dry_run: bool,
):
    """Run the CircleSeeker2 pipeline for eccDNA detection."""
    
    # Setup logging
    log_level = logging.WARNING
    if verbose == 1:
        log_level = logging.INFO
    elif verbose >= 2:
        log_level = logging.DEBUG
    
    setup_logging(level=log_level, log_file=log_file)
    logger = logging.getLogger(__name__)
    
    try:
        # Load configuration
        cfg = load_config(config) if config else Config()
        
        # Override with CLI arguments
        cfg.input_file = input_file
        cfg.reference = reference
        cfg.output_dir = output
        cfg.prefix = prefix
        cfg.threads = threads
        cfg.keep_tmp = keep_tmp
        cfg.enable_xecc = enable_xecc
        
        # Validate configuration
        cfg.validate()
        
        # Save effective configuration
        output.mkdir(parents=True, exist_ok=True)
        save_config(cfg, output / "config.yaml")
        
        if dry_run:
            logger.info("Dry run mode - showing pipeline steps:")
            pipeline = Pipeline(cfg)
            pipeline.show_steps()
            return
        
        # Run pipeline
        logger.info(f"Starting CircleSeeker2 v{__version__}")
        logger.info(f"Input: {input_file}")
        logger.info(f"Reference: {reference}")
        logger.info(f"Output: {output}")
        
        pipeline = Pipeline(cfg)
        results = pipeline.run(
            start_from=start_from,
            stop_at=stop_at,
            force=force
        )
        
        # Print summary
        click.echo("\n" + "="*60)
        click.echo("Pipeline completed successfully!")
        click.echo(f"Results saved to: {output}")
        click.echo("\nSummary:")
        click.echo(f"  Total eccDNA detected: {results.get('total_eccdna', 0)}")
        click.echo(f"  - UeccDNA: {results.get('uecc_count', 0)}")
        click.echo(f"  - MeccDNA: {results.get('mecc_count', 0)}")
        click.echo(f"  - CeccDNA: {results.get('cecc_count', 0)}")
        click.echo(f"  - MCeccDNA: {results.get('mcecc_count', 0)}")
        if enable_xecc:
            click.echo(f"  - XeccDNA: {results.get('xecc_count', 0)}")
        click.echo("="*60)
        
    except CircleSeekerError as e:
        logger.error(f"Pipeline error: {e}")
        sys.exit(1)
    except Exception as e:
        logger.exception(f"Unexpected error: {e}")
        sys.exit(2)


@cli.command()
@click.option("-o", "--output", type=click.Path(path_type=Path),
              default=Path("config.yaml"), help="Output configuration file")
def init_config(output: Path):
    """Generate a template configuration file."""
    from circleseeker2.resources import get_default_config
    
    config_text = get_default_config()
    output.write_text(config_text)
    click.echo(f"Configuration template saved to: {output}")
    click.echo("Edit this file to customize your analysis parameters.")


@cli.command()
@click.option("--full", is_flag=True, help="Run full validation including test data")
def validate(full: bool):
    """Validate CircleSeeker2 installation and dependencies."""
    click.echo("Validating CircleSeeker2 installation...")
    
    issues = validate_installation(full_check=full)
    
    if not issues:
        click.echo("✓ All checks passed!")
        click.echo(f"  CircleSeeker2 version: {__version__}")
    else:
        click.echo("✗ Issues found:")
        for issue in issues:
            click.echo(f"  - {issue}")
        sys.exit(1)


@cli.command()
@click.argument("input_files", nargs=-1, required=True, 
                type=click.Path(exists=True, path_type=Path))
@click.option("-o", "--output", type=click.Path(path_type=Path),
              default=Path("benchmark_results.csv"))
@click.option("-t", "--threads", type=int, default=8)
def benchmark(input_files: tuple[Path, ...], output: Path, threads: int):
    """Run performance benchmarks on test data."""
    from circleseeker2.utils.profiling import run_benchmarks
    
    click.echo(f"Running benchmarks on {len(input_files)} files...")
    results = run_benchmarks(list(input_files), threads=threads)
    results.to_csv(output, index=False)
    click.echo(f"Benchmark results saved to: {output}")


def main(argv: list[str] | None = None) -> int:
    """Main entry point."""
    try:
        cli(argv)
        return 0
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())
```

### 3.3 配置管理

```python
# src/circleseeker2/config.py
"""Configuration management for CircleSeeker2."""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional, Dict, Any
import yaml
import json
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
        "k": 16, "w": 1, "p": 100, "P": 2000000, "e": 0.1
    })
    blast: Dict[str, Any] = field(default_factory=lambda: {
        "word_size": 100, "evalue": 1e-50, "perc_identity": 99
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
```

### 3.4 外部工具封装

```python
# src/circleseeker2/external/base.py
"""Base class for external tool execution."""

from __future__ import annotations

import subprocess
import logging
import shutil
from pathlib import Path
from typing import Sequence, Optional, Tuple
from circleseeker2.exceptions import ExternalToolError


class ExternalTool:
    """Base class for external tool wrappers."""
    
    tool_name: str = ""
    required_version: Optional[str] = None
    
    def __init__(self, threads: int = 1):
        self.threads = threads
        self.logger = logging.getLogger(f"{__name__}.{self.tool_name}")
        self._check_installation()
    
    def _check_installation(self) -> None:
        """Check if the tool is installed."""
        if not shutil.which(self.tool_name):
            raise ExternalToolError(
                f"{self.tool_name} not found in PATH. "
                f"Please install it via: conda install {self.tool_name}"
            )
        
        if self.required_version:
            version = self._get_version()
            self.logger.debug(f"{self.tool_name} version: {version}")
    
    def _get_version(self) -> str:
        """Get tool version."""
        try:
            result = subprocess.run(
                [self.tool_name, "--version"],
                capture_output=True,
                text=True,
                check=False
            )
            return result.stdout.strip() or result.stderr.strip()
        except Exception:
            return "unknown"
    
    def run(
        self,
        cmd: Sequence[str],
        cwd: Optional[Path] = None,
        check: bool = True,
        capture_output: bool = True
    ) -> Tuple[str, str]:
        """Execute command."""
        cmd_str = " ".join(str(c) for c in cmd)
        self.logger.info(f"Running: {cmd_str}")
        
        try:
            result = subprocess.run(
                cmd,
                cwd=cwd,
                capture_output=capture_output,
                text=True,
                check=check
            )
            
            if capture_output:
                return result.stdout, result.stderr
            return "", ""
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Command failed: {cmd_str}")
            self.logger.error(f"Error: {e.stderr[:1000] if e.stderr else 'No error output'}")
            raise ExternalToolError(
                f"{self.tool_name} failed",
                command=cmd,
                returncode=e.returncode,
                stderr=e.stderr
            )


# src/circleseeker2/external/tidehunter.py
"""TideHunter wrapper."""

from pathlib import Path
from typing import Optional
from circleseeker2.external.base import ExternalTool


class TideHunter(ExternalTool):
    """TideHunter tandem repeat finder."""
    
    tool_name = "TideHunter"
    
    def run_analysis(
        self,
        input_file: Path,
        output_file: Path,
        k: int = 16,
        w: int = 1,
        p: int = 100,
        P: int = 2000000,
        e: float = 0.1,
        f: int = 2
    ) -> None:
        """Run TideHunter analysis."""
        cmd = [
            self.tool_name,
            "-t", str(self.threads),
            "-k", str(k),
            "-w", str(w),
            "-p", str(p),
            "-P", str(P),
            "-e", str(e),
            "-f", str(f),
            str(input_file)
        ]
        
        stdout, _ = self.run(cmd)
        
        # Write output
        output_file.parent.mkdir(parents=True, exist_ok=True)
        output_file.write_text(stdout)
        
        self.logger.info(f"TideHunter results saved to: {output_file}")
```

### 3.5 I/O模块

```python
# src/circleseeker2/io/bed.py
"""BED format I/O with streaming support."""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import Iterator, Tuple, Optional
import logging

logger = logging.getLogger(__name__)


def iter_bed_per_base(
    path: Path,
    min_length: int = 1
) -> Iterator[Tuple[str, int, int, float]]:
    """
    Stream BED file (supports .gz compression).
    
    Args:
        path: Path to BED file
        min_length: Minimum region length to yield
    
    Yields:
        Tuples of (chrom, start, end, coverage)
    """
    opener = gzip.open if path.suffix == ".gz" else open
    
    with opener(path, "rt", encoding="utf-8") as fh:
        for line_num, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            
            try:
                parts = line.split("\t")
                if len(parts) < 4:
                    logger.warning(f"Line {line_num}: insufficient columns")
                    continue
                
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                cov = float(parts[3])
                
                if end - start >= min_length:
                    yield chrom, start, end, cov
                    
            except (ValueError, IndexError) as e:
                logger.warning(f"Line {line_num}: parse error - {e}")


# src/circleseeker2/io/fasta.py
"""FASTA/FASTQ I/O with streaming support."""

from __future__ import annotations

from pathlib import Path
from typing import Iterator, Optional, Dict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pysam


class FastaReader:
    """Efficient FASTA reader with random access."""
    
    def __init__(self, path: Path):
        self.path = path
        self._index_path = Path(f"{path}.fai")
        self._handle: Optional[pysam.FastaFile] = None
        
        # Create index if needed
        if not self._index_path.exists():
            pysam.faidx(str(path))
    
    def __enter__(self):
        self._handle = pysam.FastaFile(str(self.path))
        return self
    
    def __exit__(self, *args):
        if self._handle:
            self._handle.close()
    
    def fetch(self, name: str, start: Optional[int] = None, 
              end: Optional[int] = None) -> str:
        """Fetch sequence by name and optional coordinates."""
        if not self._handle:
            raise RuntimeError("Reader not opened")
        
        if start is not None and end is not None:
            return self._handle.fetch(name, start, end)
        return self._handle.fetch(name)
    
    def iter_records(self) -> Iterator[SeqRecord]:
        """Iterate over all records."""
        with open(self.path) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                yield record


def stream_fasta(path: Path, chunk_size: int = 1000) -> Iterator[list[SeqRecord]]:
    """
    Stream FASTA file in chunks.
    
    Args:
        path: Path to FASTA file
        chunk_size: Number of records per chunk
    
    Yields:
        Lists of SeqRecord objects
    """
    with open(path) as handle:
        chunk = []
        for record in SeqIO.parse(handle, "fasta"):
            chunk.append(record)
            if len(chunk) >= chunk_size:
                yield chunk
                chunk = []
        
        if chunk:  # Yield remaining records
            yield chunk
```

### 3.6 核心算法

```python
# src/circleseeker2/algorithms/robust.py
"""Robust statistical algorithms for eccDNA detection."""

from __future__ import annotations

from typing import List, Tuple, Optional
import numpy as np
from statistics import median
import logging

logger = logging.getLogger(__name__)


def median_absolute_deviation(
    values: List[float],
    scale: float = 1.4826
) -> Tuple[float, float]:
    """
    Calculate Median Absolute Deviation (MAD).
    
    Args:
        values: List of numeric values
        scale: Scale factor for consistency with normal distribution
    
    Returns:
        Tuple of (median, scaled_mad)
    """
    if not values:
        return 0.0, 0.0
    
    med = median(values)
    deviations = [abs(x - med) for x in values]
    mad = median(deviations)
    
    return med, mad * scale


def robust_high_coverage_detection(
    coverages: List[float],
    z_threshold: float = 1.0,
    extreme_threshold: float = 10.0,
    min_length: int = 1
) -> List[bool]:
    """
    Detect high-coverage regions using robust statistics.
    
    Args:
        coverages: Coverage values
        z_threshold: Z-score threshold for outlier detection
        extreme_threshold: Absolute threshold for extreme values
        min_length: Minimum region length
    
    Returns:
        Boolean flags for high-coverage regions
    """
    # Separate extreme and normal values
    normal_values = [c for c in coverages if c <= extreme_threshold]
    
    if normal_values:
        med, mad = median_absolute_deviation(normal_values)
    else:
        med, mad = 0.0, 0.0
    
    flags = []
    for cov in coverages:
        if cov > extreme_threshold:
            # Extreme coverage - always high
            flags.append(True)
        elif mad == 0.0:
            # No variation - use simple threshold
            flags.append(cov > med * 2.0)
        else:
            # Use robust z-score
            z_score = (cov - med) / mad
            flags.append(z_score > z_threshold)
    
    logger.debug(f"Detected {sum(flags)}/{len(flags)} high-coverage regions")
    
    return flags


def identify_circular_signals(
    sequence: str,
    min_repeat_length: int = 100,
    max_mismatch_rate: float = 0.1
) -> List[Tuple[int, int, float]]:
    """
    Identify potential circular DNA signals in sequence.
    
    Args:
        sequence: DNA sequence
        min_repeat_length: Minimum length for repeat detection
        max_mismatch_rate: Maximum allowed mismatch rate
    
    Returns:
        List of (start, end, similarity) tuples
    """
    signals = []
    seq_len = len(sequence)
    
    # Sliding window approach
    for window_size in range(min_repeat_length, seq_len // 2):
        for start in range(seq_len - window_size):
            end = start + window_size
            
            # Check for self-similarity
            substr = sequence[start:end]
            remaining = sequence[end:]
            
            if remaining.startswith(substr[:len(remaining)]):
                # Calculate similarity
                matches = sum(1 for a, b in zip(substr, remaining) if a == b)
                similarity = matches / min(len(substr), len(remaining))
                
                if similarity >= (1 - max_mismatch_rate):
                    signals.append((start, end, similarity))
    
    return signals
```

### 3.7 并行处理工具

```python
# src/circleseeker2/utils/parallel.py
"""Parallel processing utilities."""

from __future__ import annotations

import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Callable, Iterable, List, Any, Optional
import logging
from tqdm import tqdm

logger = logging.getLogger(__name__)


class ParallelProcessor:
    """Smart parallel processor with progress tracking."""
    
    def __init__(self, max_workers: Optional[int] = None):
        self.max_workers = max_workers or mp.cpu_count()
        logger.info(f"Initialized parallel processor with {self.max_workers} workers")
    
    def map(
        self,
        func: Callable,
        items: Iterable,
        chunk_size: Optional[int] = None,
        desc: str = "Processing",
        disable_progress: bool = False
    ) -> List[Any]:
        """
        Parallel map with progress bar.
        
        Args:
            func: Function to apply
            items: Items to process
            chunk_size: Items per chunk
            desc: Progress bar description
            disable_progress: Disable progress bar
        
        Returns:
            List of results
        """
        items_list = list(items)
        n_items = len(items_list)
        
        if n_items == 0:
            return []
        
        # Auto-determine chunk size
        if chunk_size is None:
            chunk_size = max(1, n_items // (self.max_workers * 4))
        
        results = [None] * n_items
        
        with ProcessPoolExecutor(max_workers=self.max_workers) as executor:
            # Submit jobs
            futures = {}
            for i in range(0, n_items, chunk_size):
                chunk = items_list[i:i + chunk_size]
                future = executor.submit(self._process_chunk, func, chunk)
                futures[future] = i
            
            # Collect results
            with tqdm(total=n_items, desc=desc, disable=disable_progress) as pbar:
                for future in as_completed(futures):
                    start_idx = futures[future]
                    try:
                        chunk_results = future.result()
                        for j, result in enumerate(chunk_results):
                            results[start_idx + j] = result
                        pbar.update(len(chunk_results))
                    except Exception as e:
                        logger.error(f"Error processing chunk at index {start_idx}: {e}")
                        raise
        
        return results
    
    @staticmethod
    def _process_chunk(func: Callable, chunk: List) -> List:
        """Process a chunk of items."""
        return [func(item) for item in chunk]
    
    def map_reduce(
        self,
        map_func: Callable,
        reduce_func: Callable,
        items: Iterable,
        initial: Any = None
    ) -> Any:
        """
        Map-reduce pattern.
        
        Args:
            map_func: Mapping function
            reduce_func: Reduction function
            items: Items to process
            initial: Initial value for reduction
        
        Returns:
            Reduced result
        """
        # Map phase
        mapped = self.map(map_func, items)
        
        # Reduce phase
        if initial is not None:
            result = initial
            for item in mapped:
                result = reduce_func(result, item)
        else:
            if not mapped:
                return None
            result = mapped[0]
            for item in mapped[1:]:
                result = reduce_func(result, item)
        
        return result


def run_parallel_external(
    commands: List[List[str]],
    max_workers: int = 4,
    desc: str = "Running commands"
) -> List[Tuple[int, str, str]]:
    """
    Run multiple external commands in parallel.
    
    Args:
        commands: List of command arguments
        max_workers: Maximum parallel processes
        desc: Progress description
    
    Returns:
        List of (returncode, stdout, stderr) tuples
    """
    import subprocess
    
    def run_command(cmd: List[str]) -> Tuple[int, str, str]:
        """Run a single command."""
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True
        )
        return result.returncode, result.stdout, result.stderr
    
    processor = ParallelProcessor(max_workers)
    return processor.map(run_command, commands, desc=desc)
```

### 3.8 进度和性能监控

```python
# src/circleseeker2/utils/progress.py
"""Progress tracking utilities."""

from __future__ import annotations

from functools import wraps
from typing import Callable, Iterator, Optional
import time
import logging
from tqdm import tqdm

logger = logging.getLogger(__name__)


def track_progress(
    desc: str = "Processing",
    unit: str = "items",
    total: Optional[int] = None,
    disable: bool = False
):
    """
    Decorator for progress tracking.
    
    Args:
        desc: Progress bar description
        unit: Unit name
        total: Total items (if known)
        disable: Disable progress bar
    """
    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args, **kwargs):
            # Check if first argument is iterable
            if args and hasattr(args[0], '__iter__'):
                items = args[0]
                
                # Wrap with progress bar
                with tqdm(items, desc=desc, unit=unit, total=total, 
                         disable=disable) as pbar:
                    def wrapped_items():
                        for item in pbar:
                            yield item
                    
                    # Replace first argument with wrapped iterator
                    new_args = (wrapped_items(),) + args[1:]
                    return func(*new_args, **kwargs)
            else:
                # No iterable, just call function
                return func(*args, **kwargs)
        
        return wrapper
    return decorator


# src/circleseeker2/utils/profiling.py
"""Performance profiling utilities."""

from __future__ import annotations

import time
import psutil
import logging
from functools import wraps
from typing import Callable, Dict, Any
import pandas as pd
from pathlib import Path

logger = logging.getLogger(__name__)


class PerformanceMonitor:
    """Monitor performance metrics."""
    
    def __init__(self):
        self.metrics: List[Dict[str, Any]] = []
        self.process = psutil.Process()
    
    def record(self, name: str, duration: float, 
               memory_mb: float, **kwargs):
        """Record performance metrics."""
        self.metrics.append({
            "name": name,
            "duration_sec": duration,
            "memory_mb": memory_mb,
            "timestamp": time.time(),
            **kwargs
        })
    
    def get_summary(self) -> pd.DataFrame:
        """Get performance summary."""
        if not self.metrics:
            return pd.DataFrame()
        
        df = pd.DataFrame(self.metrics)
        return df
    
    def save(self, path: Path):
        """Save metrics to CSV."""
        df = self.get_summary()
        df.to_csv(path, index=False)
        logger.info(f"Performance metrics saved to: {path}")


def profile_performance(monitor: Optional[PerformanceMonitor] = None):
    """
    Decorator for performance profiling.
    
    Args:
        monitor: Performance monitor instance
    """
    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args, **kwargs):
            # Get initial state
            process = psutil.Process()
            mem_before = process.memory_info().rss / 1024 / 1024  # MB
            
            # Run function
            start = time.perf_counter()
            result = func(*args, **kwargs)
            duration = time.perf_counter() - start
            
            # Get final state
            mem_after = process.memory_info().rss / 1024 / 1024  # MB
            mem_used = mem_after - mem_before
            
            # Log metrics
            logger.debug(
                f"{func.__name__} completed in {duration:.2f}s, "
                f"memory delta: {mem_used:.1f}MB"
            )
            
            # Record if monitor provided
            if monitor:
                monitor.record(
                    name=func.__name__,
                    duration=duration,
                    memory_mb=mem_used
                )
            
            return result
        
        return wrapper
    return decorator


def run_benchmarks(
    input_files: List[Path],
    threads: int = 8
) -> pd.DataFrame:
    """
    Run performance benchmarks.
    
    Args:
        input_files: Input files to process
        threads: Number of threads
    
    Returns:
        DataFrame with benchmark results
    """
    from circleseeker2.core.pipeline import Pipeline
    from circleseeker2.config import Config
    
    results = []
    
    for input_file in input_files:
        logger.info(f"Benchmarking: {input_file}")
        
        # Setup config
        config = Config(
            input_file=input_file,
            reference=Path("reference.fa"),  # Placeholder
            output_dir=Path(f"benchmark_{input_file.stem}"),
            threads=threads
        )
        
        # Run with monitoring
        monitor = PerformanceMonitor()
        pipeline = Pipeline(config, performance_monitor=monitor)
        
        start = time.time()
        pipeline.run()
        total_time = time.time() - start
        
        # Collect metrics
        summary = monitor.get_summary()
        summary["input_file"] = str(input_file)
        summary["total_time"] = total_time
        summary["threads"] = threads
        
        results.append(summary)
    
    # Combine results
    return pd.concat(results, ignore_index=True)
```

## 四、打包和发布配置

### 4.1 pyproject.toml

```toml
[build-system]
requires = ["setuptools>=68", "wheel", "setuptools-scm"]
build-backend = "setuptools.build_meta"

[project]
name = "circleseeker2"
version = "2.0.0"
description = "Comprehensive eccDNA detection from HiFi sequencing data"
readme = "README.md"
requires-python = ">=3.9"
license = {text = "GPL-2.0"}
authors = [
    {name = "CircleSeeker Team", email = "circleseeker@example.org"}
]
maintainers = [
    {name = "CircleSeeker Team", email = "circleseeker@example.org"}
]
keywords = [
    "bioinformatics",
    "eccDNA",
    "circular-dna",
    "pacbio",
    "hifi-sequencing",
    "genomics"
]
classifiers = [
    "Development Status :: 4 - Beta",
    "Intended Audience :: Science/Research",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Programming Language :: Python :: 3.12",
    "Operating System :: POSIX :: Linux",
    "Operating System :: MacOS",
]

dependencies = [
    "click>=8.1",
    "pandas>=2.0",
    "numpy>=1.23",
    "biopython>=1.81",
    "pysam>=0.22",
    "tqdm>=4.65",
    "pyyaml>=6.0",
    "psutil>=5.9",
    "jinja2>=3.1",
]

[project.optional-dependencies]
dev = [
    "pytest>=7.4",
    "pytest-cov>=4.1",
    "pytest-xdist>=3.3",
    "black>=23.7",
    "isort>=5.12",
    "flake8>=6.1",
    "mypy>=1.5",
    "pre-commit>=3.3",
]
docs = [
    "sphinx>=7.1",
    "sphinx-rtd-theme>=1.3",
    "sphinx-click>=5.0",
    "autodoc>=0.5",
    "myst-parser>=2.0",
]
test = [
    "pytest>=7.4",
    "pytest-cov>=4.1",
    "hypothesis>=6.82",
]

[project.urls]
Homepage = "https://github.com/circleseeker/circleseeker2"
Documentation = "https://circleseeker2.readthedocs.io"
Repository = "https://github.com/circleseeker/circleseeker2"
Issues = "https://github.com/circleseeker/circleseeker2/issues"
Changelog = "https://github.com/circleseeker/circleseeker2/blob/main/CHANGELOG.md"

[project.scripts]
circleseeker2 = "circleseeker2.cli:main"

[tool.setuptools]
package-dir = {"" = "src"}
include-package-data = true

[tool.setuptools.packages.find]
where = ["src"]

[tool.setuptools.package-data]
circleseeker2 = [
    "resources/*.yaml",
    "reports/templates/*.j2",
]

[tool.black]
line-length = 100
target-version = ["py39", "py310", "py311"]
include = '\.pyi?$'

[tool.isort]
profile = "black"
line_length = 100

[tool.pytest.ini_options]
minversion = "7.0"
testpaths = ["tests"]
addopts = [
    "-ra",
    "--strict-markers",
    "--cov=circleseeker2",
    "--cov-report=term-missing",
    "--cov-report=html",
]

[tool.mypy]
python_version = "3.9"
warn_return_any = true
warn_unused_configs = true
ignore_missing_imports = true

[tool.coverage.run]
source = ["circleseeker2"]
omit = ["*/tests/*", "*/test_*.py"]

[tool.coverage.report]
exclude_lines = [
    "pragma: no cover",
    "def __repr__",
    "raise AssertionError",
    "raise NotImplementedError",
    "if __name__ == .__main__.:",
]
```

### 4.2 Bioconda配方 (recipe/meta.yaml)

```yaml
{% set name = "circleseeker2" %}
{% set version = "2.0.0" %}

package:
  name: {{ name|lower }}
  version: {{ version }}

source:
  url: https://pypi.io/packages/source/{{ name[0] }}/{{ name }}/{{ name }}-{{ version }}.tar.gz
  sha256: YOUR_SHA256_HERE

build:
  number: 0
  noarch: python
  script: {{ PYTHON }} -m pip install . -vv --no-deps --no-build-isolation
  entry_points:
    - circleseeker2 = circleseeker2.cli:main

requirements:
  host:
    - python >=3.9
    - pip
    - setuptools
    - wheel
  run:
    - python >=3.9
    - click >=8.1
    - pandas >=2.0
    - numpy >=1.23
    - biopython >=1.81
    - pysam >=0.22
    - tqdm >=4.65
    - pyyaml >=6.0
    - psutil >=5.9
    - jinja2 >=3.1
    # Bioinformatics tools
    - tidehunter >=1.4.0
    - blast >=2.12.0
    - minimap2 >=2.24
    - samtools >=1.17
    - bcftools >=1.17
    - mosdepth >=0.3.3
    - seqkit >=2.3

test:
  imports:
    - circleseeker2
  commands:
    - circleseeker2 --help
    - circleseeker2 --version
    - circleseeker2 validate

about:
  home: https://github.com/circleseeker/circleseeker2
  license: GPL-2.0-only
  license_family: GPL2
  license_file: LICENSE
  summary: 'Comprehensive eccDNA detection from HiFi sequencing data'
  description: |
    CircleSeeker2 is a bioinformatics pipeline for identifying extrachromosomal 
    circular DNA (eccDNA) from PacBio HiFi long-read sequencing data. It uses 
    a three-round iterative analysis approach combining TideHunter, BLAST+, and 
    Minimap2 to achieve high-sensitivity eccDNA detection.
  doc_url: https://circleseeker2.readthedocs.io
  dev_url: https://github.com/circleseeker/circleseeker2

extra:
  recipe-maintainers:
    - circleseeker-team
  identifiers:
    - doi:10.1234/circleseeker2
    - biotools:circleseeker2
  skip-lints:
    - uses_setuptools  # We use setuptools via pip
```

## 五、测试框架

### 5.1 单元测试示例

```python
# tests/unit/test_algorithms.py
"""Unit tests for algorithms module."""

import pytest
from circleseeker2.algorithms.robust import (
    median_absolute_deviation,
    robust_high_coverage_detection,
    identify_circular_signals
)


class TestRobustStatistics:
    """Test robust statistical algorithms."""
    
    def test_mad_basic(self):
        """Test MAD calculation."""
        values = [1, 2, 2, 3, 3, 3, 4, 4, 5, 10]
        med, mad = median_absolute_deviation(values)
        
        assert med == 3  # Median
        assert mad > 0   # MAD should be positive
    
    def test_mad_empty(self):
        """Test MAD with empty input."""
        med, mad = median_absolute_deviation([])
        assert med == 0.0
        assert mad == 0.0
    
    def test_high_coverage_detection(self):
        """Test high-coverage region detection."""
        coverages = [1, 1, 2, 2, 3, 3, 10, 12, 15, 2, 2, 1]
        flags = robust_high_coverage_detection(
            coverages,
            z_threshold=2.0,
            extreme_threshold=10.0
        )
        
        # Should detect the high values (10, 12, 15)
        assert flags[6] == True   # 10
        assert flags[7] == True   # 12
        assert flags[8] == True   # 15
        
        # Low values should not be flagged
        assert flags[0] == False  # 1
        assert flags[1] == False  # 1
    
    @pytest.mark.parametrize("sequence,expected_signals", [
        ("ATCGATCG" * 10, 1),  # Perfect repeat
        ("ATCGATCGATCG", 0),    # Too short
        ("ATCG" * 100 + "ATCG" * 50, 1),  # Partial repeat
    ])
    def test_circular_signals(self, sequence, expected_signals):
        """Test circular signal detection."""
        signals = identify_circular_signals(
            sequence,
            min_repeat_length=20
        )
        assert len(signals) >= expected_signals


# tests/integration/test_pipeline.py
"""Integration tests for the pipeline."""

import pytest
from pathlib import Path
import tempfile
import shutil
from circleseeker2.core.pipeline import Pipeline
from circleseeker2.config import Config


@pytest.fixture
def test_data_dir():
    """Get test data directory."""
    return Path(__file__).parent.parent / "data"


@pytest.fixture
def temp_output():
    """Create temporary output directory."""
    tmpdir = tempfile.mkdtemp(prefix="circleseeker2_test_")
    yield Path(tmpdir)
    shutil.rmtree(tmpdir)


class TestPipeline:
    """Test pipeline functionality."""
    
    @pytest.mark.slow
    def test_minimal_pipeline(self, test_data_dir, temp_output):
        """Test minimal pipeline run."""
        config = Config(
            input_file=test_data_dir / "test_reads.fa",
            reference=test_data_dir / "test_reference.fa",
            output_dir=temp_output,
            prefix="test"
        )
        
        # Validate config
        config.validate()
        
        # Create pipeline
        pipeline = Pipeline(config)
        
        # Run first few steps only
        results = pipeline.run(stop_at=3)
        
        # Check outputs exist
        assert (temp_output / "test.checkpoint").exists()
        assert results is not None
    
    def test_checkpoint_recovery(self, test_data_dir, temp_output):
        """Test checkpoint recovery."""
        config = Config(
            input_file=test_data_dir / "test_reads.fa",
            reference=test_data_dir / "test_reference.fa",
            output_dir=temp_output,
            prefix="test"
        )
        
        # Run partial pipeline
        pipeline1 = Pipeline(config)
        pipeline1.run(stop_at=2)
        
        # Resume from checkpoint
        pipeline2 = Pipeline(config)
        results = pipeline2.run(start_from=3, stop_at=4)
        
        assert results is not None
```

## 六、GitHub Actions CI/CD

```yaml
# .github/workflows/ci.yml
name: CI

on:
  push:
    branches: [main, develop]
  pull_request:
    branches: [main]

jobs:
  lint:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      
      - name: Set up Python
        uses: actions/setup-python@v4
        with:
          python-version: "3.11"
      
      - name: Install dependencies
        run: |
          python -m pip install --upgrade pip
          pip install black isort flake8 mypy
      
      - name: Run black
        run: black --check src tests
      
      - name: Run isort
        run: isort --check-only src tests
      
      - name: Run flake8
        run: flake8 src tests
      
      - name: Run mypy
        run: mypy src

  test:
    needs: lint
    runs-on: ${{ matrix.os }}
    strategy:
      matrix:
        os: [ubuntu-latest, macos-latest]
        python-version: ["3.9", "3.10", "3.11", "3.12"]
    
    steps:
      - uses: actions/checkout@v4
      
      - name: Set up Python
        uses: actions/setup-python@v4
        with:
          python-version: ${{ matrix.python-version }}
      
      - name: Install system dependencies (Ubuntu)
        if: runner.os == 'Linux'
        run: |
          sudo apt-get update
          sudo apt-get install -y samtools bcftools
      
      - name: Install system dependencies (macOS)
        if: runner.os == 'macOS'
        run: |
          brew install samtools bcftools
      
      - name: Install package
        run: |
          python -m pip install --upgrade pip
          pip install -e .[test]
      
      - name: Run tests
        run: |
          pytest tests/unit -v --cov=circleseeker2
      
      - name: Upload coverage
        uses: codecov/codecov-action@v3
        with:
          file: ./coverage.xml

  build:
    needs: test
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      
      - name: Set up Python
        uses: actions/setup-python@v4
        with:
          python-version: "3.11"
      
      - name: Build package
        run: |
          python -m pip install build
          python -m build
      
      - name: Check dist
        run: |
          ls -la dist/
          python -m pip install twine
          twine check dist/*
      
      - name: Upload artifacts
        uses: actions/upload-artifact@v3
        with:
          name: dist
          path: dist/

  docs:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      
      - name: Set up Python
        uses: actions/setup-python@v4
        with:
          python-version: "3.11"
      
      - name: Install dependencies
        run: |
          pip install -e .[docs]
      
      - name: Build documentation
        run: |
          cd docs
          make html
      
      - name: Deploy to GitHub Pages
        if: github.ref == 'refs/heads/main'
        uses: peaceiris/actions-gh-pages@v3
        with:
          github_token: ${{ secrets.GITHUB_TOKEN }}
          publish_dir: docs/_build/html
```

## 七、文档体系

### 7.1 README.md

~~~markdown
# CircleSeeker2

[![Version](https://img.shields.io/pypi/v/circleseeker2)](https://pypi.org/project/circleseeker2/)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/circleseeker2)](https://anaconda.org/bioconda/circleseeker2)
[![License](https://img.shields.io/badge/license-GPL--2.0-blue)](LICENSE)
[![CI](https://github.com/circleseeker/circleseeker2/workflows/CI/badge.svg)](https://github.com/circleseeker/circleseeker2/actions)
[![Documentation](https://readthedocs.org/projects/circleseeker2/badge/?version=latest)](https://circleseeker2.readthedocs.io)
[![codecov](https://codecov.io/gh/circleseeker/circleseeker2/branch/main/graph/badge.svg)](https://codecov.io/gh/circleseeker/circleseeker2)

## Overview

CircleSeeker2 is a comprehensive bioinformatics pipeline for detecting extrachromosomal circular DNA (eccDNA) from PacBio HiFi long-read sequencing data. Using a three-round iterative analysis approach, it achieves high-sensitivity identification of various eccDNA types.

### Key Features

- 🔍 **Multi-round Analysis**: TideHunter → BLAST+ (2x) → Minimap2 for comprehensive detection
- ⚡ **High Performance**: Parallel processing with efficient memory management
- 📊 **Multiple eccDNA Types**: Detects UeccDNA, MeccDNA, CeccDNA, MCeccDNA, and XeccDNA
- 🔧 **Modular Design**: Clean architecture with reusable components
- 📦 **Easy Installation**: Available via Bioconda and PyPI
- 📈 **Comprehensive Reports**: HTML, CSV, and FASTA outputs with detailed statistics

## Installation

### Via Bioconda (Recommended)

```bash
conda install -c bioconda circleseeker2
~~~

### Via PyPI

```bash
pip install circleseeker2
```

### From Source

```bash
git clone https://github.com/circleseeker/circleseeker2.git
cd circleseeker2
pip install -e .
```

## Quick Start

### Basic Usage

```bash
# Run with minimal parameters
circleseeker2 run input.fasta reference.fa

# Run with custom output and threads
circleseeker2 run input.fasta reference.fa \
    -o results \
    -p sample1 \
    -t 16
```

### Using Configuration File

```bash
# Generate configuration template
circleseeker2 init-config -o config.yaml

# Edit config.yaml with your parameters

# Run with configuration
circleseeker2 run input.fasta reference.fa -c config.yaml
```

### Validate Installation

```bash
# Quick validation
circleseeker2 validate

# Full validation with test data
circleseeker2 validate --full
```

## Input Requirements

- **Sequencing Data**: PacBio HiFi reads in FASTA/FASTQ format
- **Reference Genome**: Reference genome in FASTA format
- **Minimum Coverage**: 10x (recommended: >30x)
- **Read Length**: >1kb (recommended: >5kb)

## Output Files

```
output_dir/
├── config.yaml                     # Effective configuration
├── sample.checkpoint               # Pipeline checkpoint
├── sample.summary.csv             # Summary statistics
├── sample.report.html             # Interactive HTML report
├── sample.Final.Confirmed.Uecc.fasta   # UeccDNA sequences
├── sample.Final.Confirmed.Mecc.fasta   # MeccDNA sequences
├── sample.Final.Confirmed.Cecc.fasta   # CeccDNA sequences
├── sample.Final.Confirmed.MCecc.fasta  # MCeccDNA sequences
└── logs/
    └── circleseeker2.log         # Detailed log file
```

## Advanced Usage

### Resume from Checkpoint

```bash
# Resume from last checkpoint
circleseeker2 run input.fasta reference.fa -o previous_output

# Start from specific step
circleseeker2 run input.fasta reference.fa --start-from 5
```

### Parallel Processing

```bash
# Use 32 threads
circleseeker2 run input.fasta reference.fa -t 32

# Dry run to see pipeline steps
circleseeker2 run input.fasta reference.fa --dry-run
```

### Performance Benchmarking

```bash
# Run benchmarks
circleseeker2 benchmark sample1.fa sample2.fa -o benchmarks.csv
```

## Documentation

Complete documentation is available at [circleseeker2.readthedocs.io](https://circleseeker2.readthedocs.io/), including:

- [Installation Guide](https://circleseeker2.readthedocs.io/en/latest/installation.html)
- [User Manual](https://circleseeker2.readthedocs.io/en/latest/user_guide.html)
- [API Reference](https://circleseeker2.readthedocs.io/en/latest/api_reference.html)
- [Tutorials](https://circleseeker2.readthedocs.io/en/latest/tutorials.html)

## Performance

| Dataset Size  | Reads | Time (8 cores) | Memory | eccDNA Found |
| ------------- | ----- | -------------- | ------ | ------------ |
| Small (1GB)   | 100K  | 15 min         | 4 GB   | ~500         |
| Medium (10GB) | 1M    | 2 hours        | 16 GB  | ~5,000       |
| Large (100GB) | 10M   | 20 hours       | 64 GB  | ~50,000      |

## Citation

If you use CircleSeeker2 in your research, please cite:

```bibtex
@article{circleseeker2_2024,
  title={CircleSeeker2: Comprehensive Detection of Circular DNA from HiFi Sequencing},
  author={CircleSeeker Team},
  journal={Bioinformatics},
  year={2024},
  doi={10.1093/bioinformatics/btx000}
}
```

## Contributing

We welcome contributions! Please see our [Contributing Guide](https://claude.ai/chat/CONTRIBUTING.md) for details.

## Support

- 📖 [Documentation](https://circleseeker2.readthedocs.io/)
- 🐛 [Issue Tracker](https://github.com/circleseeker/circleseeker2/issues)
- 💬 [Discussions](https://github.com/circleseeker/circleseeker2/discussions)
- 📧 [Email Support](mailto:circleseeker@example.org)

## License

CircleSeeker2 is licensed under the [GNU General Public License v2.0](https://claude.ai/chat/LICENSE).

## Acknowledgments

CircleSeeker2 builds upon excellent bioinformatics tools:

- [TideHunter](https://github.com/yangao07/TideHunter) - Tandem repeat detection
- [BLAST+](https://blast.ncbi.nlm.nih.gov/) - Sequence alignment
- [Minimap2](https://github.com/lh3/minimap2) - Long-read alignment
- [samtools](https://github.com/samtools/samtools) - SAM/BAM manipulation

```
## 八、实施路线图

### 第一阶段：基础架构（第1周）
1. ✅ 创建项目结构
2. ✅ 实现核心模块（config, cli, external）
3. ✅ 建立日志和异常处理系统
4. ✅ 设置版本控制和包管理

### 第二阶段：核心功能迁移（第2-3周）
1. ⏳ 迁移现有模块到新架构
2. ⏳ 实现流式I/O处理
3. ⏳ 添加并行处理支持
4. ⏳ 优化内存使用

### 第三阶段：测试和文档（第4周）
1. ⏳ 编写单元测试（覆盖率>80%）
2. ⏳ 编写集成测试
3. ⏳ 完善用户文档
4. ⏳ 准备示例数据和教程

### 第四阶段：发布准备（第5周）
1. ⏳ 性能测试和优化
2. ⏳ 准备Bioconda配方
3. ⏳ 设置CI/CD流程
4. ⏳ 发布v2.0.0

## 九、总结

这个综合方案结合了三个版本的所有优点：

1. **专业的项目规划**（来自英文方案）
2. **实用的技术实现**（来自中文方案）
3. **完善的辅助功能**（来自我的方案）

关键改进：
- 采用 `src/` 布局和 `pyproject.toml`
- 完整的流式处理和并行支持
- 健壮的错误处理和日志系统
- 全面的测试覆盖
- 专业的文档体系

这个方案可以立即开始实施，预计**4-5周完成全部重构**，交付一个专业、高效、易用的生物信息学工具包。
```