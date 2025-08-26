#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_cd_hit_est.py - cd-hit-est 执行脚本

包装命令：
cd-hit-est \
  -i <input.fasta> \
  -o <output_prefix> \
  -c 0.99 -n 10 \
  -s 0.99 -aS 0.99 -G 1 \
  -d 0 \
  -M 0 \
  -T <threads>
"""

import subprocess
import shutil
import sys
import logging
import argparse
from pathlib import Path

# 配置 logging（默认 INFO，可用 -v 切换到 DEBUG）
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class CDHitEstRunner:
    """cd-hit-est 运行器"""

    def __init__(self,
                 threads=24,
                 c=0.99,
                 n=10,
                 s=0.99,
                 aS=0.99,
                 G=1,
                 d=0,
                 M=0,
                 extra=None):
        self.threads = int(threads)
        self.c = float(c)
        self.n = int(n)
        self.s = float(s)
        self.aS = float(aS)
        self.G = int(G)
        self.d = int(d)
        self.M = int(M)
        self.extra = extra or []
        self._check_installation()

    def _check_installation(self):
        """检查 cd-hit-est 是否安装"""
        if not shutil.which("cd-hit-est"):
            logger.error("cd-hit-est not found in PATH")
            logger.error("请先安装：conda install -c bioconda cd-hit")
            raise RuntimeError("cd-hit-est not installed")
        logger.debug("cd-hit-est found in PATH")

    def _build_command(self, input_fasta: Path, output_prefix: Path):
        """构建 cd-hit-est 命令行"""
        cmd = [
            "cd-hit-est",
            "-i", str(input_fasta),
            "-o", str(output_prefix),
            "-c", str(self.c),
            "-n", str(self.n),
            "-s", str(self.s),
            "-aS", str(self.aS),
            "-G", str(self.G),
            "-d", str(self.d),
            "-M", str(self.M),
            "-T", str(self.threads),
        ]
        if self.extra:
            cmd.extend(self.extra)
        return cmd

    def run(self, input_fasta: Path, output_prefix: Path):
        """
        运行 cd-hit-est

        参数：
            input_fasta: 输入 FASTA
            output_prefix: 输出前缀（与 cd-hit-est -o 含义一致）
        产出：
            - {output_prefix}         代表序列 FASTA
            - {output_prefix}.clstr  聚类信息
        """
        command = self._build_command(input_fasta, output_prefix)
        logger.info(f"Running cd-hit-est on: {input_fasta}")
        logger.debug(f"Command: {' '.join(command)}")

        try:
            result = subprocess.run(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
            # cd-hit-est 的日志多在 stdout/stderr，DEBUG 时打印
            if result.stdout:
                logger.debug("cd-hit-est stdout:\n" + result.stdout)
            if result.stderr:
                logger.debug("cd-hit-est stderr:\n" + result.stderr)

            logger.info("cd-hit-est completed successfully")
            logger.info(f"Output (representatives): {output_prefix}")
            clstr_path = Path(str(output_prefix) + ".clstr")
            if clstr_path.exists():
                logger.info(f"Cluster file: {clstr_path}")

            # 记录输出文件大小（DEBUG）
            try:
                if Path(output_prefix).exists():
                    size_main = Path(output_prefix).stat().st_size
                    logger.debug(f"{output_prefix} size: {size_main} bytes")
                if clstr_path.exists():
                    size_clstr = clstr_path.stat().st_size
                    logger.debug(f"{clstr_path} size: {size_clstr} bytes")
            except Exception:
                pass

            return True

        except subprocess.CalledProcessError as e:
            logger.error(f"cd-hit-est failed with exit code {e.returncode}")
            if e.stdout:
                logger.debug("cd-hit-est stdout (on error):\n" + e.stdout)
            if e.stderr:
                logger.error("cd-hit-est stderr:\n" + e.stderr)
            raise
        except Exception as e:
            logger.error(f"Unexpected error: {e}")
            raise


def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description="Run cd-hit-est for clustering nucleotide sequences"
    )
    parser.add_argument(
        "-i", "--input",
        type=Path,
        required=True,
        help="Input FASTA file"
    )
    parser.add_argument(
        "-o", "--output-prefix",
        type=Path,
        required=True,
        help="Output prefix (same meaning as cd-hit-est -o)"
    )
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=24,
        help="Number of threads (default: 24)"
    )
    # 与示例保持一致的默认参数，可自行覆盖
    parser.add_argument("--c", type=float, default=0.99, help="Sequence identity threshold (default: 0.99)")
    parser.add_argument("--n", type=int,   default=10,   help="Word length (default: 10)")
    parser.add_argument("--s", type=float, default=0.99, help="Minimum length of the shorter sequences / query (default: 0.99)")
    parser.add_argument("--aS", type=float, default=0.99, help="Alignment coverage for the shorter sequence (default: 0.99)")
    parser.add_argument("--G", type=int,   default=1,    help="Use global sequence identity (1) or local (0). Default: 1")
    parser.add_argument("--d", type=int,   default=0,    help="Length of description in output FASTA header (0 for full). Default: 0")
    parser.add_argument("--M", type=int,   default=0,    help="Memory limit in MB (0 for unlimited). Default: 0")
    parser.add_argument(
        "--extra",
        nargs=argparse.REMAINDER,
        help="Pass any extra args to cd-hit-est after '--extra' (e.g. --extra -bak 1)"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose (DEBUG) logging"
    )

    args = parser.parse_args()

    # 调整日志级别
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # 检查输入文件
    if not args.input.exists():
        logger.error(f"Input file not found: {args.input}")
        sys.exit(1)

    # 创建输出目录（如果需要）
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)

    try:
        runner = CDHitEstRunner(
            threads=args.threads,
            c=args.c,
            n=args.n,
            s=args.s,
            aS=args.aS,
            G=args.G,
            d=args.d,
            M=args.M,
            extra=args.extra
        )
        runner.run(args.input, args.output_prefix)
        sys.exit(0)
    except Exception as e:
        logger.error(f"Failed to run cd-hit-est: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
