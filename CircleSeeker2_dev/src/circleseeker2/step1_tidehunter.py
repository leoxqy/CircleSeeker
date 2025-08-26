#!/usr/bin/env python3
"""
run_tidehunter.py - TideHunter 执行脚本
"""

import subprocess
import shutil
import sys
import logging
import argparse
from pathlib import Path

# 配置 logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class TideHunterRunner:
    """TideHunter 运行器"""
    
    def __init__(self, num_threads=8):
        self.num_threads = num_threads
        self._check_installation()
    
    def _check_installation(self):
        """检查 TideHunter 是否安装"""
        if not shutil.which("TideHunter"):
            logger.error("TideHunter not found in PATH")
            logger.error("Please install: conda install -c bioconda tidehunter")
            raise RuntimeError("TideHunter not installed")
        logger.debug("TideHunter found in PATH")
    
    def run(self, input_fasta, output_path):
        """
        运行 TideHunter
        
        使用参数：-f 2 -t {threads} -k 16 -w 1 -p 100 -P 2000000 -e 0.1
        """
        command = [
            "TideHunter",
            "-f", "2",
            "-t", str(self.num_threads),
            "-k", "16",
            "-w", "1", 
            "-p", "100",
            "-P", "2000000",
            "-e", "0.1",
            str(input_fasta)
        ]
        
        logger.info(f"Running TideHunter on {input_fasta}")
        logger.debug(f"Command: {' '.join(command)}")
        
        try:
            with open(output_path, 'w') as out:
                result = subprocess.run(
                    command,
                    stdout=out,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True
                )
            
            logger.info(f"TideHunter completed successfully")
            logger.info(f"Output saved to: {output_path}")
            
            # 记录输出文件大小
            output_size = Path(output_path).stat().st_size
            logger.debug(f"Output file size: {output_size} bytes")
            
            return True
            
        except subprocess.CalledProcessError as e:
            logger.error(f"TideHunter failed with exit code {e.returncode}")
            logger.error(f"Error message: {e.stderr}")
            raise
        except Exception as e:
            logger.error(f"Unexpected error: {e}")
            raise


def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description="Run TideHunter for tandem repeat detection"
    )
    parser.add_argument(
        "-i", "--input",
        type=Path,
        required=True,
        help="Input FASTA file"
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        required=True,
        help="Output file path"
    )
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=8,
        help="Number of threads (default: 8)"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose logging"
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
    args.output.parent.mkdir(parents=True, exist_ok=True)
    
    # 运行 TideHunter
    try:
        runner = TideHunterRunner(num_threads=args.threads)
        runner.run(args.input, args.output)
        sys.exit(0)
    except Exception as e:
        logger.error(f"Failed to run TideHunter: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
