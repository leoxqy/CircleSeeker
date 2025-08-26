#!/usr/bin/env python3
"""
run_blast.py - BLASTN 执行脚本
用于将HeLa序列与CHM13参考基因组进行比对
"""

import subprocess
import shutil
import sys
import logging
import argparse
import time
from pathlib import Path

# 配置 logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class BlastRunner:
    """BLASTN 运行器"""
    
    def __init__(self, 
                 num_threads=8,
                 word_size=100,
                 evalue="1e-50",
                 perc_identity=99,
                 outfmt="6 std sstrand"):
        """
        初始化 BLAST 运行器
        
        Args:
            num_threads: 使用的线程数
            word_size: 初始精确匹配的长度
            evalue: E值阈值
            perc_identity: 最小序列相似度百分比
            outfmt: 输出格式
        """
        self.num_threads = num_threads
        self.word_size = word_size
        self.evalue = evalue
        self.perc_identity = perc_identity
        self.outfmt = outfmt
        self._check_installation()
    
    def _check_installation(self):
        """检查 BLAST+ 是否安装"""
        if not shutil.which("blastn"):
            logger.error("blastn not found in PATH")
            logger.error("Please install BLAST+:")
            logger.error("  Ubuntu/Debian: sudo apt-get install ncbi-blast+")
            logger.error("  macOS: brew install blast")
            logger.error("  Or download from: https://blast.ncbi.nlm.nih.gov/")
            raise RuntimeError("BLAST+ not installed")
        logger.debug("blastn found in PATH")
    
    def _check_database(self, db_name):
        """
        检查 BLAST 数据库是否存在
        
        Args:
            db_name: 数据库名称
            
        Returns:
            bool: 数据库是否存在
        """
        # Nucleotide 数据库的必需文件扩展名
        required_extensions = ['.nhr', '.nin', '.nsq']
        
        db_path = Path(db_name)
        missing_files = []
        
        for ext in required_extensions:
            db_file = Path(f"{db_name}{ext}")
            if not db_file.exists():
                missing_files.append(db_file.name)
        
        if missing_files:
            logger.error(f"BLAST database incomplete or missing: {db_name}")
            logger.error(f"Missing files: {', '.join(missing_files)}")
            logger.error("Please ensure the database is properly formatted using makeblastdb")
            return False
        
        logger.debug(f"BLAST database verified: {db_name}")
        return True
    
    def run(self, database, query_file, output_file, use_time_cmd=False):
        """
        运行 BLASTN
        
        Args:
            database: BLAST数据库路径
            query_file: 查询序列文件路径
            output_file: 输出文件路径
            use_time_cmd: 是否使用 time 命令
            
        Returns:
            bool: 运行是否成功
        """
        # 构建命令
        command = [
            "blastn",
            "-db", str(database),
            "-query", str(query_file),
            "-out", str(output_file),
            "-num_threads", str(self.num_threads),
            "-word_size", str(self.word_size),
            "-evalue", str(self.evalue),
            "-perc_identity", str(self.perc_identity),
            "-outfmt", self.outfmt
        ]
        
        # 如果使用 time 命令
        if use_time_cmd:
            command = ["time"] + command
        
        logger.info(f"Running BLASTN")
        logger.info(f"  Database: {database}")
        logger.info(f"  Query: {query_file}")
        logger.info(f"  Output: {output_file}")
        logger.debug(f"Command: {' '.join(command)}")
        
        # 记录开始时间
        start_time = time.time()
        
        try:
            # 运行 BLAST
            result = subprocess.run(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
            
            # 计算运行时间
            elapsed_time = time.time() - start_time
            
            logger.info(f"BLASTN completed successfully in {elapsed_time:.2f} seconds")
            logger.info(f"Output saved to: {output_file}")
            
            # 记录输出文件信息
            output_path = Path(output_file)
            if output_path.exists():
                output_size = output_path.stat().st_size / (1024 * 1024)  # MB
                
                # 统计结果行数
                with open(output_path, 'r') as f:
                    line_count = sum(1 for _ in f)
                
                logger.info(f"Output file size: {output_size:.2f} MB")
                logger.info(f"Number of hits: {line_count}")
            
            # 如果有 stderr 输出（通常是统计信息），记录为 debug
            if result.stderr:
                logger.debug(f"BLAST statistics:\n{result.stderr}")
            
            return True
            
        except subprocess.CalledProcessError as e:
            logger.error(f"BLASTN failed with exit code {e.returncode}")
            logger.error(f"Error message: {e.stderr}")
            raise
        except Exception as e:
            logger.error(f"Unexpected error: {e}")
            raise


def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description="Run BLASTN for sequence alignment against CHM13 reference genome",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run with default parameters
  %(prog)s -d chm13v2_db -q HeLa_rep1_10k.carousel_circular_sequences.fasta -o HeLa_vs_chm13.blast_results.tsv
  
  # Run with custom threads and E-value
  %(prog)s -d chm13v2_db -q query.fasta -o results.tsv -t 16 -e 1e-100
  
  # Run with time command
  %(prog)s -d chm13v2_db -q query.fasta -o results.tsv --time
        """
    )
    
    # 必需参数
    parser.add_argument(
        "-d", "--database",
        type=Path,
        required=True,
        default=Path("chm13v2_db"),
        help="BLAST database path (default: chm13v2_db)"
    )
    parser.add_argument(
        "-q", "--query",
        type=Path,
        required=True,
        default=Path("HeLa_rep1_10k.carousel_circular_sequences.fasta"),
        help="Query sequence file (default: HeLa_rep1_10k.carousel_circular_sequences.fasta)"
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        required=True,
        default=Path("HeLa_vs_chm13.blast_results.tsv"),
        help="Output file path (default: HeLa_vs_chm13.blast_results.tsv)"
    )
    
    # 可选参数
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=8,
        help="Number of threads (default: 8)"
    )
    parser.add_argument(
        "-w", "--word-size",
        type=int,
        default=100,
        help="Word size for initial match (default: 100)"
    )
    parser.add_argument(
        "-e", "--evalue",
        type=str,
        default="1e-50",
        help="E-value threshold (default: 1e-50)"
    )
    parser.add_argument(
        "-p", "--perc-identity",
        type=float,
        default=99,
        help="Minimum percent identity (default: 99)"
    )
    parser.add_argument(
        "-f", "--outfmt",
        type=str,
        default="6 std sstrand",
        help='Output format (default: "6 std sstrand")'
    )
    
    # 其他选项
    parser.add_argument(
        "--time",
        action="store_true",
        help="Use time command for timing"
    )
    parser.add_argument(
        "--skip-db-check",
        action="store_true",
        help="Skip database verification"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose logging"
    )
    
    # 设置默认值（用于直接运行原始命令的情况）
    parser.set_defaults(
        database=Path("chm13v2_db"),
        query=Path("HeLa_rep1_10k.carousel_circular_sequences.fasta"),
        output=Path("HeLa_vs_chm13.blast_results.tsv")
    )
    
    args = parser.parse_args()
    
    # 调整日志级别
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # 检查查询文件
    if not args.query.exists():
        logger.error(f"Query file not found: {args.query}")
        sys.exit(1)
    
    # 创建输出目录（如果需要）
    args.output.parent.mkdir(parents=True, exist_ok=True)
    
    try:
        # 创建 BLAST 运行器
        runner = BlastRunner(
            num_threads=args.threads,
            word_size=args.word_size,
            evalue=args.evalue,
            perc_identity=args.perc_identity,
            outfmt=args.outfmt
        )
        
        # 检查数据库（除非跳过）
        if not args.skip_db_check:
            if not runner._check_database(args.database):
                logger.error("Database verification failed")
                logger.info("Use --skip-db-check to skip this verification")
                sys.exit(1)
        
        # 运行 BLAST
        runner.run(
            database=args.database,
            query_file=args.query,
            output_file=args.output,
            use_time_cmd=args.time
        )
        
        logger.info("BLAST analysis completed successfully")
        sys.exit(0)
        
    except Exception as e:
        logger.error(f"Failed to run BLAST: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
