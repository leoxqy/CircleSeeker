#!/usr/bin/env python3
"""
run_minimap2_alignment.py - Minimap2 比对执行脚本

用于执行 HiFi 读长与参考基因组的比对，包括：
1. 构建参考基因组索引
2. 执行序列比对
3. 排序和索引 BAM 文件
"""

import subprocess
import shutil
import sys
import logging
import argparse
import tempfile
from pathlib import Path
from typing import Optional, List

# 配置 logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class Minimap2Aligner:
    """Minimap2 比对器"""
    
    def __init__(self, num_threads=8, preset="map-hifi", allow_secondary=True):
        """
        初始化 Minimap2 比对器
        
        Args:
            num_threads: 线程数
            preset: minimap2 预设模式 (map-hifi, map-pb, map-ont 等)
            allow_secondary: 是否允许二次比对
        """
        self.num_threads = num_threads
        self.preset = preset
        self.allow_secondary = allow_secondary
        self._check_installation()
    
    def _check_installation(self):
        """检查必需的工具是否安装"""
        required_tools = ["minimap2", "samtools"]
        missing_tools = []
        
        for tool in required_tools:
            if not shutil.which(tool):
                missing_tools.append(tool)
                logger.error(f"{tool} not found in PATH")
        
        if missing_tools:
            logger.error(f"Missing tools: {', '.join(missing_tools)}")
            logger.error("Please install via conda:")
            logger.error("conda install -c bioconda minimap2 samtools")
            raise RuntimeError(f"Required tools not installed: {', '.join(missing_tools)}")
        
        logger.debug("All required tools found in PATH")
    
    def build_index(self, reference_fasta, index_path=None):
        """
        构建参考基因组索引
        
        Args:
            reference_fasta: 参考基因组 FASTA 文件路径
            index_path: 索引文件输出路径（可选）
        
        Returns:
            索引文件路径
        """
        if index_path is None:
            # 自动生成索引文件名
            index_path = Path(str(reference_fasta).replace(".fa", ".mmi").replace(".fasta", ".mmi"))
        else:
            index_path = Path(index_path)
        
        # 如果索引已存在，跳过构建
        if index_path.exists():
            logger.info(f"Index already exists: {index_path}")
            return index_path
        
        command = [
            "minimap2",
            "-t", str(self.num_threads),
            "-d", str(index_path),
            str(reference_fasta)
        ]
        
        logger.info(f"Building minimap2 index for {reference_fasta}")
        logger.debug(f"Command: {' '.join(command)}")
        
        try:
            result = subprocess.run(
                command,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
            
            logger.info(f"Index built successfully: {index_path}")
            return index_path
            
        except subprocess.CalledProcessError as e:
            logger.error(f"minimap2 index building failed with exit code {e.returncode}")
            logger.error(f"Error message: {e.stderr}")
            raise
    
    def align(self, reference, reads, output_bam, build_index_first=True):
        """
        执行序列比对并生成排序的 BAM 文件
        
        Args:
            reference: 参考基因组 FASTA 文件或索引文件
            reads: 读长文件 (FASTA/FASTQ)
            output_bam: 输出 BAM 文件路径
            build_index_first: 是否先构建索引
        
        Returns:
            输出 BAM 文件路径
        """
        reference = Path(reference)
        reads = Path(reads)
        output_bam = Path(output_bam)
        
        # 检查输入文件
        if not reads.exists():
            raise FileNotFoundError(f"Reads file not found: {reads}")
        
        # 确定参考序列（索引或 FASTA）
        if build_index_first and reference.suffix in ['.fa', '.fasta', '.fna']:
            ref_index = self.build_index(reference)
        else:
            ref_index = reference
            if not ref_index.exists():
                raise FileNotFoundError(f"Reference file not found: {ref_index}")
        
        # 创建输出目录
        output_bam.parent.mkdir(parents=True, exist_ok=True)
        
        # 创建临时目录用于排序
        temp_dir = output_bam.parent / "tmp_sort"
        temp_dir.mkdir(parents=True, exist_ok=True)
        
        try:
            self._run_alignment_pipeline(ref_index, reads, output_bam, temp_dir)
            
            # 创建 BAM 索引
            self._index_bam(output_bam)
            
            logger.info(f"Alignment completed successfully")
            logger.info(f"Output BAM: {output_bam}")
            logger.info(f"Output BAI: {output_bam}.bai")
            
            # 输出统计信息
            self._print_stats(output_bam)
            
            return output_bam
            
        finally:
            # 清理临时目录
            if temp_dir.exists():
                shutil.rmtree(temp_dir, ignore_errors=True)
    
    def _run_alignment_pipeline(self, reference, reads, output_bam, temp_dir):
        """
        运行比对流水线：minimap2 | samtools sort
        
        Args:
            reference: 参考序列（索引或 FASTA）
            reads: 读长文件
            output_bam: 输出 BAM 文件
            temp_dir: 临时目录
        """
        # 构建 minimap2 命令
        minimap2_cmd = [
            "minimap2",
            "-t", str(self.num_threads),
            "-ax", self.preset
        ]
        
        # 添加二次比对参数
        if self.allow_secondary:
            minimap2_cmd.append("--secondary=yes")
        else:
            minimap2_cmd.append("--secondary=no")
        
        minimap2_cmd.extend([str(reference), str(reads)])
        
        # 构建 samtools sort 命令
        sort_threads = max(1, self.num_threads // 2)
        samtools_cmd = [
            "samtools", "sort",
            "-@", str(sort_threads),
            "-T", str(temp_dir / "sort_tmp"),
            "-o", str(output_bam)
        ]
        
        logger.info(f"Running alignment pipeline")
        logger.debug(f"Minimap2 command: {' '.join(minimap2_cmd)}")
        logger.debug(f"Samtools command: {' '.join(samtools_cmd)}")
        
        # 运行管道
        try:
            # 启动 minimap2 进程
            minimap2_proc = subprocess.Popen(
                minimap2_cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
            
            # 启动 samtools sort 进程
            samtools_proc = subprocess.Popen(
                samtools_cmd,
                stdin=minimap2_proc.stdout,
                stderr=subprocess.PIPE
            )
            
            # 关闭 minimap2 的 stdout（让 samtools 接管）
            minimap2_proc.stdout.close()
            
            # 等待两个进程完成
            minimap2_stderr = minimap2_proc.stderr.read().decode()
            samtools_stderr = samtools_proc.stderr.read().decode()
            
            minimap2_ret = minimap2_proc.wait()
            samtools_ret = samtools_proc.wait()
            
            # 检查返回码
            if minimap2_ret != 0:
                logger.error(f"minimap2 failed with exit code {minimap2_ret}")
                logger.error(f"minimap2 stderr: {minimap2_stderr}")
                raise subprocess.CalledProcessError(minimap2_ret, "minimap2")
            
            if samtools_ret != 0:
                logger.error(f"samtools sort failed with exit code {samtools_ret}")
                logger.error(f"samtools stderr: {samtools_stderr}")
                raise subprocess.CalledProcessError(samtools_ret, "samtools sort")
            
            # 记录日志信息
            if minimap2_stderr:
                for line in minimap2_stderr.strip().split('\n'):
                    if line:
                        logger.debug(f"minimap2: {line}")
            
        except Exception as e:
            logger.error(f"Alignment pipeline failed: {e}")
            raise
    
    def _index_bam(self, bam_file):
        """
        为 BAM 文件创建索引
        
        Args:
            bam_file: BAM 文件路径
        """
        command = ["samtools", "index", str(bam_file)]
        
        logger.info(f"Creating BAM index")
        logger.debug(f"Command: {' '.join(command)}")
        
        try:
            subprocess.run(
                command,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
            logger.info(f"BAM index created successfully")
            
        except subprocess.CalledProcessError as e:
            logger.error(f"samtools index failed with exit code {e.returncode}")
            logger.error(f"Error message: {e.stderr}")
            raise
    
    def _print_stats(self, bam_file):
        """
        输出 BAM 文件统计信息
        
        Args:
            bam_file: BAM 文件路径
        """
        try:
            # 获取比对统计
            result = subprocess.run(
                ["samtools", "flagstat", str(bam_file)],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
            
            logger.info("Alignment statistics:")
            for line in result.stdout.strip().split('\n'):
                if line:
                    logger.info(f"  {line}")
            
            # 获取文件大小
            bam_size = bam_file.stat().st_size / (1024 * 1024)  # MB
            logger.info(f"Output BAM size: {bam_size:.2f} MB")
            
        except Exception as e:
            logger.warning(f"Could not get BAM statistics: {e}")


def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description="Run minimap2 alignment for HiFi/PacBio/ONT reads",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # 输入输出参数
    parser.add_argument(
        "-r", "--reference",
        type=Path,
        required=True,
        help="Reference genome FASTA file or minimap2 index (.mmi)"
    )
    parser.add_argument(
        "-i", "--input",
        type=Path,
        required=True,
        help="Input reads file (FASTA/FASTQ)"
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        required=True,
        help="Output sorted BAM file"
    )
    
    # 比对参数
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=8,
        help="Number of threads (default: 8)"
    )
    parser.add_argument(
        "-p", "--preset",
        default="map-hifi",
        choices=["map-hifi", "map-pb", "map-ont", "asm5", "asm10", "asm20"],
        help="Minimap2 preset mode (default: map-hifi)"
    )
    parser.add_argument(
        "--no-secondary",
        action="store_true",
        help="Disable secondary alignments"
    )
    
    # 索引参数
    parser.add_argument(
        "--index-only",
        action="store_true",
        help="Only build index, don't perform alignment"
    )
    parser.add_argument(
        "--index-path",
        type=Path,
        help="Custom path for index file"
    )
    parser.add_argument(
        "--skip-index",
        action="store_true",
        help="Skip index building (assume reference is already an index)"
    )
    
    # 其他参数
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
    if not args.reference.exists():
        logger.error(f"Reference file not found: {args.reference}")
        sys.exit(1)
    
    if not args.index_only and not args.input.exists():
        logger.error(f"Input file not found: {args.input}")
        sys.exit(1)
    
    # 运行比对
    try:
        aligner = Minimap2Aligner(
            num_threads=args.threads,
            preset=args.preset,
            allow_secondary=not args.no_secondary
        )
        
        if args.index_only:
            # 只构建索引
            index_path = aligner.build_index(args.reference, args.index_path)
            logger.info(f"Index building completed: {index_path}")
        else:
            # 执行完整比对流程
            aligner.align(
                reference=args.reference,
                reads=args.input,
                output_bam=args.output,
                build_index_first=not args.skip_index
            )
        
        sys.exit(0)
        
    except Exception as e:
        logger.error(f"Failed to run alignment: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
