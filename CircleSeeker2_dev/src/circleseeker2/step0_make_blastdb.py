#!/usr/bin/env python3
"""
make_blastdb.py - BLAST数据库构建脚本
用于从参考基因组FASTA文件创建BLAST数据库并建立索引
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


class BlastDatabaseBuilder:
    """BLAST数据库构建器"""
    
    def __init__(self, dbtype='nucl', parse_seqids=False, taxid=None):
        """
        初始化数据库构建器
        
        Args:
            dbtype: 数据库类型 ('nucl' 或 'prot')
            parse_seqids: 是否解析序列ID以支持序列检索
            taxid: 分类ID（可选）
        """
        self.dbtype = dbtype
        self.parse_seqids = parse_seqids
        self.taxid = taxid
        self._check_installations()
    
    def _check_installations(self):
        """检查必需的程序是否安装"""
        missing_tools = []
        
        # 检查 makeblastdb
        if not shutil.which("makeblastdb"):
            missing_tools.append("makeblastdb (BLAST+)")
            logger.error("makeblastdb not found in PATH")
        else:
            logger.debug("makeblastdb found in PATH")
        
        # 检查 samtools
        if not shutil.which("samtools"):
            logger.warning("samtools not found in PATH")
            logger.warning("samtools is optional but recommended for FASTA indexing")
            logger.info("Install with: conda install -c bioconda samtools")
        else:
            logger.debug("samtools found in PATH")
        
        # 如果缺少必需工具，抛出错误
        if "makeblastdb (BLAST+)" in missing_tools:
            logger.error("Required tools missing:")
            logger.error("Please install BLAST+:")
            logger.error("  Ubuntu/Debian: sudo apt-get install ncbi-blast+")
            logger.error("  macOS: brew install blast")
            logger.error("  Conda: conda install -c bioconda blast")
            logger.error("  Or download from: https://blast.ncbi.nlm.nih.gov/")
            raise RuntimeError("Required tools not installed")
    
    def _check_input_file(self, input_file):
        """
        检查输入FASTA文件
        
        Args:
            input_file: 输入文件路径
            
        Returns:
            tuple: (是否有效, 文件大小MB, 序列数量估计)
        """
        input_path = Path(input_file)
        
        if not input_path.exists():
            logger.error(f"Input file not found: {input_file}")
            return False, 0, 0
        
        if not input_path.is_file():
            logger.error(f"Input path is not a file: {input_file}")
            return False, 0, 0
        
        # 获取文件大小
        file_size_mb = input_path.stat().st_size / (1024 * 1024)
        
        # 快速估计序列数量（统计'>'开头的行）
        seq_count = 0
        try:
            with open(input_path, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        seq_count += 1
        except Exception as e:
            logger.warning(f"Could not count sequences: {e}")
            seq_count = -1
        
        logger.info(f"Input file: {input_file}")
        logger.info(f"  Size: {file_size_mb:.2f} MB")
        if seq_count > 0:
            logger.info(f"  Sequences: {seq_count:,}")
        
        return True, file_size_mb, seq_count
    
    def build_database(self, input_file, output_db, title=None, max_file_size="1GB"):
        """
        构建BLAST数据库
        
        Args:
            input_file: 输入FASTA文件
            output_db: 输出数据库名称
            title: 数据库标题（可选）
            max_file_size: 最大文件大小（默认1GB）
            
        Returns:
            bool: 构建是否成功
        """
        # 构建命令
        command = [
            "makeblastdb",
            "-in", str(input_file),
            "-dbtype", self.dbtype,
            "-out", str(output_db),
            "-max_file_sz", max_file_size
        ]
        
        # 添加可选参数
        if self.parse_seqids:
            command.extend(["-parse_seqids"])
            logger.info("Sequence ID parsing enabled")
        
        if title:
            command.extend(["-title", title])
        
        if self.taxid:
            command.extend(["-taxid", str(self.taxid)])
        
        logger.info(f"Building BLAST database: {output_db}")
        logger.debug(f"Command: {' '.join(command)}")
        
        # 记录开始时间
        start_time = time.time()
        
        try:
            # 运行 makeblastdb
            result = subprocess.run(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                check=True
            )
            
            # 计算运行时间
            elapsed_time = time.time() - start_time
            
            logger.info(f"Database built successfully in {elapsed_time:.2f} seconds")
            
            # 解析输出信息
            if result.stdout:
                for line in result.stdout.split('\n'):
                    if 'sequences' in line.lower() or 'letters' in line.lower():
                        logger.info(f"  {line.strip()}")
            
            # 列出生成的文件
            self._list_database_files(output_db)
            
            return True
            
        except subprocess.CalledProcessError as e:
            logger.error(f"makeblastdb failed with exit code {e.returncode}")
            logger.error(f"Error message: {e.stdout}")
            raise
        except Exception as e:
            logger.error(f"Unexpected error: {e}")
            raise
    
    def create_fasta_index(self, fasta_file):
        """
        使用samtools创建FASTA索引
        
        Args:
            fasta_file: FASTA文件路径
            
        Returns:
            bool: 索引创建是否成功
        """
        # 检查samtools是否可用
        if not shutil.which("samtools"):
            logger.warning("Skipping FASTA indexing: samtools not installed")
            return False
        
        logger.info(f"Creating FASTA index with samtools")
        
        command = ["samtools", "faidx", str(fasta_file)]
        logger.debug(f"Command: {' '.join(command)}")
        
        try:
            # 运行 samtools faidx
            start_time = time.time()
            
            result = subprocess.run(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
            
            elapsed_time = time.time() - start_time
            
            # 检查索引文件是否创建
            index_file = Path(f"{fasta_file}.fai")
            if index_file.exists():
                index_size_kb = index_file.stat().st_size / 1024
                logger.info(f"FASTA index created successfully in {elapsed_time:.2f} seconds")
                logger.info(f"  Index file: {index_file}")
                logger.info(f"  Index size: {index_size_kb:.2f} KB")
                
                # 读取索引文件获取序列信息
                with open(index_file, 'r') as f:
                    lines = f.readlines()
                    logger.info(f"  Indexed sequences: {len(lines)}")
                
                return True
            else:
                logger.error("Index file was not created")
                return False
                
        except subprocess.CalledProcessError as e:
            logger.error(f"samtools faidx failed with exit code {e.returncode}")
            logger.error(f"Error message: {e.stderr}")
            return False
        except Exception as e:
            logger.error(f"Unexpected error during indexing: {e}")
            return False
    
    def _list_database_files(self, db_name):
        """列出生成的数据库文件"""
        db_path = Path(db_name)
        db_dir = db_path.parent if db_path.parent != Path('.') else Path('.')
        db_base = db_path.name
        
        # BLAST数据库文件扩展名
        if self.dbtype == 'nucl':
            extensions = ['.nhr', '.nin', '.nsq', '.ndb', '.njs', '.not', '.ntf', '.nto']
        else:  # prot
            extensions = ['.phr', '.pin', '.psq', '.pdb', '.pjs', '.pot', '.ptf', '.pto']
        
        logger.info("Generated database files:")
        total_size = 0
        
        for ext in extensions:
            file_path = db_dir / f"{db_base}{ext}"
            if file_path.exists():
                file_size = file_path.stat().st_size / (1024 * 1024)  # MB
                total_size += file_size
                logger.info(f"  {file_path.name}: {file_size:.2f} MB")
        
        logger.info(f"Total database size: {total_size:.2f} MB")
    
    def verify_database(self, db_name):
        """
        验证数据库完整性
        
        Args:
            db_name: 数据库名称
            
        Returns:
            bool: 数据库是否完整
        """
        # 检查必需的文件
        if self.dbtype == 'nucl':
            required_extensions = ['.nhr', '.nin', '.nsq']
        else:  # prot
            required_extensions = ['.phr', '.pin', '.psq']
        
        missing_files = []
        for ext in required_extensions:
            file_path = Path(f"{db_name}{ext}")
            if not file_path.exists():
                missing_files.append(file_path.name)
        
        if missing_files:
            logger.error(f"Database incomplete, missing: {', '.join(missing_files)}")
            return False
        
        logger.info("Database verification passed")
        return True


def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description="Build BLAST database from FASTA file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Build nucleotide database (default)
  %(prog)s -i chm13v2.0.fa -o chm13v2_db
  
  # Build protein database with sequence ID parsing
  %(prog)s -i proteins.fa -o protein_db --dbtype prot --parse-seqids
  
  # Build database with title and taxonomy ID
  %(prog)s -i genome.fa -o genome_db --title "Human Genome CHM13 v2.0" --taxid 9606
  
  # Skip samtools indexing
  %(prog)s -i genome.fa -o genome_db --no-index
        """
    )
    
    # 必需参数
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
        help="Output database name"
    )
    
    # 数据库选项
    parser.add_argument(
        "--dbtype",
        choices=['nucl', 'prot'],
        default='nucl',
        help="Database type: nucleotide or protein (default: nucl)"
    )
    parser.add_argument(
        "--title",
        type=str,
        help="Database title"
    )
    parser.add_argument(
        "--parse-seqids",
        action="store_true",
        help="Parse sequence IDs for retrieval"
    )
    parser.add_argument(
        "--taxid",
        type=int,
        help="Taxonomy ID for all sequences"
    )
    parser.add_argument(
        "--max-file-size",
        type=str,
        default="1GB",
        help="Maximum file size for database volumes (default: 1GB)"
    )
    
    # 索引选项
    parser.add_argument(
        "--no-index",
        action="store_true",
        help="Skip samtools faidx indexing"
    )
    
    # 其他选项
    parser.add_argument(
        "--no-verify",
        action="store_true",
        help="Skip database verification"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose logging"
    )
    
    # 设置默认值（用于快速运行）
    parser.set_defaults(
        input=Path("chm13v2.0.fa"),
        output=Path("chm13v2_db")
    )
    
    args = parser.parse_args()
    
    # 调整日志级别
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    try:
        # 创建数据库构建器
        builder = BlastDatabaseBuilder(
            dbtype=args.dbtype,
            parse_seqids=args.parse_seqids,
            taxid=args.taxid
        )
        
        # 检查输入文件
        valid, size_mb, seq_count = builder._check_input_file(args.input)
        if not valid:
            logger.error("Input file validation failed")
            sys.exit(1)
        
        # 警告大文件
        if size_mb > 1000:  # > 1GB
            logger.warning(f"Large input file ({size_mb:.2f} MB), this may take some time...")
        
        # 构建数据库
        logger.info("=" * 60)
        logger.info("Starting database construction")
        logger.info("=" * 60)
        
        builder.build_database(
            input_file=args.input,
            output_db=args.output,
            title=args.title,
            max_file_size=args.max_file_size
        )
        
        # 验证数据库（除非跳过）
        if not args.no_verify:
            if not builder.verify_database(args.output):
                logger.error("Database verification failed")
                sys.exit(1)
        
        # 创建FASTA索引（除非跳过）
        if not args.no_index:
            logger.info("-" * 60)
            builder.create_fasta_index(args.input)
        
        logger.info("=" * 60)
        logger.info("Database construction completed successfully")
        logger.info(f"Database ready for use: {args.output}")
        logger.info("=" * 60)
        
        sys.exit(0)
        
    except Exception as e:
        logger.error(f"Failed to build database: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
