import argparse
import subprocess
import os

def run_command(command):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        raise Exception(f"Command failed: {command}\nError: {stderr.decode()}")
    return stdout.decode()

def extract_unmatched_reads(input_fastq, read_list_csv, output_fasta):
    # 从CSV文件中提取匹配的read名称
    command = f"cut -d',' -f1 {read_list_csv} > matched_reads.txt"
    run_command(command)

    # 使用seqkit提取未匹配的reads
    command = f"seqkit grep -v -f matched_reads.txt {input_fastq} > {output_fasta}"
    run_command(command)

    # 清理临时文件
    os.remove("matched_reads.txt")

def count_reads(input_fasta):
    command = f"seqkit stats {input_fasta} | awk 'NR>1 {{print $4}}'"
    count = run_command(command).strip()
    return count

def main(args):
    # 1. 提取未匹配的reads
    print("正在提取未匹配的reads...")
    extract_unmatched_reads(args.input, args.csv, args.unmatched)

    # 2. 计算未匹配reads的数量
    print("正在计算未匹配reads的数量...")
    read_count = count_reads(args.unmatched)

    print("处理完成。结果：")
    print(f"1. 未匹配的reads文件: {args.unmatched}")
    print(f"2. 未匹配reads的数量: {read_count}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="提取未匹配reads并计算数量")
    parser.add_argument("-i", "--input", required=True, help="输入的FASTQ文件")
    parser.add_argument("-c", "--csv", required=True, help="包含匹配reads的CSV文件")
    parser.add_argument("-u", "--unmatched", default="unSel.fasta", help="输出的未匹配reads文件名 (默认: unSel.fasta)")
    args = parser.parse_args()

    main(args)
