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

def run_hifiasm_assembly(input_fasta, output_prefix, threads):
    command = f"hifiasm -o {output_prefix} -t {threads} {input_fasta}"
    run_command(command)

def get_unused_reads(input_fasta, assembly_gfa, output_fasta):
    # 从GFA文件中提取所有使用的read ID
    command = f"grep '^S' {assembly_gfa} | cut -f2 > used_reads.txt"
    run_command(command)

    # 使用seqkit获取未使用的reads
    command = f"seqkit grep -v -f used_reads.txt {input_fasta} > {output_fasta}"
    run_command(command)

    # 清理临时文件
    os.remove("used_reads.txt")

def calculate_length_stats(input_fasta, output_file):
    command = f"seqkit stats {input_fasta} | awk 'NR>1 {{print $5}}'"
    stats = run_command(command)
    with open(output_file, 'a') as f:
        f.write(f"{os.path.basename(input_fasta)} 长度: {stats.strip()}\n")

def main(args):
    # 1. 提取未匹配的reads
    print("正在提取未匹配的reads...")
    extract_unmatched_reads(args.input, args.csv, args.unmatched)

    # 2. 使用Hifiasm进行拼装
    print(f"正在使用Hifiasm进行拼装 (使用 {args.threads} 个线程)...")
    run_hifiasm_assembly(args.unmatched, args.assembly_prefix, args.threads)

    # 3. 提取拼装后的contigs
    # Hifiasm输出GFA格式，需要转换为FASTA
    print("正在转换GFA到FASTA...")
    command = f"awk '/^S/{{print \">\"$2\"\\n\"$3}}' {args.assembly_prefix}.bp.p_ctg.gfa > {args.contigs}"
    run_command(command)

    # 4. 提取未使用的reads
    print("正在提取未使用的reads...")
    get_unused_reads(args.unmatched, f"{args.assembly_prefix}.bp.p_ctg.gfa", args.unused)

    # 5. 计算contigs和未使用reads的长度
    print("正在计算长度统计...")
    with open(args.stats, 'w') as f:
        f.write("长度统计:\n")
    calculate_length_stats(args.contigs, args.stats)
    calculate_length_stats(args.unused, args.stats)

    print("处理完成。结果文件：")
    print(f"1. 未匹配的reads: {args.unmatched}")
    print(f"2. 拼装后的contigs: {args.contigs}")
    print(f"3. 未使用的reads: {args.unused}")
    print(f"4. 长度统计: {args.stats}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="提取未匹配reads并进行Hifiasm拼装")
    parser.add_argument("-i", "--input", required=True, help="输入的FASTQ文件")
    parser.add_argument("-c", "--csv", required=True, help="包含匹配reads的CSV文件")
    parser.add_argument("-u", "--unmatched", default="unSel.fasta", help="输出的未匹配reads文件名 (默认: unSel.fasta)")
    parser.add_argument("-p", "--assembly-prefix", default="assembly", help="Hifiasm拼装输出前缀 (默认: assembly)")
    parser.add_argument("-g", "--contigs", default="assembled_contigs.fa", help="输出的拼装contigs文件名 (默认: assembled_contigs.fa)")
    parser.add_argument("-r", "--unused", default="unused_reads.fa", help="输出的未使用reads文件名 (默认: unused_reads.fa)")
    parser.add_argument("-s", "--stats", default="length_stats.txt", help="输出的长度统计文件名 (默认: length_stats.txt)")
    parser.add_argument("-t", "--threads", type=int, default=16, help="Hifiasm使用的线程数 (默认: 16)")
    args = parser.parse_args()

    main(args)
