import pandas as pd
from itertools import combinations
import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def read_and_preprocess_data(file_path):
    columns = [
        'readName', 'repN', 'copyNum', 'readLen', 'start', 'end', 
        'consLen', 'aveMatch', 'fullLen', 'subPos', 'consSeq'
    ]
    
    df = pd.read_csv(file_path, sep='\t', header=None, names=columns)
    df['subPos'] = df['subPos'].str.replace(',', '|')
    df['Effective_Length'] = ((df['end'] - df['start'] + 1) / df['readLen'] * 100).round(2)
    df = df[df['aveMatch'] > 99]
    df = df.drop(columns=['fullLen'])
    
    return df

def combine_rows_data(rows):
    new_row = rows.iloc[0].copy()
    new_row['repN'] = 'repM0'
    new_row['copyNum'] = rows['copyNum'].sum()
    new_row['start'] = rows['start'].min()
    new_row['end'] = rows['end'].max()
    min_start_row = rows.loc[rows['start'].idxmin()]
    new_row['consLen'] = min_start_row['consLen']
    new_row['aveMatch'] = min_start_row['aveMatch']
    new_row['subPos'] = '|'.join(rows['subPos'])
    new_row['consSeq'] = min_start_row['consSeq']
    new_row['Effective_Length'] = min(100, rows['Effective_Length'].sum())
    return new_row

def combine_rows(group):
    group = group.reset_index(drop=True)
    
    if len(group) == 2:
        cons_len_diff = abs(group['consLen'].iloc[0] - group['consLen'].iloc[1])
        if cons_len_diff <= 15:
            return pd.DataFrame([combine_rows_data(group)])
        else:
            total_effective_length = group['Effective_Length'].sum()
            max_single_effective_length = group['Effective_Length'].max()
            
            if abs(total_effective_length - 100) <= abs(max_single_effective_length - 100):
                return group
            else:
                return group[group['Effective_Length'] == max_single_effective_length]
    
    if group['Effective_Length'].sum() <= 102:
        return pd.DataFrame([combine_rows_data(group)])
    
    best_combination = None
    best_total = 0
    for r in range(1, len(group) + 1):
        for combo in combinations(range(len(group)), r):
            subset = group.iloc[list(combo)]
            subset_total = subset['Effective_Length'].sum()
            if 98 <= subset_total <= 102:
                cons_len_range = subset['consLen'].max() - subset['consLen'].min()
                if cons_len_range <= 15:
                    return pd.DataFrame([combine_rows_data(subset)])
            elif subset_total < 98 and subset_total > best_total:
                best_combination = subset
                best_total = subset_total
    
    if best_combination is not None:
        return pd.DataFrame([combine_rows_data(best_combination)])
    return group

def process_repeated_readnames(df):
    repeated_readnames = df['readName'].value_counts()[df['readName'].value_counts() > 1].index
    repeated_df = df[df['readName'].isin(repeated_readnames)]
    
    processed_repeated = []
    for name, group in repeated_df.groupby('readName'):
        processed_group = combine_rows(group)
        processed_repeated.append(processed_group)
    
    processed_repeated = pd.concat(processed_repeated, ignore_index=True)
    
    single_occurrence = df[~df['readName'].isin(repeated_readnames)]
    result_df = pd.concat([single_occurrence, processed_repeated], ignore_index=True)
    
    return result_df

def classify_readnames(df):
    read_list_df = pd.DataFrame(columns=['readName', 'readClass'])
    
    repeated_readnames = df['readName'].value_counts()[df['readName'].value_counts() > 1].index
    read_list_df = pd.concat([read_list_df, 
                              pd.DataFrame({'readName': repeated_readnames, 'readClass': 'CTCR-M'})], 
                             ignore_index=True)
    
    single_occurrence = df[~df['readName'].isin(repeated_readnames)]
    single_repM = single_occurrence[single_occurrence['repN'].str.startswith('repM', na=False)]
    read_list_df = pd.concat([read_list_df, 
                              pd.DataFrame({'readName': single_repM['readName'], 'readClass': 'CTCR-R'})], 
                             ignore_index=True)
    
    remaining = df[~df['readName'].isin(read_list_df['readName'])]
    
    ctcr_p = remaining[remaining['Effective_Length'] >= 99]
    read_list_df = pd.concat([read_list_df, 
                              pd.DataFrame({'readName': ctcr_p['readName'], 'readClass': 'CTCR-P'})], 
                             ignore_index=True)
    
    ctcr_i = remaining[(remaining['Effective_Length'] >= 70) & (remaining['Effective_Length'] < 99)]
    read_list_df = pd.concat([read_list_df, 
                              pd.DataFrame({'readName': ctcr_i['readName'], 'readClass': 'CTCR-I'})], 
                             ignore_index=True)
    
    return read_list_df

def read_fasta(file_path):
    sequences = list(SeqIO.parse(file_path, "fasta"))
    print(f"读取了 {len(sequences)} 个序列")
    return sequences

def circularize_sequences(sequences):
    circular_sequences = []
    for seq in sequences:
        circular_seq = seq.seq + seq.seq
        circular_record = SeqRecord(circular_seq, 
                                    id=seq.id + "|circular", 
                                    description="")
        circular_sequences.append(circular_record)
    print(f"环化了 {len(circular_sequences)} 个序列")
    return circular_sequences

def write_fasta(sequences, output_file):
    SeqIO.write(sequences, output_file, "fasta")
    print(f"写入了 {len(sequences)} 个序列到 {output_file}")

def main(input_file, output_file, read_list_file, circular_fasta_file):
    # 读取和预处理数据
    df = read_and_preprocess_data(input_file)
    
    # 处理重复的readNames
    result_df = process_repeated_readnames(df)
    
    # 将copyNum四舍五入到一位小数
    result_df['copyNum'] = result_df['copyNum'].round(1)
    
    # 添加新列并将其调整为第一列，包含copyNum
    result_df['readName|repN|consLen|copyNum'] = (
        result_df['readName'] + '|' + 
        result_df['repN'] + '|' + 
        result_df['consLen'].astype(str) + '|' + 
        result_df['copyNum'].astype(str)
    )
    columns = ['readName|repN|consLen|copyNum'] + [col for col in result_df.columns if col != 'readName|repN|consLen|copyNum']
    result_df = result_df[columns]
    
    # 保存处理后的结果
    result_df.to_csv(output_file, index=False)
    
    # 分类readNames
    read_list_df = classify_readnames(result_df)
    
    # 保存read_list结果
    read_list_df.to_csv(read_list_file, index=False)
    
    # 生成环化FASTA文件
    sequences = [SeqRecord(Seq(row['consSeq']), id=row['readName|repN|consLen|copyNum'], description="") for _, row in result_df.iterrows()]
    circular_sequences = circularize_sequences(sequences)
    
    # 将环化序列分成两组
    short_sequences = [seq for seq in circular_sequences if len(seq.seq) < 800]
    long_sequences = [seq for seq in circular_sequences if len(seq.seq) >= 800]
    
    # 自动生成短序列和长序列的文件名
    circular_fasta_file_short = f"{os.path.splitext(circular_fasta_file)[0]}_short.fasta"
    circular_fasta_file_long = f"{os.path.splitext(circular_fasta_file)[0]}_long.fasta"
    
    # 写入短序列和长序列
    write_fasta(short_sequences, circular_fasta_file_short)
    write_fasta(long_sequences, circular_fasta_file_long)
    
    print(f"处理完成。结果已保存到 {output_file} 和 {read_list_file}")
    print(f"总行数: {len(result_df)}")
    print(f"分类后的read数量: {len(read_list_df)}")
    print(f"短环化序列 (<800bp) 已保存到 {circular_fasta_file_short}")
    print(f"长环化序列 (>=800bp) 已保存到 {circular_fasta_file_long}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="处理CTCR读取结果")
    parser.add_argument("-i", "--input", required=True, help="输入文件路径")
    parser.add_argument("-o", "--output", required=True, help="处理后的结果输出文件路径")
    parser.add_argument("-r", "--read_list", required=True, help="读取列表输出文件路径")
    parser.add_argument("-c", "--circular_fasta", required=True, help="输出环化FASTA文件路径（将自动添加_short和_long后缀）")
    args = parser.parse_args()

    main(args.input, args.output, args.read_list, args.circular_fasta)

print("处理完成。")
