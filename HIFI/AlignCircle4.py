import pandas as pd
import numpy as np
import argparse
import csv
from Bio import SeqIO
import os

def read_blast_results(file_path):
    columns = ['query_id', 'subject_id', 'identity', 'alignment_length', 'mismatches', 
               'gap_opens', 'q_start', 'q_end', 's_start', 's_end', 'evalue', 'bit_score']
    df = pd.read_csv(file_path, sep='\t', header=None, names=columns)
    
    df['strand'] = np.where(df['s_start'] < df['s_end'], '+', '-')
    df.loc[df['strand'] == '-', ['s_start', 's_end']] = df.loc[df['strand'] == '-', ['s_end', 's_start']].values
    
    return df

def process_query_id(df):
    split_cols = df['query_id'].str.split('|', expand=True)
    df['readName'] = split_cols[0]
    df['repN'] = split_cols[1]
    df['consLen'] = split_cols[2].astype(int)
    df['copyNum'] = split_cols[3].astype(float)
    return df

def add_calculated_columns(df):
    df['Rlength'] = df['s_end'] - df['s_start'] + 1
    df['gap_Length'] = df['consLen'] - df['Rlength']
    df['Gap_Percentage'] = ((df['gap_Length'].abs() / df['consLen']) * 100).round(2)
    return df

def filter_low_gap_percentage(df):
    return df[(df['Gap_Percentage'] <= 1) | (df['gap_Length'].abs() <= 30)].copy()

def process_group(group):
    if len(group) == 1:
        return 'Fecc', group.index[0]
    
    start_diffs = np.abs(group['s_start'].values[:, None] - group['s_start'].values)
    end_diffs = np.abs(group['s_end'].values[:, None] - group['s_end'].values)
    
    close_pairs = (start_diffs <= 5) | (end_diffs <= 5)
    np.fill_diagonal(close_pairs, False)
    
    if close_pairs.any():
        min_gap_index = group['Gap_Percentage'].idxmin()
        return 'Fecc', min_gap_index
    else:
        return 'Mecc', None

def process_fecc(df):
    fecc_df = df[df['eClass'] == 'Fecc'].copy()
    columns = ['query_id', 'subject_id', 's_start', 's_end', 'strand', 'consLen', 'copyNum', 'Gap_Percentage', 'eClass', 'gap_Length', 'readName']
    fecc_df = fecc_df[columns].copy()

    column_mapping = {
        'subject_id': 'eChr',
        'consLen': 'eLength',
        'copyNum': 'eRepeatNum',
        'strand': 'eStrand',
        'readName': 'eReads'
    }
    fecc_df.rename(columns=column_mapping, inplace=True)

    fecc_df['eStart'] = fecc_df['s_start'] - fecc_df['gap_Length']
    fecc_df['eEnd'] = fecc_df['s_end']

    mask = fecc_df['eStart'] < 1
    fecc_df.loc[mask, 'eStart'] = fecc_df.loc[mask, 's_start']
    fecc_df.loc[mask, 'eEnd'] = fecc_df.loc[mask, 's_end'] + fecc_df.loc[mask, 'gap_Length']

    fecc_df['MatDegree'] = 100 - fecc_df['Gap_Percentage']
    fecc_df['eStart'] = fecc_df['eStart'].astype(int)
    fecc_df['eEnd'] = fecc_df['eEnd'].astype(int)
    fecc_df['eName'] = fecc_df['eChr'] + "-" + fecc_df['eStart'].astype(str) + "-" + fecc_df['eEnd'].astype(str)

    result_columns = ['eName', 'eChr', 'eStart', 'eEnd', 'eStrand', 'eLength', 'eRepeatNum', 'MatDegree', 'eClass', 'eReads', 'query_id']
    return fecc_df[result_columns].copy()

def generate_align_region(row):
    if row['s_start'] - row['gap_Length'] < 1:
        return f"{row['identity']}|{row['subject_id']}-{row['s_start']}-{row['s_end'] + row['gap_Length']}"
    else:
        return f"{row['identity']}|{row['subject_id']}-{row['s_start'] - row['gap_Length']}-{row['s_end']}"

def process_mecc_group(group):
    q_start_mode = group['q_start'].mode().iloc[0]
    q_start_row = group[group['q_start'] == q_start_mode].iloc[0]
    
    q_end_mode = group['q_end'].mode().iloc[0]
    q_end_row = group[group['q_end'] == q_end_mode].iloc[0]
    
    align_regions = group.apply(generate_align_region, axis=1).tolist()
    
    return pd.Series({
        'query_id': group['query_id'].iloc[0],
        'q_start': q_start_row['q_start'],
        'q_end': q_end_row['q_end'],
        'copyNum': group['copyNum'].iloc[0],
        'consLen': group['consLen'].iloc[0],
        'readName': group['readName'].iloc[0],
        'eClass': group['eClass'].iloc[0],
        'AlignRegion': ';'.join(align_regions),
        'MatchNum': len(group)
    })

def process_mecc(df):
    mecc_df = df[df['eClass'] == 'Mecc'].copy()
    result_df = mecc_df.groupby('query_id').apply(process_mecc_group).reset_index(drop=True)

    result_df['gap_Length_2'] = result_df['consLen'] - (result_df['q_end'] - result_df['q_start'] + 1)
    result_df['tStart'] = result_df['q_start'] - result_df['gap_Length_2']
    result_df['tEnd'] = result_df['q_end']

    mask = result_df['tStart'] < 1
    result_df.loc[mask, 'tStart'] = result_df.loc[mask, 'q_start']
    result_df.loc[mask, 'tEnd'] = result_df.loc[mask, 'q_end'] + result_df.loc[mask, 'gap_Length_2']

    result_df['tStart'] = result_df['tStart'].astype(int)
    result_df['tEnd'] = result_df['tEnd'].astype(int)

    return result_df

def extract_query_ids(fecc_df, mecc_df):
    fecc_ids = set(fecc_df['query_id'])
    mecc_ids = set(mecc_df['query_id'])
    return fecc_ids, mecc_ids

def process_sequences(input_fasta, fecc_ids, mecc_ids, mecc_df):
    tpm_fecc = []
    tpm_mecc = []
    uncircle = []
    multi_mapped = {}
    
    # 从mecc_df获取q_start和q_end信息
    mecc_info = dict(zip(mecc_df['query_id'], zip(mecc_df['q_start'], mecc_df['q_end'])))
    
    for record in SeqIO.parse(input_fasta, "fasta"):
        if record.id in fecc_ids:
            tpm_fecc.append(record)
        elif record.id in mecc_ids:
            tpm_mecc.append(record)
            if record.id in mecc_info:
                start, end = mecc_info[record.id]
                multi_mapped[record.id] = record[start-1:end]  # 调整为0基索引
        else:
            # 只保留前半序列
            half_length = len(record.seq) // 2
            uncircle.append(record[:half_length])
    
    # 写入序列到文件
    SeqIO.write(tpm_fecc, "tpm.fecc.fa", "fasta")
    SeqIO.write(tpm_mecc, "tpm.mecc.fa", "fasta")
    SeqIO.write(uncircle, "UnCircle.fa", "fasta")
    SeqIO.write(multi_mapped.values(), "Multi-mapped.eccDNA.fa", "fasta")
    
    # 删除临时文件
    os.remove("tpm.fecc.fa")
    os.remove("tpm.mecc.fa")

def main(input_file, fecc_output, mecc_output, input_fasta):
    # 读取和处理数据
    df = read_blast_results(input_file)
    df = process_query_id(df)
    df = add_calculated_columns(df)
    
    # 过滤低Gap百分比的数据
    low_gap_df = filter_low_gap_percentage(df)
    
    # 添加新列
    low_gap_df['occurrence_count'] = low_gap_df.groupby('query_id')['query_id'].transform('count')
    low_gap_df['eClass'] = ''
    low_gap_df.loc[low_gap_df['occurrence_count'] == 1, 'eClass'] = 'Fecc'

    # 处理重复组
    duplicate_mask = low_gap_df['occurrence_count'] > 1
    duplicate_groups = low_gap_df[duplicate_mask].groupby('query_id')

    results = []
    for name, group in duplicate_groups:
        eclass, keep_index = process_group(group)
        if eclass == 'Fecc':
            results.append((name, eclass, [keep_index]))
        else:
            results.append((name, eclass, group.index.tolist()))

    # 更新eClass并过滤数据
    for name, eclass, indices in results:
        low_gap_df.loc[low_gap_df['query_id'] == name, 'eClass'] = eclass
        if eclass == 'Fecc':
            low_gap_df = low_gap_df[~((low_gap_df['query_id'] == name) & (~low_gap_df.index.isin(indices)))]

    # 删除不需要的列
    low_gap_df = low_gap_df.drop('occurrence_count', axis=1)

    # 处理Fecc和Mecc数据
    fecc_results = process_fecc(low_gap_df)
    fecc_results.to_csv(fecc_output, index=False)

    mecc_results = process_mecc(low_gap_df)
    mecc_results.to_csv(mecc_output, index=False)

    # 新增：处理序列
    fecc_ids, mecc_ids = extract_query_ids(fecc_results, mecc_results)
    process_sequences(input_fasta, fecc_ids, mecc_ids, mecc_results)

    print("处理完成。输出文件：Multi-mapped.eccDNA.fa 和 UnCircle.fa")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='处理BLAST结果并生成Fecc和Mecc输出文件')
    parser.add_argument('-i', '--input', required=True, help='输入BLAST结果文件路径')
    parser.add_argument('-f', '--fecc_output', default='processed_Fecc_results.csv', help='Fecc结果输出文件路径 (默认: processed_Fecc_results.csv)')
    parser.add_argument('-m', '--mecc_output', default='processed_Mecc_results.csv', help='Mecc结果输出文件路径 (默认: processed_Mecc_results.csv)')
    parser.add_argument('-fa', '--input_fasta', required=True, help='输入FASTA文件路径')
    
    args = parser.parse_args()
    
    main(args.input, args.fecc_output, args.mecc_output, args.input_fasta)
