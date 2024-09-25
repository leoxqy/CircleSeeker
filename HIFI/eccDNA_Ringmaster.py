#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import argparse
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
    return df[(df['Gap_Percentage'] <= 10) | (df['gap_Length'].abs() <= 50)].copy()

def process_group(group):
    if len(group) == 1:
        return 'Uecc', group.index[0]

    start_diffs = np.abs(group['s_start'].values[:, None] - group['s_start'].values)
    end_diffs = np.abs(group['s_end'].values[:, None] - group['s_end'].values)

    close_pairs = (start_diffs <= 5) | (end_diffs <= 5)
    np.fill_diagonal(close_pairs, False)

    if close_pairs.any():
        min_gap_index = group['Gap_Percentage'].idxmin()
        return 'Uecc', min_gap_index
    else:
        return 'Mecc', None

def process_Uecc(df):
    Uecc_df = df[df['eClass'] == 'Uecc'].copy()
    columns = ['query_id', 'subject_id', 's_start', 's_end', 'strand', 'consLen', 'copyNum', 'Gap_Percentage', 'eClass', 'gap_Length', 'readName']
    Uecc_df = Uecc_df[columns].copy()

    column_mapping = {
        'subject_id': 'eChr',
        'consLen': 'eLength',
        'copyNum': 'eRepeatNum',
        'strand': 'eStrand',
        'readName': 'eReads'
    }
    Uecc_df.rename(columns=column_mapping, inplace=True)

    Uecc_df['eStart'] = Uecc_df['s_start'] - Uecc_df['gap_Length']
    Uecc_df['eEnd'] = Uecc_df['s_end']

    mask = Uecc_df['eStart'] < 1
    Uecc_df.loc[mask, 'eStart'] = Uecc_df.loc[mask, 's_start']
    Uecc_df.loc[mask, 'eEnd'] = Uecc_df.loc[mask, 's_end'] + Uecc_df.loc[mask, 'gap_Length']

    Uecc_df['MatDegree'] = 100 - Uecc_df['Gap_Percentage']
    Uecc_df['eStart'] = Uecc_df['eStart'].astype(int)
    Uecc_df['eEnd'] = Uecc_df['eEnd'].astype(int)
    Uecc_df['eName'] = Uecc_df['eChr'] + "-" + Uecc_df['eStart'].astype(str) + "-" + Uecc_df['eEnd'].astype(str)

    result_columns = ['eName', 'eChr', 'eStart', 'eEnd', 'eStrand', 'eLength', 'eRepeatNum', 'MatDegree', 'eClass', 'eReads', 'query_id']
    return Uecc_df[result_columns].copy()

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

def extract_query_ids(Uecc_df, mecc_df, Cecc_df):
    Uecc_ids = set(Uecc_df['query_id'])
    mecc_ids = set(mecc_df['query_id'])
    Cecc_ids = set(Cecc_df['query_id'])
    return Uecc_ids, mecc_ids, Cecc_ids

def process_sequences(input_fasta, Uecc_ids, mecc_ids, Cecc_ids, mecc_df, Cecc_df, mecc_fa='MeccDNA.fa', cecc_fa='CeccDNA.fa', xecc_fa='XeccDNA.fa'):
    tpm_Uecc = []
    tpm_mecc = []
    tpm_Cecc = []
    xecc_dna = []
    mecc_mapped = []
    Cecc_mapped = []

    # 从mecc_df获取q_start和q_end信息
    mecc_info = dict(zip(mecc_df['query_id'], zip(mecc_df['q_start'], mecc_df['q_end'])))
    Cecc_info = dict(zip(Cecc_df['query_id'], zip(Cecc_df['q_start'], Cecc_df['q_end'])))

    for record in SeqIO.parse(input_fasta, "fasta"):
        if record.id in Uecc_ids:
            tpm_Uecc.append(record)
        elif record.id in mecc_ids:
            if record.id in mecc_info:
                start, end = mecc_info[record.id]
                sub_record = record[start-1:end]  # 调整为0基索引
                sub_record.id = record.id
                sub_record.description = record.description
                mecc_mapped.append(sub_record)
        elif record.id in Cecc_ids:
            if record.id in Cecc_info:
                start, end = Cecc_info[record.id]
                sub_record = record[start-1:end]
                sub_record.id = record.id
                sub_record.description = record.description
                tpm_Cecc.append(sub_record)
        else:
            # 只保留序列的前一半
            half_length = len(record.seq) // 2
            sub_record = record[:half_length]
            sub_record.id = record.id
            sub_record.description = record.description
            xecc_dna.append(sub_record)

    # 写入序列到文件
    SeqIO.write(mecc_mapped, mecc_fa, "fasta")
    SeqIO.write(tpm_Cecc, cecc_fa, "fasta")
    SeqIO.write(xecc_dna, xecc_fa, "fasta")

def process_other_blast_rows(df):
    # 筛选既不是Uecc也不是Mecc的行
    other_df = df[~df['eClass'].isin(['Uecc', 'Mecc'])].copy()

    # 按照 query_id 分组，删除行数超过100的组
    group_sizes = other_df.groupby('query_id').size()
    query_ids_to_remove = group_sizes[group_sizes > 100].index

    # 删除这些 query_id
    other_df = other_df[~other_df['query_id'].isin(query_ids_to_remove)]

    return other_df

def calculate_score(row1, row2):
    q_diff = min(abs(row1['q_start'] - row2['q_start']), abs(row1['q_end'] - row2['q_end']))
    subject_score = 100 if row1['subject_id'] == row2['subject_id'] else 0
    identity_score = min(row1['identity'], row2['identity'])
    s_diff = min(abs(row1['s_start'] - row2['s_start']), abs(row1['s_end'] - row2['s_end']))

    s_diff_score = 1 / (s_diff + 1)

    return q_diff * 0.4 + subject_score * 0.3 + identity_score * 0.2 + s_diff_score * 0.1

def analyze_group(group):
    group = group.sort_values('q_start')
    Prelists = []
    Prelist_summaries = []
    processed = set()

    query_id = group['query_id'].iloc[0]
    cons_len = group['consLen'].iloc[0]  # Assuming consLen is the same for all rows in a group

    for i, row in group.iterrows():
        if i in processed:
            continue

        current_chain = [row]
        processed.add(i)

        while True:
            potential_next = group[(group['q_start'] >= current_chain[-1]['q_end'] - 10) &
                                   (group['q_start'] <= current_chain[-1]['q_end'] + 10) &
                                   (~group.index.isin(processed))]

            if potential_next.empty:
                break

            if len(potential_next) > 1:
                scores = potential_next.apply(lambda x: calculate_score(current_chain[-1], x), axis=1)
                next_row = potential_next.loc[scores.idxmax()]
            else:
                next_row = potential_next.iloc[0]

            current_chain.append(next_row)
            processed.add(next_row.name)

        Prelist = pd.DataFrame(current_chain)
        start_q_start = Prelist['q_start'].min()
        end_q_end = Prelist['q_end'].max()
        pre_length = end_q_end - start_q_start + 1

        # 只保留符合条件的Prelist
        if pre_length > cons_len and len(Prelist) > 1:
            # 将第一行移到最后
            Prelist = pd.concat([Prelist.iloc[1:], Prelist.iloc[:1]]).reset_index(drop=True)

            Prelist['Prelist_ID'] = len(Prelists) + 1
            Prelists.append(Prelist)

            # 创建Prelist摘要
            summary = {
                'query_id': query_id,
                'Prelist_ID': len(Prelists),
                'start_q_start': start_q_start,
                'end_q_end': end_q_end,
                'PreLength': pre_length,
                'consLen': cons_len,
                'num_rows': len(Prelist)
            }
            Prelist_summaries.append(summary)

    if Prelists:
        return pd.concat(Prelists, ignore_index=True), pd.DataFrame(Prelist_summaries)
    else:
        return pd.DataFrame(), pd.DataFrame()

def process_detailed_results(df):
    # 根据 query_id 和 Prelist_ID 分组
    grouped = df.groupby(['query_id', 'Prelist_ID'])

    results = []

    for (query_id, prelist_id), group in grouped:
        # 标记 s_start 重复的行
        group['Repeat1'] = group.duplicated(subset=['s_start'], keep=False)

        # 标记 s_end 重复的行
        group['Repeat2'] = group.duplicated(subset=['s_end'], keep=False)

        # 提取被标记为 Repeat1 或 Repeat2 的行
        repeated_rows = group[(group['Repeat1'] | group['Repeat2'])].copy()

        if not repeated_rows.empty:
            # 为每个组创建 LocNum
            repeated_rows.loc[:, 'LocNum'] = range(1, len(repeated_rows) + 1)

            # 计算 s_length
            repeated_rows.loc[:, 's_length'] = repeated_rows['s_end'] - repeated_rows['s_start'] + 1

            # 选择需要的列
            selected_columns = ['Prelist_ID', 'query_id', 'q_start', 'q_end', 'consLen', 'LocNum',
                                'identity', 'alignment_length', 'subject_id', 's_start', 's_end',
                                's_length', 'strand', 'readName', 'repN', 'copyNum']

            results.append(repeated_rows[selected_columns])

    if results:
        final_df = pd.concat(results, ignore_index=True)
        return final_df
    else:
        return pd.DataFrame()

def process_detailed_results_further(df):
    # 按 query_id 和 Prelist_ID 分组
    grouped = df.groupby(['query_id', 'Prelist_ID'])

    results = []

    for (query_id, prelist_id), group in grouped:
        # 如果组内行数小于等于2，直接跳过这个组
        if len(group) <= 2:
            continue

        # 删除最后两行
        group = group.iloc[:-2]

        # 按 s_start 和 s_end 去重
        group = group.drop_duplicates(subset=['s_start', 's_end'])

        # 计算 s_length 总和
        total_s_length = group['s_length'].sum()

        # 获取 consLen（假设每个组内 consLen 相同）
        cons_len = group['consLen'].iloc[0]

        # 计算比率
        ratio = total_s_length / cons_len

        # 如果比率在 0.8-1.2 之间，保留这个组
        if 0.8 <= ratio <= 1.2:
            # 添加比率列
            group['ratio'] = ratio
            results.append(group)

    if results:
        final_df = pd.concat(results, ignore_index=True)

        # 计算每个 query_id 对应的 Prelist_ID 数量
        prelist_counts = final_df.groupby('query_id')['Prelist_ID'].nunique()

        # 更新 eClass 列
        final_df['eClass'] = final_df.apply(lambda row: 'Cecc' if prelist_counts[row['query_id']] == 1 else 'Ceccm', axis=1)

        return final_df
    else:
        return pd.DataFrame()

def main(input_file, Uecc_output_csv, mecc_output_csv, Cecc_output_csv, input_fasta, mecc_fa, cecc_fa, xecc_fa):
    # 读取和处理数据
    df = read_blast_results(input_file)
    df = process_query_id(df)
    df = add_calculated_columns(df)

    # 过滤低Gap百分比的数据
    low_gap_df = filter_low_gap_percentage(df)

    # 添加新列
    low_gap_df['occurrence_count'] = low_gap_df.groupby('query_id')['query_id'].transform('count')
    low_gap_df['eClass'] = ''
    low_gap_df.loc[low_gap_df['occurrence_count'] == 1, 'eClass'] = 'Uecc'

    # 处理重复组
    duplicate_mask = low_gap_df['occurrence_count'] > 1
    duplicate_groups = low_gap_df[duplicate_mask].groupby('query_id')

    results = []
    for name, group in duplicate_groups:
        eclass, keep_index = process_group(group)
        if eclass == 'Uecc':
            results.append((name, eclass, [keep_index]))
        else:
            results.append((name, eclass, group.index.tolist()))

    # 更新eClass并过滤数据
    for name, eclass, indices in results:
        low_gap_df.loc[low_gap_df['query_id'] == name, 'eClass'] = eclass
        if eclass == 'Uecc':
            low_gap_df = low_gap_df[~((low_gap_df['query_id'] == name) & (~low_gap_df.index.isin(indices)))]

    # 删除不需要的列
    low_gap_df = low_gap_df.drop('occurrence_count', axis=1)

    # 将 eClass 信息合并回原始的 df
    eClass_mapping = low_gap_df[['query_id', 'eClass']].drop_duplicates()
    df = df.merge(eClass_mapping, on='query_id', how='left')

    # 处理Uecc和Mecc数据
    Uecc_results = process_Uecc(low_gap_df)
    Uecc_results.to_csv(Uecc_output_csv, index=False)

    mecc_results = process_mecc(low_gap_df)
    mecc_results.to_csv(mecc_output_csv, index=False)

    # 处理其他blast行，使用合并了 eClass 的原始 df
    other_results = process_other_blast_rows(df)

    # 如果other_results为空，保留所有不是Uecc和Mecc的行
    if other_results.empty:
        print("其他BLAST结果为空。")
        Cecc_results = pd.DataFrame()
    else:
        # 分析other_results，生成CeccDNA数据
        detailed_results_list = []
        summary_results_list = []

        grouped = other_results.groupby('query_id')
        for name, group in grouped:
            detailed, summary = analyze_group(group)
            if not detailed.empty:
                detailed_results_list.append(detailed)
                summary_results_list.append(summary)

        if detailed_results_list:
            detailed_results = pd.concat(detailed_results_list, ignore_index=True)
            processed_detailed_results = process_detailed_results(detailed_results)
            if not processed_detailed_results.empty:
                Cecc_results = process_detailed_results_further(processed_detailed_results)
                if not Cecc_results.empty:
                    Cecc_results.to_csv(Cecc_output_csv, index=False)
                else:
                    Cecc_results = pd.DataFrame()
            else:
                Cecc_results = pd.DataFrame()
        else:
            Cecc_results = pd.DataFrame()

    # 处理序列
    if not Cecc_results.empty:
        Cecc_ids = set(Cecc_results['query_id'])
    else:
        Cecc_ids = set()

    Uecc_ids, mecc_ids, Cecc_ids = extract_query_ids(Uecc_results, mecc_results, Cecc_results)
    process_sequences(input_fasta, Uecc_ids, mecc_ids, Cecc_ids, mecc_results, Cecc_results, mecc_fa, cecc_fa, xecc_fa)

    print(f"处理完成。输出文件：{Uecc_output_csv}、{mecc_output_csv}、{Cecc_output_csv}、{mecc_fa}、{cecc_fa} 和 {xecc_fa}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='处理BLAST结果并生成输出文件')
    parser.add_argument('-i', '--input', required=True, help='输入BLAST结果文件路径')
    parser.add_argument('-fa', '--input_fasta', required=True, help='输入FASTA文件路径')
    parser.add_argument('-u', '--Uecc_output_csv', default='UeccDNA.csv', help='Uecc输出CSV文件路径')
    parser.add_argument('-m', '--Mecc_output_csv', default='MeccDNA.csv', help='Mecc输出CSV文件路径')
    parser.add_argument('-c', '--Cecc_output_csv', default='CeccDNA.csv', help='Cecc输出CSV文件路径')
    parser.add_argument('--mecc_fa', default='MeccDNA.fa', help='Mecc输出FASTA文件路径')
    parser.add_argument('--cecc_fa', default='CeccDNA.fa', help='Cecc输出FASTA文件路径')
    parser.add_argument('--xecc_fa', default='XeccDNA.fa', help='Xecc输出FASTA文件路径')

    args = parser.parse_args()

    main(args.input, args.Uecc_output_csv, args.Mecc_output_csv, args.Cecc_output_csv, args.input_fasta, args.mecc_fa, args.cecc_fa, args.xecc_fa)
