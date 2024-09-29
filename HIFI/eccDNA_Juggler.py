#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process BLAST results and generate analysis outputs.')
    parser.add_argument('-i', '--input_blast', required=True, help='Input BLAST results file')
    parser.add_argument('-f', '--input_fai', required=True, help='Input FASTA index file')
    parser.add_argument('-n', '--output_num_linr', required=True, help='Output Num_LinR CSV file')
    parser.add_argument('-o', '--output_tecc', required=True, help='Output tecc_analysis_results CSV file')
    parser.add_argument('-r', '--output_read_class', required=True, help='Output read classification CSV file')
    args = parser.parse_args()
    return args

def read_and_merge_data(input_blast, input_fai):
    # 读取 BLAST 结果文件
    df1 = pd.read_csv(input_blast, sep='\s+', header=None,
                      names=['qname', 'tname', 'identity', 'alnlen', 'mismatch', 'gapopen',
                             'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits'])

    # 读取 FASTA 索引文件
    df2 = pd.read_csv(input_fai, sep='\s+', header=None,
                      names=['qname', 'length', 'offset', 'linebases', 'linewidth'])

    # 合并两个 DataFrame，使用 'qname' 作为合并键
    merged_df = pd.merge(df1, df2, on='qname', how='left')

    # 添加 Strand 列
    merged_df['Strand'] = np.where(merged_df['tstart'] < merged_df['tend'], '+', '-')

    # 确保 tstart < tend
    merged_df.loc[merged_df['tstart'] > merged_df['tend'], ['tstart', 'tend']] = \
        merged_df.loc[merged_df['tstart'] > merged_df['tend'], ['tend', 'tstart']].values

    # 计算 eff_Ratio 并保留两位小数
    merged_df['eff_Ratio'] = (merged_df['alnlen'] / merged_df['length'] * 100).round(2)

    # 统计 qname 的频数
    qname_counts = merged_df['qname'].value_counts().reset_index()
    qname_counts.columns = ['qname', 'frequency']

    # 将频数信息合并回原 DataFrame
    merged_df = pd.merge(merged_df, qname_counts, on='qname', how='left')

    return merged_df

def process_single_occurrence(merged_df, output_num_linr):
    # 筛选出 qname 仅出现一次的行
    single_occurrence = merged_df[merged_df['frequency'] == 1]

    # 在单次出现的行中，保留 eff_Ratio 大于 98 的行
    filtered_single_df = single_occurrence[single_occurrence['eff_Ratio'] > 98]

    # 计算被保留的行数
    kept_rows = len(filtered_single_df)

    # 将 Num_LinR 输出到 CSV 文件
    num_linr_df = pd.DataFrame({'Num_LinR': [kept_rows]})
    num_linr_df.to_csv(output_num_linr, index=False)

    return filtered_single_df

def check_overlap(group):
    if len(group) != 2:
        return None

    # 检查 tname 是否一致
    if group['tname'].nunique() != 1:
        return None

    tstart_min = group['tstart'].min()
    tend_max = group['tend'].max()

    # 检查是否有重叠
    if max(group['tstart']) <= min(group['tend']):
        overlap_length = min(group['tend']) - max(group['tstart'])
        length = tend_max - tstart_min + 1

        # 检查 overlap_length 是否小于 length
        if overlap_length >= length:
            return None

        strand = '+' if group['Strand'].nunique() == 1 else '+'

        return pd.Series({
            'qname': group['qname'].iloc[0],
            'PreCtcRlist_number': 0,
            'eChr': group['tname'].iloc[0],
            'eStart': tstart_min,
            'eEnd': tend_max,
            'eStrand': strand,
            'eLength': length,
            'eRepeatNum': 1,
            'MatDegree': round(group['identity'].mean(), 2),
            'eClass': 'Tecc',
            'eReads': group['qname'].iloc[0],
            'query_id': 'blast2',
            'inversion': 0  # 对于两行的 Blast 结果，inversion 列填充为 0
        })
    return None

def process_two_line_blast(merged_df):
    # 筛选出 qname 出现 2 次的行
    df_twice = merged_df[merged_df['frequency'] == 2]

    # 按 qname 分组并应用 check_overlap 函数
    two_line_blast_df = df_twice.groupby('qname').apply(check_overlap).reset_index(drop=True)

    # 删除 None 结果（即没有重叠的组）
    two_line_blast_df = two_line_blast_df.dropna()

    return two_line_blast_df

def calculate_score(row1, row2):
    """计算两行之间的得分"""
    q_diff = min(abs(row1['qstart'] - row2['qstart']), abs(row1['qend'] - row2['qend']))
    subject_score = 100 if row1['tname'] == row2['tname'] else 0
    identity_score = min(row1['identity'], row2['identity'])
    s_diff = min(abs(row1['tstart'] - row2['tstart']), abs(row1['tend'] - row2['tend']))
    s_diff_score = 1 / (s_diff + 1)  # 避免除以零

    # 综合得分，权重可以根据需求调整
    return q_diff * 0.4 + subject_score * 0.3 + identity_score * 0.2 + s_diff_score * 0.1

def analyze_group(group):
    """分析每个 qname 组，构建比对链"""
    group = group.sort_values('qstart')
    PreCtcRlists = []
    processed = set()
    list_number = 1

    for i, row in group.iterrows():
        if i in processed:
            continue

        current_chain = [row]
        processed.add(i)

        while True:
            potential_next = group[
                (group['qstart'] >= current_chain[-1]['qend'] - 10) &
                (group['qstart'] <= current_chain[-1]['qend'] + 10) &
                (~group.index.isin(processed))
            ]

            if potential_next.empty:
                break

            if len(potential_next) > 1:
                scores = potential_next.apply(lambda x: calculate_score(current_chain[-1], x), axis=1)
                next_row = potential_next.loc[scores.idxmax()]
            else:
                next_row = potential_next.iloc[0]

            current_chain.append(next_row)
            processed.add(next_row.name)

        # 为当前链中的每一行添加 PreCtcRlist 编号
        chain_df = pd.DataFrame(current_chain)
        chain_df['PreCtcRlist_number'] = list_number
        PreCtcRlists.append(chain_df)
        list_number += 1

    return pd.concat(PreCtcRlists, ignore_index=True)

def process_multi_occurrence(merged_df):
    # 筛选出 qname 出现大于 2 次小于 100 次的行
    multi_occurrence_df = merged_df[(merged_df['frequency'] > 2) & (merged_df['frequency'] < 100)]

    # 按 qname 分组并应用 analyze_group 函数
    prectcrlist_df = multi_occurrence_df.groupby('qname', group_keys=False).apply(analyze_group).reset_index(drop=True)

    return prectcrlist_df

def merge_overlapping_lists(prectcrlist_df):
    def merge_lists(group):
        prectcrlists = group.groupby('PreCtcRlist_number')
        merged = []
        merged_numbers = set()

        for num1, list1 in prectcrlists:
            if num1 in merged_numbers or len(list1) < 3:
                continue

            merged_list = list1.copy()
            merged_numbers.add(num1)

            for num2, list2 in prectcrlists:
                if num2 in merged_numbers or len(list2) < 3:
                    continue

                if set(list1['tstart']).intersection(set(list2['tstart'])) or \
                   set(list1['tend']).intersection(set(list2['tend'])):
                    merged_list = pd.concat([merged_list, list2])
                    merged_numbers.add(num2)

            if len(merged_list) >= 3:
                merged_list['PreCtcRlist_number'] = f"{num1}m" if len(merged_numbers) > 1 else num1
                merged.append(merged_list)

        return pd.concat(merged) if merged else pd.DataFrame()

    # 首先过滤掉所有短于 3 行的 PreCtcRlist
    prectcrlist_df = prectcrlist_df.groupby(['qname', 'PreCtcRlist_number']).filter(lambda x: len(x) >= 3)

    # 筛选出有多个 PreCtcRlist_number 的 qname
    multi_list_qnames = prectcrlist_df.groupby('qname')['PreCtcRlist_number'].nunique()
    multi_list_qnames = multi_list_qnames[multi_list_qnames > 1].index

    # 对这些 qname 进行处理
    processed_df = prectcrlist_df[prectcrlist_df['qname'].isin(multi_list_qnames)].groupby('qname').apply(merge_lists)
    processed_df = processed_df.reset_index(drop=True)

    # 将未处理的行（只有一个 PreCtcRlist_number 的 qname）添加回结果
    single_list_df = prectcrlist_df[~prectcrlist_df['qname'].isin(multi_list_qnames)]
    final_prectcrlist_df = pd.concat([processed_df, single_list_df]).reset_index(drop=True)

    return final_prectcrlist_df

def process_group(group):
    # 确保 tname 一致
    if len(group['tname'].unique()) > 1:
        return None

    # 检查 tstart 和 tend 列中的重复
    def check_repeats(series, tolerance=30):
        sorted_vals = sorted(series)
        repeats = []
        for i in range(len(sorted_vals) - 1):
            if abs(sorted_vals[i] - sorted_vals[i + 1]) <= tolerance:
                repeats.extend([sorted_vals[i], sorted_vals[i + 1]])
        return list(set(repeats))

    tstart_repeats = check_repeats(group['tstart'])
    tend_repeats = check_repeats(group['tend'])

    # 标记重复的行
    group['Repeat1'] = group['tstart'].apply(lambda x: x in tstart_repeats)
    group['Repeat2'] = group['tend'].apply(lambda x: x in tend_repeats)

    # 检查是否有同时被标记为 Repeat1 和 Repeat2 的行
    repeat_rows = group[group['Repeat1'] & group['Repeat2']]
    if not repeat_rows.empty:
        # 计算出现次数最多的 Strand 值
        most_common_strand = repeat_rows['Strand'].mode().iloc[0]

        # 检查 Strand 是否一致
        inversion_flag = 1 if group['Strand'].nunique() > 1 else 0

        # 返回一个 DataFrame，注意将标量值包装在列表中
        return pd.DataFrame({
            'qname': [group['qname'].iloc[0]],
            'PreCtcRlist_number': [group['PreCtcRlist_number'].iloc[0]],
            'eChr': [group['tname'].iloc[0]],
            'eStart': [min(repeat_rows['tstart'])],
            'eEnd': [max(repeat_rows['tend'])],
            'eStrand': [most_common_strand],
            'eLength': [max(repeat_rows['tend']) - min(repeat_rows['tstart']) + 1],
            'eRepeatNum': [len(repeat_rows)],
            'MatDegree': [round(repeat_rows['identity'].mean(), 2)],
            'eClass': ['Tecc'],
            'eReads': [group['qname'].iloc[0]],
            'query_id': ['blast'],
            'inversion': [inversion_flag]
        })
    return None

def process_final_prectcrlist(final_prectcrlist_df):
    # 按照 qname 和 PreCtcRlist_number 分组并应用处理函数
    result_dfs = []
    for (qname, prectc), group in final_prectcrlist_df.groupby(['qname', 'PreCtcRlist_number']):
        result = process_group(group)
        if result is not None:
            result_dfs.append(result)

    # 合并所有结果
    if result_dfs:
        output_df = pd.concat(result_dfs, ignore_index=True)
    else:
        output_df = pd.DataFrame()

    return output_df

def generate_read_classification(output_df, output_read_class):
    # 初始化一个空的列表来存储 Read_df 的数据
    read_classification = []

    # 按 qname 分组
    grouped = output_df.groupby('qname')

    for qname, group in grouped:
        prectc_numbers = group['PreCtcRlist_number'].unique()
        inversions = group['inversion'].unique()
        eLengths = group['eLength'].unique()

        # 检查 PreCtcRlist_number 是否包含 'm'
        if any('m' in str(num) for num in prectc_numbers):
            read_classification.append({'readName': qname, 'readClass': 'CtcR-hybrid'})
            output_df.loc[output_df['qname'] == qname, 'eClass'] = 'Uecc'

        # 如果只有一个 PreCtcRlist_number
        elif len(prectc_numbers) == 1:
            if inversions[0] == 0:
                read_classification.append({'readName': qname, 'readClass': 'CtcR-perfect'})
                output_df.loc[output_df['qname'] == qname, 'eClass'] = 'Uecc'
            else:
                read_classification.append({'readName': qname, 'readClass': 'CtcR-inversion'})
                output_df.loc[output_df['qname'] == qname, 'eClass'] = 'Uecc'

        # 如果有多个 PreCtcRlist_number
        else:
            # 检查 eLength 是否在 50bp 内相等
            max_length = max(eLengths)
            min_length = min(eLengths)
            if max_length - min_length <= 50:
                read_classification.append({'readName': qname, 'readClass': 'CtcR-perfect'})
                output_df.loc[output_df['qname'] == qname, 'eClass'] = 'Mecc'
            else:
                read_classification.append({'readName': qname, 'readClass': 'CtcR-multiple'})
                output_df.loc[output_df['qname'] == qname, 'eClass'] = 'Uecc'

    # 将 read_classification 转换为 DataFrame
    read_df = pd.DataFrame(read_classification)

    # 保存 Read_df 到 CSV 文件，使用逗号分隔
    read_df.to_csv(output_read_class, index=False)

    return output_df

def main():
    args = parse_arguments()

    # 读取并合并数据
    merged_df = read_and_merge_data(args.input_blast, args.input_fai)

    # 处理单次出现的 qname
    filtered_single_df = process_single_occurrence(merged_df, args.output_num_linr)

    # 处理出现两次的 qname
    two_line_blast_df = process_two_line_blast(merged_df)

    # 处理多次出现的 qname
    prectcrlist_df = process_multi_occurrence(merged_df)

    # 合并重叠的列表
    final_prectcrlist_df = merge_overlapping_lists(prectcrlist_df)

    # 处理最终的 PreCtcRlist 数据
    output_df = process_final_prectcrlist(final_prectcrlist_df)

    # 将 two_line_blast_df 追加到 output_df
    if not two_line_blast_df.empty:
        output_df = pd.concat([output_df, two_line_blast_df], ignore_index=True)

    # 生成 Read_df 并更新 output_df
    output_df = generate_read_classification(output_df, args.output_read_class)

    # 保存最终的 output_df 到 CSV 文件
    output_df.to_csv(args.output_tecc, index=False)

if __name__ == '__main__':
    main()
