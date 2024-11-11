import argparse
import logging
import os
import shutil
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

class MergeMecc:
    def __init__(self, meccdna_part1, tecc_analysis, mecc_part1, mecc_part2, output_csv, output_fasta):
        self.meccdna_part1 = meccdna_part1
        self.tecc_analysis = tecc_analysis
        self.mecc_part1 = mecc_part1
        self.mecc_part2 = mecc_part2
        self.output_csv = output_csv
        self.output_fasta = output_fasta
        
        # 设置日志
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )

    def run_merge_mecc(self):
        """执行完整的合并流程"""
        # 读取输入文件
        logging.info("Reading input files...")
        df1 = self.read_csv(self.meccdna_part1)
        df2 = self.read_csv(self.tecc_analysis)
        
        # 处理数据
        logging.info("Processing MeccDNA part1 data...")
        df1_processed = self.process_meccdna_part1(df1)
        
        logging.info("Processing Tecc analysis data...")
        df2_processed = self.process_tecc_analysis(df2)
        
        # 合并数据框
        logging.info("Combining processed data...")
        combined_df = pd.concat([df1_processed, df2_processed], ignore_index=True)
        
        # 处理查询组
        logging.info("Processing query groups...")
        query_groups = self.process_query_groups(combined_df)
        
        logging.info("Merging similar feature groups...")
        merged_query_groups = self.merge_same_feature_groups(query_groups)
        
        logging.info("Integrating data...")
        integrated_df = self.integrate_data(merged_query_groups, combined_df)
        
        logging.info("Adding names and states...")
        final_df = self.add_name_and_state(integrated_df)
        
        # 处理FASTA文件
        logging.info("Processing FASTA files...")
        self.merge_fasta_files(self.mecc_part1, self.mecc_part2, "temp_merged.fa")
        
        if not final_df.empty:
            name_map = dict(zip(final_df['query_id'].unique(), 
                              final_df.groupby('query_id')['eName'].first()))
            self.update_fasta_names("temp_merged.fa", name_map, self.output_fasta)
        
        # 清理临时文件
        if os.path.exists("temp_merged.fa"):
            os.remove("temp_merged.fa")
        
        # 准备并保存最终结果
        logging.info("Preparing final results...")
        final_results = self.prepare_final_dataframe(final_df)
        final_results.to_csv(self.output_csv, index=False)
        
        logging.info("Processing completed successfully.")

    @staticmethod
    def is_file_empty(file_path):
        """检查文件是否为空或不存在"""
        return not os.path.exists(file_path) or os.path.getsize(file_path) == 0

    def read_csv(self, file_path):
        """读取CSV文件并处理空文件情况"""
        if self.is_file_empty(file_path):
            logging.warning(f"File {file_path} is empty or doesn't exist. Creating empty DataFrame.")
            return pd.DataFrame()
        return pd.read_csv(file_path)

    def process_meccdna_part1(self, df):
        """处理meccdna_part1数据"""
        if df.empty:
            logging.warning("MeccDNA part1 DataFrame is empty. Returning empty DataFrame.")
            return pd.DataFrame(columns=['eChr', 'eStart', 'eEnd', 'eLength', 'eRepeatNum', 
                                      'MatDegree', 'eClass', 'eReads', 'query_id'])
        
        rows = []
        for _, row in df.iterrows():
            align_regions = row['AlignRegion'].split(';')
            for region in align_regions:
                mat_degree, location = region.split('|')
                chr, start_end = location.split('-', 1)
                start, end = map(int, start_end.split('-'))
                rows.append({
                    'eChr': chr,
                    'eStart': start,
                    'eEnd': end,
                    'eLength': row['consLen'],
                    'eRepeatNum': row['copyNum'],
                    'MatDegree': float(mat_degree),
                    'eClass': 'Mecc',
                    'eReads': row['readName'],
                    'query_id': row['query_id']
                })
        return pd.DataFrame(rows)

    def process_tecc_analysis(self, df):
        """处理tecc_analysis数据"""
        if df.empty:
            logging.warning("Tecc analysis DataFrame is empty. Returning empty DataFrame.")
            return pd.DataFrame(columns=['eChr', 'eStart', 'eEnd', 'eLength', 'eRepeatNum', 
                                      'MatDegree', 'eClass', 'eReads', 'query_id'])

        mecc_df = df[df['eClass'] == 'Mecc'].copy()
        if mecc_df.empty:
            logging.warning("No Mecc records found in tecc analysis. Returning empty DataFrame.")
            return mecc_df

        group_stats = mecc_df.groupby('eReads').agg({
            'eLength': 'mean',
            'eRepeatNum': 'mean'
        }).round().astype(int)

        group_stats['query_id'] = group_stats.apply(
            lambda x: f"{x.name}|Second|{x['eLength']}|{x['eRepeatNum']}|circular", 
            axis=1
        )

        mecc_df['query_id'] = mecc_df['eReads'].map(group_stats['query_id'])
        
        return mecc_df[['eChr', 'eStart', 'eEnd', 'eLength', 'eRepeatNum', 'MatDegree', 
                        'eClass', 'eReads', 'query_id']]

    def process_query_groups(self, df):
        """处理Query组和特征值"""
        df_processed = df.copy()
        df_processed['location_feature'] = df_processed.apply(
            lambda row: f"{row['eChr']}-{row['eStart']}-{row['eEnd']}", 
            axis=1
        )
        
        query_groups = df_processed.groupby('query_id').agg({
            'location_feature': list,
            'eRepeatNum': 'first',
            'eReads': 'first',
        }).reset_index()
        
        query_groups['feature_count'] = query_groups['location_feature'].apply(len)
        valid_groups = query_groups[query_groups['feature_count'] >= 1].copy()
        valid_groups['location_feature'] = valid_groups['location_feature'].apply(sorted)
        
        return valid_groups

    def merge_same_feature_groups(self, query_groups):
        """合并具有相同特征值的query组"""
        query_groups['feature_str'] = query_groups['location_feature'].apply(lambda x: ';'.join(sorted(x)))
        feature_groups = query_groups.groupby('feature_str').agg(list).reset_index()
        
        merged_groups = []
        for _, group in feature_groups.iterrows():
            merged_groups.append({
                'query_id': group['query_id'][0],
                'location_feature': query_groups[query_groups['query_id'] == group['query_id'][0]]['location_feature'].iloc[0],
                'eRepeatNum': sum(group['eRepeatNum']),
                'eReads': ';'.join(group['eReads'])
            })
        
        return pd.DataFrame(merged_groups)

    def integrate_data(self, merged_df, raw_combined_df):
        """整合merged_df和原始数据 - 优化版本"""
        # 1. 创建 query_id 到 merged info 的映射
        merged_info_map = merged_df.set_index('query_id')[['eRepeatNum', 'eReads']].to_dict('index')
        
        # 2. 过滤有效的 raw 数据
        valid_query_ids = set(merged_df['query_id'])
        filtered_raw = raw_combined_df[raw_combined_df['query_id'].isin(valid_query_ids)].copy()
        
        # 3. 使用矢量化操作更新数据
        filtered_raw['eRepeatNum'] = filtered_raw['query_id'].map(
            lambda x: merged_info_map[x]['eRepeatNum'])
        filtered_raw['eReads'] = filtered_raw['query_id'].map(
            lambda x: merged_info_map[x]['eReads'])
        
        # 4. 选择需要的列
        result = filtered_raw[['eChr', 'eStart', 'eEnd', 'eLength', 'MatDegree', 
                            'eClass', 'eRepeatNum', 'eReads', 'query_id']]
        
        return result

    def add_name_and_state(self, df):
        """添加eName和eState"""
        unique_query_ids = sorted(df['query_id'].unique())
        query_id_to_num = {qid: idx+1 for idx, qid in enumerate(unique_query_ids)}
        
        df['eName'] = df['query_id'].map(lambda x: f"Mecc|Confirmed|{query_id_to_num[x]}")
        df['eState'] = 'Confirmed-eccDNA'
        
        df['eName_num'] = df['eName'].str.extract(r'(\d+)').astype(int)
        df = df.sort_values('eName_num')
        df = df.drop('eName_num', axis=1)
        
        return df

    def merge_fasta_files(self, file1, file2, output_file):
        """合并FASTA文件"""
        if self.is_file_empty(file1) and self.is_file_empty(file2):
            logging.warning("Both FASTA files are empty or don't exist. Creating empty output file.")
            with open(output_file, 'w') as f:
                pass
            return
        
        if self.is_file_empty(file1):
            shutil.copy(file2, output_file)
            return
        
        if self.is_file_empty(file2):
            shutil.copy(file1, output_file)
            return

        try:
            records = list(SeqIO.parse(file1, "fasta")) + list(SeqIO.parse(file2, "fasta"))
            SeqIO.write(records, output_file, "fasta")
        except Exception as e:
            logging.error(f"Error merging FASTA files: {str(e)}")
            with open(output_file, 'w') as f:
                pass

    def update_fasta_names(self, fasta_file, name_map, output_file):
        """更新FASTA文件中的序列名称"""
        if self.is_file_empty(fasta_file):
            logging.warning(f"Input FASTA file {fasta_file} is empty. Creating empty output file.")
            with open(output_file, 'w') as f:
                pass
            return

        try:
            records = []
            for record in SeqIO.parse(fasta_file, "fasta"):
                base_id = record.id.split('|Second')[0]
                matched_query_id = None
                for query_id in name_map.keys():
                    if query_id.startswith(base_id):
                        matched_query_id = query_id
                        break
                
                if matched_query_id:
                    record.id = name_map[matched_query_id]
                    record.description = ""
                    records.append(record)
            
            SeqIO.write(records, output_file, "fasta")
            logging.info(f"Processed FASTA records: {len(records)}")
            
        except Exception as e:
            logging.error(f"Error updating FASTA names: {str(e)}")
            with open(output_file, 'w') as f:
                pass

    def prepare_final_dataframe(self, df):
        """整理最终数据框，添加MapNum列"""
        final_df = df.copy()
        
        # 计算eReadNum
        final_df['eReadNum'] = final_df['eReads'].apply(lambda x: len(str(x).split(';')))
        
        # 计算MapNum（每个query_id的组内行数）
        map_num_series = df.groupby('query_id').size()
        final_df['MapNum'] = final_df['query_id'].map(map_num_series)
        
        # 设置列顺序
        columns_order = [
            'eName','MapNum', 'eRepeatNum', 'eState', 'eChr', 'eStart', 'eEnd',
            'eLength', 'MatDegree', 'eClass', 'eReadNum', 'eReads'
        ]
        
        # 排序
        final_df['sort_key'] = final_df['eName'].str.extract(r'(\d+)').astype(int)
        final_df = final_df.sort_values('sort_key')
        final_df = final_df.drop('sort_key', axis=1)
        
        # 选择最终列
        return final_df[columns_order]

def main():
    parser = argparse.ArgumentParser(description="Process Mecc DNA data, generate output files, and remove duplicates.")
    parser.add_argument("--meccdna_part1", required=True, help="Input MeccDNA_part1.csv file")
    parser.add_argument("--tecc_analysis", required=True, help="Input tecc_analysis_results.csv file")
    parser.add_argument("--mecc_part1", required=True, help="Input mecc_part1.fa file")
    parser.add_argument("--mecc_part2", required=True, help="Input mecc_part2.fa file")
    parser.add_argument("--output_csv", required=True, help="Output CSV file")
    parser.add_argument("--output_fasta", required=True, help="Output FASTA file (deduplicated)")
    args = parser.parse_args()

    merger = MergeMecc(args.meccdna_part1, args.tecc_analysis, args.mecc_part1, 
                      args.mecc_part2, args.output_csv, args.output_fasta)
    merger.run_merge_mecc()

if __name__ == "__main__":
    main()
