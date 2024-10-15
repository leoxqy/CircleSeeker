# CircleSeeker: 完整工作流程教程

本教程描述了使用CircleSeeker工具套件进行环状DNA (eccDNA) 识别和分析的完整流程。

## 0. 准备工作

确保您已经安装了所有必要的工具和依赖项，包括samtools、TideHunter、BLAST、minimap2等。

### 0.1 构建BLAST数据库

在开始分析之前，我们需要使用参考基因组构建BLAST数据库。这一步骤只需执行一次，除非您更换了参考基因组。

```bash
# 假设您的参考基因组文件名为 hg38.fa
makeblastdb -in hg38.fa -dbtype nucl -out hg38_db
```

这个命令会创建一个名为 `hg38_db` 的BLAST数据库，我们将在后续的BLAST比对步骤中使用它。

## 1. CircleSeeker_Hunter: 使用TideHunter进行初始识别

首先，我们使用TideHunter对HIFI读段进行初始识别。

```bash
TideHunter -f 2 -t 18 -k 16 -w 1 -p 100 -P 2000000 -e 0.1 All_HIFI_reads.fasta > eccDNA_candidates_TH.txt
```

## 2. CircleSeeker_Carousel: 处理TideHunter的结果

处理TideHunter的输出，生成初步的候选序列。

```bash
python eccDNA_Carousel.py -i eccDNA_candidates_TH.txt -o processed_candidates.csv -r CtcR.List_patr1.csv -c DL_circular_sequences.fasta
```

## 3. CircleSeeker_Aligner: 将候选序列比对到参考基因组

使用BLAST将候选序列比对到人类参考基因组（hg38）。

```bash
blastn -db hg38_db -query DL_circular_sequences.fasta -out aligned_candidates.txt -num_threads 18 -word_size 100 -evalue 1e-50 -perc_identity 99 -outfmt 6
```

## 4. CircleSeeker_Ringmaster: 分类和处理比对结果

根据比对结果对候选序列进行分类。

```bash
python eccDNA_Ringmaster_1.py -i aligned_candidates.txt -fa DL_circular_sequences.fasta --uecc_fa uecc_part1.fa --mecc_fa mecc_part1.fa --cecc_fa cecc.fa --xecc_fa xecc.fa -u UeccDNA_part1.csv -m MeccDNA_part1.csv
```

## 5. CircleSeeker_Sieve_1: 筛选未被选中的读段

筛选出未在初始识别中被选中的读段。

```bash
python eccDNA_Sieve_1.py -i All_HIFI_reads.fasta -l eccDNA_candidates_TH.txt -u unSelected_reads.fa -d t
```

## 6. CircleSeeker_RealignUnselected: 重新比对未被选中的读段

将未被选中的读段重新比对到参考基因组。

```bash
blastn -db hg38_db -query unSelected_reads.fa -out unSelected.aligned_results.txt -num_threads 18 -word_size 100 -evalue 1e-50 -perc_identity 99 -outfmt 6
```

## 7. CircleSeeker_IndexUnselected: 为未被选中的读段创建索引

为未被选中的读段创建索引，以便后续处理。

```bash
samtools faidx unSelected_reads.fa
```

## 8. CircleSeeker_Juggler: 分析未被选中的读段的比对结果

分析未被选中读段的比对结果，生成第二部分的候选序列。

```bash
python eccDNA_Juggler_1.py -i unSelected.aligned_results.txt -f unSelected_reads.fa.fai -n Num_LinR.csv -o tecc_analysis_results.csv -r CtcR.List_patr2.csv
```

## 9. CircleSeeker_Tamer: 筛选出part2的UeccDNA和MeccDNA

从第二部分的分析结果中筛选出UeccDNA和MeccDNA。

```bash
python eccDNA_Tamer.py -i tecc_analysis_results.csv -f unSelected_reads.fa -u uecc_part2.fa -m mecc_part2.fa
```

## 10. CircleSeeker_Sieve_2: 进一步筛选未分类的读段

进一步筛选出仍未分类的读段。

```bash
python eccDNA_Sieve_1.py -i unSelected_reads.fa -l CtcR.List_patr2.csv -u Final_unclassified.fa -d c
```

## 11. CircleSeeker_FinalAlign: 最终未分类读段的比对

使用minimap2对最终未分类的读段进行比对。

```bash
minimap2 -t 18 -x map-hifi -a hg38.fa Final_unclassified.fa > Final_unclassified.sam
```

## 12. CircleSeeker_SortAndIndex: 对最终比对结果进行排序和索引

对最终比对结果进行排序和索引，以便后续分析。

```bash
samtools sort -o Final_unclassified_sorted.bam Final_unclassified.sam
samtools index Final_unclassified_sorted.bam
```

## 13. CircleSeeker_Astrologer: 最终分析和评分

对最终的比对结果进行分析和评分。

```bash
python eccDNA_Astrologer.py -i Final_unclassified_sorted.bam -o circular_dna_results_scored.csv
```

## 14. 合并UeccDNA结果

合并两部分的UeccDNA结果。

```bash
python eccDNA_Merge.Uecc.py --ueccdna_part1 UeccDNA_part1.csv --tecc_analysis tecc_analysis_results.csv --uecc_part1 uecc_part1.fa --uecc_part2 uecc_part2.fa --output_csv Uecc.Confirmed.csv --output_fasta UeccDNA.Confirmed.fa
```

## 15. 合并MeccDNA结果

合并两部分的MeccDNA结果。

```bash
python eccDNA_Merge.Mecc.py --meccdna_part1 MeccDNA_part1.csv --tecc_analysis tecc_analysis_results.csv --mecc_part1 mecc_part1.fa --mecc_part2 mecc_part2.fa --output_csv Mecc.Confirmed.csv --output_fasta MeccDNA.Confirmed.fa
```

## 16. 处理其他类型的eccDNA

处理其他类型的eccDNA，包括XeccDNA和CeccDNA。

```bash
samtools faidx xecc.fa
python eccDNA_Treat.Other.py --xecc_fai xecc.fa.fai --xecc_fa xecc.fa --cecc_csv CeccDNA.csv --cecc_fa cecc.fa --uecc_csv circular_dna_results_scored.csv
```

## 17. 生成最终分析报告

使用 `CircleSeeker_ReportGenerator.py` 脚本整合所有分析结果，生成最终的综合报告。

```bash
samtools faidx All_HIFI_reads.fasta
python eccDNA_ReportGenerator.py \
  --fai All_HIFI_reads.fasta.fai \
  --ctcr1 CtcR.List_patr1.csv \
  --ctcr2 CtcR.List_patr2.csv \
  --linr Num_LinR.csv \
  --uecc Uecc.Confirmed.csv \
  --mecc Mecc.Confirmed.csv \
  --xecc Xecc.Confirmed.csv \
  --cecc Cecc.Confirmed.csv \
  --uecc_inferred Uecc.Inferred.csv \
  --output eccDNA_analysis_report.txt
```

这个脚本将整合以下信息：
- 原始HIFI读段的信息（来自 `All_HIFI_reads.fasta.fai`）
- 两部分的环状读段列表（`CtcR.List_patr1.csv` 和 `CtcR.List_patr2.csv`）
- 线性读段数量信息（`Num_LinR.csv`）
- 各类型eccDNA的分析结果（UeccDNA、MeccDNA、XeccDNA、CeccDNA）
- 推断的UeccDNA信息

最终报告将保存在 `eccDNA_analysis_report.txt` 文件中，提供整个分析过程的综合概览。

## 18. 清理中间文件

为了节省存储空间并保持工作目录整洁，我们可以删除不再需要的中间文件。请注意，在删除之前，确保你已经成功生成了最终报告，并备份了所有重要的结果文件。

```bash
# 删除中间文件
rm eccDNA_candidates_TH.txt
rm processed_candidates.csv
rm aligned_candidates.txt
rm unSelected_reads.fa
rm unSelected.aligned_results.txt
rm tecc_analysis_results.csv
rm Final_unclassified.fa
rm Final_unclassified.sam
rm Final_unclassified_sorted.bam
rm Final_unclassified_sorted.bam.bai

# 删除临时FASTA文件
rm uecc_part1.fa mecc_part1.fa cecc.fa xecc.fa
rm uecc_part2.fa mecc_part2.fa

# 删除临时CSV文件
rm UeccDNA_part1.csv MeccDNA_part1.csv
rm CtcR.List_patr1.csv CtcR.List_patr2.csv
rm Num_LinR.csv

# 保留最终结果文件
# UeccDNA.Confirmed.fa, MeccDNA.Confirmed.fa, CeccDNA.Confirmed.fa, XeccDNA.Confirmed.fa
# Merge.Uecc.csv, Merge.Mecc.csv, Xecc.Confirmed.csv, Cecc.Confirmed.csv
# Uecc.Inferred.csv
# eccDNA_analysis_report.txt
```

注意：在运行此清理脚本之前，请确保：
1. 您已经成功完成了所有分析步骤。
2. 最终报告（eccDNA_analysis_report.txt）已经生成。
3. 您已经备份了所有重要的结果文件。



完成这一步后，您的工作目录将只包含最终的分析结果和报告，使得结果更易于管理和解释。

完成以上所有步骤后，您将获得全面的eccDNA分析结果，包括不同类型的eccDNA、它们的特征，以及一份详细的分析报告。这份报告可以帮助您更好地理解和解释分析结果，为后续的研究提供重要参考。

