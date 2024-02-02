#!/bin/bash
Dir=/share/home/zhangkx/uploadFromHuicao/20peaks_overlap/Intersect_S1ENDSeq/withoutpromoter/intersect_mergedHDNA/jaccard/
# 存储 BED 文件列表
files=("CHD4.withoutPromoter.mergedHDNA.srt.bed" "DDX17.withoutPromoter.mergedHDNA.srt.bed" "DDX21.withoutPromoter.mergedHDNA.srt.bed" "DDX5.withoutPromoter.mergedHDNA.srt.bed" "DNMT1.withoutPromoter.mergedHDNA.srt.bed" "HNRNPC.withoutPromoter.mergedHDNA.srt.bed" "HNRNPK.withoutPromoter.mergedHDNA.srt.bed" "MATR3.withoutPromoter.mergedHDNA.srt.bed" "NUP93.withoutPromoter.mergedHDNA.srt.bed" "PCBP1.withoutPromoter.mergedHDNA.srt.bed" "PCBP2.withoutPromoter.mergedHDNA.srt.bed" "PDS5A.withoutPromoter.mergedHDNA.srt.bed" "PHB2.withoutPromoter.mergedHDNA.srt.bed" "PHB2.withoutPromoter.mergedHDNA.srt.bed" "PSIP1.withoutPromoter.mergedHDNA.srt.bed" "SFPQ.withoutPromoter.mergedHDNA.srt.bed" "SMARCA5.withoutPromoter.mergedHDNA.srt.bed" "TIF1B.withoutPromoter.mergedHDNA.srt.bed" "TOP2A.withoutPromoter.mergedHDNA.srt.bed" "TOP2B.withoutPromoter.mergedHDNA.srt.bed" "XRCC5.withoutPromoter.mergedHDNA.srt.bed")  # 替换为实际的文件名

# 创建一个空白的矩阵文件
matrix_file="jaccard_matrix.txt"
truncate -s 0 "$matrix_file"

# 循环遍历计算 Jaccard Index
for ((i = 0; i < ${#files[@]}; i++)); do
    for ((j = i + 1; j < ${#files[@]}; j++)); do
        file1="${files[$i]}"
        file2="${files[$j]}"

        # 提取文件名前缀作为输出文件的标识
        prefix=$(echo "$file1" | cut -d'.' -f1)_$(echo "$file2" | cut -d'.' -f1)

        # 使用 bedtools 计算 Jaccard Index
        bedtools jaccard -a "$file1" -b "$file2" > "${prefix}_jaccard.txt"
        # 提取 Jaccard Index 的值
        jaccard=$(awk 'NR==2 {print $3}' "${file1%%.*}_${file2%%.*}_jaccard.txt")

        # 写入矩阵文件
        echo -n "$jaccard " >> "$matrix_file"
    done
    echo >> "$matrix_file"
done
