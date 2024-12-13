#!/bin/bash
# bash convert_standard_fa.sh PUlearning_pep.fa  PUlearning_pep_mod.fa 
# 用于处理非标准格式的 FASTA 文件。它会将输入文件中的每个序列提取出来，并将序列 ID 和序列内容写入为标准的fasta文件（第一行为序列id，第二行为序列）中。

# 检查参数个数
if [ $# -ne 2 ]; then
    echo "Usage: $0 <input_file> <output_file>"
    exit 1
fi

# 输入文件
input_file=$1
# 输出文件
output_file=$2

# 清空输出文件
echo -n "" > $output_file

# 初始化变量
sequence_id=""
sequence=""

# 读取非标准的fasta文件
while IFS= read -r line
do
    # 检查是否为序列id
    if [[ $line == \>* ]]; then
        # 如果sequence变量不为空，则输出上一个序列
        if [[ -n $sequence ]]; then
            echo -e "$sequence_id\n$sequence" >> $output_file
        fi

        # 更新序列id
        sequence_id=$line
        # 清空sequence变量
        sequence=""
    else
        # 将行添加到sequence变量
        sequence+=$line
    fi
done < "$input_file"

# 输出最后一个序列
if [[ -n $sequence ]]; then
    echo -e "$sequence_id\n$sequence" >> $output_file
fi