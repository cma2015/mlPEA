#!/bin/bash
index=$2
thread=$3
dir=$4
flag=0
if [ "$1" == "PE" ]; then
    IFS=',' read -r -a R1_list <<< "$5"
    IFS=',' read -r -a R2_list <<< "$6"
    length=${#R1_list[@]}
    for i in $(seq 0 $((length - 1)));do
        /home/galaxy/miniconda3/envs/bio/bin/bowtie2 -p $thread -x $index \
            -k 20 -1 ${R1_list[$i]} \
            -2 ${R2_list[$i]} 2> $dir/tmp_R${flag}_mapping_info |
        /home/galaxy/miniconda3/envs/bio/bin/samtools view -@ 10 -F 4 -Sb - |
        /home/galaxy/miniconda3/envs/bio/bin/samtools sort -@ 10 -o $dir/R${i}.tmp.bam  - 
        /home/galaxy/miniconda3/envs/bio/bin/samtools index $dir/R${i}.tmp.bam
        /home/galaxy/miniconda3/envs/bio/bin/samtools view \
            -bS -@ 10 -h -f 3 -F 12 \
            -o $dir/R${i}.clean.bam  $dir/R${i}.tmp.bam
        /home/galaxy/miniconda3/envs/bio/bin/samtools index $dir/R${i}.clean.bam
        rm -rf $dir/*tmp*.bam*
    done
else
    IFS=',' read -r -a SE_list <<< "$5"
    for SE in "${SE_list[@]}"
    do
        flag=$((flag + 1))
        /home/galaxy/miniconda3/envs/bio/bin/bowtie2 -p $thread -x $index \
            -k 20 -U $SE 2> $dir/tmp_R${flag}_mapping_info |
        /home/galaxy/miniconda3/envs/bio/bin/samtools view -@ $thread -F 4 -Sb - |
        /home/galaxy/miniconda3/envs/bio/bin/samtools sort -@ $thread -o $dir/R${flag}.clean.bam -
        /home/galaxy/miniconda3/envs/bio/bin/samtools index $dir/R${flag}.clean.bam
    done
fi 
