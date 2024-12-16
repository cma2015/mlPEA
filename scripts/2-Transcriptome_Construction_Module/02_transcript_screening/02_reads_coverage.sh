#!/bin/bash

fasta_dir=$1
output=$2
thread=$3
data_type=$4
/home/galaxy/miniconda3/envs/bio/bin/samtools merge -@ $thread ${output}/merged.input.clean.bam ${output}/R*.clean.bam 
            
/home/galaxy/miniconda3/envs/bio/bin/samtools faidx ${fasta_dir}
cut -f 1-2 ${fasta_dir}.fai > ${fasta_dir}.genome
if [ "$data_type" == "PE" ]; then
    /home/galaxy/miniconda3/envs/bio/bin/bedtools genomecov -dz -pc -ibam ${output}/merged.input.clean.bam -g ${fasta_dir}.genome > ${output}/merged.input.txt
else
    /home/galaxy/miniconda3/envs/bio/bin/bedtools genomecov -dz -ibam ${output}/merged.input.clean.bam -g ${fasta_dir}.genome > ${output}/merged.input.txt
fi