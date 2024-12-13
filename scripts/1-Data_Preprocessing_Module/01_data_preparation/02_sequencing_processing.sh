#!/bin/bash
echo "Script:$0";
echo "SRA ID:$1";
echo "minReadLen:$2";



dirpath=/galaxy/server/tools/01_data_preparation/$3/

/home/galaxy/miniconda3/envs/bio/bin/fastq-dump  $1 -M $2 --split-3 -O $dirpath ;


count=$(ls ${dirpath} | wc -l)
echo $count
#ls -l ${dirpath} | wc -l > $6



name=$3
if [ ${count} == "2" ] ; then
    mv ${dirpath}*_1.fastq $4 ; #cat ${dirpath}*_1.fastq > $4 ;
    mv ${dirpath}*_2.fastq $5 ;
else
    mv ${dirpath}*.fastq  $4 ;
fi