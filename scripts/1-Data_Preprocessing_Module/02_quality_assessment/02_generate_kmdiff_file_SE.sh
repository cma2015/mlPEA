#!/bin/bash

thread=$4
kmer_len=$5
hard_min=$6
echo "input: $1" >> $3/fof.txt
echo "ip: $2" >> $3/fof.txt

/home/galaxy/miniconda3/envs/bio/bin/kmdiff count --file $1/fof.txt \
    -d $3/tmp \
    --kmer-size $kmer_len --hard-min $hard_min --threads $thread &&

/home/galaxy/miniconda3/envs/bio/bin/kmdiff diff --km-run $3/tmp \
    -1 1 -2 1 --output-dir $1 -s 0.05 \
    --threads $thread