#!/bin/bash

thread=$6
kmer_len=$7
hard_min=$8
echo "input: $2;$3" >> $1/fof.txt
echo "ip: $4;$5" >> $1/fof.txt

/home/galaxy/miniconda3/envs/bio/bin/kmdiff count --file $1/fof.txt \
    -d $1/tmp \
    --kmer-size $kmer_len --hard-min $hard_min --threads $thread &&

/home/galaxy/miniconda3/envs/bio/bin/kmdiff diff --km-run $1/tmp \
    -1 1 -2 1 --output-dir $1 -s 0.05 \
    --threads $thread