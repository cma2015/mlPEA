echo "genome: ${1}"

if [ ! -f "${1}.1.bt2" ]; then 
    /home/galaxy/miniconda3/envs/bio/bin/bowtie2-build ${1} ${1}
fi 