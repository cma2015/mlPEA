echo "genome: ${1}" 

if [ ! -f "${1}.bwt" ]; then 
    /home/galaxy/miniconda3/envs/bio/bin/bwa index ${1}
fi 