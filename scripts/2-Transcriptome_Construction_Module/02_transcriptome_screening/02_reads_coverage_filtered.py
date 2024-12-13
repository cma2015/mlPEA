# python /home/malab29/Documents/disk4/mlPEA/scripts/read_coverage_filtered.py /home/malab29/Documents/disk8/Data   Ath_sample_8 "WT_0h_leaf WT_3h_leaf WT_0h_flower WT_3h_flower MUT_0h_leaf MUT_3h_leaf MUT_0h_flower MUT_3h_flower" sample_merge_tools_merge 
import os
import sys
from collections import defaultdict
import subprocess

def calculate_seq_lengths(fasta_file, length_file):
    seq_lengths = {}
    with open(fasta_file) as f, open(length_file, "w") as lf:
        seq_id = None
        seq = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq_id:
                    length = len("".join(seq))
                    seq_lengths[seq_id] = length
                    lf.write(f"{seq_id}\t{length}\n")
                seq_id = line[1:]
                seq = []
            else:
                seq.append(line)
        if seq_id:
            length = len("".join(seq))
            seq_lengths[seq_id] = length
            lf.write(f"{seq_id}\t{length}\n")
    return seq_lengths

def filter_coverage(coverage_file, threshold=3):
    seq_counts = defaultdict(int)
    with open(coverage_file) as f:
        for line in f:
            seq, pos, count = line.strip().split()
            count = int(float(count))  # Convert count to float first, then to int
            if count >= threshold:
                seq_counts[seq] += 1
    return seq_counts

def filter_sequences(seq_lengths, seq_counts, output_file):
    with open(output_file, "w") as of:
        for seq, length in seq_lengths.items():
            count = seq_counts.get(seq, 0)
            if count >= length / 2:
                of.write(f"{seq}\n")

def run_seqkit_commands(fasta_file, output_file, output_dir, results):
    filtered_fasta = f"{output_dir}/{results}_filtered.fa"
    seqkit_grep = f"/home/galaxy/miniconda3/envs/bio/bin/seqkit grep -f {output_file} {fasta_file} > {filtered_fasta}"
    seqkit_fx2tab = f"/home/galaxy/miniconda3/envs/bio/bin/seqkit fx2tab --length --name {fasta_file} > {output_dir}/seq_lengths.txt"
    seqkit_grep_v = f"/home/galaxy/miniconda3/envs/bio/bin/seqkit grep -v -f {output_dir}/long_seq_ids.txt {fasta_file} > {output_dir}/temp.fasta"

    # Execute commands and check results
    subprocess.run(seqkit_grep, shell=True, check=True)
    subprocess.run(seqkit_fx2tab, shell=True, check=True)

    # # Debugging: check seq_lengths.txt content
    # with open(f"{output_dir}/seq_lengths.txt") as lf:
    #     for line in lf:
    #         print(line.strip())  # Print lengths for debugging

    with open(f"{output_dir}/seq_lengths.txt") as lf, open(f"{output_dir}/long_seq_ids.txt", "w") as lf_out:
        for line in lf:
            seq, length = line.strip().split()
            length = int(length)
            if length > 30000:
                lf_out.write(f"{seq}\n")

    subprocess.run(seqkit_grep_v, shell=True, check=True)
    os.replace(f"{output_dir}/temp.fasta", fasta_file)

    # Debugging: Check final fasta file
    grep_cmd = f"grep '>' {fasta_file} | wc -l"
    subprocess.run(grep_cmd, shell=True, check=True)

    grep_cmd = f"/home/galaxy/miniconda3/envs/bio/bin/seqkit fx2tab --length --name {fasta_file} | awk '$2 > 30000'"
    subprocess.run(grep_cmd, shell=True, check=True)

    os.remove(f"{output_dir}/seq_lengths.txt")
    os.remove(f"{output_dir}/long_seq_ids.txt")


def process_sample(fasta_file, coverage_file, output_dir, results):
    os.makedirs(output_dir, exist_ok=True)
    output_file = f"{output_dir}/{results}_filtered_ids.txt"
    length_file = f"{output_dir}/{results}.lengths.txt"

    print(f"Processing sample: {results}")
    print(f"Fasta file: {fasta_file}")
    print(f"Coverage file: {coverage_file}")
    print(f"Output dir: {output_dir}")
    print(f"Output file: {output_file}")
    print(f"Length file: {length_file}")

    seq_lengths = calculate_seq_lengths(fasta_file, length_file)
    seq_counts = filter_coverage(coverage_file)
    filter_sequences(seq_lengths, seq_counts, output_file)
    run_seqkit_commands(fasta_file, output_file, output_dir, results)

def main(fasta, input_file, output_dir, threshold=3):

    fasta_file = fasta
    coverage_file = input_file
    output_dir = output_dir
    process_sample(fasta_file, coverage_file, output_dir, "fasta")

if __name__ == "__main__":
    fasta = sys.argv[1]
    input_file = sys.argv[2]
    output_dir = sys.argv[3]
    threshold = sys.argv[4]
    main(fasta, input_file, output_dir, threshold)