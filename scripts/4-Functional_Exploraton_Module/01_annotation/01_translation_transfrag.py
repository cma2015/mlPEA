import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# python script.py input.fasta input.txt output.fasta

# 输入标准fa文件，1-based 编码区域文件
# 读取fasta文件
def read_fasta(file):
    fasta_sequences = SeqIO.to_dict(SeqIO.parse(file, "fasta"))
    return fasta_sequences

# 读取txt文件
def read_txt(file):
    with open(file, 'r') as f:
        lines = f.readlines()
    transcript_info = {line.split()[0]: (int(line.split()[1]), int(line.split()[2])) for line in lines}
    return transcript_info

# 翻译序列
def translate_sequences(fasta_sequences, transcript_info):
    protein_sequences = {}
    for transcript, (start, end) in transcript_info.items():
        if transcript in fasta_sequences:
            coding_rna = fasta_sequences[transcript].seq[start-1:end]  # 注意这里是1-based，所以要减1
            # 检查序列长度是否是3的倍数
            if len(coding_rna) % 3 != 0:
                print(f"序列 {transcript} 的长度不是3的倍数，正在补全N。")
                while len(coding_rna) % 3 != 0:
                    coding_rna += 'N'
            # 直接翻译RNA序列
            protein_seq = coding_rna.translate(to_stop=True)
            protein_sequences[transcript] = protein_seq
    return protein_sequences

# 保存为fasta文件
def save_as_fasta(protein_sequences, output_file):
    records = [SeqRecord(seq, id=id, description="") for id, seq in protein_sequences.items()]
    SeqIO.write(records, output_file, "fasta")

# 处理命令行参数
parser = argparse.ArgumentParser(description='Translate coding regions to protein sequences.')
parser.add_argument('fasta_file', type=str, help='The fasta file containing transcript sequences.')
parser.add_argument('txt_file', type=str, help='The txt file containing start and end positions of coding regions.')
parser.add_argument('output_file', type=str, help='The output fasta file to save protein sequences.')
args = parser.parse_args()

fasta_sequences = read_fasta(args.fasta_file)
transcript_info = read_txt(args.txt_file)
protein_sequences = translate_sequences(fasta_sequences, transcript_info)
save_as_fasta(protein_sequences, args.output_file)