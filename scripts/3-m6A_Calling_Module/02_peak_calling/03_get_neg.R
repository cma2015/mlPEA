library(argparse)

parser <- ArgumentParser()
parser$add_argument("-dir", help="the path of output dir", type="character", dest="dir")
parser$add_argument("-fasta", help="Fasta file", type="character", dest="fasta")
parser$add_argument("-peak", help="peak file", type="character", default=2, dest="peak")

bed_path <- "/home/galaxy/miniconda3/envs/bio/bin/bedtools"

args <- parser$parse_args()
output_dir <- args$dir
fasta_file <- args$fasta
peak_file <- args$peak

suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(future.apply))
suppressMessages(library(tidyr))

transfrags <- seqinr::read.fasta(fasta_file, seqtype = "DNA", as.string = T, forceDNAtolower = F)
seq_length <- sapply(transfrags, nchar) 
pos_peak <- read.table(peak_file, sep = "\t")

peakseq <- seqinr::read.fasta(sprintf("%s/pos.fa", output_dir), seqtype = "DNA", as.string = T, forceDNAtolower = F)

write.table(data.frame("name"=names(seq_length), "len"=seq_length), 
            file = sprintf("%s/transcript.size", output_dir), quote = F, row.names = F, col.names = F, sep = "\t")

system(command = sprintf("%s shuffle -i %s -g %s/transcript.size -chrom -excl > %s/neg.bed", bed_path, peak_file, output_dir, output_dir))
