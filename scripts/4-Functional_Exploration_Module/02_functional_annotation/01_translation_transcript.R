library(argparse)

parser <- ArgumentParser()
parser$add_argument("-fasta_file", help="Fasta file", type="character", dest="fasta_file")
parser$add_argument("-transcript_anno", help="transcript annotation", type="character", dest="trans_anno")
parser$add_argument("-dir", help="dir", type="character", dest="dir")

args <- parser$parse_args()

fasta_file_name <- args$fasta_file
dir <- args$dir

suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(seqinr))

pu_translationai <- read.table(sprintf("%s" , args$trans_anno),
                                header = FALSE, sep = "\t", fill = TRUE)
pu_translationai <- pu_translationai[pu_translationai$V2 != "", ]

transcript_info <- strsplit(pu_translationai$V1, "[:,(),]")

transfrags <- sapply(transcript_info, function(x) x[1])
pu_translationai_df <- data.frame(transfrags)
pre_info <- strsplit(as.character(pu_translationai$V2), ",")
pre_info <- do.call(rbind, pre_info)

pu_translationai_df$pre_start <- as.integer(pre_info[, 1])
pu_translationai_df$pre_end <- as.integer(pre_info[, 2])
pu_translationai_df$pre_start_score <- pre_info[, 3]
pu_translationai_df$pre_end_score <- pre_info[, 4]

pu_translationai_df$pre_start <- pu_translationai_df$pre_start + 1 
pu_translationai_df$pre_end <- pu_translationai_df$pre_end  + 2 + 1 
pu_translationai_df$transfrags <- gsub(">", "", pu_translationai_df$transfrags)

fasta_file <- sprintf(fasta_file_name)
f_transfrags <- seqinr::read.fasta(fasta_file, seqtype = "DNA", as.string = T,
                                    forceDNAtolower = F)
pu_translationai_df$length <- sapply(pu_translationai_df$transfrags, function(x) nchar(f_transfrags[[x]]))
pu_translationai_df$pre_end <- ifelse(pu_translationai_df$pre_end > pu_translationai_df$length, pu_translationai_df$length, pu_translationai_df$pre_end)

write.table(pu_translationai_df,sprintf("%s/tmp_pu_translationai_df.txt", dir),sep = "\t",col.names = T,row.names = F,quote = F)


write.table(pu_translationai_df[,c("transfrags","pre_start","pre_end")],sprintf("%s/tmp_PUlearning_ORF.txt", dir),
            sep = "\t",col.names = F,row.names = F,quote = F)
input <- sprintf("%s/tmp_PUlearning_ORF.txt", dir)
output <- sprintf("%s/tmp_PUlearning_pep.fa_tmp", dir)
translation_command <- sprintf("/home/galaxy/miniconda3/envs/bio/bin/python /galaxy/server/tools/07_annotation/01_translation_transcript.py %s %s %s",
                                fasta_file, input, output)
system(command = translation_command)
system(command = sprintf("/home/galaxy/miniconda3/envs/bio/bin/seqkit seq %s -w 0 > %s/tmp_PUlearning_pep.fa", output, dir))