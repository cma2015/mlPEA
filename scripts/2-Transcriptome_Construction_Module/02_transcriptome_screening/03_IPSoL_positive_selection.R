library(getopt)

spec <- matrix(
  c(
    "tmp_fasta", "a", 1, "character", "This is data directory!",
    "assem_length", "b", 1, "character", "This is results directory!",
    "combine_assembly_df", "c", 1, "character", "This is identity and coverage!",
    "transfrags_fasta", "d", 1, "character", "This is identity and coverage!",
    "output_set", "e", 1, "character", "This is identity and coverage!",
    "output_label", "f", 1, "character", "This is identity and coverage!"
  ),
  byrow = TRUE,
  ncol = 5
)
# 使用 getopt 解析命令行参数
opt <- getopt(spec=spec)

tmp_fasta <- opt$tmp_fasta
assem_length <- opt$assem_length
combine_assembly_df <- opt$combine_assembly_df
transfrags_fasta <- opt$transfrags_fasta
combine_assembly_df_path <- opt$combine_assembly_df
output_set <- opt$output_set
output_label <- opt$output_label

suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
library(Biostrings)

# mmseqs_dir聚类结果存放文件路径
system(command = sprintf('awk \'/^>/ { if (seq) { print length(seq), id } seq=""; id=substr($0, 2); next } { seq = seq $0 } END { if (seq) { print length(seq), id } }\'  %s > %s',
                         tmp_fasta,assem_length))
system(command = sprintf("sed -i 's/description:.*//g' %s",assem_length))
system(command = sprintf("sed -i 's/gene_symbol:.*//g' %s",assem_length))

plants_assembly_info <- fread(assem_length,
                              sep = " ",header = F, fill = TRUE)
plants_assembly_info <- plants_assembly_info[,-c(3,4,5,6,7)]
colnames(plants_assembly_info) <- c("transcript_length","transcript_id")

combine_assembly_df <- fread(combine_assembly_df,
                             sep = "\t",header = T,fill = T)

combine_assembly_df <- merge(combine_assembly_df,plants_assembly_info,by.x = "member_id",by.y = "transcript_id")
rm(plants_assembly_info)

grouped_df <- combine_assembly_df %>% group_by(cluster_id)
rm(combine_assembly_df)

# 组装得到的fa文件
fasta_file <- readDNAStringSet(transfrags_fasta)
assembly_id <- names(fasta_file)

assembly_grouped_df <- subset(grouped_df,grouped_df$cluster_id %in% assembly_id)
positive_set_1 <- unique(assembly_grouped_df$cluster_id)

non_assembly_grouped_df <- subset(grouped_df,!(grouped_df$cluster_id %in% assembly_grouped_df$cluster_id))
best_non_assembly_grouped_df <- non_assembly_grouped_df %>%
  group_by(cluster_id) %>%
  slice_max(transcript_length) %>%
  slice_head(n=1)
positive_set_2 <- unique(best_non_assembly_grouped_df$member_id)
positive_set <- c(positive_set_1,positive_set_2)
write.table(positive_set, sprintf("%s",output_set), 
            sep = "\t", col.names = F, row.names = FALSE, quote = FALSE)

seq_label <- data.frame(Seqname = assembly_id,
                        Label = rep(0,length(assembly_id))
                      )
seq_label$Label <- ifelse(seq_label$Seqname %in% positive_set,"1","0")
seq_label$Seqname <- paste0(">", seq_label$Seqname)
write.csv(seq_label,file = sprintf("%s",output_label),quote = F,row.names = F)