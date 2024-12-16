library(getopt)

spec <- matrix(
  c(
    "clu_rep_fasta", "a", 1, "character", "This is data directory!",
    "clu_rep_txt", "b", 1, "character", "This is results directory!",
    "clustered_tsv", "c", 1, "character", "This is identity and coverage!",
    "transfrags_fasta", "d", 1, "character", "This is identity and coverage!",
    "combine_assembly_df", "e", 1, "character", "This is identity and coverage!",
    "assembly_id_out", "f", 1, "character", "This is identity and coverage!"
  ),
  byrow = TRUE,
  ncol = 5
)
opt <- getopt(spec=spec)

clu_rep_fasta <- opt$clu_rep_fasta
clu_rep_txt <- opt$clu_rep_txt
clustered_tsv <- opt$clustered_tsv
transfrags_fasta <- opt$transfrags_fasta
combine_assembly_df_path <- opt$combine_assembly_df
assembly_id_out <- opt$assembly_id_out

suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
library(Biostrings)

cluster <- fread(clustered_tsv, sep = "\t", col.names =c("cluster_id", "member_id"))

cluster_result <- cluster %>%
  group_by(cluster_id) %>%
  summarize(member_count = n_distinct(member_id))

pre_cluster <- cluster %>%
  left_join(cluster_result, by = "cluster_id")
rm(cluster)
rm(cluster_result)
system(command = sprintf('awk \'/^>/ { if (seq) { print length(seq), id } seq=""; id=substr($0, 2); next } { seq = seq $0 } END { if (seq) { print length(seq), id } }\'  %s > %s',clu_rep_fasta, clu_rep_txt))
system(command = sprintf("sed -i 's/gene_biotype:.*//g' %s ",clu_rep_txt))

cluster_info_file = sprintf("%s", clu_rep_txt)

cluster_info <- fread(cluster_info_file, sep = " ", header = F, fill = TRUE,
                      col.names = c("cluster_length","cluster_id","cluster_type","else","cluster_geneid"))
cluster_info_df <- separate(data = cluster_info, col = "else", into = c("chr_type","cluster_genome","cluster_chr",      "cluster_start",        "cluster_end","cluster_strand"), sep = ":")
rm(cluster_info)

cluster_info_df <- merge(cluster_info_df,pre_cluster,by = "cluster_id")

fasta_file <- readDNAStringSet(transfrags_fasta)
assembly_id <- names(fasta_file)
write.table(assembly_id, sprintf("%s",assembly_id_out), 
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
assembly_info_df <- subset(cluster_info_df,cluster_info_df$member_id %in% assembly_id)
rm(cluster_info_df)

assembly_info_df_list <- split(assembly_info_df,assembly_info_df$cluster_id)
only_assembly_list <- list()
combine_assembly_list <- list()

for (i in 1:length(assembly_info_df_list)) {
  cluster_ids <- unique(assembly_info_df_list[[i]][["cluster_id"]])
  member_ids <- assembly_info_df_list[[i]][["member_id"]]
  member_count <- unique(assembly_info_df_list[[i]][["member_count"]])
  if(member_count==1){
    only_assembly_list <- append(only_assembly_list, list(assembly_info_df_list[[i]]))
  } else{
    if (all(cluster_ids %in% assembly_id) && length(member_ids) == member_count ) {
      only_assembly_list <- append(only_assembly_list, list(assembly_info_df_list[[i]]))
    } else {
      combine_assembly_list <- append(combine_assembly_list, list(assembly_info_df_list[[i]]))
    }
  }
}

combine_assembly_df <- do.call(rbind, lapply(combine_assembly_list, data.frame))
write.table(combine_assembly_df, sprintf("%s",combine_assembly_df_path), 
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)