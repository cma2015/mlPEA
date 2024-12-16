library(argparse)
library(utils)
library(QNB)
source("/galaxy/server/tools/07_annotation/00_Reference_free.CountRead.R")
source("/galaxy/server/tools/07_annotation/00_Reference_free.callPeakFisher.R")
source("/galaxy/server/tools/07_annotation/00_Reference_free.MergeRep.Peakcalling.R")
source("/galaxy/server/tools/07_annotation/00_Reference_free.MergeRep.Calculate_logOR.R")

parser <- ArgumentParser()
parser$add_argument("-fasta_file", help="Fasta file", type="character", dest="fasta_file")
parser$add_argument("-case_ip_bam", help="transcript annotation", type="character", dest="case_ip_bam")
parser$add_argument("-control_ip_bam", help="transcript annotation", type="character", dest="control_ip_bam")
parser$add_argument("-case_input_bam", help="transcript annotation", type="character", dest="case_input_bam")
parser$add_argument("-control_input_bam", help="transcript annotation", type="character", dest="control_input_bam")
parser$add_argument("-dir", help="dir", type="character", dest="dir")
parser$add_argument("-threads", help="thread", type="numeric", dest="threads")

args <- parser$parse_args()

fasta_file_name <- args$fasta_file
dir <- args$dir
case_input_bam <- args$case_input_bam
control_input_bam <- args$control_input_bam
case_ip_bam <- args$case_ip_bam
control_ip_bam <- args$control_ip_bam
threads <- args$threads

case_ip_path_c <- strsplit(case_ip_bam, ",")[[1]]
control_ip_path_c <- strsplit(control_ip_bam, ",")[[1]]
case_input_path_c <- strsplit(case_input_bam, ",")[[1]]
control_input_path_c <- strsplit(control_input_bam, ",")[[1]]
case_nums <- length(case_ip_path_c)
control_nums <- length(control_ip_path_c)

case_path <- dirname(case_ip_path_c)[1]

case_ip_file_name <- basename(case_ip_path_c)
control_ip_file_name <- basename(control_ip_path_c)
case_input_file_name <- basename(case_input_path_c)
control_input_file_name <- basename(control_input_path_c)

total_ip_file_name <- c(case_ip_file_name, control_ip_file_name)
total_input_file_name <- c(case_input_file_name, control_input_file_name)

case_ip_colnames <- paste0("case_rep", 1:case_nums,"_ip")
control_ip_colnames <- paste0("control_rep", 1:control_nums,"_ip")
case_input_colnames <- paste0("case_rep", 1:case_nums,"_input")
control_input_colnames <- paste0("control_rep", 1:control_nums,"_input")

samples <- c("case", "control")
samplenames_case <- paste0("case_rep", 1:case_nums)
samplenames_control <- paste0("control_rep", 1:control_nums)
samplenames <- c(samplenames_case, samplenames_control)

total_ip_colnames <- c(case_ip_colnames, control_ip_colnames)
total_input_colnames <- c(case_input_colnames, control_input_colnames)

Reference_free_readCounts <- countReads(IP_filenames = total_ip_file_name,# file name of ip
           Input_filenames=total_input_file_name,# file name of input
           IP_colnames=total_ip_colnames, # name of ip output column
           Input_colnames=total_input_colnames, # name of input output column
           fasta_fliename=fasta_file_name, # fasta file used for peak calling
           samplenames=samplenames, # name of sample
           fragmentLength = 150,
           bamFolder=case_path,
           outputDir=dir,
           binSize = 50,
           strandToKeep = "opposite",
           paired = T,
           threads = threads,
           saveOutput = T)

MeRIPdata <- callPeakFisher(Reference_free_readCounts,threads = 10)

MergeRep_PeakCallResult <- MergePeakCallResult(MeRIP_peakCallResult = MeRIPdata@peakCallResult,
                                               samples = samples,rep_number = case_nums,rep_filter_threshold = 1)


MergeRep_JointPeak <- reportJointPeak_MergeRep(joint_threshold = 1,MeRIP_peakCallResult = MeRIPdata,
                                               MergeRep_PeakCallResult = MergeRep_PeakCallResult,threads=15)
MeRIP_peakCallResult = MeRIPdata
MergeRep_jointPeakCount <- jointPeakCount(MergeRep_JointPeak=MergeRep_JointPeak)

MergeRep_ORlist <- list()
MergeRep_JointPeaklist <- list()
transcriptBins_list <- list()

collect_list <- MergeRep_Calculate_OR(MergeRep_jointPeakCount=MergeRep_jointPeakCount,
                                      samples = samples,rep_number = case_nums,rep_filter_threshold = 1)

MergeRep_ORlist <- collect_list[[1]]
MergeRep_JointPeaklist <- collect_list[[2]]
transcriptBins_list <- collect_list[[3]]
transcriptBins_pairs_list <- list(collect_list[[4]])

if (is.double(MergeRep_ORlist)) {
  MergeRep_OR <- data.frame(MergeRep_ORlist)
} else if (is.list(MergeRep_ORlist)) {
  MergeRep_OR <- do.call(rbind, MergeRep_ORlist)
} else {
  stop("MergeRep_ORlist is not in a recognized format.")
}


transcript_info <- rownames(MergeRep_OR)
transcript <- sapply(transcript_info, function(x) strsplit(x, ",")[[1]][1])
start <- sapply(transcript_info, function(x) strsplit(x, ",")[[1]][2])


MergeRep_OR <- data.frame(transcript = transcript, start = start, MergeRep_OR, row.names = NULL)


MergeRep_OR <- MergeRep_OR[, c("transcript", "start", names(MergeRep_OR)[-c(1, 2)])]

write.table(MergeRep_OR, file = sprintf("%s/MergeRep_OR.txt", dir),
            sep = "\t", quote = F, row.names = F)

MergeRep_JointPeak <-as.data.frame(t(do.call(rbind,MergeRep_JointPeaklist)))
MergeRep_transcriptBins <- as.data.frame(t(do.call(rbind, transcriptBins_list)))
MergeRep_jointPeak_id_pairs<- as.data.frame(do.call(rbind,transcriptBins_pairs_list))

MergeRep_JointPeak_df <- merge(MergeRep_JointPeak,MergeRep_OR,by = c("transcript","start"))

MergeRep_JointPeak_df$peak_id <- paste("jointpeak", seq_len(nrow(MergeRep_JointPeak_df)), sep = "_")

write.table(MergeRep_JointPeak_df, file = sprintf("%s/joint_peak_quantification.txt", dir),
            sep = "\t", quote = F, row.names = F)

jointPeak <- as.data.frame(MergeRep_jointPeakCount@jointPeaks)
quantification_input <- as.data.frame(MergeRep_jointPeakCount@jointPeak_input)
quantification_ip <-  as.data.frame(MergeRep_jointPeakCount@jointPeak_ip)
merged_quantification <- merge(quantification_input, quantification_ip, by = "row.names")
quant_matrix <- merged_quantification



write.table(quant_matrix, file = sprintf("%s/quant_matrix.txt", dir),
            sep = "\t", quote = F, row.names = F)
# Running QNB quantification
group_id_1 <- samples[1]
group_id_2 <- samples[2]
# Run the QNB by using the average counts of peaks
meth1 <- quant_matrix[,control_ip_colnames]
meth2 <-  quant_matrix[,case_ip_colnames]
unmeth1 <- quant_matrix[,control_input_colnames]
unmeth2 <- quant_matrix[,case_input_colnames]

# Run QNB test
result <- qnbtest(meth1, meth2, unmeth1, unmeth2, mode = "auto", output.dir = dir)
colnames(result) <- c("p.treated", "p.control", "log2FC", "log2.OR", "pvalue", "qvalue", "padj")
  
result$peak_id <- quant_matrix$Row.names

write.table(result, file = sprintf("%s/diff_peaks.txt",dir),
              sep = "\t", quote = FALSE)

transcript_info <- result$peak_id
transcript <- sapply(transcript_info, function(x) strsplit(x, ",")[[1]][1])
start <- sapply(transcript_info, function(x) strsplit(x, ",")[[1]][2])
result$transcript <- transcript
result$start <- start

merge_results <- merge(result, MergeRep_JointPeak_df, by = c("transcript","start"), all.x = TRUE)

write.table(merge_results, file = sprintf("%s/all_diffm6A.txt",dir),
            sep = "\t", quote = FALSE,row.names = FALSE)

diff_meth <- merge_results[merge_results$padj < 0.05 & abs(merge_results$log2FC) >=1,]

diff_meth <- diff_meth[,c("transcript","start","end","name","score","strand","thickStart","thickEnd",
                          "itemRgb","blockCount", "blockSizes",  "blockStarts", "log2FC","padj")]

write.table(diff_meth, file = sprintf("%s/diffm6A.txt",dir),
            sep = "\t", quote = FALSE,row.names = FALSE)

awk_command <- sprintf(
  'awk \'BEGIN {OFS="\\t"} NR > 1 {print $1, $2-1, $3, "diff_peak" (NR-1)}\' %s/diffm6A.txt > %s/diffm6A_peak.bed',
  dir,
  dir
)
system(awk_command)