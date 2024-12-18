library(argparse)

parser <- ArgumentParser()
parser$add_argument("-inputfile", help="Input file", type="character", dest="inputfile")
parser$add_argument("-peak", help="Peak file", type="character", dest="peak")
parser$add_argument("-threshold", help="threshold", type="numeric", default=0.5, dest="threshold")
parser$add_argument("-output_dir", help="Output dir", type="character", dest="output_dir")

args <- parser$parse_args()
inputfile <- args$inputfile
peak <- args$peak
threshold <- args$threshold
output_dir <- args$output_dir


prediction_peak <- read.table(peak, sep = "\t")
prediction_score <- read.table(inputfile)
prediction_peak_score <- cbind(prediction_peak, prediction_score)
colnames(prediction_peak_score) <- c("transfrags", "start", "end", "peak_id", "score")

weakRM_remain_peak <- prediction_peak_score[prediction_peak_score$score > threshold, ]
write.table(weakRM_remain_peak[,1:4], file = sprintf("%s/remain_peak.bed", output_dir), sep = "\t", col.names = F, row.names = F, quote = F)

weakRM_remove_peak <- prediction_peak_score[prediction_peak_score$score <= threshold, ]
write.table(weakRM_remove_peak[,1:4], file = sprintf("%s/remove_peak.bed", output_dir), sep = "\t", col.names = F, row.names = F, quote = F)