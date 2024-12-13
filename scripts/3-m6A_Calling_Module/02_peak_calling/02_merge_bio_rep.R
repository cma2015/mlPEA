library(argparse)

parser <- ArgumentParser()
parser$add_argument("-input", help="Input file", type="character", dest="input")
parser$add_argument("-output", help="Output file", type="character", dest="output")
parser$add_argument("-num", help="Number of replicates", type="integer", default=2, dest="num")

args <- parser$parse_args()

input <- args$input
output <- args$output
num <- args$num

merge_peaks_txt <- function(peak_files, cover_rep_num=2) {
  library(GenomicRanges)
  library(IRanges)
  
  peaks_list <- lapply(peak_files, function(file) {
    peaks_data <- read.table(file, header=TRUE)
    
    peaks_data$Start <- as.numeric(peaks_data$start)
    peaks_data$End <- as.numeric(peaks_data$end)
    
    gr <- GRanges(seqnames = peaks_data$seqnames, 
                  ranges = IRanges(start = peaks_data$start, end = peaks_data$end),
                  mean_FDR = peaks_data$mean_FDR,
                  max_FDR = peaks_data$max_FDR,
                  min_FDR = peaks_data$min_FDR,
                  mean_ratio = peaks_data$mean_ratio,
                  max_ratio = peaks_data$max_ratio,
                  min_ratio = peaks_data$min_ratio)
    return(gr)
  })
  peaks_granges_list <- GRangesList(peaks_list)
  peak_coverage <- coverage(peaks_granges_list)
  
  covered_ranges <- IRanges::slice(peak_coverage, lower=cover_rep_num, rangesOnly=TRUE)
  merged_peaks <- GRanges(covered_ranges)
  return(merged_peaks)
}

peak_files <- unlist(strsplit(input, ","))
merged_peaks <- merge_peaks_txt(peak_filesk, num)

free_peak <- as.data.frame(merged_peaks)[,1:4]
free_peak <- free_peak[free_peak$width >= 100,]
free_peak$start <- free_peak$start - 1
free_peak$peak_id <- sprintf("peak_%s", 1:nrow(free_peak))

write.table(free_peak, sprintf("%s", output), 
            col.names = F, row.names = F, quote = F, sep = "\t")