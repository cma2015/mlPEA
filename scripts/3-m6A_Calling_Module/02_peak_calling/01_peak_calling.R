library(argparse)

parser <- ArgumentParser()
parser$add_argument("-fasta_file", help="Fasta file", type="character", dest="fasta_file")
parser$add_argument("-inputfile", help="Input file", type="character", dest="inputfile")
parser$add_argument("-ipfile", help="IP file", type="character", dest="ipfile")
parser$add_argument("-outputfile", help="Output file", type="character", dest="outputfile")
parser$add_argument("-fc", help="Fold change", type="numeric", default=1.5, dest="fc")
parser$add_argument("-fdr", help="FDR", type="numeric", default=0.01, dest="fdr")
parser$add_argument("-readsCount", help="Reads count", type="numeric", default=10, dest="readsCount")
parser$add_argument("-sliding_step", help="Sliding step", type="numeric", default=25, dest="sliding_step")
parser$add_argument("-cpus", help="CPUs", type="numeric", default=5, dest="cpus")

args <- parser$parse_args()
# fc: log2(fold change)的阈值
# fdr: FDR的阈值
# readsCount: 至少支持的reads数
# 表示至少有多少个重复覆盖的peak才会被保留

fasta_file <- args$fasta_file
inputfile <- args$inputfile
ipfile <- args$ipfile
outputfile <- args$outputfile
fc <- args$fc
fdr <- args$fdr
readsCount <- args$readsCount
sliding_step <- args$sliding_step
cpus <- args$cpus

suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(future.apply))
suppressMessages(library(tidyr))
#improved_PEA for free peak

# call peak


fIPEATransciptPeak <- function(fasta_file, inputfile, ipfile, outputfile,fc=1.5, fdr=0.01, readsCount=10, sliding_step=25, cpus=15){

  library(GenomicRanges)
  library(future)
  library(future.apply)
  library(progress)
  library(rlang)
  library(data.table)
  
  ratio <- -fc
  level <- fdr
  readsCount <- readsCount
  concatenate = 4
  sliding_step <- sliding_step
  
  .winCount <- function(x, y, z){
    data.frame(
      "Position" = subset(z, !duplicated(x)),
      "window" = subset(x, !duplicated(x)),
      "curWave" = as.numeric(by(y, factor(x, levels = unique(x)), mean))
    )
  }
  
  .intervalCov <- function(dzFile){
    baseCov <- fread(file = dzFile, sep = "\t", header = F, stringsAsFactors = F)
    winDow <- rep(NA, nrow(baseCov))
    baseCov <- cbind(baseCov, winDow)
    colnames(baseCov) <- c("Chr", "Position", "readsNumber", "window")
    baseCov$Position <- baseCov$Position + 1
    baseCov$window <- ceiling(baseCov$Position/as.numeric(sliding_step))
    
    baseCovOutput <- baseCov %>% 
      group_by(Chr) %>% 
      do(.winCount(x=.$window, y=.$readsNumber, z=.$Position))
    
    baseCovOutput <- split(baseCovOutput, factor(baseCovOutput$Chr, levels = unique(baseCovOutput$Chr)))
    baseCovOutput
  }
  
  .getPvalue <- function(inputVec, mappedInput, mappedRIP){
    testMat <- matrix(c(as.numeric(inputVec[8]),
                        mappedInput,
                        as.numeric(inputVec[9]),
                        mappedRIP), nrow = 2, ncol = 2)
    p.value <- fisher.test(x = testMat)$p
    ratio <- log2(((as.numeric(inputVec[8]) + 1)*mappedRIP)/((as.numeric(inputVec[9]) + 1)*mappedInput))
    res <- c(ratio, p.value)
    res
  }
  
  .findContinuous <- function(inputVec){
    Breaks <- c(0, which(diff(inputVec) != 1), length(inputVec))
    res <- lapply(seq(length(Breaks) - 1),
                  function(i) inputVec[(Breaks[i] + 1):Breaks[i+1]])
    res
  }
  
  # sequences and its lengths of transfrags
  transfrags <- seqinr::read.fasta(fasta_file, seqtype = "DNA", as.string = T, forceDNAtolower = F)
  seq_length <- sapply(transfrags, nchar)
  
  input <- .intervalCov(dzFile = inputfile)
  RIP <- .intervalCov(dzFile = ipfile)
  
  if(!all(names(input) == names(RIP))){
    cat("Note: The chromosomes in the input and RIP are not consistent!\n",
        "the interactions will be used!")
    interNames <- intersect(names(input), names(RIP))
    input <- input[interNames]
    RIP <- RIP[interNames]
  }
  
  resList <- vector("list", length = length(RIP))
  names(resList) <- names(RIP)
  
  cl <- parallel::makeCluster(cpus)
  plan(cluster, workers = cl)
  system.time(resList <- future_lapply(names(input), function(x){
    curMat <- merge(input[[x]], RIP[[x]], by = 'window', all=TRUE)
    curMat$curWave.x[is.na(curMat$curWave.x)] <- 0
    curMat$curWave.y[is.na(curMat$curWave.y)] <- 0
    curMat <- curMat %>%
      mutate(windowave.input = round(curWave.x),
             windowave.RIP = round(curWave.y)) %>%
      filter(windowave.RIP != 0, windowave.input != 0) # %>%
    
    if (nrow(curMat) == 0) return(NULL)
    
    ## 对剩余行求各自的中值，每个窗口中的 read count通过该基因的所有窗口的中位数计数进行标准化
    median_input <- median(curMat$windowave.input)
    median_ip <- median(curMat$windowave.RIP)
    # 检查中值是否为 NA，如果是，则跳过当前序列
    if (is.na(median_input) || is.na(median_ip)) {
      return(NULL)
    }
    
    curPvalue <- rep(NA, nrow(curMat))
    curRatio <- rep(NA, nrow(curMat))
    curFDR <- rep(NA, nrow(curMat))
    curMat <- cbind(curMat, curPvalue, curRatio, curFDR)
    pvalue <- t(apply(curMat, 1, .getPvalue, mappedInput = median_input, mappedRIP = median_ip))
    curMat$curPvalue <- as.numeric(pvalue[,2])
    curMat$curRatio <- as.numeric(pvalue[,1])
    curMat$curFDR <- p.adjust(curMat$curPvalue, "fdr")
    curMat
  }))
  # future:::ClusterRegistry("stop")
  parallel::stopCluster(cl)
  
  resMat <- do.call(what = rbind, args = resList)
  resMat <- resMat %>%
    filter(curRatio < ratio, curFDR < level, curWave.y >= readsCount)
  
  tt <- by(resMat$window, factor(x = resMat$Chr.y, levels = unique(resMat$Chr.y)), .findContinuous)
  resPeaks <- NULL
  
  pb <- progress_bar$new(total = length(tt), clear = FALSE)
  
  for(i in 1:length(tt)){
    pb$tick()
    curList <- tt[[i]]
    curLen <- unlist(lapply(curList, length))
    curList <- curList[which(curLen >= concatenate)]
    if(length(curList) == 0){
      next
    }
    Start <- unlist(lapply(curList, function(x) (x[1]-1)*as.numeric(sliding_step)+1))
    End <- unlist(lapply(curList, function(x) (x[length(x)]*as.numeric(sliding_step))))
    curMat <- subset(resMat, resMat$Chr.y == names(tt)[i])
    curFDR <- lapply(curList, function(x)  curMat$curFDR[match(x, curMat$window)])
    meanFDR <- unlist(lapply(curFDR, mean))
    maxFDR <- unlist(lapply(curFDR, max))
    minFDR <- unlist(lapply(curFDR, min))
    curRatio <- lapply(curList, function(x)  curMat$curRatio[match(x, curMat$window)])
    meanRatio <- unlist(lapply(curRatio, mean))
    maxRatio <- unlist(lapply(curRatio, max))
    minRatio <- unlist(lapply(curRatio, min))
    windowNumber <- unlist(lapply(curList, length))
    curPeaks <- cbind(names(tt)[i], Start, End, windowNumber,
                      meanFDR, maxFDR, minFDR,
                      meanRatio, maxRatio, minRatio)
    resPeaks <- rbind(resPeaks, curPeaks)
  }
  colnames(resPeaks) <- c("Chromosome", "Start(1-based)", "End", "Bin number",
                          "Mean FDR", "Max FDR", "Minimum FDR",
                          "Mean Ratio", "Max Ratio", "Minimum Ratio")
  ## peaks
  peaks_granges <- GRanges(seqnames = resPeaks[,1], 
                           ranges = IRanges(start = as.numeric(resPeaks[,2]), 
                                            end = as.numeric(resPeaks[,3])),
                           mean_FDR = resPeaks[,5],
                           max_FDR = resPeaks[,6],
                           min_FDR = resPeaks[,7],
                           mean_ratio = resPeaks[,8],
                           max_ratio = resPeaks[,9],
                           min_ratio = resPeaks[,10])
  
  peak <- data.frame(seqnames = seqnames(peaks_granges),
                     start = start(peaks_granges),
                     end = end(peaks_granges),
                     mean_FDR = peaks_granges$mean_FDR,
                     max_FDR = peaks_granges$max_FDR,
                     min_FDR = peaks_granges$min_FDR,
                     mean_ratio = peaks_granges$mean_ratio,
                     max_ratio = peaks_granges$max_ratio,
                     min_ratio = peaks_granges$min_ratio)
  
  write.table(peak, file = sprintf("%s", outputfile), sep = "\t", row.names = FALSE, quote = FALSE)
  return(peaks_granges)
}


peaks_granges <- fIPEATransciptPeak(fasta_file, inputfile, ipfile, outputfile, fc,fdr,readsCount,sliding_step, cpus)