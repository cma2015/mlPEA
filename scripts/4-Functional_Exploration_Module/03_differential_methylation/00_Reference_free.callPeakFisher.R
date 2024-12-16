library(dplyr)
library(tidyr)
library(foreach)
library(doParallel)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(UpSetR)
library(ggplotify)
library(ggpubr)

setGeneric("transcriptBins", function(object) {
  standardGeneric("transcriptBins")
})

setMethod("transcriptBins", "MeRIP", function(object) {
  object@transcriptBins
})


# helper function 
.fisher_exact_test <- function(IP, input, IP_overall, input_overall, pseudo_count=1){
  
  test.m <- matrix(c(IP, input, IP_overall, input_overall), nrow = 2, byrow = FALSE,
                   dimnames = list(c("IP", "input"), c("bin", "overall")))
  
  # add a pseudo_count (1) to avoid zeros in denominators when calculating the odds ratio
  test.m <- test.m + pseudo_count
  
  fisher_test <- fisher.test(test.m, alternative="greater")
  
  fisher_result <- data.frame(pvalue = NA, odds_ratio = NA)
  
  fisher_result$pvalue <- fisher_test$p.value
  # fisher_result$odds_ratio <- fisher_test$estimate
  
  # use sample odds ratio (unconditional MLE) rather than the conditional Maximum Likelihood Estimate from Fisher's exact test
  a <- test.m[1,1] # IP in bin
  b <- test.m[1,2] # IP overall
  c <- test.m[2,1] # input in bin
  d <- test.m[2,2] # input overall
  
  # fisher_result$odds_bin <- a/c
  
  fisher_result$odds_ratio <- (a*d)/(b*c)
  
  return(fisher_result)
}

callPeakFisher <- function (MeRIP, min_counts = 15, peak_cutoff_fdr = 0.05, peak_cutoff_oddRatio = 1, 
          threads = 1) 
{
  if (!is(MeRIP, "MeRIP")) {
    stop("The input MeRIP must be a MeRIP dataset!")
  }else if (is(MeRIP, "MeRIP.Peak")) {
    cat(paste0("Input is an object of MeRIP.Peak, will override the peakCallResult by current call of fisher exact test!\n"))
  }else {
    cat("Performing fisher exact test on MeRIP dataset...\n")
  }
  input <- as.matrix(MeRIP@reads[, 1:length(MeRIP@samplenames)])
  ip <- as.matrix(MeRIP@reads[, (1 + length(MeRIP@samplenames)):(2 * length(MeRIP@samplenames))])
  colnames(input) <- colnames(ip) <- MeRIP@samplenames

  tmp <- as.data.frame(input)
  tmp$transcript <- rownames(tmp)
  transcriptBins <- tmp %>%
    separate(transcript, into = c("transcript", "bin"), sep = ",") %>%
    mutate(bin = as.integer(bin)) %>%
    dplyr::select(transcript, bin)
  rm(tmp)
  
  batch_id_list <- unique(transcriptBins$transcript)
  num_batch_ids <- length(batch_id_list)
  cat("Calling peaks for ", num_batch_ids, "transcripts... \n")
  registerDoParallel(cores = threads)
  start_time <- Sys.time()
  cat(paste("Hyper-thread registered:", getDoParRegistered(), 
            "\n"))
  cat(paste("Using", getDoParWorkers(), "thread(s) to call peaks in continuous bins...\n"))
  peak_call_batches <- foreach(i = 1:num_batch_ids, .combine = rbind) %dopar% 
    {
      idx_batch <- which(transcriptBins$transcript == batch_id_list[i])
      batch_input <- input[idx_batch, ]
      batch_ip <- ip[idx_batch, ]
      overall_input <- round(apply(batch_input, 2, median, na.rm = TRUE))
      overall_ip <- round(apply(batch_ip, 2, median, na.rm = TRUE))
      fisher_exact_test_p <- NULL
      fisher_exact_test_oddRatio <- NULL
      for (j in 1:length(overall_input)) {
        fisher_result <- t(mapply(.fisher_exact_test, batch_ip[, j], batch_input[, j], overall_ip[j], overall_input[j]))
        fisher_exact_test_p <- cbind(fisher_exact_test_p, fisher_result[, 1])
        fisher_exact_test_oddRatio <- cbind(fisher_exact_test_oddRatio, fisher_result[, 2])
      }
      above_thresh_counts <- ((batch_input + batch_ip) >= min_counts)
      fisher_exact_test_fdr <- matrix(1, nrow = nrow(fisher_exact_test_p), ncol = ncol(fisher_exact_test_p))
      if (sum(rowSums(above_thresh_counts) > (length(overall_input)/2)) > 1) {
        fisher_exact_test_fdr[rowSums(above_thresh_counts) > (length(overall_input)/2), ] <- 
          apply(fisher_exact_test_p[which(rowSums(above_thresh_counts) > (length(overall_input)/2)), ], 2, p.adjust, method = "fdr")
      }
      fisher_exact_test_peak <- (fisher_exact_test_fdr < peak_cutoff_fdr & 
                                   fisher_exact_test_oddRatio > peak_cutoff_oddRatio & 
                                   above_thresh_counts)
      fisher_exact_test_peak
    }
  rm(list = ls(name = foreach:::.foreachGlobals), pos = foreach:::.foreachGlobals)
  end_time <- Sys.time()
  cat(paste("Time used to call peaks:", difftime(end_time, 
                                                 start_time, units = "mins"), "mins... \n"))
  colnames(peak_call_batches) <- MeRIP@samplenames
  data.out <- as(MeRIP, "MeRIP.Peak")
  data.out@transcriptBins <- transcriptBins
  data.out@peakCallResult <- peak_call_batches
  data.out@peakCalling <- "fisher's exact test"
  return(data.out)
}


PlotBinFisherHeatmap <- function(MeRIP_peakCallResult,transcript){
  # plot heatmap for Merip fisher test result of all gene bin
  single_transcript_data <- MeRIP_peakCallResult[grep(transcript,rownames(MeRIP_peakCallResult)),]
  library(pheatmap)
  single_transcript_data[which(single_transcript_data==T)]<-1
  single_transcript_data[which(single_transcript_data==F)]<-0
  pheatmap::pheatmap(t(single_transcript_data),cluster_rows = F,cluster_cols = F,show_colnames = F,height = 5,width = 10,main = transcript)
}

