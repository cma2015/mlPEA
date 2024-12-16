library(MeRIPtools)

.noZero <- function(x){sapply(x,max,1)}

MergeRep_Calculate_OR <- function(MergeRep_jointPeakCount,samples,rep_number,rep_filter_threshold){
  LCLs <- MergeRep_jointPeakCount
  newLCLs <- MeRIPtools::filter(LCLs, !apply(extractInput(LCLs), 1, function(x) any(x == 0 )))
  T0 <- colSums(counts(newLCLs)[,1:length(newLCLs@samplenames)] )
  T1 <- colSums(counts(newLCLs)[,(length(newLCLs@samplenames)+1) : (2*length(newLCLs@samplenames)) ] )
  test_OR <- t(t(extractIP(newLCLs))/T1 )/t(t(extractInput(newLCLs) )/T0)
  test_OR_merge <- sapply(seq(1, ncol(test_OR), rep_number), function(j){
    q <- rowSums(test_OR[,j+(0:(rep_number-1))]>1)
    q[which(q<rep_filter_threshold)] <- FALSE
    q[which(q>=rep_filter_threshold)] <- TRUE
    return(as.logical(q))
  })
  rownames(test_OR_merge) <- rownames(test_OR)
  colnames(test_OR_merge) <- samples

  enrichFlag <- apply(test_OR_merge,1,function(x){sum(x)>= newLCLs@jointPeak_threshold})
  newLCLs <-  MeRIPtools::filter(newLCLs, enrichFlag)

  OR <- t(apply(extractIP(newLCLs),1,.noZero)/T1)/ t(t(extractInput(newLCLs))/T0 )
  OR <- sapply(seq(1, ncol(OR), rep_number), function(j) round(rowMeans(OR[, j + (0:(rep_number-1))]), 6))

  colnames(OR) <- samples

  newLCLs@transcriptBins$bin_num <- 1:nrow(newLCLs@transcriptBins)

  return(list(OR,newLCLs@jointPeaks,newLCLs@transcriptBins,newLCLs@jointPeak_id_pairs))

}


MergeRep_OR_Heantmap <- function(MergeRep_OR,MergeRep_transcriptBins,MergeRep_jointPeak_id_pairs,transcriptid){
  transcript_bin <- MergeRep_transcriptBins[which(MergeRep_transcriptBins$transcript==transcriptid),]
  
  transcript_OR <- MergeRep_OR[MergeRep_OR$transcript %in% transcript_bin$transcript,]
  transcript_OR <- transcript_OR[, -c(1,2)]
  heatmap_matrix <- matrix(0,nrow = dim(transcript_OR)[2],ncol = dim(transcript_bin)[1])
  
  rownames(heatmap_matrix) <- colnames(transcript_OR)
  colnames(heatmap_matrix) <- transcript_bin$bin
  
  
  transcript_jointPeak <- MergeRep_jointPeak_id_pairs[grepl(unique(transcript_bin$transcript),rownames(MergeRep_jointPeak_id_pairs)),]
  bin_idx <- lapply(1:nrow(transcript_jointPeak), function(i) {
    which(transcript_bin$bin_num %in% transcript_jointPeak[i, 1]:transcript_jointPeak[i, 2])
  })
  
  for (j in 1:nrow(transcript_OR)) {
    heatmap_matrix[,bin_idx[[j]]] <- as.numeric(transcript_OR[j,])
  }
  
  return(heatmap_matrix)
}


