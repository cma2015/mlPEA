.getPeakBins <- function(transcriptGRList,transcriptName,slidingStarts,binSize){
  
  transcriptModel =reduce( transcriptGRList[transcriptName][[1]] )## merge overlapping exons
  dna.range = as.data.frame(range(transcriptModel) )
  
  exon.length = width(transcriptModel)
  #creat a corresponding map from RNA to DNA
  RNA2DNA <- start(transcriptModel):(start(transcriptModel) + exon.length - 1)
  
  #create center points of continuous window
  # Generate part of a mapping of DNA positions to RNA positions
  if(exon.length <= binSize | (dna.range$strand == "+" & as.numeric(slidingStarts[2]) == ( exon.length - binSize - exon.length %% binSize + 1) ) ){ # for genes that is shorter than bin size || if it is the rightmost bin (buffer bin) to be retrieved (positive strand only)
    mapping = data.frame(start = RNA2DNA[ slidingStarts[1] ], end = RNA2DNA[exon.length]  )
  }else if( dna.range$strand == "-" & as.numeric(slidingStarts[2]) == 1 ){ # when retrieve the leftmost (buffer) bin on reverse strand
    mapping = data.frame(start = RNA2DNA[ slidingStarts[1] ], end = RNA2DNA[slidingStarts[2] + binSize + exon.length %% binSize - 1 ]  )
  }else{ # retrieve regular bins
    mapping = data.frame(start = RNA2DNA[ slidingStarts[1] ], end = RNA2DNA[slidingStarts[2] + binSize - 1 ] )
  }
  
  mapping$transcript = as.character(dna.range$seqnames)
  return(mapping[,c("transcript","start","end")])
  
}


.peakExons <- function(peak,y){
  exonID <- peak$start <= y$end & peak$end >= y$start
  if(sum(exonID) == 1){
    return(data.frame(start = peak$start, end = peak$end, width = peak$end - peak$start + 1))
  }else if(sum(exonID) > 1){
    peakexon <- y[exonID,]
    peakexon[1,"start"] <- peak$start
    peakexon[sum(exonID),"end"] <- peak$end
    return(data.frame(start = peakexon$start, end = peakexon$end, width = peakexon$end - peakexon$start + 1))
  }
}

## function
# 
MergePeakCallResult <- function(MeRIP_peakCallResult,samples,rep_number=3,rep_filter_threshold=2){
  # Keep at least {rep_filter_threshold} significant bin of {rep_number} repeats
  peakCallResult_merge <- sapply(seq(1, ncol(MeRIP_peakCallResult), rep_number), function(j){
    q <- rowSums(MeRIP_peakCallResult[,j+(0:(rep_number-1))])
    q[which(q<rep_filter_threshold)] <- FALSE
    q[which(q>=rep_filter_threshold)] <- TRUE
    return(as.logical(q))
  })
  
  rownames(peakCallResult_merge) <- rownames(MeRIP_peakCallResult)
  colnames(peakCallResult_merge) <- samples
  
  return(peakCallResult_merge)
}

reportJointPeak_MergeRep <- function(joint_threshold=1,MeRIP_peakCallResult,MergeRep_PeakCallResult,threads){
  cat(paste0("Reporting joint peak at joint threshold "),joint_threshold, "\n")
  
  transcriptBins <- MeRIP_peakCallResult@transcriptBins
  ID <- (rowSums(MergeRep_PeakCallResult) >= joint_threshold)
  num_lines <- length(ID)
  # Find the starting position of the joint peak that meets certain conditions
  start_id <- which((ID[2:num_lines] - ID[1:num_lines - 1] == 1) | ((transcriptBins$transcript[2:num_lines] != transcriptBins$transcript[1:num_lines - 1]) & (ID[2:num_lines] == TRUE)))
  start_id <- start_id + 1
  if (ID[1] == TRUE) {
    start_id <- c(1, start_id)
  }
  # Find the ending position of the joint peak that meets certain conditions
  end_id <- which((ID[1:num_lines - 1] - ID[2:num_lines] == 1) | ((transcriptBins$transcript[1:num_lines - 1] != transcriptBins$transcript[2:num_lines]) & (ID[1:num_lines - 1] == TRUE)))
  if (ID[num_lines] == TRUE) {
    end_id <- c(end_id, num_lines)
  }
  peak_id_pairs <- cbind(start_id, end_id)
  num_peaks <- nrow(peak_id_pairs)
  transcriptGRList <- MeRIP_peakCallResult@transcriptModel
  peakGenes <- as.character(transcriptBins[peak_id_pairs[, 1],"transcript"])
  if (num_peaks == 0) {
    return(data.frame())
  }else {
    start_time <- Sys.time()
    registerDoParallel(cores = threads)
    cat(paste("Hyper-thread registered:", getDoParRegistered(),"\n"))
    cat(paste("Using", getDoParWorkers(), "thread(s) to report merged report...\n"))
    merged.report <- foreach(p = 1:num_peaks, .combine = rbind) %dopar%{
      peak_row_id <- peak_id_pairs[p, ]
      geneExons <- reduce(transcriptGRList[peakGenes[p]][[1]])
      ## Obtain the starting and ending locations of peak on the genome
      peak <- .getPeakBins(transcriptGRList, peakGenes[p],c(transcriptBins$bin[peak_row_id[1]], transcriptBins$bin[peak_row_id[2]]),MeRIP_peakCallResult@binSize)
      ## The location of the exon with peak on the genome
      peakE <- .peakExons(peak, as.data.frame(geneExons))
      data.frame(transcript = peak$transcript, 
                 start = peak$start,
                 end = peak$end, 
                 name = peakGenes[p], 
                 score = 0,
                 strand = as.character(strand(geneExons))[1],
                 thickStart = peak$start, 
                 thickEnd = peak$end,
                 itemRgb = 0, 
                 blockCount = nrow(peakE), 
                 blockSizes = paste(peakE$width,collapse = ","), 
                 blockStarts = paste(peakE$start -replicate(nrow(peakE), peakE$start[1]),collapse = ","))
    }
    rm(list = ls(name = foreach:::.foreachGlobals),
       pos = foreach:::.foreachGlobals)
    end_time <- Sys.time()
    cat(paste("Time used to report peaks:", difftime(end_time,
                                                     start_time, units = "mins"), "mins... \n"))
  }
  data.out <- MeRIP_peakCallResult
  data.out@jointPeaks <- merged.report
  data.out@jointPeak_id_pairs <- peak_id_pairs
  data.out@jointPeak_threshold <- joint_threshold
  if (MeRIP_peakCallResult@jointPeak_threshold != joint_threshold & MeRIP_peakCallResult@jointPeak_threshold > 0) {
    cat(paste0("joint threshold was previous set at ",MeRIP_peakCallResult@jointPeak_threshold, ".\nWill remove joint-peak read count,test statistics etc. to avoid inconsistency with new joint peaks. Please re-fetch joint peak counts!\n"))
    data.out@jointPeak_ip <- new("matrix")
    data.out@jointPeak_input <- new("matrix")
    data.out@norm.jointPeak_ip <- new("matrix")
    data.out@jointPeak_adjExpr <- new("matrix")
    data.out@test.est <- new("matrix")
    data.out@test.method <- "none"
  }
  
  # jointPeaks <- data.out@jointPeaks
  # # 对 jointPeaks 进行去重处理
  # unique_jointPeaks <- jointPeaks[!duplicated(paste(jointPeaks$chr,
  #                                                   jointPeaks$start, jointPeaks$end,
  #                                                   jointPeaks$strand, sep = ":")), ]
  # # 将去重后的结果重新赋值回 data.out@jointPeaks
  # data.out@jointPeaks <- unique_jointPeaks
  
  ##Filter out duplicated peaks due to duplicated gene names
  data.out <- MeRIPtools::filter(data.out, !duplicated(paste(data.out@jointPeaks$transcript,
                                                             data.out@jointPeaks$start, data.out@jointPeaks$end,
                                                             data.out@jointPeaks$strand, sep = ":")))
  return(data.out)
}

jointPeakCount <- function(transcriptBins,MergeRep_JointPeak){
  
  #transcriptBins <- MeRIP_peakCallResult@transcriptBins
  
  peak_id_pairs <- MergeRep_JointPeak@jointPeak_id_pairs
  
  
  ip <- MergeRep_JointPeak@reads[,(1+length(MergeRep_JointPeak@samplenames)):(2*length(MergeRep_JointPeak@samplenames))]
  input <- MergeRep_JointPeak@reads[,1:length(MergeRep_JointPeak@samplenames)]
  
  joint_peak_ip <- t( apply(peak_id_pairs,1,function(x,y){
    if(x[1]==x[2]){
      return(y[x[1]:x[2],])
    }else{
      colSums(y[x[1]:x[2],])
    }
  },y = ip) )
  
  joint_peak_input <- t( apply(peak_id_pairs,1,function(x,y){
    if(x[1]==x[2]){
      return(y[x[1]:x[2],])
    }else{
      colSums(y[x[1]:x[2],])
    }
  },y = input) )
  
  data.out <- MergeRep_JointPeak
  data.out@jointPeak_ip <- joint_peak_ip
  data.out@jointPeak_input <- joint_peak_input
  
  return(data.out)
  
}
