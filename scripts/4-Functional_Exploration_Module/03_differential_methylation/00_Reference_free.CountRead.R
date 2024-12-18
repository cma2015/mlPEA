library(GenomicFeatures)
library(Rsamtools)
library(doParallel)

#' @export
MeRIP <- setClass("MeRIP",
                  representation( reads = "matrix",
                                  binSize = "numeric",
                                  transcriptModel = "GRangesList",
                                  bamPath.input = "character",
                                  bamPath.ip = "character",
                                  samplenames = "character",
                                  transcriptBins = "data.frame",
                                  transcriptSum = "matrix",
                                  mode = "character"
                  ),
                  validity = function(object){
                    errors <- character()
                    ## Check basic element
                    if( all(dim(object@reads) <= 1) ){
                      errors <- c(errors, paste0("read count table is a required data for MeRIP instance!"))
                    }
                    if( length(object@binSize)==0){
                      errors <- c(errors, paste0("bin size to slice the transcript is a required parameter for MeRIP instance!"))
                    }
                    if( length(object@transcriptModel)<=0 ){
                      errors <- c(errors, paste0("transcriptModel is a required data for MeRIP instance!"))
                    }
                    if( length(object@bamPath.input)<=0){
                      errors <- c(errors, paste0("path to the input BAM files are the required data for MeRIP instance!"))
                    }
                    if( length(object@bamPath.ip)<=0){
                      errors <- c(errors, paste0("path to the IP BAM files are the required data for MeRIP instance!"))
                    }
                    if( length(object@samplenames)<=0){
                      errors <- c(errors, paste0("Sample names are the required data for MeRIP instance!"))
                    }
                    ## match data dimension with num of samples
                    if(ncol(object@reads) !=  2 * length(object@samplenames) ){
                      errors <- c(errors, paste0("The number of colnumns of read count is ",ncol(object@reads),". The number of samples is ",length(object@samplenames),". The number of colnumns of read count should be 2x the number of samples!"))
                    }
                    if(length(object@bamPath.input) != length(object@samplenames) ){
                      errors <- c(errors, paste0("The number of input bam file path(es) should equal to the number of samplenames!"))
                    }
                    if(length(object@bamPath.ip) != length(object@samplenames) ){
                      errors <- c(errors, paste0("The number of IP bam file path(es) should equal to the number of samplenames!"))
                    }
                    
                    if (length(errors) == 0) TRUE else errors
                  },
                  prototype(transcriptBins = data.frame(transcript = character(), bin = character() ), transcriptSum = matrix(), mode = "mRNA" )
)

#' @export
MeRIP.Peak <- setClass("MeRIP.Peak",representation( peakCallResult = "matrix",
                                                    jointPeak_id_pairs = "matrix",
                                                    jointPeaks = "data.frame",
                                                    jointPeak_ip = "matrix",
                                                    jointPeak_input = "matrix",
                                                    norm.jointPeak_ip = "matrix",
                                                    sizeFactor = "data.frame",
                                                    variate = "data.frame",
                                                    jointPeak_adjExpr = "matrix",
                                                    test.est = "matrix",
                                                    peakCalling = "character",
                                                    jointPeak_threshold = "numeric",
                                                    test.method = "character"),
                       contains = "MeRIP",
                       prototype(peakCalling = "none", jointPeak_threshold = 0, test.method = "none")
)



## a helper function to count reads from bam files
.countReadFromBam <-
  function(bam,which,DNA2RNA,reads.strand,fragmentLength,left,sliding, binSize ,paired){
    
    if(paired){  # pair-end sequence
      ba = scanBam(bam, param=ScanBamParam( which=which, what =c("pos","strand","qwidth","isize"), flag = scanBamFlag(isProperPair = T,isSecondMateRead = F) ) )
      ba = data.frame( pos = ba[[1]]$pos, strand = ba[[1]]$strand ,readLength = ba[[1]]$qwidth, isize = ba[[1]]$isize )
      ba = ba[ ba$pos > left, ]
      ba = ba[ba$strand == reads.strand | reads.strand == "*", ] ## filter for strand
      ba$pos = DNA2RNA[ba$pos - left] ## convert mapped read pos into RNA position
      ba = ba[ ba$pos > 0, ] ## drop intron reads.
      ##shift the read pos to the center of the reads
      if(reads.strand == "+"){ba$pos = ba$pos + round(abs(ba$isize)/2) }else{ba$pos = ba$pos + ba$readLength - round(abs(ba$isize)/2)  }
      ##count the reads in the sliding windows
      no.window = length(sliding)
      windowCounts = vector(length = no.window)
      for(j in 1:no.window){
        windowCounts[j]= sum( ba$pos >= sliding[j] &  ba$pos < (sliding[j] + binSize)  ) 
      } # count regular bins
      
      if( as.character( strand(which) ) == "+" ){
        windowCounts[no.window] = windowCounts[no.window] + sum( ba$pos >=  (sliding[no.window] + binSize) & ba$pos <= max(DNA2RNA) ) #count the extra part of the last bin on positive strand
      }else{
        windowCounts[1] = windowCounts[1] + sum( ba$pos >=  (sliding[1] + binSize) & ba$pos <= (sliding[1] + binSize + max(DNA2RNA) %% binSize ) ) #count the extra part of the first bin on negative strand
      }
      
      
    }else{ # single-end sequence
      ba = scanBam(bam, param=ScanBamParam( which=which, what =c("pos","strand","qwidth") ) )
      ba = data.frame( pos = ba[[1]]$pos, strand = ba[[1]]$strand, readLength = ba[[1]]$qwidth )
      ba = ba[ ba$pos > left, ]
      ba = ba[ba$strand == reads.strand | reads.strand == "*", ] ## filter for strand
      ba$pos = DNA2RNA[ba$pos - left] ## convert mapped read pos into RNA position
      ba = ba[ ba$pos > 0, ] ## drop intron reads.
      ##shift the read pos to the center of the reads
      if(reads.strand == "+"){ba$pos = ba$pos + round(fragmentLength/2) }else{ba$pos = ba$pos + ba$readLength - round(fragmentLength/2)  }
      ##count the reads in the sliding windows
      no.window = length(sliding)
      windowCounts = vector(length = no.window)
      for(j in 1:no.window){
        windowCounts[j]= sum( ba$pos >= sliding[j] &  ba$pos < (sliding[j] + binSize)  ) 
      } # count regular bins
      
      if( as.character( strand(which) ) == "+" ){
        windowCounts[no.window] = windowCounts[no.window] + sum( ba$pos >=  (sliding[no.window] + binSize) & ba$pos <= max(DNA2RNA) ) #count the extra part of the last bin on positive strand
      }else{
        windowCounts[1] = windowCounts[1] + sum( ba$pos >=  (sliding[1] + binSize) & ba$pos <= (sliding[1] + binSize + max(DNA2RNA) %% binSize ) ) #count the extra part of the first bin on negative strand
      }
      
      
    }
    
    return(windowCounts)
  }


## main function that takes a fasta file and bam files to count read for continuous windows
countReads <- function(
    IP_filenames,# file name of ip
    Input_filenames,# file name of input
    IP_colnames, # name of ip output column
    Input_colnames, # name of input output column
    fasta_fliename, # fasta file used for peak calling
    samplenames, # name of sample
    fragmentLength = 150,
    bamFolder,
    outputDir=NA,
    binSize = 50,
    strandToKeep = "opposite",
    paired = FALSE,
    threads = 1,
    saveOutput = FALSE
){
  
  ##read bam files
  bamPath.input = paste0(bamFolder,"/",Input_filenames)
  bamPath.IP = paste0(bamFolder,"/",IP_filenames)
  
  ## Check for missing files and index bam files
  if( !all(file.exists(bamPath.input)) ) stop( "input bam file missing!!!" )
  if( !all(file.exists(bamPath.IP)) ) stop( "IP bam file missing!!!" )
  num_bam_files <- length(bamPath.input)
  for (ibam in 1:num_bam_files) {
    inputfile = bamPath.input[ibam]
    IPfile = bamPath.IP[ibam]
    if (! file.exists(paste(inputfile,'.bai',sep=""))) {
      print(paste("Stage: index bam file", inputfile))
      indexBam(inputfile)
    }
    if (! file.exists(paste(IPfile,'.bai',sep=""))) {
      print(paste("Stage: index bam file", IPfile))
      indexBam(IPfile)
    }
  }
  
  
  
  ## This step removes ambiguous annotations and returns transcript model
  cat("Reading fasta file to obtain transcript model\nFilter out ambiguous model...\n")
  transcript_fasta <- seqinr::read.fasta(fasta_fliename)
  transcript_fasta_df <- as.data.frame(sapply(transcript_fasta,length))
  colnames(transcript_fasta_df) <- "transcript_length"
  transcript_fasta_df$transcript_name <- rownames(transcript_fasta_df)
  transcript_fasta_gr <- GRanges(seqnames = transcript_fasta_df$transcript_name, ranges = IRanges(start = 1, end = transcript_fasta_df$transcript_length), strand = "+")
  transcriptGRList <- split(transcript_fasta_gr, seqnames(transcript_fasta_gr))
  cat("transcript model obtained from fasta file...\n")
  
  ## Check BAM headers and remove chr in transcriptModel that is not in BAM file. 
  bamHeader <- scanBamHeader(bamPath.input, what=c("targets") )
  seqLevels <- unique( unlist( lapply( bamHeader, function(x) names( x$targets) ) ) )
  transcriptGRList <- transcriptGRList[ unlist( runValue( seqnames( transcriptGRList ) ) ) %in% seqLevels ]
  
  
  no.transcripts=length(transcriptGRList)## define number of transcripts
  
  cat("counting reads for each transcripts, this step may takes a few hours....\n")
  start_time <- Sys.time()
  registerDoParallel( cores = threads)
  cat(paste("Hyper-thread registered:",getDoParRegistered(),"\n"))
  cat(paste("Using",getDoParWorkers(),"thread(s) to count reads in continuous bins...\n"))
  reads <- foreach(i = 1:no.transcripts, .combine = rbind) %dopar%{
    
    transcriptName = names(transcriptGRList)[i]
    print(paste0("Counting Reads for: ",transcriptName))
    transcriptModel =IRanges::reduce( transcriptGRList[transcriptName][[1]] )## merge overlapping exons
    
    # DNA location to transcript location conversion
    df.transcriptModel= as.data.frame(transcriptModel) ##data frame of transcript model
    dna.range = as.data.frame(range(transcriptModel))
    
    DNA2RNA = rep(1,dna.range$end - dna.range$start +1)
    exon.length = sum(DNA2RNA)
    DNA2RNA=cumsum(DNA2RNA)*DNA2RNA
    
    ## skip any transcript with length less than 200
    if(exon.length < 200) {return(NULL)}
    
    ## switch strand because stranded RNA library protocol sequence reverse strand
    if(strandToKeep == "opposite"){
      reads.strand = character()
      if(dna.range$strand == "+"){reads.strand = "-"}else{reads.strand = "+"} ## switch strand on RNA reads for Truseq protocol
    }else if(strandToKeep == "same"){
      reads.strand = as.character(dna.range$strand)
    }else{
      cat("Currently m6Amonter only support strand specific RNA-seq data.\nCounting reads at opposite strand by defalt...\n")
      reads.strand = character()
      if(dna.range$strand == "+"){reads.strand = "-"}else{reads.strand = "+"}
    }
    
    #create start points of continuous window
    if(exon.length <= binSize){
      slidingStart = 1
    }else{
      ## use the 3' end terminal bin as a elastic-size bin
      if(dna.range$strand == "+"){
        slidingStart = seq(from = 1, to = ( exon.length - binSize - exon.length %% binSize + 1) , length.out = floor(exon.length/binSize) ) 
      }else{ # make the first bin elastic bin if a transcript is on reverse strand
        slidingStart = c(1, seq(from = binSize + exon.length %% binSize + 1, to = ( exon.length - binSize + 1) , length.out = floor(exon.length/binSize) - 1 )  )
      }
    }
    
    
    #count reads in all samples
    ba.IP = sapply(bamPath.IP,.countReadFromBam,which = range(transcriptModel),reads.strand = reads.strand,DNA2RNA = DNA2RNA,fragmentLength=fragmentLength,left=dna.range$start,sliding = slidingStart, binSize = binSize, paired = paired)
    ba.input = sapply(bamPath.input,.countReadFromBam,which = range(transcriptModel),reads.strand = reads.strand,DNA2RNA = DNA2RNA,fragmentLength=fragmentLength,left=dna.range$start,sliding = slidingStart, binSize = binSize, paired = paired)
    
    if(is.vector(ba.IP) ){# if there is only one window for this transcript, make it a matrix to avoid bug
      ba.IP = matrix(ba.IP, nrow = 1 )
      ba.input = matrix( ba.input, nrow = 1 )
    }
    ba.counts <- cbind(ba.input,ba.IP)
    rownames(ba.counts) <-  paste(transcriptName,slidingStart,sep = ",")
    
    ba.counts
  }
  rm(list=ls(name=foreach:::.foreachGlobals), pos=foreach:::.foreachGlobals)
  end_time <- Sys.time()
  cat(paste("Time used to count reads:",difftime(end_time, start_time, units = "mins"),"mins... \n"))
  
  colnames(reads) <- c(Input_colnames,IP_colnames)
  
  
  data.out <- MeRIP(reads = reads, binSize = binSize, transcriptModel = transcriptGRList, bamPath.input = bamPath.input, bamPath.ip = bamPath.IP, samplenames = samplenames)
  if(saveOutput){
    ## create output directory
    dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)
    saveRDS(data.out,paste0(outputDir,"/Reference_free_readCounts.RDS"))
  }
  
  
  return(data.out)
}