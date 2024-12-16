############################################################
##general functions in ML
fPUbaggingPSOL <- function(transfrags_id_path, model_dir, positive_sample, feature, iterator, cpus, bagging_model = "randomForest") {
  transfrags_id <- scan(file = transfrags_id_path ,what = character())
  psolResDic <- model_dir
  
  feature <- feature
  
  # Execute PU bagging
  pu_python_input <- sprintf("%s/PU_input.txt",model_dir)

  write.table(feature, file = pu_python_input, quote = F, sep = "\t")

  
  pu_python_output <- sprintf("%s/PU_output.txt",model_dir)
  
  if (bagging_model == "randomForest") {
      PU_command <- sprintf("/home/galaxy/miniconda3/envs/tf26/bin/python %s/03_PUlearning_RF.py %s %s",
                            script_dir, pu_python_input, pu_python_output)
   }else if (bagging_model == "svm") {
      PU_command <- sprintf("/home/galaxy/miniconda3/envs/tf26/bin/python %s/03_PUlearning_SVM.py %s %s",
                            script_dir, pu_python_input, pu_python_output)
  }else if (bagging_model == "XGBoost") {
      PU_command <- sprintf("/home/galaxy/miniconda3/envs/tf26/bin/python %s/03_PUlearning_XGBoost.py %s %s",
                            script_dir, pu_python_input, pu_python_output)
  }else if (bagging_model == "decision_tree") {
      PU_command <- sprintf("/home/galaxy/miniconda3/envs/tf26/bin/python %s/03_PUlearning_DT.py %s %s",
                            script_dir, pu_python_input, pu_python_output)
  } else if (bagging_model == "Logistic_regression") {
      PU_command <- sprintf("/home/galaxy/miniconda3/envs/tf26/bin/python %s/03_PUlearning_LR.py %s %s",
                            script_dir, pu_python_input, pu_python_output)
  } else {
      stop("Error: bagging_model should be one of Random Forest, Support Vector Machine, XGBoost, Decision Tree, and Logistic Regression.")
  }
  system(command = PU_command)
  
  positiveSamples <- positive_sample
  unlabelSamples <- setdiff(transfrags_id, positiveSamples)
  
  featureMat <- feature[,-ncol(feature)]
  positiveSamples <- positiveSamples[positiveSamples %in% rownames(featureMat)]
  unlabelSamples <- setdiff(rownames(featureMat), positiveSamples)
  
  pu_result <-  read.table(sprintf("%s/PU_output.txt",model_dir), sep = ",",row.names = 1)
  pre_neg <- dplyr::arrange(pu_result, V2)
  
  set.seed(9729)
  col_conditions <- colSums(pu_result == 0) <= length(positiveSamples)
  
  if (all(col_conditions)) {
    initial_negatives <- rownames(pre_neg)[1:length(positiveSamples)]
  } else {
    eligible_negatives <- rownames(pre_neg)[table(pre_neg$V2 == 0)[2]]
    initial_negatives <- sample(eligible_negatives, length(positiveSamples))
  }
  
  unlabel_new <- setdiff(unlabelSamples,initial_negatives)
  
  res <- list(positives = positiveSamples, negatives = initial_negatives,
              unlabels = unlabel_new)
  save(res, file = sprintf("%s/PSOL_InitialNegativeSelection_Res.RData",psolResDic))
  
  ####PSOL Negative Exapansion
  #' @export
  .PSOL_NegativeExpansion <- function (featureMat, positives, negatives, unlabels, cpus = 1, 
                                       iterator = 50, cross = 5, TPR = 0.98, method = "randomForest", plot = TRUE, trace = TRUE, PSOLResDic, ...) 
  {
    call <- match.call()
    if (length(method) > 1) {
      method <- method[1]
    }
    dir.create(path = PSOLResDic, showWarnings = FALSE)
    if (TPR > 1 | TPR < 0) 
      stop("Error: TPR is a value ranged from 0 to 1.")
    if (is.null(rownames(featureMat))) 
      stop("Error: rownames should be given to featureMat")
    thresholdIdx <- floor((1 - TPR) * length(positives))
    if (thresholdIdx <= 0) 
      thresholdIdx <- 1
    finalUnlabels <- unlabels
    finalNegatives <- negatives
    negCount <- length(finalNegatives)
    numMat <- matrix(0, nrow = iterator, ncol = 5)
    rownames(numMat) <- paste("Iter", 1:iterator, sep = "")
    colnames(numMat) <- c("IterNo", "AUC_On_TrainingDataSet", 
                          "AUC_On_TestingDataSet", "Negative_Sample_Num", "Unlabeled_Sample_Num")
    numMat[, 1] <- 1:iterator
    if (length(setdiff(c(positives, negatives, unlabels), rownames(featureMat))) > 
        0) {
      stop("Error: some samples not included in the featureMat.\n")
    }
    featureMat <- featureMat[c(positives, negatives, unlabels), 
    ]
    zeroNumCount = 0
    iter <- 0
    
    while (length(finalUnlabels) > 0) {
      iter <- iter + 1
      if (iter > iterator) {
        iter <- iter - 1
        break
      }
      
      permutRes <- cross_validation(seed = .randomSeed(), 
                                    method = method, featureMat = featureMat, positives = positives, 
                                    negatives = finalNegatives, cross = cross, cpus = cpus, 
                                    ...)
      maxAUC_Classifer <- .find_ClassifierWithMaxAUC(permutRes)
      prediction.score <- .predictor(method = method, classifier = maxAUC_Classifer$classifier, 
                                     featureMat = featureMat)
      positives.score <- sort(prediction.score[positives], decreasing = FALSE)
      negatives.score <- prediction.score[finalNegatives]
      unlabels.score <- prediction.score[finalUnlabels]
      positive.Threshold <- as.numeric(positives.score[thresholdIdx])
      
      num <- length(which(unlabels.score < positive.Threshold))
      
      # Update finalUnlabels and finalNegatives
      finalUnlabels <- names(unlabels.score[which(unlabels.score > positive.Threshold)])
      finalNegatives <- unique(c(names(unlabels.score[which(unlabels.score <= positive.Threshold)]), finalNegatives))
      
      # If no more unlabeled samples are left, break the loop
      if (length(finalUnlabels) == 0) {
        break
      }
      
      AUCMat <- .obtain_CV_AUCMat(permutRes)
      numMat[iter, 2] <- mean(AUCMat[, 1])
      numMat[iter, 3] <- mean(AUCMat[, 2])
      numMat[iter, 4] <- length(negatives.score)
      numMat[iter, 5] <- length(unlabels.score)
      
      if (trace == TRUE) {
        resultDir <- sprintf("%s/PSOL_Iteration_%s.RData",PSOLResDic, iter)
        
        iterRes <- list(permutRes = permutRes, method = method, 
                        classifier = maxAUC_Classifer, predictionScores = prediction.score, 
                        negativeScores = negatives.score, unlabelScores = unlabels.score, 
                        threshold = positive.Threshold, positives = positives, 
                        negatives = names(positives.score), unlabels = names(unlabels.score), 
                        finalNegatives = finalNegatives, finalUnlabels = finalUnlabels)
        save(iterRes, file = resultDir)
      }
      
      cat("\nPSOL_Iteration: ", iter, "\tAUC: ", numMat[iter, 3], "\tCurrentPosNum:", length(positives), "\tCurrentNegNum: ", numMat[iter, 4], "\tCurrentUnlabelNum: ", numMat[iter, 5], "\tIncreased negatives Num: ", num, "\n")
    }
    
    
    numMat <- numMat[1:iter, ]
    write.table(numMat, sprintf("%s/PSOL_NegativeIncreasement.txt",PSOLResDic), sep = "\t", quote = F)
    if (plot == TRUE) {
      pdf(sprintf("%s/PSOL_NegativeIncreasement.pdf",PSOLResDic), height = 10, width = 10)
      par(mar = c(5, 12, 4, 4) + 0.1)
      plot(numMat[, 1], numMat[, 3], axes = F, ylim = c(0, 
                                                        1), xlab = "", ylab = "", type = "l", col = "red", 
           main = "")
      points(numMat[, 1], numMat[, 3], pch = 20, col = "red", 
             cex = 0.8)
      axis(2, ylim = c(0, 1), col = "red", lwd = 2)
      mtext(2, text = "AUC", line = 2)
      par(new = T)
      plot(numMat[, 1], numMat[, 4], axes = F, ylim = c(0, 
                                                        max(numMat[, 4])), xlab = "", ylab = "", type = "l", 
           col = "black", lty = 2, main = "", lwd = 2)
      axis(2, ylim = c(0, max(numMat[, 4])), lwd = 2, line = 3.5, 
           col = "black")
      points(numMat[, 1], numMat[, 4], pch = 20, col = "black", 
             cex = 0.8)
      mtext(2, text = "Number of \"filtered-out\" genes ", 
            line = 5.5)
      axis(1, numMat[, 1])
      mtext("Iteration Number", side = 1, col = "black", line = 2)
      dev.off()
    }
  }
  
  
  PSOL <- function (featureMatrix = NULL, positives, balanced = FALSE, 
                    ratio = 10, unlabels, cpus = 20, PSOLResDic = NULL, iterator = iterator, 
                    TPR = 0.995, cross = 5, method = c("randomForest", "svm"), 
                    ...) 
  {
    if (is.null(featureMatrix)) {
      stop("Parameter featureMatrix is necessary for psol.", 
           "\n")
    }
    cat("Start using PSOL to select negative samples......", 
        "\n")
    if (is.null(PSOLResDic)) {
      PSOLResDic <- getwd()
      PSOLResDic <- paste0(PSOLResDic, "/")
    }
    intialNegatives <- res
    if (cpus >= cross) {
      cpus <- cross
    }
    negativeExpand <- .PSOL_NegativeExpansion(featureMat = featureMatrix, 
                                              positives = positives, negatives = intialNegatives$negatives, 
                                              unlabels = intialNegatives$unlabels, cpus = cpus, iterator = iterator, 
                                              cross = cross, PSOLResDic = PSOLResDic, method = method, 
    )
    load(sprintf("%s/PSOL_Iteration_%s.RData",PSOLResDic, iterator))
    finalNegatives <- iterRes$finalNegatives
    featureMat <- featureMatrix[c(positives, finalNegatives), 
    ]
    posLen <- length(positives)
    negLen <- length(finalNegatives)
    label <- c(rep(1, posLen), rep(0, negLen))
    fmat <- data.frame(featureMat)
    model <- randomForest(x = fmat, y = factor(label))
    iterRes[["model"]] <- model
    if (balanced) {
      finalNegatives <- sample(finalNegatives, size = length(positives))
    }
    iterRes[["finalNegatives"]] <- finalNegatives
    iterRes
  }
  
  # Starting running PSOL  
  psolRes <- PSOL(featureMatrix = featureMat, positives = positiveSamples, balanced = FALSE,
                  unlabels = unlabelSamples, PSOLResDic = psolResDic, cpus = cpus, iterator = iterator)
  
  return(psolRes)
}

#Function: get sample index for cv cross validation
.cvSampleIndex <- function( len, cross = 5, seed = 1 ) {
  
  cv <- cross
  sample_matrix <- matrix(0, nrow = len, ncol = cv)
  colnames(sample_matrix) <- paste("cv", c(1:cv), sep = "" )
  
  #random samples 
  set.seed(seed)
  index <- sample(1:len, len, replace = FALSE )
  step = floor( len/cv )
  
  start <- NULL
  end <- NULL
  train_lens <- rep(0, cv)
  for( i in c(1:cv) ) {
    start <- step*(i-1) + 1
    end <- start + step - 1
    if( i == cv ) 
      end <- len
    
    train <- index[-c(start:end)]
    test <- index[start:end]
    train_lens[i] <- length(train)
    
    sample_matrix[,i] <- c(train, test)
  }#end for i
  
  return( list( train_lens = train_lens, sample_matrix = sample_matrix))
}

#' @export
classifier <- function( method = c("randomForest", "svm"), featureMat, positiveSamples, 
                        negativeSamples, ...) {
  call <- match.call()
  
  if( length(method) > 1){
    method <- method[1]
  } 
  
  
  if( is.null(rownames(featureMat) ) )
    stop("Error: no row names (i.e., sample IDs) were assigned for featureMat." )
  if( is.null(colnames(featureMat) ) )
    stop("Error: no colnames were defined for featureMat." )
  
  positiveSamples <- intersect( rownames(featureMat), positiveSamples )
  negativeSamples <- intersect( rownames(featureMat), negativeSamples )
  posLen <- length(positiveSamples)
  negLen <- length(negativeSamples)
  if( posLen == 0 )
    stop("Error: no positive samples included in featureMat." )
  if( negLen == 0 )
    stop("Error: no negative samples were included in featureMat." )
  
  label <- c( rep(1, posLen), rep(0, negLen) )
  fmat <- data.frame( featureMat[c(positiveSamples, negativeSamples), ] )
  tmpData <- cbind( fmat, label )
  colnames(tmpData) <- c(colnames(fmat), "Class")
  if( method == "randomForest" ) {
    obj <- randomForest(x = fmat, y = factor(label), ... )
  }else{
    obj <- svm(x = fmat, y = factor(label), ... )
  }
  obj
}

.predictor <- function( method = c("randomForest", "svm"), classifier, featureMat ) {
  
  if(length(method) > 1){
    method <- method[1]
  }
  
  if( method == "randomForest") {
    res <- predict(classifier, data.frame(featureMat), type= "vote" )[,"1"]
  }else {
    res <- predict( classifier, data.frame(featureMat), type = "raw") 
  }
  names(res) <- rownames(featureMat)
  res
}

.find_ClassifierWithMaxAUC <- function( cvRes ) {
  
  classifier <- NA
  maxAUC <- 0
  for( i in 1:length(cvRes) ) {
    res <- cvRes[[i]]
    if( res$test.AUC > maxAUC) {
      maxAUC <- res$test.AUC
      classifier <- res$classifier
    }
  }#end for i
  
  return( list(maxAUC = maxAUC, classifier = classifier))
}

.obtain_CV_AUCMat <- function( cvRes ) {
  cv <- length(cvRes)
  AUCMat <- matrix(0, nrow = cv, ncol = 2 )
  rownames(AUCMat) <- paste( "cv", 1:cv, sep = "" )
  colnames(AUCMat) <- c("trainingSet", "testingSet")
  
  for( i in 1:cv ) {
    res <- cvRes[[i]]
    AUCMat[i,2] <- res$test.AUC
    AUCMat[i,1] <- res$train.AUC
  }#end for i
  
  AUCMat
}

##get system time for seed and then generate random index
.randomSeed <- function() {
  curtime <- format(Sys.time(), "%H:%M:%OS4")
  XXX <- unlist(strsplit(curtime, ":"))
  curtimeidx <- (as.numeric(XXX[1])*3600 + as.numeric(XXX[2])*60 + as.numeric(XXX[3]))*10000
  curtimeidx
}

.one_cross_validation <- function( cv, method, featureMat, positives, negatives, posSample_cv, negSample_cv, balanced = TRUE, ratio = 10, ... ) {
  call <- match.call()
  j <- cv
  
  #for train samples
  train_genes_p <- positives[ (posSample_cv$sample_matrix[,j][1:posSample_cv$train_lens[j]] ) ]
  test_genes_p <- positives[ (posSample_cv$sample_matrix[,j][-c(1:posSample_cv$train_lens[j])]) ]
  
  #trained negatives randomly selected, and tested on all negatives
  train_genes_n <- negatives[(negSample_cv$sample_matrix[,j][1:negSample_cv$train_lens[j]] ) ]
  test_genes_n <- negatives[ (negSample_cv$sample_matrix[,j][-c(1:negSample_cv$train_lens[j])]) ]
  
  #select part of train_genes_n
  if( balanced == TRUE ) {
    if( length(train_genes_n) > ratio*length(train_genes_p) ) {
      train_genes_n <- train_genes_n[sample(1:length(train_genes_n), replace=FALSE)[1:(ratio*length(train_genes_p))]]
    }
  }
  
  
  
  obj <- classifier( method = method, featureMat = featureMat, positiveSamples = train_genes_p, negativeSamples = train_genes_n, ... )
  bestmodel <- obj
  
  positives.train.score <- .predictor( method = method, classifier = bestmodel, featureMat = featureMat[train_genes_p,])
  negatives.train.score <- .predictor( method = method, classifier = bestmodel, featureMat = featureMat[train_genes_n,])
  positives.test.score <- .predictor( method = method, classifier = bestmodel, featureMat = featureMat[test_genes_p,])
  negatives.test.score <- .predictor( method = method, classifier = bestmodel, featureMat = featureMat[test_genes_n,])
  
  
  
  train.AUC <- roc( c(rep(1, length(train_genes_p)), rep(0, length(train_genes_n))), 
                    c(positives.train.score, negatives.train.score) )$auc[1]
  test.AUC <- roc( c(rep(1, length(test_genes_p)), rep(0, length(test_genes_n))), 
                   c(positives.test.score, negatives.test.score) )$auc[1]
  
  res <- ( list( positves.train = train_genes_p, negatives.train = train_genes_n, 
                 positives.test = test_genes_p, negatives.test = test_genes_n, 
                 ml = method, classifier = bestmodel, 
                 positives.train.score = positives.train.score,
                 negatives.train.score = negatives.train.score,
                 positives.test.score = positives.test.score,
                 negatives.test.score = negatives.test.score,
                 train.AUC = train.AUC,
                 test.AUC = test.AUC) )
  
  res
}

#' @export
cross_validation <- function( seed = 1, method = c("randomForest", "svm"), 
                              featureMat, positives, negatives, cross = 5, 
                              cpus = 1, ... ){
  
  call <- match.call()
  
  if( length(method) > 1){
    method <- method[1]
  } 
  
  #sample index for cv
  posSample_cv <- .cvSampleIndex(length(positives), cross = cross, seed = seed)
  negSample_cv <- .cvSampleIndex(length(negatives), cross = cross, seed = seed)
  
  cvRes <- list()
  if( cpus > 1 ) {
    #require(snowfall)
    sfInit(parallel = TRUE, cpus = cpus)
    sfExport("classifier", namespace = "PEA")
    sfExport(".predictor", namespace = "PEA")
    sfExport(".one_cross_validation", namespace = "PEA")
    sfLibrary("pROC", character.only = TRUE)
    sfLibrary("e1071", character.only = TRUE)
    sfLibrary("randomForest", character.only = TRUE)
    
    cvRes <- sfApply( matrix(1:cross, ncol = 1), 1,  .one_cross_validation, method = method, featureMat = featureMat, positives = positives, negatives = negatives, posSample_cv = posSample_cv, negSample_cv = negSample_cv, ...)
    sfStop()
  }else {
    for( j in 1:cross ) {
      cvRes[[j]] <- .one_cross_validation( cv = j, method = method, featureMat = featureMat, positives = positives, negatives = negatives, posSample_cv = posSample_cv, negSample_cv = negSample_cv, ... )
    }
  }
  cvRes
}

classifier <- function( method = c("randomForest", "svm"), featureMat, positiveSamples, 
                        negativeSamples, ...) {
  call <- match.call()
  
  if( length(method) > 1){
    method <- method[1]
  } 
  
  
  if( is.null(rownames(featureMat) ) )
    stop("Error: no row names (i.e., sample IDs) were assigned for featureMat." )
  if( is.null(colnames(featureMat) ) )
    stop("Error: no colnames were defined for featureMat." )
  
  positiveSamples <- intersect( rownames(featureMat), positiveSamples )
  negativeSamples <- intersect( rownames(featureMat), negativeSamples )
  posLen <- length(positiveSamples)
  negLen <- length(negativeSamples)
  if( posLen == 0 )
    stop("Error: no positive samples included in featureMat." )
  if( negLen == 0 )
    stop("Error: no negative samples were included in featureMat." )
  
  label <- c( rep(1, posLen), rep(0, negLen) )
  fmat <- data.frame( featureMat[c(positiveSamples, negativeSamples), ] )
  tmpData <- cbind( fmat, label )
  colnames(tmpData) <- c(colnames(fmat), "Class")
  if( method == "randomForest" ) {
    obj <- randomForest(x = fmat, y = factor(label), ... )
  }else{
    obj <- svm(x = fmat, y = factor(label), ... )
  }
  obj
}

plotROC <- function(cvRes) {
  
  
  cvListPredictions <- list()
  cvListLabels <- list()
  AUCVec <- rep(0, length(cvRes) )
  for( i in 1:length(cvRes) ) {
    curCV <- cvRes[[i]]
    cvListPredictions[[i]] <- c( curCV$positives.test.score, curCV$negatives.test.score )
    cvListLabels[[i]] <- c( rep(1, length(curCV$positives.test.score)), rep(0, length(curCV$negatives.test.score) ) )
    AUCVec[i] <- curCV$test.AUC
  }
  mAUC <- format( mean(AUCVec), digits= 3)
  
  pred <- prediction(cvListPredictions, cvListLabels)
  perf <- performance(pred,"tpr","fpr")
  
  
  par(mar=c(5,6,4,2))   
  plot(perf, col= "gray", lty=3, main = paste( "AUC = ", mAUC, sep = ""), cex.lab = 2.5, cex.axis = 2, cex.main = 3, mgp = c(4,1.8,0) )
  plot(perf, col = "black",  lwd= 3, avg="vertical", spread.estimate="none", add=TRUE)  
  
}

prediction <- function(predictions, labels, label.ordering=NULL) {
  
  ## bring 'predictions' and 'labels' into list format,
  ## each list entry representing one x-validation run
  
  ## convert predictions into canonical list format
  if (is.data.frame(predictions)) {
    names(predictions) <- c()
    predictions <- as.list(predictions)
  } else if (is.matrix(predictions)) {
    predictions <- as.list(data.frame(predictions))
    names(predictions) <- c()
  } else if (is.vector(predictions) && !is.list(predictions)) {
    predictions <- list(predictions)
  } else if (!is.list(predictions)) {
    stop("Format of predictions is invalid. It couldn't be coerced to a list.",
         call. = FALSE)
  }
  ## if predictions is a list -> keep unaltered
  if(any(vapply(predictions,anyNA,logical(1)))){
    stop("'predictions' contains NA.", call. = FALSE)
  }
  
  ## convert labels into canonical list format
  if (is.data.frame(labels)) {
    names(labels) <- c()
    labels <- as.list( labels)
  } else if (is.matrix(labels)) {
    labels <- as.list( data.frame( labels))
    names(labels) <- c()
  } else if ((is.vector(labels) ||
              is.ordered(labels) ||
              is.factor(labels)) &&
             !is.list(labels)) {
    labels <- list( labels)
  } else if (!is.list(labels)) {
    stop("Format of labels is invalid. It couldn't be coerced to a list.",
         call. = FALSE)
  }
  ## if labels is a list -> keep unaltered
  
  ## Length consistency checks
  if (length(predictions) != length(labels))
    stop(paste("Number of cross-validation runs must be equal",
               "for predictions and labels."))
  if (! all(sapply(predictions, length) == sapply(labels, length)))
    stop(paste("Number of predictions in each run must be equal",
               "to the number of labels for each run."))
  
  ## only keep prediction/label pairs that are finite numbers
  for (i in 1:length(predictions)) {
    finite.bool <- is.finite( predictions[[i]] )
    predictions[[i]] <- predictions[[i]][ finite.bool ]
    labels[[i]] <- labels[[i]][ finite.bool ]
  }
  
  ## abort if 'labels' format is inconsistent across
  ## different cross-validation runs
  label.format=""  ## one of 'normal','factor','ordered'
  if (all(sapply( labels, is.factor)) &&
      !any(sapply(labels, is.ordered))) {
    label.format <- "factor"
  } else if (all(sapply( labels, is.ordered))) {
    label.format <- "ordered"
  } else if (all(sapply( labels, is.character)) ||
             all(sapply( labels, is.numeric)) ||
             all(sapply( labels, is.logical))) {
    label.format <- "normal"
  } else {
    stop(paste("Inconsistent label data type across different",
               "cross-validation runs."))
  }
  
  ## abort if levels are not consistent across different
  ## cross-validation runs
  if (! all(sapply(labels, levels)==levels(labels[[1]])) ) {
    stop(paste("Inconsistent factor levels across different",
               "cross-validation runs."))
  }
  
  ## convert 'labels' into ordered factors, aborting if the number
  ## of classes is not equal to 2.
  levels <- c()
  if ( label.format == "ordered" ) {
    if (!is.null(label.ordering)) {
      stop(paste("'labels' is already ordered. No additional",
                 "'label.ordering' must be supplied."))
    } else {
      levels <- levels(labels[[1]])
    }
  } else {
    if ( is.null( label.ordering )) {
      if ( label.format == "factor" ) levels <- sort(levels(labels[[1]]))
      else levels <- sort( unique( unlist( labels)))
    } else {
      ## if (!setequal( levels, label.ordering)) {
      if (!setequal( unique(unlist(labels)), label.ordering )) {
        stop("Label ordering does not match class labels.")
      }
      levels <- label.ordering
    }
    for (i in 1:length(labels)) {
      if (is.factor(labels))
        labels[[i]] <- ordered(as.character(labels[[i]]),
                               levels=levels)
      else labels[[i]] <- ordered( labels[[i]], levels=levels)
    }
    
  }
  
  if (length(levels) != 2) {
    message <- paste("Number of classes is not equal to 2.\n",
                     "ROCR currently supports only evaluation of ",
                     "binary classification tasks.",sep="")
    stop(message)
  }
  
  ## determine whether predictions are continuous or categorical
  ## (in the latter case stop; scheduled for the next ROCR version)
  if (!is.numeric( unlist( predictions ))) {
    stop("Currently, only continuous predictions are supported by ROCR.")
  }
  
  ## compute cutoff/fp/tp data
  
  cutoffs <- list()
  fp <- list()
  tp <- list()
  fn <- list()
  tn <- list()
  n.pos <- list()
  n.neg <- list()
  n.pos.pred <- list()
  n.neg.pred <- list()
  for (i in 1:length(predictions)) {
    n.pos <- c( n.pos, sum( labels[[i]] == levels[2] ))
    n.neg <- c( n.neg, sum( labels[[i]] == levels[1] ))
    ans <- .compute.unnormalized.roc.curve( predictions[[i]], labels[[i]] )
    cutoffs <- c( cutoffs, list( ans$cutoffs ))
    fp <- c( fp, list( ans$fp ))
    tp <- c( tp, list( ans$tp ))
    fn <- c( fn, list( n.pos[[i]] - tp[[i]] ))
    tn <- c( tn, list( n.neg[[i]] - fp[[i]] ))
    n.pos.pred <- c(n.pos.pred, list(tp[[i]] + fp[[i]]) )
    n.neg.pred <- c(n.neg.pred, list(tn[[i]] + fn[[i]]) )
  }
  
  
  return( new("prediction", predictions=predictions,
              labels=labels,
              cutoffs=cutoffs,
              fp=fp,
              tp=tp,
              fn=fn,
              tn=tn,
              n.pos=n.pos,
              n.neg=n.neg,
              n.pos.pred=n.pos.pred,
              n.neg.pred=n.neg.pred))
}


## fast fp/tp computation based on cumulative summing
.compute.unnormalized.roc.curve <- function( predictions, labels ) {
  ## determine the labels that are used for the pos. resp. neg. class :
  pos.label <- levels(labels)[2]
  neg.label <- levels(labels)[1]
  
  pred.order <- order(predictions, decreasing=TRUE)
  predictions.sorted <- predictions[pred.order]
  tp <- cumsum(labels[pred.order]==pos.label)
  fp <- cumsum(labels[pred.order]==neg.label)
  
  ## remove fp & tp for duplicated predictions
  ## as duplicated keeps the first occurrence, but we want the last, two
  ## rev are used.
  ## Highest cutoff (Infinity) corresponds to tp=0, fp=0
  dups <- rev(duplicated(rev(predictions.sorted)))
  tp <- c(0, tp[!dups])
  fp <- c(0, fp[!dups])
  cutoffs <- c(Inf, predictions.sorted[!dups])
  
  return(list( cutoffs=cutoffs, fp=fp, tp=tp ))
}