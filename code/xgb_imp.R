library(xgboost)
#library(caret)
#library(tidyverse)
library(doParallel)
library(doRNG)
library(itertools)
numCores <- parallel::detectCores()
#cl <- makeCluster(5)
registerDoParallel(numCores - 1)

#source('./SimulationDesign1/constants.R')
#source('./SimulationDesign1/simulateData.R')
#Y_sim <- simulateData()$Y
#Y_mis <- generate_missing(Y_sim)

xgboost_imputation <- function(data, maxiter = 10, verbose = TRUE, max.depth = 2, nrounds = 50, decreasing = FALSE, parallelize = 'variables') {
  xmis <- data
  #xmis <- Y_mis
  n <- nrow(xmis)
  p <- ncol(xmis)
  
  # Remove complete missing columns
  if (any(apply(is.na(xmis), 2, sum) == n)){
    indCmis <- which(apply(is.na(xmis), 2, sum) == n)
    xmis <- xmis[,-indCmis]
    p <- ncol(xmis)
    cat('  removed variable(s)', indCmis,
        'due to the missingness of all entries\n')
  } 
  
  # Initial mean imputation
  ximp <- xmis
  varType <- character(p)
  for (t.co in 1:p) {
    ximp[is.na(xmis[,t.co]),t.co] <- colMeans(xmis, na.rm = TRUE)[t.co] #mean(xmis[,1], na.rm = TRUE)
    varType[t.co] <- 'numeric'
    #if (is.numeric(xmis[[t.co]])) {
    # varType[t.co] <- 'numeric'
    #ximp[is.na(xmis[,t.co]),t.co] <- mean(xmis[,t.co], na.rm = TRUE)
    #next()
    #}
  }
  
  # Extract missing location
  NAloc <- is.na(xmis)            # where are missings
  noNAvar <- apply(NAloc, 2, sum) # how many are missing in the vars
  sort.j <- order(noNAvar)        # indices of increasing amount of NA in vars
  #decreasing = FALSE
  if (decreasing)
    sort.j <- rev(sort.j)
  sort.noNAvar <- noNAvar[sort.j]
  
  # Extract the columns to be imputed 
  nzsort.j <- sort.j[sort.noNAvar > 0]
  #parallelize <- 'variables'
  if (parallelize == 'variables') {
    '%cols%' <- get('%dorng%')
    idxList <- as.list(isplitVector(nzsort.j, chunkSize = getDoParWorkers()))
  } 
  
  ## initialize parameters of interest
  iter <- 0
  k <- length(unique(varType))
  convNew <- rep(0, k)
  names(convNew) <- c('numeric')
  convOld <- rep(Inf, k)
  convergence <- c()
  #OOBerror <- numeric(p)
  #names(OOBerror) <- varType
  
  #maxiter <- 10
  Ximp <- vector('list', maxiter)
  
  # Stop criteria
  stopCriterion <- function(varType, convNew, convOld, iter, maxiter){
    k <- length(unique(varType))
    if (k == 1){
      (convNew < convOld) & (iter < maxiter)
    } else {
      ((convNew[1] < convOld[1]) | (convNew[2] < convOld[2])) & (iter < maxiter)
    }
  }
  
  # Main Imputation Function
  while (stopCriterion(varType, convNew, convOld, iter, maxiter)) {
    if (iter != 0) {
      convOld <- convNew
      # OOBerrOld <- OOBerr
    }
    if (verbose) {
      cat("  XGBoost iteration", iter+1, "in progress...")
    }
    t.start <- proc.time()
    ximp.old <- ximp
    
    for (idx in idxList) {
      results <- foreach(varInd = idx, .packages = 'xgboost') %cols% {
        obsi <- !NAloc[, varInd] # which i's are observed
        misi <- NAloc[, varInd] # which i's are missing
        obsY <- ximp[obsi, varInd] # training response
        misY <- ximp[misi, varInd] # testing response
        obsX <- ximp[obsi, seq(1, p)[-varInd]] # training variables
        misX <- ximp[misi, seq(1, p)[-varInd]] # prediction variables
        xgb_obsX = xgb.DMatrix(data = as.matrix(obsX), label = as.matrix(obsY)) # xgboost version of train
        xgb_misX = xgb.DMatrix(data = as.matrix(misX), label = as.matrix(misY))
        #typeY <- varType[varInd]
        
        # xgboost on train
        xgb <- xgboost(data = xgb_obsX, max.depth = max.depth, nrounds = nrounds)
        
        misY <- predict(xgb, xgb_misX) ## predict missing values in column varInd
        
        list(varInd = varInd, misY = misY) #oerr = oerr
      }
      
      for (res in results) {
        misi <- NAloc[,res$varInd]
        ximp[misi, res$varInd] <- res$misY
        #OOBerror[res$varInd] <- res$oerr
      }
    }
    
    if (verbose){
      cat('done!\n')
    }
    
    iter <- iter + 1
    Ximp[[iter]] <- ximp
    
    t.co2 <- 1
    ## check the difference between iteration steps
    for (t.type in names(convNew)){
      t.ind <- which(varType == t.type)
      convNew[t.co2] <- sum((ximp[, t.ind] - ximp.old[, t.ind])^2) / sum(ximp[, t.ind]^2)
      t.co2 <- t.co2 + 1
    }
  }  
  
  # Extrtact the output
  if (iter == maxiter) {
    out <- Ximp[[iter]]
  } else {
    out <- Ximp[[iter - 1]]
  }
  
  return(out)
}

