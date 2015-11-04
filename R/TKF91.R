TKF91LikelihoodFunctionWrapperR <- function(x, seq1Int, seq2Int, 
                                            mu, expectedLength,
                                            substModel, substModelBF){
  ansTemp <- .Call("TKF91LikelihoodFunctionWrapper",
                   seq1Int, seq2Int, x, mu,
                   expectedLength, substModel,
                   substModelBF)
  return(ansTemp["negLogLikelihood"])
}

smartOptimBrent <- function(fn, par, lower, upper, ...){
  message("The range is ", lower, " to ", upper)
  res <- optim(par, TKF91LikelihoodFunctionWrapperR,
               gr=NULL, ..., 
               method="Brent",
               lower=lower, upper=upper, hessian=TRUE)
  if(isTRUE(all.equal(res$par[1], lower, tolerance=1e-4))){
    lower <- lower * 1.1
    res <- smartOptimBrent(fn, par, lower=lower, upper=upper, ...)
  }
  if(isTRUE(all.equal(res$par[1], upper, tolerance=1e-4))){
    upper <- upper * 0.9
    res <- smartOptimBrent(fn, par, lower=lower, upper=upper, ...)
  }
  return(res)
}


TKF91Pair <- function(seq1, seq2, mu=NULL, distance=NULL,
                      ## mu: by default is 0.001 from median of mu values 
                      ## from Fungi dataset.
                      method=c("gsl", "nlopt", "NM", "Sbplx", "COBYLA", 
                               "BOBYQA", "PRAXIS"),
                      expectedLength=362, 
                      substModel, substModelBF){
  if(!all(seq1 %in% AACharacterSet) || !all(seq2 %in% AACharacterSet)){
    stop("This implementation currently only supports 20 AA characters ", 
         paste(AACharacterSet, collapse=" "))
  }
  method <- match.arg(method)
  methodsOpt <- c("NM", "Sbplx", "COBYLA", "BOBYQA", "PRAXIS")
  seq1Int <- AAToInt(seq1)
  seq2Int <- AAToInt(seq2)
  ## for the C matrix index
  seq1Int <- seq1Int - 1L
  seq2Int <- seq2Int - 1L

  expectedLength <- as.numeric(expectedLength)
  
  if(is.null(mu) && is.null(distance)){ 
    ## Do the 2D optimisation
    if(method == "nlopt"){
      ## We try all the optimisation methods and select the best one
      ans_all <- lapply(methodsOpt, 
                        function(x){.Call("TKF91LikelihoodFunction2DMain_nlopt",
                                          seq1Int, seq2Int, expectedLength, 
                                          substModel, substModelBF, x)}
                        )
      ans <- ans_all[[which.min(sapply(ans_all, "[", "negLogLikelihood"))]]
    }else if(method == "gsl"){
      ans <- .Call("TKF91LikelihoodFunction2DMainNM", seq1Int, seq2Int,
                   expectedLength, substModel, substModelBF)
    }else{
      ans <- .Call("TKF91LikelihoodFunction2DMain_nlopt", seq1Int, seq2Int,
                   expectedLength, substModel, substModelBF, method)
    }
    ansHessian <- hessian(function(x, seq1Int, seq2Int, expectedLength, 
                                   substModel, substModelBF){
                          ansTemp <- .Call("TKF91LikelihoodFunctionWrapper", 
                                           seq1Int, seq2Int, x[1], x[2], 
                                           expectedLength, substModel, 
                                           substModelBF)
                          return(ansTemp["negLogLikelihood"])
                 }, c(ans["PAM"], ans["Mu"]),
                 seq1Int=seq1Int, seq2Int=seq2Int,
                 expectedLength=expectedLength, substModel=substModel,
                 substModelBF=substModelBF)
    if(any(is.nan(ansHessian))){
      message("Hessian matrix calculation failed on current optimal points! Use the same values as variance.")
      invHessian <- matrix(c(ans["PAM"], NaN, NaN, ans["Mu"]), ncol=2)
    }else{
      invHessian <- solve(ansHessian)
    }
    return(c(ans, "PAMVariance"=invHessian[1,1],
             "MuVariance"=invHessian[2,2],
             "coVariance"=invHessian[1,2]))
  }else if(!is.null(mu) && is.null(distance)){
    ## Do the 1D distance optimisation
    #if(method1D == "brent"){
    #  message("Using method: ", method1D)
    #ans <- .Call("TKF91LikelihoodFunction1DMain", seq1Int, seq2Int, mu,
    #             expectedLength, substModel, substModelBF)
    #ansHessian <- hessian(TKF91LikelihoodFunctionWrapperR,
    #             ans["PAM"], 
    #             seq1Int=seq1Int, seq2Int=seq2Int,
    #             mu=mu, expectedLength=expectedLength,
    #             substModel=substModel, 
    #             substModelBF=substModelBF)
    #}else if(method1D == "optimise"){
    #  message("Using method: ", method1D)
    #  distanceMin <- optimise(TKF91LikelihoodFunctionWrapperR,
    #             interval=c(0.0494497, 1000), 
    #             seq1Int=seq1Int, seq2Int=seq2Int,
    #             mu=mu, expectedLength=expectedLength,
    #             substModel=substModel,
    #             substModelBF=substModelBF)
    #  ansHessian <- hessian(TKF91LikelihoodFunctionWrapperR,
    #                        distanceMin$minimum,
    #                        seq1Int=seq1Int, seq2Int=seq2Int,
    #                        mu=mu, expectedLength=expectedLength,
    #                        substModel=substModel,
    #                        substModelBF=substModelBF)
    #  ans <- c("PAM"=distanceMin$minimum, "Mu"=mu, 
    #           "negLogLikelihood"=distanceMin$objective)
    #}else if(method1D == "optim"){
    #  message("Using method: ", method1D)
      res <- smartOptimBrent(TKF91LikelihoodFunctionWrapperR, par=100,
                             lower=0.0494497, upper=2000,
                   seq1Int=seq1Int, seq2Int=seq2Int,
                   mu=mu, expectedLength=expectedLength,
                   substModel=substModel, substModelBF=substModelBF)
      ansHessian <- res$hessian
      ans <- c("PAM"=res$par[1], "Mu"=mu, "negLogLikelihood"=res$value)
    invHessian <- solve(ansHessian)
    return(c(ans, "PAMVariance"=invHessian[1,1]))
  }else if(!is.null(mu) && !is.null(distance)){
    ## Just calculate the likelihood, given mu and distance
    ans <- .Call("TKF91LikelihoodFunctionWrapper", seq1Int, seq2Int, 
                 distance, mu, expectedLength, substModel, substModelBF)
    ansHessian <- hessian(function(x, seq1Int, seq2Int, expectedLength, 
                                   substModel, substModelBF){
                          ansTemp <- .Call("TKF91LikelihoodFunctionWrapper", 
                                           seq1Int, seq2Int, x[1], x[2], 
                                           expectedLength, substModel, 
                                           substModelBF)
                          return(ansTemp["negLogLikelihood"])
                 }, c(ans["PAM"], ans["Mu"]),
                 seq1Int=seq1Int, seq2Int=seq2Int,
                 expectedLength=expectedLength, substModel=substModel,
                 substModelBF=substModelBF)
    #invHessian <- chol2inv(chol(ansHessian))
    invHessian <- solve(ansHessian)
    return(c(ans, "PAMVariance"=invHessian[1,1],
             "MuVariance"=invHessian[2,2],
             "coVariance"=invHessian[1,2]))
  }else{
    stop("You cannot estimate mu alone!")
  }
}


TKF91 <- function(fasta, mu=NULL, 
                  method=c("gsl", "nlopt", "NM", "Sbplx", "COBYLA", 
                           "BOBYQA", "PRAXIS"),
                  expectedLength=362, 
                  substModel, substModelBF,
                  skipFailure=FALSE){
  method <- match.arg(method)
  seqnames <- names(fasta)
  nSeqs <- length(fasta)
  distanceMatrix <- matrix(NA, ncol=nSeqs, nrow=nSeqs,
                           dimnames=list(seqnames, seqnames))
  diag(distanceMatrix) <- 0
  varianceMatrix <- distanceMatrix
  if(is.null(mu)){
    muMatrix <- distanceMatrix
  }
  negLoglikelihoodMatrix <- distanceMatrix  
  for(i in 1:(nSeqs-1L)){
    for(j in (i+1L):nSeqs){
      message(seqnames[i], " vs ", seqnames[j])
      if(skipFailure){
        ans <- try(TKF91Pair(fasta[[i]], fasta[[j]],
                             mu=mu, method=method,
                             expectedLength=expectedLength,
                             substModel=substModel, substModelBF=substModelBF))
        if(class(ans) == "try-error")
          next
      }else{
      ans <- TKF91Pair(fasta[[i]], fasta[[j]], 
                       mu=mu, method=method,
                       expectedLength=expectedLength,
                       substModel=substModel, substModelBF=substModelBF)
      }
      distanceMatrix[i,j] <- distanceMatrix[j,i] <- ans["PAM"]
      varianceMatrix[i,j] <- varianceMatrix[j,i] <- ans["PAMVariance"]
      negLoglikelihoodMatrix[i,j] <- negLoglikelihoodMatrix[j,i] <- 
        ans["negLogLikelihood"]
      if(is.null(mu)){
        muMatrix[i,j] <- muMatrix[j,i] <- ans["Mu"]
      }
    }
  }
  if(is.null(mu)){
    return(list(distanceMatrix=distanceMatrix,
                varianceMatrix=varianceMatrix,
                muMatrix=muMatrix,
                negLoglikelihoodMatrix=negLoglikelihoodMatrix
                ))
  }else{
    return(list(distanceMatrix=distanceMatrix,
                varianceMatrix=varianceMatrix,
                negLoglikelihoodMatrix=negLoglikelihoodMatrix
                ))
  }
}



