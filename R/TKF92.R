TKF92LikelihoodFunctionWrapperR <- function(x, seq1Int, seq2Int, 
                                            mu, r, expectedLength,
                                            substModel, substModelBF){
  ansTemp <- .Call("TKF92LikelihoodFunctionWrapper",
                   seq1Int, seq2Int, x, mu, r, 
                   expectedLength, substModel,
                   substModelBF)
  return(ansTemp["negLogLikelihood"])
}

TKF923DLikelihoodFunctionWrapperR <- function(x, seq1Int, seq2Int, 
                                              expectedLength,
                                              substModel, substModelBF){
  ansTemp <- .Call("TKF92LikelihoodFunctionWrapper",
                   seq1Int, seq2Int, x[1], x[2], x[3],
                   expectedLength, substModel,
                   substModelBF)
  return(ansTemp["negLogLikelihood"])
}

smartOptimBrentTKF92 <- function(fn, par, lower, upper, ...){
  message("The range is ", lower, " to ", upper)
  res <- optim(par, TKF92LikelihoodFunctionWrapperR,
               gr=NULL, ..., 
               method="Brent",
               lower=lower, upper=upper, hessian=TRUE)
  if(isTRUE(all.equal(res$par[1], lower, tolerance=1e-4))){
    lower <- lower * 1.1
    res <- smartOptimBrentTKF92(fn, par, lower=lower, upper=upper, ...)
  }
  if(isTRUE(all.equal(res$par[1], upper, tolerance=1e-4))){
    upper <- upper * 0.9
    res <- smartOptimBrentTKF92(fn, par, lower=lower, upper=upper, ...)
  }
  return(res)
}


TKF92Pair <- function(seq1, seq2, mu=NULL, r=NULL, distance=NULL,
                      ## mu: by default is 0.001 from median of mu values 
                      ## of Fungi dataset.
                      ## r: the probability in the geometric distribution 
                      ## of fragment length. By default, it should be 0.8480
                      ## from median of r values of Fungi dataset.
                      method=c("NM", "constrOptim"),
                      expectedLength=362, 
                      substModel, substModelBF){
  if(!all(seq1 %in% AACharacterSet) || !all(seq2 %in% AACharacterSet)){
    stop("This implementation currently only supports 20 AA characters ",
         paste(AACharacterSet, collapse=" "))
  }
  method <- match.arg(method)
  #methodsOpt <- c("NM", "Sbplx", "COBYLA", "BOBYQA", "PRAXIS")

  seq1Int <- AAToInt(seq1)
  seq2Int <- AAToInt(seq2)
  ## for the C matrix index
  seq1Int <- seq1Int - 1L
  seq2Int <- seq2Int - 1L

  expectedLength <- as.numeric(expectedLength)
  
  if(is.null(mu) && is.null(distance) && is.null(r)){ 
    ## Do the 3D optimisation
#     if(method == "nlopt"){
#       ## We try all the optimisation methods and select the best one
#       ans_all <- lapply(methodsOpt,
#                         function(x){.Call("TKF92LikelihoodFunction3DMain_nlopt",
#                                           seq1Int, seq2Int, 
#                                           expectedLength, 
#                                           substModel, substModelBF,
#                                           x)}
#                         )
#       ans <- ans_all[[which.min(sapply(ans_all, "[", "negLogLikelihood"))]]
#     }else if(method == "gsl"){
    if(method == "NM"){
      ans <- .Call("TKF92LikelihoodFunction3DMainNM", seq1Int, seq2Int,
                   expectedLength, substModel, substModelBF)
    }else{
      #ans <- .Call("TKF92LikelihoodFunction3DMain_nlopt", seq1Int, seq2Int,
      #             expectedLength, substModel, substModelBF, method)
      res <- constrOptim(theta=c(100, exp(-3), 0.5),
                         f=TKF923DLikelihoodFunctionWrapperR,
                         grad=NULL,
                         ui=matrix(c(1,0,-1,0,0, 0, 
                                     0,1,0,-1, 0, 0,
                                     0,0,0,0,1,-1), ncol=3),
                         ci=c(0.0494497, 1e-20, -2000, -1+1e-20, 
                              1e-20, -1+1e-20),
                         seq1Int=seq1Int, seq2Int=seq2Int,
                         expectedLength=expectedLength,
                         substModel=substModel, substModelBF=substModelBF
      )
      ans <- c("PAM"=res$par[1], "Mu"=res$par[2], "r"=res$par[3], 
               "negLogLikelihood"=unname(res$value))
    }
    ansHessian <- hessian(function(x, seq1Int, seq2Int, expectedLength, 
                                   substModel, substModelBF){
                          ansTemp <- .Call("TKF92LikelihoodFunctionWrapper", 
                                           seq1Int, seq2Int, x[1], x[2], x[3], 
                                           expectedLength, substModel, 
                                           substModelBF)
                          return(ansTemp["negLogLikelihood"])
                 }, c(ans["PAM"], ans["Mu"], ans["r"]),
                 seq1Int=seq1Int, seq2Int=seq2Int,
                 expectedLength=expectedLength, substModel=substModel,
                 substModelBF=substModelBF)
    if(any(is.nan(ansHessian))){
      message("Hessian matrix calculation failed on current optimal points! 
              Use the same values as variance.")
      invHessian <- matrix(c(ans["PAM"], NaN, NaN, 
                             NaN, ans["Mu"], NaN,
                             NaN, NaN, ans["r"]), ncol=3)
    }else{
      invHessian <- solve(ansHessian)
    }
    return(c(ans, "PAMVariance"=invHessian[1,1],
             "MuVariance"=invHessian[2,2],
             "rVariance"=invHessian[3,3],
             "coVariancePAMMu"=invHessian[1,2],
             "coVariancePAMr"=invHessian[1,3],
             "coVarianceMur"=invHessian[2,3]
             )
           )
  }else if(!is.null(mu) && is.null(distance) && !is.null(r)){
    ## Do the 1D distance optimisation
    message("Using method: brent")
    res <- smartOptimBrentTKF92(TKF92LikelihoodFunctionWrapperR, par=100,
                           lower=0.0494497, upper=2000,
                           seq1Int=seq1Int, seq2Int=seq2Int,
                           mu=mu, r=r, expectedLength=expectedLength,
                           substModel=substModel, substModelBF=substModelBF)
    ansHessian <- res$hessian
    ans <- c("PAM"=res$par[1], "Mu"=mu, "r"=r, "negLogLikelihood"=res$value)
    invHessian <- solve(ansHessian)
    return(c(ans, "PAMVariance"=invHessian[1,1]))
  }else if(!is.null(mu) && !is.null(distance) && !is.null(r)){
    ## Just calculate the likelihood, given mu and distance
    ans <- .Call("TKF92LikelihoodFunctionWrapper", seq1Int, seq2Int, 
                 distance, mu, r, expectedLength, substModel, substModelBF)
    ansHessian <- hessian(function(x, seq1Int, seq2Int, expectedLength, 
                                   substModel, substModelBF){
                          ansTemp <- .Call("TKF92LikelihoodFunctionWrapper", 
                                           seq1Int, seq2Int, x[1], x[2], x[3], 
                                           expectedLength, substModel, 
                                           substModelBF)
                          return(ansTemp["negLogLikelihood"])
                 }, c(ans["PAM"], ans["Mu"], ans["r"]),
                 seq1Int=seq1Int, seq2Int=seq2Int,
                 expectedLength=expectedLength, substModel=substModel,
                 substModelBF=substModelBF)
    invHessian <- solve(ansHessian)
    return(c(ans, "PAMVariance"=invHessian[1,1],
             "MuVariance"=invHessian[2,2],
             "rVariance"=invHessian[3,3],
             "coVariancePAMMu"=invHessian[1,2],
             "coVariancePAMr"=invHessian[1,3],
             "coVarianceMur"=invHessian[2,3]
             )
           )
  }else{
    stop("You cannot estimate mu or r alone!")
  }
}

TKF92 <- function(fasta, mu=NULL, r=NULL, 
                  method=c("NM", "constrOptim"),
                  expectedLength=362,
                  substModel, substModelBF,
                  skipFailure=FALSE){
  method <- match.arg(method)
  seqnames <- names(fasta)
  nSeqs <- length(fasta)
  distanceMatrix <- matrix(0, ncol=nSeqs, nrow=nSeqs,
                           dimnames=list(seqnames, seqnames))
  varianceMatrix <- distanceMatrix
  negLoglikelihoodMatrix <- distanceMatrix
  if(is.null(mu) && is.null(r)){
    muMatrix <- distanceMatrix
    rMatrix <- distanceMatrix
  }
  for(i in 1:(nSeqs-1L)){
    for(j in (i+1L):nSeqs){
      message(seqnames[i], " vs ", seqnames[j])
      if(skipFailure){
        ans <- try(TKF92Pair(fasta[[i]], fasta[[j]],
                             mu=mu, r=r, method=method, 
                             expectedLength=expectedLength,
                             substModel=substModel, substModelBF=substModelBF))
        if(class(ans) == "try-error")
          next
      }else{
        ans <- TKF92Pair(fasta[[i]], fasta[[j]],
                         mu=mu, r=r, method=method,
                         expectedLength=expectedLength,
                         substModel=substModel, substModelBF=substModelBF)
      }
      distanceMatrix[i,j] <- distanceMatrix[j,i] <- ans["PAM"]
      varianceMatrix[i,j] <- varianceMatrix[j,i] <- ans["PAMVariance"]
      negLoglikelihoodMatrix[i,j] <- negLoglikelihoodMatrix[j,i] <-
        ans["negLogLikelihood"]
      if(is.null(mu) && is.null(r)){
        muMatrix[i,j] <- muMatrix[j,i] <- ans["Mu"]
        rMatrix[i,j] <- rMatrix[j,i] <- ans["r"]
      }
    }
  }
  if(is.null(mu) && is.null(r)){
    return(list(distanceMatrix=distanceMatrix,
                varianceMatrix=varianceMatrix,
                muMatrix=muMatrix,
                rMatrix=rMatrix,
                negLoglikelihoodMatrix=negLoglikelihoodMatrix))
  }else{
    return(list(distanceMatrix=distanceMatrix,
                varianceMatrix=varianceMatrix,
                negLoglikelihoodMatrix=negLoglikelihoodMatrix))
  }
}
