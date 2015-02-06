#########################################################################
# File Name: ../R/TKF91.R
# Author: Ge Tan
# mail: gtan@me.com
# Created Time: Sun Sep  7 14:23:14 2014
#########################################################################

TKF91Pair <- function(seq1, seq2, mu=NULL, distance=NULL,
                  ## mu: by default is 0.001 from median of mu values from Fungi dataset.
                  expectedLength=362, 
                  substModel=GONNET, substModelBF=GONNETBF){
  seq1Int <- AAToInt(seq1)
  seq2Int <- AAToInt(seq2)
  ## for the C matrix index
  seq1Int <- seq1Int - 1L
  seq2Int <- seq2Int - 1L

  expectedLength <- as.numeric(expectedLength)
  
  if(is.null(mu) && is.null(distance)){ 
    ## Do the 2D optimisation

  }else if(!is.null(mu) && is.null(distance)){
    ## Do the 1D distance optimisation
    ans <- .Call("TKF91LikelihoodFunction1DMain", seq1Int, seq2Int, mu,
                 expectedLength, substModel, substModelBF)
    ansHessian <- hessian(function(x, seq1Int, seq2Int, mu, expectedLength, substModel, substModelBF){
                          ansTemp <- .Call("TKF91LikelihoodFunctionWrapper", seq1Int, seq2Int, x, mu, expectedLength, substModel, substModelBF)
                          return(ansTemp["negLogLikelihood"])
                 }, ans["PAM"], 
                 seq1Int=seq1Int, seq2Int=seq2Int,
                 mu=mu, expectedLength=expectedLength,
                 substModel=substModel, 
                 substModelBF=substModelBF)
    return(c(ans, "PAMVariance"=ginv(ansHessian)[1,1]))
  }else if(!is.null(mu) && !is.null(distance)){
    ## Just calculate the likelihood, given mu and distance
    ans <- .Call("TKF91LikelihoodFunctionWrapper", seq1Int, seq2Int, 
                 distance, mu, expectedLength, substModel, substModelBF)
  }else{
    stop("You cannot estimate mu alone!")
  }
  return(ans)
}



