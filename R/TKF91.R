#########################################################################
# File Name: ../R/TKF91.R
# Author: Ge Tan
# mail: gtan@me.com
# Created Time: Sun Sep  7 14:23:14 2014
#########################################################################

TKF91 <- function(seq1, seq2, mu=NULL, distance=NULL,
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
  }else if(!is.null(mu) && !is.null(distance)){
    ## Just calculate the likelihood, given mu and distance

  }

}
