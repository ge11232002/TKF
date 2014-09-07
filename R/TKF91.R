#########################################################################
# File Name: ../R/TKF91.R
# Author: Ge Tan
# mail: gtan@me.com
# Created Time: Sun Sep  7 14:23:14 2014
#########################################################################




TKF911D <- function(seq1, seq2, mu=0.001, expectedLength=362,
                    substModel=GONNET){
  seq1Int <- AAToInt(seq1)
  seq2Int <- AAToInt(seq2)
  ## for the C matrix index
  seq1Int <- seq1Int - 1L
  seq2Int <- seq2Int - 1L
  .Call("TKF91LikelihoodFunction1DMain", seq1Int, seq2Int, mu,
        expectedLength, substModel)
}
