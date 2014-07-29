#########################################################################
# File Name: ScoringMatrix.R
# Author: Ge Tan
# mail: gtan@me.com
# Created Time: Tue Jul 29 20:48:49 2014
#########################################################################

### -----------------------------------------------------------------
### Calculate the n-PAM matrices from PAM1 mutation matrix and n.
### To compute n-PAM matrices, we multiply the PAM1 matrix through itself 
### N times, which is most efficiently achieved through 
### n additions in log space.
###
PAMn <- function(PAM1, n){
  ## Validated by Darwin
  require(expm)
  ans <- expm.Higham08(n * logm(PAM1))
  dimnames(ans) <- dimnames(PAM1)
  return(ans)
}
### PAM250 <- PAMn(GONNET, 250)

### -----------------------------------------------------------------
### Computing Dayhoff matrices from PAM mutation matrices and AA frequency.
Dayhoffn <- function(PAM1, BF, n){
  ## Validated by Darwin
  require(expm)
  pamn <- PAMn(PAM1, n)
  ans <- 10 * log10(sweep(pamn, MARGIN=2, BF, FUN="/"))
  return(ans)
}
### Dayhoff250 <- Dayhoffn(GONNET, GONNETBF, 250)



