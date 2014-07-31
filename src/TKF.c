/*************************************************************************
    > File Name: TKF.c
    > Author: Ge Tan
    > Mail: gtan@me.com 
    > Created Time: Wed Jul 30 17:55:28 2014
 ************************************************************************/

#include<stdio.h>
#include<gsl/gsl_matrix.h>
#include "Rdefines.h"
#include "matrix.h"

SEXP TKF91LikelihoodFunction1D(SEXP distanceR, SEXP probMatR){
  int ncol, nrow;
  ncol = INTEGER(GET_DIM(probMatR))[1];
  nrow = INTEGER(GET_DIM(probMatR))[0];
  double distance;
  distance = REAL(distanceR)[0];
  int i, j; 
  // matrix allocation and setting
  gsl_matrix *probMat = gsl_matrix_alloc(nrow, ncol);
  gsl_matrix *probMatN = gsl_matrix_alloc(nrow, ncol);
  for(i = 0; i < nrow; i++) 
    for(j = 0; j < ncol; j++)
      gsl_matrix_set(probMat, i, j, REAL(probMatR)[i+j*ncol]);
  
  // print the GSL matrix
  printGSLMatrix(probMat);
  printf("\n");
  PAMn(probMat, distance, probMatN); 
  printGSLMatrix(probMat);
  printf("\n");
  printGSLMatrix(probMatN);
  // free the GSL matrix
  gsl_matrix_free(probMat);
  gsl_matrix_free(probMatN);

  Rprintf("the ncol is %d\n", ncol);
  Rprintf("the nrow is %d\n", nrow);
  return R_NilValue;
}




