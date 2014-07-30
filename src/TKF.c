/*************************************************************************
    > File Name: TKF.c
    > Author: Ge Tan
    > Mail: gtan@me.com 
    > Created Time: Wed Jul 30 17:55:28 2014
 ************************************************************************/

#include<stdio.h>
#include<gsl/gsl_matrix.h>
#include "Rdefines.h"

/*-------------------------------------------------------------------
 * Print a GSL matrix in a proper shape
 * ----------------------------------------------------------------*/
int printGSLMatrix(const gsl_matrix *m){
  int status, n = 0;
  for(size_t i = 0; i < m->size1; i++){
    for(size_t j = 0; j < m->size2; j++){
      if((status = printf("%g ", gsl_matrix_get(m, i, j))) < 0)
        return -1;
       n += status;
    }
    if((status = printf("\n")) < 0)
      return -1;
    n += status;
  }
  return n;
}



SEXP TKF91(SEXP probMat){
  int ncol, nrow;
  ncol = INTEGER(GET_DIM(probMat))[1];
  nrow = INTEGER(GET_DIM(probMat))[0];
  int i, j; 
  // matrix allocation and setting
  gsl_matrix *probPAMMat = gsl_matrix_alloc(nrow, ncol);
  for(i = 0; i < nrow; i++) 
    for(j = 0; j < ncol; j++)
      gsl_matrix_set(probPAMMat, i, j, REAL(probMat)[i+j*ncol]);
  
  // print the matrix
  printGSLMatrix(probPAMMat);
  
  gsl_matrix_free(probPAMMat);

  Rprintf("the ncol is %d\n", ncol);
  Rprintf("the nrow is %d\n", nrow);
  return R_NilValue;
}

