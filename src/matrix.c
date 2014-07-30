/*************************************************************************
    > File Name: matrix.c
    > Author: Ge Tan
    > Mail: gtan@me.com 
    > Created Time: Wed Jul 30 22:55:24 2014
 ************************************************************************/

#include<stdio.h>
#include<gsl/gsl_matrix.h>
#include "matrix.h"

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

//gsl_matrix *PAMnC(gsl_matrix *PAM1, const int n){

//}
