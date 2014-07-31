/*************************************************************************
    > File Name: matrix.c
    > Author: Ge Tan
    > Mail: gtan@me.com 
    > Created Time: Wed Jul 30 22:55:24 2014
 ************************************************************************/

#include<stdio.h>
#include<gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
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

/********************************************************************
 * Make diagonal matrix from vector
 * *****************************************************************/

//void *my_diag_alloc(const gsl_vector_complex *X, gsl_matrix_complex *mat)
//{
//    gsl_vector_complex_view diag = gsl_matrix_diagonal(mat);
//    gsl_matrix_set_all(mat, 0.0); //or whatever number you like
//    gsl_vector_memcpy(&diag.vector, X);
//}

/********************************************************************
 * Matrix logarithm only for diagonalizable matrix.
 * *****************************************************************/

void logm(gsl_matrix *m, gsl_matrix *mlog){
  int nrow;
  //int i;
  nrow = m->size1;
  gsl_vector_complex *eval = gsl_vector_complex_alloc(nrow);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc (nrow, nrow);
  gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc(nrow);
  gsl_eigen_nonsymmv(m, eval, evec, w);
  
  // free the eigen nonsymmv workspace
  gsl_eigen_nonsymmv_free(w);
  gsl_eigen_nonsymmv_sort(eval, evec,
      GSL_EIGEN_SORT_ABS_DESC);
 
  {
    int i, j;

    for (i = 0; i < 20; i++)
      {
        gsl_complex eval_i 
           = gsl_vector_complex_get (eval, i);
        gsl_vector_complex_view evec_i 
           = gsl_matrix_complex_column (evec, i);

        printf ("eigenvalue = %g + %gi\n",
                GSL_REAL(eval_i), GSL_IMAG(eval_i));
        printf ("eigenvector = \n");
        for (j = 0; j < 20; ++j)
          {
            gsl_complex z = 
              gsl_vector_complex_get(&evec_i.vector, j);
            printf("%g + %gi\n", GSL_REAL(z), GSL_IMAG(z));
          }
      }
  }

  gsl_vector_complex *evalLog = gsl_vector_complex_alloc(nrow);
  //for(i = 0; i < nrow; i++){
  //  gsl_vector_set(evalLog, i, gsl_vector_get ());
  //}
  // free the allocated eigen values and vectors
  gsl_vector_complex_free(eval);
  gsl_matrix_complex_free(evec);
  gsl_vector_complex_free(evalLog);
}

