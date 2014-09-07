/*************************************************************************
    > File Name: TKF.c
    > Author: Ge Tan
    > Mail: gtan@me.com 
    > Created Time: Wed Jul 30 17:55:28 2014
 ************************************************************************/

#include<stdio.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include "Rdefines.h"
#include "matrix.h"
#include "MathFunctions.h"

struct TKF91LikelihoodFunction1D_params
{
      double len, mu;
      gsl_matrix *substModel;
      gsl_vector *seq1Int, *seq2Int;
};

double TKF91LikelihoodFunction1D(double distance, void *params){
  struct TKF91LikelihoodFunction1D_params *p = (struct TKF91LikelihoodFunction1D_params *) params;
  double len = p->len;
  double mu = p->mu;
  gsl_matrix *substModel = p->substModel;
  gsl_vector *seq1Int = p->seq1Int;
  gsl_vector *seq2Int = p->seq2Int;
  
  double lambda = len / (len + 1) * mu;
  double alpha = -mu * distance;
  double lmt = exp((lambda-mu)*distance);
  double lbeta = log1x(-exp((lambda-mu)*distance)) - (log(mu) + log1x(-lambda/mu * exp((lambda-mu)*distance)));
  double beta = exp(lbeta); // beta is  not a very small number.
  double lP1t = alpha + log1x(-lambda * beta);
  double lP11t = log1x(-exp(alpha) - mu * beta) + log1x(-lambda * beta);
  double lP12t = log1x(-lambda * beta); 
  double lP01t = log(mu) + lbeta;
  
  // test
  int i;
  for(i = 0; i < seq1Int->size; i++){
    printf ("v_%d = %g\n", i, gsl_vector_get (seq1Int, i));
  }

  return 0.1;
}

SEXP TKF91LikelihoodFunction1DMain(SEXP seq1Int, SEXP seq2Int, SEXP mu,
    SEXP expectedLength, SEXP probMatR){
  int ncol, nrow;
  ncol = INTEGER(GET_DIM(probMatR))[1];
  nrow = INTEGER(GET_DIM(probMatR))[0];
  //double distance;
  //distance = REAL(distanceR)[0];
  int i, j; 
  // matrix allocation and setting
  gsl_matrix *probMat = gsl_matrix_alloc(nrow, ncol);
  //gsl_matrix *probMatN = gsl_matrix_alloc(nrow, ncol);
  for(i = 0; i < nrow; i++) 
    for(j = 0; j < ncol; j++)
      gsl_matrix_set(probMat, i, j, REAL(probMatR)[i+j*ncol]);
 
  // seq preparation
  gsl_vector *seq1IntGSL = gsl_vector_alloc(GET_LENGTH(seq1Int));
  gsl_vector *seq2IntGSL = gsl_vector_alloc(GET_LENGTH(seq2Int));
  for(i = 0; i < GET_LENGTH(seq1Int); i++){
    gsl_vector_set(seq1IntGSL, i, INTEGER(seq1Int)[i]);
  }
  for(i = 0; i < GET_LENGTH(seq2Int); i++){
    gsl_vector_set(seq2IntGSL, i, INTEGER(seq2Int)[i]);
  }
  for(i = 0; i < seq1IntGSL->size; i++){
    printf ("v_%d = %g\n", i, gsl_vector_get (seq1IntGSL, i));
  }
  // GSL minimizer 
  /*int status;
  int iter = 0, max_iter = 100;
  const gsl_min_fminimizer_type *T;
  gsl_min_fminimizer *s;
  gsl_function F;
  struct TKF91LikelihoodFunction1D_params *params;
  params->len = REAL(expectedLength)[0];
  params->mu = REAL(mu)[0];
  params->substModel = probMat;
  params->seq1Int = seq1IntGSL;
  params->seq2Int = seq2IntGSL;
  F.function = &TKF91LikelihoodFunction1D;
  F.params = params;
  double x_lo = 0, x_hi = 100000; 
  double x = 0.11;
  double mEps = 0.001;
  T = gsl_min_fminimizer_brent;
  s = gsl_min_fminimizer_alloc (T);
  gsl_min_fminimizer_set (s, &F, x, x_lo, x_hi);
  
  printf("using %s method\n", 
      gsl_min_fminimizer_name (s));
  printf("%5s [%9s, %9s] %9s %10s %9s\n",
      "iter", "lower", "upper", "min", 
      "err", "err(est)");
  do
    {
      iter++;
      status = gsl_min_fminimizer_iterate (s);
      x = gsl_min_fminimizer_x_minimum (s);
      x_lo = gsl_min_fminimizer_x_lower (s);
      x_hi = gsl_min_fminimizer_x_upper (s);
      status = gsl_min_test_interval(x_lo, x_hi,
                                     mEps*mEps, mEps);
      if (status == GSL_SUCCESS)
        printf ("Converged:\n");
      printf ("%5d [%.7f, %.7f] "
              "%.7f %.7f\n",
              iter, x_lo,  x_hi,
              x, b - a);
    }
  while (status == GSL_CONTINUE && iter < max_iter);
*/
  //PAMn(probMat, distance, probMatN); 
  //printGSLMatrix(probMatN);
  
  // free the GSL matrix
  gsl_matrix_free(probMat);
  //gsl_min_fminimizer_free (s);
  gsl_vector_free(seq1IntGSL);
  gsl_vector_free(seq2IntGSL);
  //gsl_matrix_free(probMatN);

  return R_NilValue;
}

