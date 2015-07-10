#include "TKF91_nlopt.h"

double TKF91LikelihoodFunction2D_nlopt(unsigned n, const double* x,
    double* grad, void* params){
  gsl_vector *x_opt = gsl_vector_alloc(2);
  gsl_vector_set(x_opt, 0, x[0]); // The distance
  gsl_vector_set(x_opt, 1, x[1]); // The mu
  double likelihood = TKF91LikelihoodFunction2D(x_opt, params);
  gsl_vector_free(x_opt);
  return likelihood;
}

/********************************************************************
 * nlopt main functions
 * *****************************************************************/
SEXP TKF91LikelihoodFunction2DMain_nlopt(SEXP seq1IntR, SEXP seq2IntR,
    SEXP expectedLength, SEXP probMatR, SEXP eqFrequenciesR){
  int ncol, nrow;
  ncol = INTEGER(GET_DIM(probMatR))[1];
  nrow = INTEGER(GET_DIM(probMatR))[0];
  int i, j;
  
  // probMat
  gsl_matrix *probMat = gsl_matrix_alloc(nrow, ncol);
  for(i = 0; i < nrow; i++)
    for(j = 0; j < ncol; j++)
      gsl_matrix_set(probMat, i, j, REAL(probMatR)[i+j*ncol]);

  // eqFrequenciesR
  gsl_vector *eqFrequencies = gsl_vector_alloc(GET_LENGTH(eqFrequenciesR));
  for(i = 0; i < GET_LENGTH(eqFrequenciesR); i++){
    gsl_vector_set(eqFrequencies, i, REAL(eqFrequenciesR)[i]);
  }

  // seqInt preparation
  int *seq1Int, *seq2Int;
  seq1Int = (int *) R_alloc(GET_LENGTH(seq1IntR), sizeof(int));
  seq2Int = (int *) R_alloc(GET_LENGTH(seq2IntR), sizeof(int));
  for(i = 0; i < GET_LENGTH(seq1IntR); i++){
    seq1Int[i] = INTEGER(seq1IntR)[i];
  }
  for(i = 0; i < GET_LENGTH(seq2IntR); i++){
    seq2Int[i] = INTEGER(seq2IntR)[i];
  }

  // construct params
  struct TKF91LikelihoodFunction2D_params params;
  params.len = REAL(expectedLength)[0];
  params.substModel = probMat;
  params.eqFrequencies = eqFrequencies;
  params.seq1Int = seq1Int;
  params.seq2Int = seq2Int;
  params.SA = GET_LENGTH(seq1IntR);
  params.SB = GET_LENGTH(seq2IntR);

  // nlopt main procedure
  double lb[2] = {0.0494497, 1e-20}; // lower bounds
  double ub[2] = {2000, 1-1e-20};    // upper bounds

  nlopt_opt opt;
  opt = nlopt_create(NLOPT_LN_NELDERMEAD, 2); /* algorithm and dimensionality */
  nlopt_set_lower_bounds(opt, lb);
  nlopt_set_upper_bounds(opt, ub);
  
  nlopt_set_min_objective(opt, TKF91LikelihoodFunction2D_nlopt, &params);
  nlopt_set_ftol_rel(opt, 1e-9); // stopping criteria

  double x[2] = {100, exp(-3)};  /* some initial guess */
  double minf; /* the minimum objective value, upon return */
  if (nlopt_optimize(opt, x, &minf) < 0) {
    printf("nlopt failed!\n");
  }else{
    printf("found minimum at f(%g,%g) = %0.10g\n", x[0], x[1], minf);
  }

  // free everything
  nlopt_destroy(opt);
  gsl_vector_free(eqFrequencies);
  gsl_matrix_free(probMat);

  return R_NilValue;
}

