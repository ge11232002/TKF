#include "TKF92HG_nlopt.h"

/********************************************************************
 * Constants regarding the optimisation
 * *****************************************************************/
#define MAX_ITER 30000
#define F_TOL 1e-5   // Relative tolerance on function value


double TKF92HGLikelihoodFunction5D_nlopt(unsigned n, const double* x,
    double* grad, void* params){
  R_CheckUserInterrupt();
  gsl_vector *x_opt = gsl_vector_alloc(5);
  gsl_vector_set(x_opt, 0, x[0]); // The distance
  gsl_vector_set(x_opt, 1, x[1]); // The mu
  gsl_vector_set(x_opt, 2, x[2]); // The r
  gsl_vector_set(x_opt, 3, x[3]); // The ps
  gsl_vector_set(x_opt, 4, x[4]); // The kf
  double likelihood = TKF92HGLikelihoodFunction5D(x_opt, params);
  gsl_vector_free(x_opt);
  return likelihood;
}

/********************************************************************
 * nlopt main function
 * *****************************************************************/
SEXP TKF92HGLikelihoodFunction5DMain_nlopt(SEXP seq1IntR, SEXP seq2IntR,
    SEXP expectedLength, SEXP probMatR, SEXP eqFrequenciesR,
    SEXP method){
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
  struct TKF92HGLikelihoodFunction5D_params params;
  params.len = REAL(expectedLength)[0];
  params.substModel = probMat;
  params.eqFrequencies = eqFrequencies;
  params.seq1Int = seq1Int;
  params.seq2Int = seq2Int;
  params.SA = GET_LENGTH(seq1IntR);
  params.SB = GET_LENGTH(seq2IntR);

  // nlopt main procedure
  double lb[5] = {0.0494497, 1e-20, 1e-20, 1e-20, 1}; // lower bounds
  double ub[5] = {2000, 1-1e-20, 1-1e-20, 1-1e-20, 20};  // upper bounds
  //double dx[3] = {20, 0.01, 0.1}; // The initial step size

  nlopt_opt opt;
  if(strcmp(CHAR(STRING_ELT(method, 0)), "NM") == 0){
    opt = nlopt_create(NLOPT_LN_NELDERMEAD, 5); /* algorithm and dimensionality */
  }else if(strcmp(CHAR(STRING_ELT(method, 0)), "Sbplx") == 0){
    opt = nlopt_create(NLOPT_LN_SBPLX, 5);
  }else if(strcmp(CHAR(STRING_ELT(method, 0)), "COBYLA") == 0){
    opt = nlopt_create(NLOPT_LN_COBYLA, 5);
  }else if(strcmp(CHAR(STRING_ELT(method, 0)), "BOBYQA") == 0){
    opt = nlopt_create(NLOPT_LN_BOBYQA, 5);
  }else if(strcmp(CHAR(STRING_ELT(method, 0)), "PRAXIS") == 0){
    opt = nlopt_create(NLOPT_LN_PRAXIS, 5);
  }else{
    error("Wrong optimisation method!");
  }

  nlopt_set_lower_bounds(opt, lb);
  nlopt_set_upper_bounds(opt, ub);

  nlopt_set_min_objective(opt, TKF92HGLikelihoodFunction5D_nlopt, &params);
  nlopt_set_ftol_rel(opt, F_TOL); // stopping criteria
  //nlopt_set_initial_step(opt, dx); // initial step size
  nlopt_set_maxeval(opt, MAX_ITER);

  double x[5] = {100, exp(-3), 0.5, 0.5, 1.5};  /* some initial guess */
  double minf; /* the minimum objective value, upon return */
  if (nlopt_optimize(opt, x, &minf) < 0) {
    Rprintf("nlopt failed!\n");
  }else{
    Rprintf("found minimum at f(%g,%g,%g,%g,%g) = %0.10g using %s algorithm\n",
        x[0], x[1], x[2], x[3], x[4], minf, CHAR(STRING_ELT(method, 0)));
  }

  SEXP ans, ansNames;
  PROTECT(ans = NEW_NUMERIC(6)); // a vector of distance, mu, r, ps, kf and the negative log likelihood
  PROTECT(ansNames = NEW_CHARACTER(6));
  REAL(ans)[0] = x[0];
  REAL(ans)[1] = x[1];
  REAL(ans)[2] = x[2];
  REAL(ans)[3] = x[3];
  REAL(ans)[4] = x[4];
  REAL(ans)[5] = minf;
  SET_STRING_ELT(ansNames, 0, mkChar("PAM"));
  SET_STRING_ELT(ansNames, 1, mkChar("Mu"));
  SET_STRING_ELT(ansNames, 2, mkChar("r"));
  SET_STRING_ELT(ansNames, 3, mkChar("Ps"));
  SET_STRING_ELT(ansNames, 4, mkChar("Kf"));
  SET_STRING_ELT(ansNames, 5, mkChar("negLogLikelihood"));
  SET_NAMES(ans, ansNames);

  // free everything
  nlopt_destroy(opt);
  gsl_vector_free(eqFrequencies);
  gsl_matrix_free(probMat);

  UNPROTECT(2);
  return ans;


}

