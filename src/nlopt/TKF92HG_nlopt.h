#include <nlopt.h>
#include "TKF92HG.h"

/********************************************************************
 * TKF92HG nlopt objective function
 * *****************************************************************/
double TKF92HGLikelihoodFunction5D_nlopt(unsigned n, const double* x,
    double* grad, void* params);

/**** TKF92HG with nlopt implementation ****/
SEXP TKF92HGLikelihoodFunction5DMain_nlopt(SEXP seq1IntR, SEXP seq2IntR,
    SEXP expectedLength, SEXP probMatR, SEXP eqFrequenciesR,
    SEXP method);

