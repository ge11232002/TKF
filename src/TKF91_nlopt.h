#include <nlopt.h>
#include "TKF91.h"

/********************************************************************
 * TKF91 nlopt objective function
 * *****************************************************************/
double TKF91LikelihoodFunction2D_nlopt(unsigned n, const double* x, 
    double* grad, void* params);

/**** TKF91 with nlopt implementation ****/
SEXP TKF91LikelihoodFunction2DMain_nlopt(SEXP seq1IntR, SEXP seq2IntR,
    SEXP expectedLength, SEXP probMatR, SEXP eqFrequenciesR,
    SEXP method);


