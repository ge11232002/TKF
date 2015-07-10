#include <nlopt.h>
#include "TKF92.h"

/********************************************************************
 * TKF92 nlopt objective function
 * *****************************************************************/
double TKF92LikelihoodFunction3D_nlopt(unsigned n, const double* x, 
    double* grad, void* params);

/**** TKF91 with nlopt implementation ****/
SEXP TKF92LikelihoodFunction3DMain_nlopt(SEXP seq1IntR, SEXP seq2IntR,
    SEXP expectedLength, SEXP probMatR, SEXP eqFrequenciesR,
    SEXP method);


