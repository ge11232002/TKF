/*************************************************************************
    > File Name: .matrix.h
    > Author: Ge Tan
    > Mail: gtan@me.com 
    > Created Time: Wed Jul 30 22:57:00 2014
 ************************************************************************/

#include<stdio.h>
#include <gsl/gsl_matrix.h>
#define Eps 1e-12
#include "Rdefines.h"
#include "Rinternals.h"

int printGSLMatrix(const gsl_matrix *m);
void PAMn(const gsl_matrix *m, const double distance, gsl_matrix *mPAM);

