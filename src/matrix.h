/*************************************************************************
    > File Name: .matrix.h
    > Author: Ge Tan
    > Mail: gtan@me.com 
    > Created Time: Wed Jul 30 22:57:00 2014
 ************************************************************************/

#include<stdio.h>
#include <gsl/gsl_matrix.h>


int printGSLMatrix(const gsl_matrix *m);
void PAMn(const gsl_matrix *m, const double distance, gsl_matrix *mPAM);

//gsl_matrix *PAMnC(gsl_matrix *PAM1, const int n);
