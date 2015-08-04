#include "TKF91_nlopt.h"
#include "TKF92_nlopt.h"
#include "TKF92HG_nlopt.h"

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

static const R_CallMethodDef callMethods[] = {
  /* TKF91.c */
  CALLMETHOD_DEF(TKF91LikelihoodFunctionWrapper, 7),
  CALLMETHOD_DEF(TKF91LikelihoodFunction1DMain, 6),
  CALLMETHOD_DEF(TKF91LikelihoodFunction2DMainNM, 5),
  CALLMETHOD_DEF(TKF91LikelihoodFunction2DMain_nlopt, 6),

  /* TKF92.c */
  CALLMETHOD_DEF(TKF92LikelihoodFunctionWrapper, 8),
  CALLMETHOD_DEF(TKF92LikelihoodFunction1DMain, 7),
  CALLMETHOD_DEF(TKF92LikelihoodFunction3DMainNM, 5),
  CALLMETHOD_DEF(TKF92LikelihoodFunction3DMain_nlopt, 6),

  /* TKF92HG.c */
  CALLMETHOD_DEF(TKF92HGLikelihoodFunctionWrapper, 10),
  CALLMETHOD_DEF(TKF92HGLikelihoodFunction1DMain, 9),
  CALLMETHOD_DEF(TKF92HGLikelihoodFunction5DMainNM, 5),
  CALLMETHOD_DEF(TKF92HGLikelihoodFunction5DMain_nlopt, 6),

  {NULL, NULL, 0}
};

void R_init_TKF(DllInfo *info)
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  return;
}

