/*! \file LLR_Mp_Av.c
\brief MaxPhase fitness function.

The function LLR_Mp_Av() is a wrapper around the function
LogLikelihoodRatioMP5(). Any fitness function sent to the
PSO function ptapso() should accept standardized particle positions
(where each coordinate lies in [0,1]). LLR_Mp_Av() accepts
standardized particle coordinates and converts then to real
coordinates before calling LogLikelihoodRatioMP5().
 */
#include "LLR_Mp_Av.h"
#include "gslhdf5_io.h"
#include "hdf5_hl.h"
#include "loadRAAPTR.h"
#include "maxphase.h"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_roots.h> // root finding functions
#include <gsl/gsl_sort.h>
#include <stdio.h>
#include <stdlib.h>

/*! Deep copy of fitness function parameter struct. Copies are required
 in an OpenMP code since sharing this struct will lead to overwriting
of critical fields. */
struct fitFuncParams *ffparams_clone(struct fitFuncParams *srcffp) {

  struct fitFuncParams *dstffp = ffparam_alloc(srcffp->nDim);

  struct llr_pso_params *srcllp = (struct llr_pso_params *)srcffp->splParams;

  struct llr_pso_params *llp = llrparam_alloc(srcllp->N, srcllp->Np);

  gsl_vector_memcpy(dstffp->rmin, srcffp->rmin);
  gsl_vector_memcpy(dstffp->rangeVec, srcffp->rangeVec);
  gsl_vector_memcpy(dstffp->realCoord, srcffp->realCoord);
  dstffp->fitEvalFlag = srcffp->fitEvalFlag;

  size_t lpc1, lpc2;
  for (lpc1 = 0; lpc1 < srcllp->Np; lpc1++) {
    llp->sd[lpc1] = srcllp->sd[lpc1];
    llp->alphaP[lpc1] = srcllp->alphaP[lpc1];
    llp->deltaP[lpc1] = srcllp->deltaP[lpc1];
    llp->phiI[lpc1] = srcllp->phiI[lpc1];
    for (lpc2 = 0; lpc2 < srcllp->N; lpc2++) {
      llp->s[lpc1][lpc2] = srcllp->s[lpc1][lpc2];
    }
  }
  for (lpc2 = 0; lpc2 < srcllp->N; lpc2++) {
    llp->yr[lpc2] = srcllp->yr[lpc2];
  }

  dstffp->splParams = llp;

  return dstffp;
}

struct fitFuncParams *RAAPTR_clone(struct fitFuncParams *srcffp) {

  struct fitFuncParams *dstffp = ffparam_alloc(srcffp->nDim);

  struct RAAPTR_data *srcllp = (struct RAAPTR_data *)srcffp->splParams;

  struct RAAPTR_data *llp = llrparam_alloc_RAAPTR(srcllp->Np);
  // dynamic alloc memory for llp
  for (size_t psr = 0; psr < srcllp->Np; psr++) {
    llp->s[psr] = (double *)malloc(sizeof(double) * srcllp->N[psr]);
    llp->sd[psr] = (double *)malloc(sizeof(double) * srcllp->N[psr]);
    llp->yr[psr] = (double *)malloc(sizeof(double) * srcllp->N[psr]);
  }

  gsl_vector_memcpy(dstffp->rmin, srcffp->rmin);
  gsl_vector_memcpy(dstffp->rangeVec, srcffp->rangeVec);
  gsl_vector_memcpy(dstffp->realCoord, srcffp->realCoord);
  dstffp->fitEvalFlag = srcffp->fitEvalFlag;

  size_t lpc1, lpc2;
  for (lpc1 = 0; lpc1 < srcllp->Np; lpc1++) {
    llp->alphaP[lpc1] = srcllp->alphaP[lpc1];
    llp->deltaP[lpc1] = srcllp->deltaP[lpc1];
    llp->phiI[lpc1] = srcllp->phiI[lpc1];
    llp->N[lpc1] = srcllp->N[lpc1];
    for (lpc2 = 0; lpc2 < srcllp->N[lpc1]; lpc2++) {
      llp->s[lpc1][lpc2] = srcllp->s[lpc1][lpc2];
      llp->sd[lpc1][lpc2] = srcllp->sd[lpc1][lpc2];
      llp->yr[lpc1][lpc2] = srcllp->yr[lpc1][lpc2];
    }
  }

  dstffp->splParams = llp;

  return dstffp;
}

/*!
fitFuncVal = LLR_Mp_Av(xVec,P)
The fitness value is returned in fitFuncVal.
xVec is standardized, that is 0<=xVec(i,j)<=1.
The values used to convert xVec(i,j)
internally before computing fitness are given in P.rmin and
P.rangeVec:
xVec(j) -> xVec(j)*rangevec(j)+rmin(j).
fitFuncVal = infty if the point xVec falls
outside the hypercube defined by 0<=xVec(j)<=1.
The real coordinates are returned in P.realCoord.

\author Y. Wang
\date Jan 2016
\remark Derived from ptapsotestfunc.c.
*/
double
LLR_av(gsl_vector *xVec,     /*!< Standardized Particle Coordinates*/
       void *inParamsPointer /*!< Fitness function parameter structure
                                                    containing information for
                                conversion of standardized to real coordinates*/
) {

  // Set to 0 if angular variables do not have a periodic boundary conditions
  size_t wrapAngles = 1;

  size_t validPt;
  /* Cast from void * to known structure pointer before any
           of the fields can be accessed.
         */
  struct fitFuncParams *inParams = (struct fitFuncParams *)inParamsPointer;

  double fitFuncVal;

  validPt = llrpsochkcoord(wrapAngles, inParams->rmin, inParams->rangeVec, xVec,
                           inParams->realCoord);

  if (validPt) {
    inParams->fitEvalFlag = 1;
    fitFuncVal = AvPhaseLLR(inParams);
  } else {
    inParams->fitEvalFlag = 0;
    fitFuncVal = GSL_POSINF;
  }
  return fitFuncVal;
}

double f(double x, void *p) {
  // double alpha0 = *(double *) params[0];
  struct avPhase_param *params = (struct avPhase_param *)p;
  double b[5], norm;
  b[0] = (params->p0);
  b[1] = (params->p1);
  b[2] = (params->p2);
  b[3] = (params->p3);
  b[4] = (params->p4);
  norm = (params->p5); // p5 is 'norm'
  // printf("f: b = %f\n", b[0]);
  // printf("f: norm = %f\n", norm);
  // printf("f: p->b = %f\n", (*params).b[0]);
  // double *norm = params->norm;
  // double norm1 = norm;
  // printf("f: norm = %f\n", norm1);
  // printf("f: b = %f\n", *(b+2);
  // b[1]=bb[1];
  // b[2]=bb[2];
  // b[3]=bb[3];
  // b[4]=bb[4];
  double f = exp(b[1] * cos(x) + b[2] * sin(x) + b[3] * sin(2.0 * x) +
                 b[4] * pow(cos(x), 2.0) - norm);
  return f;
}

double
LLR_mp(gsl_vector *xVec,     /*!< Standardized Particle Coordinates*/
       void *inParamsPointer /*!< Fitness function parameter structure
                                                    containing information for
                                conversion of standardized to real coordinates*/
) {

  // Set to 0 if angular variables do not have a periodic boundary conditions
  size_t wrapAngles = 1;

  size_t validPt;
  /* Cast from void * to known structure pointer before any
           of the fields can be accessed.
         */
  struct fitFuncParams *inParams = (struct fitFuncParams *)inParamsPointer;

  double fitFuncVal;

  validPt = llrpsochkcoord(wrapAngles, inParams->rmin, inParams->rangeVec, xVec,
                           inParams->realCoord);

  if (validPt) {
    inParams->fitEvalFlag = 1;
    fitFuncVal = LogLikelihoodRatioMP5(inParams);
  } else {
    inParams->fitEvalFlag = 0;
    fitFuncVal = GSL_POSINF;
  }

  return fitFuncVal;
}
// EOF LLR_Mp_Av

double LLR_av_RAAPTR(
    gsl_vector *xVec,     /*!< Standardized Particle Coordinates*/
    void *inParamsPointer /*!< Fitness function parameter structure
                                                 containing information for
                             conversion of standardized to real coordinates*/
) {

  // Set to 0 if angular variables do not have a periodic boundary conditions
  size_t wrapAngles = 1;

  size_t validPt;
  /* Cast from void * to known structure pointer before any
           of the fields can be accessed.
         */
  struct fitFuncParams *inParams = (struct fitFuncParams *)inParamsPointer;

  double fitFuncVal;

  validPt = llrpsochkcoord(wrapAngles, inParams->rmin, inParams->rangeVec, xVec,
                           inParams->realCoord);

  if (validPt) {
    inParams->fitEvalFlag = 1;
    fitFuncVal = AvPhaseLLR_RAAPTR(inParams);
  } else {
    inParams->fitEvalFlag = 0;
    fitFuncVal = GSL_POSINF;
  }
  return fitFuncVal;
}

// average/marginalize LLR over pulsar phases, Jan. 2017
double AvPhaseLLR(struct fitFuncParams *inParams) {

  struct llr_pso_params *splParams =
      (struct llr_pso_params *)inParams->splParams;

  struct cfunc_OUTPUT *output;
  output = (struct cfunc_OUTPUT *)malloc(1 * sizeof(struct cfunc_OUTPUT));
  output->c = (double *)malloc(4 * sizeof(double));
  output->v = (double *)malloc(9 * sizeof(double));

  // struct avPhase_param * avPhaseParam;
  // avPhaseParam = (struct avPhase_param *)malloc(1 * sizeof(struct
  // avPhase_param)); avPhaseParam->b = (double *)malloc(4 * sizeof(double));
  // avPhaseParam->norm = (double *)malloc(1 * sizeof(double));

  double tmp;
  double res;
  unsigned int lpr, i, j, j1, j2, j3, jj;
  //  size_t pp[6];   not used.
  const size_t kk = 1, stride = 1, nn = 6; // gsl_sort_largest_index
  // double src[6] = {1.0, 3.0, 6.0, 4.0, 5.0, 2.0};  //
  //  src = (double *)malloc(6 * sizeof(double));
  //  *(src+0)=1.0;
  //  *(src+1)=3.0;
  //  *(src+2)=6.0;
  //  *(src+3)=4.0;
  //  *(src+4)=5.0;
  //  *(src+5)=2.0;

  unsigned int Np = splParams->Np; // number of pulsars in PTA
  unsigned int N = splParams->N;   // number of samples

  // transfer parameters from structure inParams
  // printf("MP5: Np = %d\n", Np);
  // printf("MP5: N = %d\n", N);
  double *yr;
  // yr = (double *)malloc(N * sizeof(double));
  double *sd;
  // sd = malloc(Np * sizeof(double));
  double *alphaP, *deltaP;
  // alphaP = (double *)malloc(Np * sizeof(double));
  // deltaP = (double *)malloc(Np * sizeof(double));
  sd = splParams->sd;
  alphaP = splParams->alphaP;
  deltaP = splParams->deltaP;
  yr = splParams->yr;

  double *Phi;
  Phi = (double *)malloc(N * sizeof(double));
  double theta;
  double alpha, delta, omega, phi0, Amp, iota, thetaN;
  double *phiI;

  unsigned int nDim = inParams->nDim; // dimension of fitness function
  // printf("MP5: nDim = %d\n", nDim);

  double b[5]; // quartic equation coefficients, closed form solution
  // b[0]=1.0;
  // b[1]=2.0;
  // b[2]=1.2;
  // b[3]=2.5;
  // b[4]=0.8;

  // avPhaseParam.b[0] = b[0];
  // printf("LLR_PSOav: avPhaseParam.b = %f\n",);
  // double z[8];  // 4 (real) + 4 (complex) solutions from
  // gsl_poly_complex_solve

  // unsigned int nr;  // number of EFFECTIVE roots (real number && abs(r)<1)

  double LLR = 0.0; // log likelihood ratio
  double **s;
  double C = 0.0;

  gsl_vector *skyLocSrc =
      gsl_vector_calloc(3); // sky location of source in Cartesian coordinate
  gsl_vector *skyLocPulsar =
      gsl_vector_calloc(3); // sky location of pulsars in Cartesian coordinate
  // gsl_vector * realCoord = gsl_vector_calloc(nDim);
  // gsl_vector_memcpy(realCoord,inParams->realCoord);

  alpha = gsl_vector_get(inParams->realCoord, 0);
  delta = gsl_vector_get(inParams->realCoord, 1);
  omega = gsl_vector_get(inParams->realCoord, 2);
  phi0 = gsl_vector_get(inParams->realCoord, 3);
  Amp = gsl_vector_get(inParams->realCoord, 4);
  Amp = pow(10, Amp); // physical amplitude
  gsl_vector_set(inParams->realCoord, 4, Amp);
  iota = gsl_vector_get(inParams->realCoord, 5);
  thetaN = gsl_vector_get(inParams->realCoord, 6);

  s = (double **)malloc(Np * sizeof(double));
  for (i = 0; i < Np; i++) {
    //*(s+i) = (double *)malloc(N * sizeof(double));  // not needed!
    s[i] = splParams->s[i];
  }

  phiI = malloc(Np * sizeof(double));
  // printf("MP5: s[0][1] = %f\n", *(*(s+0)+1));
  // printf("MP5: s[1][2] = %f\n", *(*(s+1)+5));

  double *norm, *LRn, M, *norm2;
  norm = malloc(4 * sizeof(double));
  // norm1 = malloc(4 * sizeof(double));
  gsl_vector *norm1 = gsl_vector_alloc(4);
  norm2 = malloc(4 * sizeof(double));
  LRn = malloc(4 * sizeof(double));
  double sign[4][5] = {{1.0, 1.0, 1.0, 1.0, 1.0},
                       {1.0, -1.0, 1.0, -1.0, 1.0},
                       {1.0, -1.0, -1.0, 1.0, 1.0},
                       {1.0, 1.0, -1.0, -1.0, 1.0}};
  double intup[4] = {M_PI_2, M_PI, 3 * M_PI_2,
                     2 * M_PI};                     // upper limits of quadrants
  double intlow[4] = {0, M_PI_2, M_PI, 3 * M_PI_2}; // lower limits of quadrants

  gsl_vector_set(skyLocSrc, 0, cos(delta) * cos(alpha));
  gsl_vector_set(skyLocSrc, 1, cos(delta) * sin(alpha));
  gsl_vector_set(skyLocSrc, 2, sin(delta));

  for (i = 0; i < N; i++) {
    Phi[i] = yr[i] * omega;
    // printf("MP5: *(Phi+i) = %e, *(yr+i) = %e\n", *(Phi+i), *(yr+i));
  }

  double bs[5], tmp0, NN, result, error;

  gsl_integration_workspace *w =
      gsl_integration_workspace_alloc(1000); // # double precision intervals

  // double bx[2];
  // bx[0]= 1.0;
  // bx[1]= 2.0;
  // double test;
  // test=9.999;
  //
  // norm[0] = 1.0;
  // norm[1] = 1.1;
  // norm[2] = 1.2;
  // norm[3] = 1.3;

  // gsl_function F;
  // struct avPhase_param avPhaseParam = {b[0],b[1],b[2],b[3],b[4],test};
  // //avPhaseParam->b = b;
  // //avPhaseParam->norm = &test;
  // //printf("LLR_PSOav: avPhaseParam.b = %f\n",(*avPhaseParam).b[1]);
  // //printf("LLR_PSOav: test = %f\n",*test);
  // //printf("LLR_PSOav: avPhaseParam.norm = %f\n",(*avPhaseParam).norm);
  // F.function = &f;
  // F.params = &avPhaseParam;
  // //F.norm = 2.22;
  // gsl_integration_qags(&F,0,1.57,0,1e-5,1000,w,&result,&error);
  // printf("LLR_PSOav: result = %.18f\n",result );
  // printf("AvPhaseLLR: Np = %d\n", Np);

  for (i = 0; i < Np; i++) {
    // printf("AvPhaseLLR: i = %d\n", i);
    M = 0.0;
    // printf("MP5: i = %d\n", i);
    // printf("alphaP = %f, deltaP = %f \n", *(alphaP+i), *(deltaP+i));
    gsl_vector_set(skyLocPulsar, 0, cos(deltaP[i]) * cos(alphaP[i]));
    gsl_vector_set(skyLocPulsar, 1, cos(deltaP[i]) * sin(alphaP[i]));
    gsl_vector_set(skyLocPulsar, 2, sin(deltaP[i]));

    gsl_blas_ddot(skyLocSrc, skyLocPulsar, &res);
    // printf("MP5: kp*k = %f\n", res);

    theta = acos(res);
    // printf("theta = %f\n", theta);

    cfunc(N, alpha, delta, alphaP[i], deltaP[i], theta, Amp, omega, iota,
          thetaN, phi0, Phi, s[i], (sd + i), output);
    tmp = (*output).v[2] - 0.5 * ((*output).v[6] + (*output).v[8]);
    b[0] = tmp;
    tmp = (*output).v[0] - (*output).v[5];
    b[1] = tmp;
    tmp = (*output).v[1] - (*output).v[7];
    b[2] = tmp;
    tmp = -0.5 * (*output).v[4];
    b[3] = tmp;
    tmp = 0.5 * ((*output).v[6] - (*output).v[3]);
    b[4] = tmp;

    // quadrant I-IV
    for (j = 0; j < 4; j++) {

      for (j1 = 0; j1 < 5; j1++) {
        bs[j1] = b[j1] * sign[j][j1];
      }

      tmp0 = 0;
      for (jj = 1; jj < 5; jj++) {
        if (bs[jj] > 0) {
          tmp0 = tmp0 + bs[jj];
        }
      }
      norm[j] = tmp0;

      // printf("AvPhaseLLR: PSR i = %d, norm = %f\n",i,norm[j]);

      gsl_function F;
      struct avPhase_param avPhaseParam = {b[0], b[1], b[2],
                                           b[3], b[4], norm[j]};
      F.function = &f;
      F.params = &avPhaseParam;

      gsl_integration_qags(&F, intlow[j], intup[j], 0, 1e-5, 1000, w, &result,
                           &error);
      // printf("result = %.18f\n",result );

      LRn[j] = result;

      if (gsl_isinf(LRn[j])) {
        printf("AvPhaseLLR: Inf for PSR i = %d\n", i);
        printf("LRn at j = %d\n", j);
      } else if (gsl_isnan(LRn[j])) {
        printf("AvPhaseLLR: NAN for PSR i = %d\n", i);
        printf("LRn at j = %d\n", j);
      }

      // norm1[j] = norm[j] + log(LRn[j] + b[1]);
      gsl_vector_set(norm1, j,
                     norm[j] + log(LRn[j]) + b[0]); // b[0] = b[1] in Matlab
    }

    // printf("AvPhaseLLR: = %f\n",norm1[j]);

    // NN = gsl_max_dbl(norm1);
    NN = gsl_vector_max(norm1);

    for (j2 = 0; j2 < 4; j2++) {
      norm2[j2] = gsl_vector_get(norm1, j2) - NN;
    }

    for (j3 = 0; j3 < 4; j3++) {
      M += exp(norm2[j3]);
    }

    LLR = LLR + NN + log(M);
  }

  gsl_integration_workspace_free(w);
  free(output->c);
  free(output->v);
  free(output);
  free(Phi);
  free(s);
  // for (i = 0; i < Np; i++) {
  //   free(phiItmp[i]);
  //   free(lh[i]);
  // }
  // free(phiItmp);
  // free(lh);
  free(phiI);
  gsl_vector_free(skyLocSrc);
  gsl_vector_free(skyLocPulsar);
  gsl_vector_free(norm1);

  LLR = -LLR;
  return LLR / M_PI;
}

double AvPhaseLLR_RAAPTR(struct fitFuncParams *inParams) {

  struct RAAPTR_data *splParams = (struct RAAPTR_data *)inParams->splParams;

  struct cfunc_OUTPUT *output;
  output = (struct cfunc_OUTPUT *)malloc(1 * sizeof(struct cfunc_OUTPUT));
  output->c = (double *)malloc(4 * sizeof(double));
  output->v = (double *)malloc(9 * sizeof(double));

  // struct avPhase_param * avPhaseParam;
  // avPhaseParam = (struct avPhase_param *)malloc(1 * sizeof(struct
  // avPhase_param)); avPhaseParam->b = (double *)malloc(4 * sizeof(double));
  // avPhaseParam->norm = (double *)malloc(1 * sizeof(double));

  double tmp;
  double res;
  unsigned int lpr, i, j, j1, j2, j3, jj;
  //  size_t pp[6];   not used.
  const size_t kk = 1, stride = 1, nn = 6; // gsl_sort_largest_index
  // double src[6] = {1.0, 3.0, 6.0, 4.0, 5.0, 2.0};  //
  //  src = (double *)malloc(6 * sizeof(double));
  //  *(src+0)=1.0;
  //  *(src+1)=3.0;
  //  *(src+2)=6.0;
  //  *(src+3)=4.0;
  //  *(src+4)=5.0;
  //  *(src+5)=2.0;

  unsigned int Np = splParams->Np; // number of pulsars in PTA
  // unsigned int N = splParams->N;   // number of samples

  // transfer parameters from structure inParams
  // printf("MP5: Np = %d\n", Np);
  // printf("MP5: N = %d\n", N);
  double *alphaP, *deltaP;
  // alphaP = (double *)malloc(Np * sizeof(double));
  // deltaP = (double *)malloc(Np * sizeof(double));
  alphaP = splParams->alphaP;
  deltaP = splParams->deltaP;

  double theta;
  double alpha, delta, omega, phi0, Amp, iota, thetaN;
  double *phiI;

  unsigned int nDim = inParams->nDim; // dimension of fitness function
  // printf("MP5: nDim = %d\n", nDim);

  double b[5]; // quartic equation coefficients, closed form solution
  // b[0]=1.0;
  // b[1]=2.0;
  // b[2]=1.2;
  // b[3]=2.5;
  // b[4]=0.8;

  // avPhaseParam.b[0] = b[0];
  // printf("LLR_PSOav: avPhaseParam.b = %f\n",);
  // double z[8];  // 4 (real) + 4 (complex) solutions from
  // gsl_poly_complex_solve

  // unsigned int nr;  // number of EFFECTIVE roots (real number && abs(r)<1)

  double LLR = 0.0; // log likelihood ratio
  double C = 0.0;

  gsl_vector *skyLocSrc =
      gsl_vector_calloc(3); // sky location of source in Cartesian coordinate
  gsl_vector *skyLocPulsar =
      gsl_vector_calloc(3); // sky location of pulsars in Cartesian coordinate
  // gsl_vector * realCoord = gsl_vector_calloc(nDim);
  // gsl_vector_memcpy(realCoord,inParams->realCoord);

  alpha = gsl_vector_get(inParams->realCoord, 0);
  delta = gsl_vector_get(inParams->realCoord, 1);
  omega = gsl_vector_get(inParams->realCoord, 2);
  phi0 = gsl_vector_get(inParams->realCoord, 3);
  Amp = gsl_vector_get(inParams->realCoord, 4);
  Amp = pow(10, Amp); // physical amplitude
  gsl_vector_set(inParams->realCoord, 4, Amp);
  iota = gsl_vector_get(inParams->realCoord, 5);
  thetaN = gsl_vector_get(inParams->realCoord, 6);

  double **yr, **s, **sd,
      **Phi; // asign for jagged array, for different pulsars
  // sign value from splParams
  yr = splParams->yr;
  s = splParams->s;
  sd = splParams->sd;

  // Phi is not malloced yet
  Phi = (double **)malloc(Np * sizeof(double *));

  phiI = (double *)malloc(Np * sizeof(double));
  // printf("MP5: s[0][1] = %f\n", *(*(s+0)+1));
  // printf("MP5: s[1][2] = %f\n", *(*(s+1)+5));

  double *norm, *LRn, M, *norm2;
  norm = (double *)malloc(4 * sizeof(double));
  // norm1 = malloc(4 * sizeof(double));
  gsl_vector *norm1 = gsl_vector_alloc(4);
  norm2 = (double *)malloc(4 * sizeof(double));
  LRn = (double *)malloc(4 * sizeof(double));
  double sign[4][5] = {{1.0, 1.0, 1.0, 1.0, 1.0},
                       {1.0, -1.0, 1.0, -1.0, 1.0},
                       {1.0, -1.0, -1.0, 1.0, 1.0},
                       {1.0, 1.0, -1.0, -1.0, 1.0}};
  double intup[4] = {M_PI_2, M_PI, 3 * M_PI_2,
                     2 * M_PI};                     // upper limits of quadrants
  double intlow[4] = {0, M_PI_2, M_PI, 3 * M_PI_2}; // lower limits of quadrants

  gsl_vector_set(skyLocSrc, 0, cos(delta) * cos(alpha));
  gsl_vector_set(skyLocSrc, 1, cos(delta) * sin(alpha));
  gsl_vector_set(skyLocSrc, 2, sin(delta));
  for (i = 0; i < Np; i++) {
    Phi[i] = (double *)malloc(splParams->N[i] * sizeof(double));
    for (j = 0; j < splParams->N[i]; j++) {
      Phi[i][j] = yr[i][j] * omega;
      // printf("MP5: *(Phi+i) = %e, *(yr+i) = %e\n", *(Phi+i), *(yr+i));
    }
  }
  double bs[5], tmp0, NN, result, error;

  gsl_integration_workspace *w =
      gsl_integration_workspace_alloc(1000); // # double precision intervals

  // double bx[2];
  // bx[0]= 1.0;
  // bx[1]= 2.0;
  // double test;
  // test=9.999;
  //
  // norm[0] = 1.0;
  // norm[1] = 1.1;
  // norm[2] = 1.2;
  // norm[3] = 1.3;

  // gsl_function F;
  // struct avPhase_param avPhaseParam = {b[0],b[1],b[2],b[3],b[4],test};
  // //avPhaseParam->b = b;
  // //avPhaseParam->norm = &test;
  // //printf("LLR_PSOav: avPhaseParam.b = %f\n",(*avPhaseParam).b[1]);
  // //printf("LLR_PSOav: test = %f\n",*test);
  // //printf("LLR_PSOav: avPhaseParam.norm = %f\n",(*avPhaseParam).norm);
  // F.function = &f;
  // F.params = &avPhaseParam;
  // //F.norm = 2.22;
  // gsl_integration_qags(&F,0,1.57,0,1e-5,1000,w,&result,&error);
  // printf("LLR_PSOav: result = %.18f\n",result );
  // printf("AvPhaseLLR: Np = %d\n", Np);

  for (i = 0; i < Np; i++) {
    // printf("AvPhaseLLR: i = %d\n", i);
    M = 0.0;
    // printf("MP5: i = %d\n", i);
    // printf("alphaP = %f, deltaP = %f \n", *(alphaP+i), *(deltaP+i));
    gsl_vector_set(skyLocPulsar, 0, cos(deltaP[i]) * cos(alphaP[i]));
    gsl_vector_set(skyLocPulsar, 1, cos(deltaP[i]) * sin(alphaP[i]));
    gsl_vector_set(skyLocPulsar, 2, sin(deltaP[i]));

    gsl_blas_ddot(skyLocSrc, skyLocPulsar, &res);
    // printf("MP5: kp*k = %f\n", res);

    theta = acos(res);
    // printf("theta = %f\n", theta);

    cfunc_raaptr(splParams->N[i], alpha, delta, alphaP[i], deltaP[i], theta,
                 Amp, omega, iota, thetaN, phi0, Phi[i], s[i], sd[i], output);
    tmp = (*output).v[2] - 0.5 * ((*output).v[6] + (*output).v[8]);
    b[0] = tmp;
    tmp = (*output).v[0] - (*output).v[5];
    b[1] = tmp;
    tmp = (*output).v[1] - (*output).v[7];
    b[2] = tmp;
    tmp = -0.5 * (*output).v[4];
    b[3] = tmp;
    tmp = 0.5 * ((*output).v[6] - (*output).v[3]);
    b[4] = tmp;

    // quadrant I-IV
    for (j = 0; j < 4; j++) {

      for (j1 = 0; j1 < 5; j1++) {
        bs[j1] = b[j1] * sign[j][j1];
      }

      tmp0 = 0;
      for (jj = 1; jj < 5; jj++) {
        if (bs[jj] > 0) {
          tmp0 = tmp0 + bs[jj];
        }
      }
      norm[j] = tmp0;

      // printf("AvPhaseLLR: PSR i = %d, norm = %f\n",i,norm[j]);

      gsl_function F;
      struct avPhase_param avPhaseParam = {b[0], b[1], b[2],
                                           b[3], b[4], norm[j]};
      F.function = &f;
      F.params = &avPhaseParam;

      gsl_integration_qags(&F, intlow[j], intup[j], 0, 1e-5, 1000, w, &result,
                           &error);
      // printf("result = %.18f\n",result );

      LRn[j] = result;

      if (gsl_isinf(LRn[j])) {
        printf("AvPhaseLLR: Inf for PSR i = %d\n", i);
        printf("LRn at j = %d\n", j);
      } else if (gsl_isnan(LRn[j])) {
        printf("AvPhaseLLR: NAN for PSR i = %d\n", i);
        printf("LRn at j = %d\n", j);
      }

      // norm1[j] = norm[j] + log(LRn[j] + b[1]);
      gsl_vector_set(norm1, j,
                     norm[j] + log(LRn[j]) + b[0]); // b[0] = b[1] in Matlab
    }

    // printf("AvPhaseLLR: = %f\n",norm1[j]);

    // NN = gsl_max_dbl(norm1);
    NN = gsl_vector_max(norm1);

    for (j2 = 0; j2 < 4; j2++) {
      norm2[j2] = gsl_vector_get(norm1, j2) - NN;
    }

    for (j3 = 0; j3 < 4; j3++) {
      M += exp(norm2[j3]);
    }

    LLR = LLR + NN + log(M);
  }

  // free memory
  gsl_integration_workspace_free(w);
  free(output->c);
  free(output->v);
  free(output);
  for (i = 0; i < Np; i++)
    free(Phi[i]);

  free(Phi);
  free(phiI);
  free(norm);
  free(norm2);
  free(LRn);
  gsl_vector_free(skyLocSrc);
  gsl_vector_free(skyLocPulsar);
  gsl_vector_free(norm1);

  LLR = -LLR;
  return LLR / M_PI;
}

// double f (double x, void * params) {
//   double alpha0 = *(double *) params;
//   //double b = *(double *)bb;
//   // printf("b = %f\n", b);
//   //b[1]=bb[1];
//   //b[2]=bb[2];
//   //b[3]=bb[3];
//   //b[4]=bb[4];
//   double f = exp(alpha0) * cos(x); // + b[2] * sin(x) + b[3] * sin(2.0 * x) +
//               //  b[4] * pow(cos(x),2.0); // - norm;
//   return f;
// }

/*! Check coordinate validity after wrapping angles (if requested). Search range
for angles need not cover the whole sphere or circle. So, after wrapping angles,
check against actual search range.

\author Yan Wang
\author Soumya D. Mohanty (angle wrapping)

*/
size_t llrpsochkcoord(const size_t wrapAngles, const gsl_vector *rmin,
                      const gsl_vector *rangeVec, gsl_vector *xVec,
                      gsl_vector *realCoord) {

  // Default is to assume coordinates are valid
  size_t validPt = 1;

  /* Convert from standardized to real coordinates. Return real coordinates
           through realCoord. */
  s2rvector(xVec, rmin, rangeVec, realCoord);

  /* if angle wrapping not needed, just check standardized coordinates only */
  if (!wrapAngles) {
    validPt = chkstdsrchrng(xVec);
    return validPt;
  }

  // Indices of coordinates in realCoord
  /*
           coord #0: range of alpha is [0, 2*pi]
                 #1: range of delta is [pi/2,-pi/2]
           Coord #2: Frequency
                 #3: range of phi0 is [0,pi]
           Coord #4: Amplitude
                 #5: range of iota is [0,pi]
                 #6: range of thetaN is [0,pi]
        */
  size_t alphaIndx = 0;
  size_t deltaIndx = 1;
  size_t omgIndx = 2;
  size_t phi0Indx = 3;
  size_t ampIndx = 4;
  size_t iotaIndx = 5;
  size_t thetaNIndx = 6;

  /* First check non-angular parameters. Standardized coordinates are enough. */
  if (gsl_vector_get(xVec, omgIndx) < 0 || gsl_vector_get(xVec, omgIndx) > 1) {
    validPt = 0;
    return validPt;
  }
  if (gsl_vector_get(xVec, ampIndx) < 0 || gsl_vector_get(xVec, ampIndx) > 1) {
    validPt = 0;
    return validPt;
  }

  // Check alpha
  double twoPi = 2 * M_PI;
  double mn_alpha = gsl_vector_get(rmin, alphaIndx);
  double mx_alpha = gsl_vector_get(rangeVec, alphaIndx) + mn_alpha;
  double rng_alpha = mx_alpha - mn_alpha;
  if (mn_alpha < 0 || mx_alpha > twoPi) {
    printf("Warning: Check the limits on alpha\n");
  }
  double alpha = gsl_vector_get(realCoord, alphaIndx);
  if (alpha < mn_alpha || alpha > mx_alpha) {
    if (alpha < 0) {
      alpha = fmod(alpha, -twoPi);
      alpha = twoPi + alpha;
    } else if (alpha > twoPi) {
      alpha = fmod(alpha, twoPi);
    }
    // Reset real coordinate values
    gsl_vector_set(realCoord, alphaIndx, alpha);
    // Reset standardized coordinate value
    gsl_vector_set(xVec, alphaIndx, (alpha - mn_alpha) / rng_alpha);
    // check with boundaries in case search reagion does not cover the whole
    // sphere
    if (alpha < mn_alpha || alpha > mx_alpha) {
      validPt = 0;
      return validPt;
    }
  }

  // check delta
  double mn_delta = gsl_vector_get(rmin, deltaIndx);
  double mx_delta = gsl_vector_get(rangeVec, deltaIndx) + mn_delta;
  if (mn_delta < -M_PI_2 || mx_delta > M_PI_2) {
    printf("Warning: Check the limits on delta\n");
  }
  double delta = gsl_vector_get(realCoord, deltaIndx);
  if (delta < mn_delta || delta > mx_delta) {
    double polTheta = M_PI_2 - delta;
    if (polTheta < 0) {
      polTheta = fmod(polTheta, -twoPi);
      if (polTheta > -M_PI) {
        alpha = alpha + M_PI;
        polTheta = -polTheta;
        alpha = fmod(alpha, twoPi);
        // reset alpha
        gsl_vector_set(realCoord, alphaIndx, alpha);
        // Reset standardized value of alpha
        gsl_vector_set(xVec, alphaIndx, (alpha - mn_alpha) / rng_alpha);
      } else if (polTheta <= -M_PI) {
        polTheta = twoPi + polTheta;
      }
    } else if (polTheta > M_PI) {
      polTheta = fmod(polTheta, twoPi);
      if (polTheta < M_PI) {
        // Do nothing
      } else if (polTheta >= M_PI) {
        polTheta = twoPi - polTheta;
        alpha = alpha + M_PI;
        // Reset alpha
        gsl_vector_set(realCoord, alphaIndx, alpha);
        // Reset standardized value of alpha
        gsl_vector_set(xVec, alphaIndx, (alpha - mn_alpha) / rng_alpha);
      }
    }
    delta = M_PI_2 - polTheta;
    // Reset real coordinate
    gsl_vector_set(realCoord, deltaIndx, delta);
    // Reset standardized coordinate value
    gsl_vector_set(xVec, deltaIndx, (delta - mn_delta) / (mx_delta - mn_delta));
    // check with boundaries in case search reagion does not cover the whole
    // sphere
    if (alpha < mn_alpha || alpha > mx_alpha || delta < mn_delta ||
        delta > mx_delta) {
      validPt = 0;
      return validPt;
    }
  }

  // Check phi0
  double mn_phi0 = gsl_vector_get(rmin, phi0Indx);
  double mx_phi0 = gsl_vector_get(rangeVec, phi0Indx) + mn_phi0;
  if (mn_phi0 < 0 || mx_phi0 > M_PI) {
    printf("Warning: Check the limits for phi0\n");
  }
  double phi0 = gsl_vector_get(realCoord, phi0Indx);
  if (phi0 < mn_phi0 || phi0 > mx_phi0) {
    phi0 = wraphalfcircangle(phi0);
    // Reset real coordinate
    gsl_vector_set(realCoord, phi0Indx, phi0);
    // Reset standardized coordinate
    gsl_vector_set(xVec, phi0Indx, (phi0 - mn_phi0) / (mx_phi0 - mn_phi0));
    if (phi0 < mn_phi0 || phi0 > mx_phi0) {
      validPt = 0;
      return validPt;
    }
  }

  // Check iota
  double mn_iota = gsl_vector_get(rmin, iotaIndx);
  double mx_iota = gsl_vector_get(rangeVec, iotaIndx) + mn_iota;
  if (mn_iota < 0 || mx_iota > M_PI) {
    printf("Warning: Check the limits for iota\n");
    abort();
  }
  double iota = gsl_vector_get(realCoord, iotaIndx);
  if (iota < mn_iota || iota > mx_iota) {
    iota = wraphalfcircangle(iota);
    // Reset real coordinate
    gsl_vector_set(realCoord, iotaIndx, iota);
    // Reset standardized coordinate
    gsl_vector_set(xVec, iotaIndx, (iota - mn_iota) / (mx_iota - mn_iota));
    if (iota < mn_iota || iota > mx_iota) {
      validPt = 0;
      return validPt;
    }
  }

  // Check thetaN
  double mn_thetaN = gsl_vector_get(rmin, thetaNIndx);
  double mx_thetaN = gsl_vector_get(rangeVec, thetaNIndx) + mn_thetaN;
  if (mn_thetaN < 0 || mx_thetaN > M_PI) {
    printf("Warning: Check the limits for thetaN\n");
  }
  double thetaN = gsl_vector_get(realCoord, thetaNIndx);
  if (thetaN < mn_thetaN || thetaN > mx_thetaN) {
    thetaN = wraphalfcircangle(thetaN);
    // Reset real coordinate
    gsl_vector_set(realCoord, thetaNIndx, thetaN);
    // Reset standardized coordinate
    gsl_vector_set(xVec, thetaNIndx,
                   (thetaN - mn_thetaN) / (mx_thetaN - mn_thetaN));
    if (thetaN < mn_thetaN || thetaN > mx_thetaN) {
      validPt = 0;
      return validPt;
    }
  }

  // Wrap up: No violations were found, so validPt = 1.
  return validPt;
}

/*! Wrap angles into [0, pi].*/
double wraphalfcircangle(double phi) {
  double twoPi = 2 * M_PI;
  double twoPhi = 2 * phi;
  if (twoPhi < 0) {
    twoPhi = twoPi + fmod(twoPhi, -twoPi);
  } else if (twoPhi > twoPi) {
    twoPhi = fmod(twoPhi, twoPi);
  }
  return twoPhi / 2;
}

/*! log likelihood ratio function which maximizing over pulsar phase-MP,
 this function should be implemented as efficient as possible, since it
 would be called for a large amount of times, single source
\date 01/21/14, Yan Wang: adopt and modify from LogLikelihoodRatio()
\date 11/11/14, add fzero() after roots(), defaul TolX for fzero() is
eps=2.2204e-16 \date 11/13/14, compare fitness values at boundaries and at the
stationary points \date 12/19/14, cancel the changes made on 11/11/14, but still
check boundary \date Converted to C, Dec. 29, 2015.
*/
double LogLikelihoodRatioMP5(struct fitFuncParams *inParams) {

  struct llr_pso_params *splParams =
      (struct llr_pso_params *)inParams->splParams;

  struct cfunc_OUTPUT *output;
  output = (struct cfunc_OUTPUT *)malloc(1 * sizeof(struct cfunc_OUTPUT));
  output->c = (double *)malloc(4 * sizeof(double));
  output->v = (double *)malloc(9 * sizeof(double));

  struct lh_OUTPUT *lhoutput;
  lhoutput = (struct lh_OUTPUT *)malloc(1 * sizeof(struct lh_OUTPUT));
  lhoutput->phiI = (double *)malloc(1 * sizeof(double));
  lhoutput->lhI = (double *)malloc(1 * sizeof(double));

  double tmp;
  double res;
  unsigned int lpr, i, j;
  size_t pp[6];
  const size_t kk = 1, stride = 1, nn = 6; // gsl_sort_largest_index
  // double src[6] = {1.0, 3.0, 6.0, 4.0, 5.0, 2.0};  //
  //  src = (double *)malloc(6 * sizeof(double));
  //  *(src+0)=1.0;
  //  *(src+1)=3.0;
  //  *(src+2)=6.0;
  //  *(src+3)=4.0;
  //  *(src+4)=5.0;
  //  *(src+5)=2.0;

  unsigned int Np = splParams->Np; // number of pulsars in PTA
  unsigned int N = splParams->N;   // number of samples

  // transfer parameters from structure inParams
  // printf("MP5: Np = %d\n", Np);
  // printf("MP5: N = %d\n", N);
  double *yr;
  // yr = (double *)malloc(N * sizeof(double));
  double *sd;
  // sd = malloc(Np * sizeof(double));
  double *alphaP, *deltaP;
  // alphaP = (double *)malloc(Np * sizeof(double));
  // deltaP = (double *)malloc(Np * sizeof(double));
  sd = splParams->sd;
  alphaP = splParams->alphaP;
  deltaP = splParams->deltaP;
  yr = splParams->yr;

  double *Phi;
  Phi = (double *)malloc(N * sizeof(double));
  double theta;
  double alpha, delta, omega, phi0, Amp, iota, thetaN;
  double **phiItmp, **lh, *phiI;

  unsigned int nDim = inParams->nDim; // dimension of fitness function
  // printf("MP5: nDim = %d\n", nDim);

  double e[5]; // quartic equation coefficients, closed form solution
  double z[8]; // 4 (real) + 4 (complex) solutions from gsl_poly_complex_solve

  unsigned int nr; // number of EFFECTIVE roots (real number && abs(r)<1)

  double LLR = 0.0; // log likelihood ratio
  double **s;
  double C = 0.0;

  gsl_vector *skyLocSrc =
      gsl_vector_calloc(3); // sky location of source in Cartesian coordinate
  gsl_vector *skyLocPulsar =
      gsl_vector_calloc(3); // sky location of pulsars in Cartesian coordinate
  // gsl_vector * realCoord = gsl_vector_calloc(nDim);
  // gsl_vector_memcpy(realCoord,inParams->realCoord);

  alpha = gsl_vector_get(inParams->realCoord, 0);
  delta = gsl_vector_get(inParams->realCoord, 1);
  omega = gsl_vector_get(inParams->realCoord, 2);
  phi0 = gsl_vector_get(inParams->realCoord, 3);
  Amp = gsl_vector_get(inParams->realCoord, 4);
  Amp = pow(10, Amp); // physical amplitude
  gsl_vector_set(inParams->realCoord, 4, Amp);
  iota = gsl_vector_get(inParams->realCoord, 5);
  thetaN = gsl_vector_get(inParams->realCoord, 6);
  // printf("MP5: alpha = %f, delta = %f\n", alpha, delta);
  // printf("MP5: Amp = %f\n", Amp);

  // for (lpr=0; lpr<Np; lpr++){
  //   //sd = &inParams->sd[lpr];
  //   printf("MP5 test pointer: sd[%d] = %f, alphaP = %f, deltaP = %f\n",lpr,
  //   sd[lpr], alphaP[lpr], deltaP[lpr]);
  //   //sd++;
  // }

  // for (size_t lpr=0; lpr<Np; lpr++){
  //   s[lpr] = (*inParams).s[lpr];
  // }
  // sd = (*inParams).sd;
  // alphaP = M_PI/3.0;  //(*inParams).alphaP;
  // deltaP = M_PI/6.0; //(*inParams).deltaP;
  // kp = (*inParams).kp;
  // yr = inParams.yr;  //yr = (*inParams).yr;
  // xmaxmin = (*inParams).xmaxmin;

  // c = (double **)malloc(Np * sizeof(double));
  // for (size_t i = 0; i < Np; i++) {
  //    *(c+i) = malloc(N * sizeof(double));
  // }

  s = (double **)malloc(Np * sizeof(double));
  for (i = 0; i < Np; i++) {
    //*(s+i) = (double *)malloc(N * sizeof(double));  // not needed!
    s[i] = splParams->s[i];
  }

  // printf("MP5: s[0][1] = %f\n", *(*(s+0)+1));
  // printf("MP5: s[1][2] = %f\n", *(*(s+1)+5));
  phiI = malloc(Np * sizeof(double));
  phiItmp = (double **)malloc(Np * sizeof(double));
  lh = (double **)malloc(Np * sizeof(double));
  for (i = 0; i < Np; i++) {
    phiItmp[i] = (double *)malloc(6 * sizeof(double));
    lh[i] = (double *)malloc(6 * sizeof(double));
  }

  // gsl_matrix * c = gsl_matrix_calloc(Np,4);
  // gsl_matrix * e = gsl_matrix_calloc(Np,5);  // quartic equation
  // coefficients, closed form solution
  // gsl_matrix * lh = gsl_matrix_calloc(Np,6);  // likelihood for each pulsar,
  // two for boundaries gsl_matrix * phiItmp = gsl_matrix_calloc(Np,6);
  // gsl_vector * phiI = gsl_vector_alloc(6);
  //
  // gsl_vector * alphaP = gsl_vector_alloc(Np);
  // gsl_vector * deltaP = gsl_vector_alloc(Np);

  // calculate c for each pulsar
  // Phi=omega*yr;  // Phi=omega*t, N by 1 matrix
  // //c=zeros(Np,4);  // Eq. 25
  //
  // k[0]=cos(delta)*cos(alpha);
  // k[1]=cos(delta)*sin(alpha);
  // k[2]=sin(delta);
  // k[0]=cos(M_PI/4.0)*cos(M_PI/6.0);
  // k[1]=cos(M_PI/4.0)*sin(M_PI/6.0);
  // k[2]=sin(M_PI/4.0);

  // gsl_vector_set(k,0, cos(M_PI/4.0)*cos(M_PI/6.0));
  // gsl_vector_set(k,1, cos(M_PI/4.0)*sin(M_PI/6.0));
  // gsl_vector_set(k,2, sin(M_PI/4.0));

  gsl_vector_set(skyLocSrc, 0, cos(delta) * cos(alpha));
  gsl_vector_set(skyLocSrc, 1, cos(delta) * sin(alpha));
  gsl_vector_set(skyLocSrc, 2, sin(delta));

  // printf("MP5: omega = %f\n", omega);
  for (i = 0; i < N; i++) {
    Phi[i] = yr[i] * omega;
    // printf("MP5: *(Phi+i) = %e, *(yr+i) = %e\n", *(Phi+i), *(yr+i));
  }

  gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc(5);

  for (i = 0; i < Np; i++) {
    // printf("MP5: i = %d\n", i);
    // printf("alphaP = %f, deltaP = %f \n", *(alphaP+i), *(deltaP+i));
    gsl_vector_set(skyLocPulsar, 0, cos(deltaP[i]) * cos(alphaP[i]));
    gsl_vector_set(skyLocPulsar, 1, cos(deltaP[i]) * sin(alphaP[i]));
    gsl_vector_set(skyLocPulsar, 2, sin(deltaP[i]));

    gsl_blas_ddot(skyLocSrc, skyLocPulsar, &res);
    // printf("MP5: kp*k = %f\n", res);

    theta = acos(res);
    // printf("theta = %f\n", theta);

    cfunc(N, alpha, delta, alphaP[i], deltaP[i], theta, Amp, omega, iota,
          thetaN, phi0, Phi, s[i], (sd + i), output);

    tmp = 4 * (pow((*output).c[2], 2.0) + pow((*output).c[3], 2.0));
    // printf("MP5: tmp = %f\n", tmp);
    e[4] = tmp; // gsl_matrix_set(e,i,0,tmp);  // y^4
    // printf("MP5: e[4] = %f\n", e[4]);
    tmp =
        4 * ((*output).c[0] * (*output).c[2] + (*output).c[1] * (*output).c[3]);
    e[3] = tmp; // gsl_matrix_set(e,i,1,tmp);  // y^3
    // printf("MP5: e[3] = %f\n", e[3]);
    tmp = pow((*output).c[0], 2.0) + pow((*output).c[1], 2.0) -
          4.0 * (pow((*output).c[2], 2.0) + pow((*output).c[3], 2.0));
    e[2] = tmp; // gsl_matrix_set(e,i,2,tmp);  // y^2
    // printf("MP5: e[2] = %f\n", e[2]);
    tmp = -2.0 * (*output).c[1] * (*output).c[3] -
          4.0 * (*output).c[0] * (*output).c[2];
    e[1] = tmp; // gsl_matrix_set(e,i,3,tmp);  // y^1
    // printf("MP5: e[1] = %f\n", e[1]);
    tmp = pow((*output).c[3], 2.0) - pow((*output).c[0], 2.0);
    e[0] = tmp; // gsl_matrix_set(e,i,4,tmp);  // y^0
    // printf("MP5: e[0] = %f\n", e[0]);

    // calculate phiI(i) which maximize the likelihood ratio function
    // calculate roots of a quartic equation (there may be complex solutions)
    gsl_poly_complex_solve(e, 5, w, z);

    // for (j = 0; j < 4; j++) {
    //    printf ("MP5: z%d = %+.18f %+.18f\n", j, z[2*j], z[2*j+1]);
    // }

    nr = 0;

    for (j = 0; j < 4; j++) {

      if (z[2 * j + 1] == 0.0 &&
          fabs(z[2 * j]) <= 1.0) // real solution for y=cos(2\phi)
      {
        likelihood((z + 2 * j), output->v, lhoutput);
        //("MP5: z[2*j] = %f\n", *(z+2*j));
        phiItmp[i][j] = *(*lhoutput).phiI;
        lh[i][j] = *(*lhoutput).lhI;
        // printf("MP5: after i = %d\n", i);
        // printf("MP5: phiItmp[%d][%d]= %f, lh[i][j] = %f\n",
        // i,j,*(*(phiItmp+i)+j),*(*(lh+i)+j));

        nr = nr + 1;
      } else {
        // printf("MP5: There is a problem\n");
        phiItmp[i][j] = GSL_NAN; // NAN; //*(*lhoutput).phiI;
        lh[i][j] = GSL_NEGINF;   // INFINITY; //gsl_neginf; //*(*lhoutput).lhI;
                                 // printf("MP5: phiItmp[%d][%d] = %f\n", i, j,
                                 // phiItmp[i][j]);
      }
    }

    // likelihood((sd+1), output->c, lhoutput);

    // find fitness values at boundaries even with nr!=0
    tmp = -1.0;
    likelihood(&tmp, output->v, lhoutput);
    phiItmp[i][4] = M_PI / 2.0;
    lh[i][4] = *(*lhoutput).lhI;

    tmp = 1.0;
    likelihood(&tmp, output->v, lhoutput);
    phiItmp[i][5] = 0.0;
    lh[i][5] = *(*lhoutput).lhI;

    if (nr > 0) {

      // printf("Going to break here \n");
      // printf("%f, %f, %f, %f\n", *(src+0), *(src+1),*(src+4), *(src+5));
      // gsl_sort_largest_index(pp, kk, src, stride, nn);
      //*(*(lh+i)+0) = 5000.0;
      gsl_sort_largest_index(pp, kk, lh[i], stride, nn);
      //  printf("MP5: largest index = %zu  %zu  %zu %zu  %zu  %zu\n",
      //         pp[0],pp[1],pp[2],pp[3],pp[4],pp[5]);

      // C = *(src+pp[0]);
      // printf("MP5: C = %f\n", C);

      if (pp[0] == 4 || pp[0] == 5) {
        printf("Using fitness at boundary for PSR: %d. \n", i);
      }

      // phiI = *(*(phiItmp+i)+pp[0])
      //(*lhoutput).phiI[i] = phiItmp[i][pp[0]];
      phiI[i] = phiItmp[i][pp[0]];
      LLR = LLR + lh[i][pp[0]];
    }

    else if (nr == 0) {

      printf(
          "MP5: NO effective root (nr=0) for PSR: %d, use boundary fitness.\n",
          i);

      if (lh[i][4] > lh[i][5]) {
        C = lh[i][4];
        //(*lhoutput).phiI[i] = M_PI/2.0;
        phiI[i] = M_PI / 2.0;
      } else {
        C = lh[i][5];
        //(*lhoutput).phiI[i] = 0.0;
        phiI[i] = 0.0;
        // printf("There is a problem \n");
      }

      LLR = LLR + C;
    }

    // printf("phiI[%d] = %f\n", i, phiI[i]);
    // printf("phiI[0] = %f\n", phiI[0]);
  }

  /* Copy out additional special output from the fitness function
through the special parameter structure
*/
  size_t lpPhi;
  for (lpPhi = 0; lpPhi < Np; lpPhi++) {
    splParams->phiI[lpPhi] = phiI[lpPhi];
  }

  gsl_poly_complex_workspace_free(w);

  // printf("MP5: cfunc return = %f\n", (*output).fitVal);

  /* Free allocated memory */
  // free(yr);  // no need for deallocation of memo for them
  // free(alphaP);
  // free(deltaP);
  //  free(s);
  //  free(output);
  //  free(lhoutput);
  //  free(Phi);
  //  free(phiItmp);
  //  free(lh);

  free(output->c);
  free(output->v);
  free(output);
  free(lhoutput->phiI);
  free(lhoutput->lhI);
  free(lhoutput);
  free(Phi);
  free(s);
  for (i = 0; i < Np; i++) {
    free(phiItmp[i]);
    free(lh[i]);
  }
  free(phiItmp);
  free(lh);
  free(phiI);
  gsl_vector_free(skyLocSrc);
  gsl_vector_free(skyLocPulsar);

  LLR = -LLR; // return -log likelihood ratio, pso search for min
  return LLR;
}
// EOF: LogLikelihoodRatioMP5

/*! Core function called by \ref LogLikelihoodRatioMP5. */
void likelihood(double *y, double *inn, struct lh_OUTPUT *output) {
  unsigned int i;
  double x[2];
  double lh[2];
  double yy = *y;
  double tmp, lhI, phiI;

  // struct lh_OUTPUT * output = (struct lh_OUTPUT *)outputPtr;

  // printf("likelihood: y = %f\n", *y);

  //  there are two solution (+-x,y) for each y=cos(2\phi)
  x[0] = sqrt(1.0 - pow(yy, 2.0));
  x[1] = -sqrt(1.0 - pow(yy, 2.0));

  for (i = 0; i < 2; i++) {
    lh[i] = inn[0] * yy + inn[1] * x[i] + inn[2] -
            0.5 * (inn[3] * (yy * yy) + inn[6] * (x[i] * x[i]) +
                   2.0 * inn[4] * x[i] * yy + 2.0 * inn[5] * yy +
                   2.0 * inn[7] * x[i] + inn[8]); //  Eq 22
  }

  if (lh[0] > lh[1]) {
    tmp = x[0];
    lhI = lh[0];
  } else {
    tmp = x[1];
    lhI = lh[1];
  }

  // if y>=0
  if (tmp >= 0.0) {
    phiI = atan2(tmp, yy) / 2.0;
    //  y=cos(2*phiI)==X in Matlab atan2 func
  } else {
    phiI = (atan2(tmp, yy) + 2.0 * M_PI) / 2.0;
    // phiI = (rt_atan2d_snf(tmp, y) + 6.2831853071795862) / 2.0;
  }

  // return lhI, phiI;
  *(*output).phiI = phiI; // 2.46;
  *(*output).lhI = lhI;   // 3.69;

  // return output;
}
// EOF: likelihood

/*! \brief Function for quartic solver.
 */
void cfunc(unsigned int N, double alpha, double delta, double alphaP,
           double deltaP, double theta, double Amp, double omega, double iota,
           double thetaN, double phi0, double *Phi, double *s, double *sd,
           struct cfunc_OUTPUT *varargout) {
  unsigned int i;
  double alphatilde;
  double a, b, c, d, e, f;
  double Pp, Pc, Fp, Fc;
  double A, psi;
  double *x, *y, *z;
  double sx, sy, sz, xx, xy, xz, yy, yz, zz;

  // double res = 9.876;
  //(*varargout).fitVal = res;

  // printf("cfunc: alpha = %f, alphaP = %f, Phi[5] = %f, sd = %f\n", alpha,
  // alphaP, *(Phi+4), *sd);

  // for (i = 0; i < 6; i++) {
  //   printf("cfunc: s[%d] = %f\n", i, *(s+i));
  // }

  x = (double *)malloc(N * sizeof(double));
  y = (double *)malloc(N * sizeof(double));
  z = (double *)malloc(N * sizeof(double));

  // varargout->c = (double *)malloc(4 * sizeof(double));
  // varargout->v = (double *)malloc(9 * sizeof(double));

  // printf("cfunc: varargout.fitVal = %f\n", (*varargout).fitVal);
  alphatilde = alpha - alphaP;
  a = cos(deltaP);
  b = sin(deltaP);
  c = cos(alphatilde);
  d = sin(alphatilde);
  e = cos(delta);
  f = sin(delta);

  Pp = -pow(a, 2.0) * (1.0 - 2.0 * pow(c, 2.0) + pow(c, 2.0) * pow(e, 2.0)) +
       pow(b, 2.0) * pow(e, 2.0) -
       0.5 * sin(2.0 * deltaP) * c * sin(2.0 * delta);

  Pc = 2.0 * a * d * (a * c * f - b * e);
  // printf("cfunc: Pp = %f\n", Pp);
  Fp = Pp / (1.0 - cos(theta));
  Fc = Pc / (1.0 - cos(theta));
  // printf("cfunc: Fp = %f\n", Fp);
  // printf("cfunc: Amp = %f\n", Amp);
  A = 2.0 * Amp *
      sqrt(pow(1.0 + pow(cos(iota), 2.0), 2.0) *
               pow(Fp * cos(2.0 * thetaN) - Fc * sin(2.0 * thetaN), 2.0) +
           4.0 * pow(cos(iota), 2.0) *
               pow((Fp * sin(2.0 * thetaN) + Fc * cos(2.0 * thetaN)), 2.0));
  // printf("cfunc: A = %f\n", A);
  //  tmp=-2*cos(iota)/(1+cos(iota)^2)*(Fp*sin(2*thetaN)+Fc*cos(2*thetaN))/(Fp*cos(2*thetaN)-Fc*sin(2*thetaN));
  //  solve psi, atan or atan2 ?
  //  psi=atan(tmp);

  // psi = rt_atan2d_snf( -2.0 * cos(iota) * (Fp * sin(2.0 * thetaN) + Fc *
  // cos(2.0 * thetaN)),
  //       (1.0 + pow(cos(iota),2.0)) * (Fp * cos(2.0 * thetaN) - Fc * sin(2.0 *
  //       thetaN)) );
  //
  psi = atan2(-2.0 * cos(iota) *
                  (Fp * sin(2.0 * thetaN) + Fc * cos(2.0 * thetaN)),
              (1.0 + pow(cos(iota), 2.0)) *
                  (Fp * cos(2.0 * thetaN) - Fc * sin(2.0 * thetaN)));
  //("cfunc: psi = %f\n", psi);

  for (i = 0; i < N; i++) {
    x[i] = 0.5 * A * cos(psi + Phi[i]);
    y[i] = -0.5 * A * sin(psi + Phi[i]);
    z[i] = -0.5 * A * cos(2 * phi0 + psi + Phi[i]);
  }
  // printf("cfunc: x[0] = %f, x[1] = %f\n", *(x+0), *(x+1));

  // c are combination of inner weighted product of s with X,Y,Z
  // printf("cfunc: N= %d, s[5] = %f\n", N, *(s+5) );
  sx = InnProduct(N, s, x, *sd); //  scalar
  sy = InnProduct(N, s, y, *sd);
  sz = InnProduct(N, s, z, *sd);
  xx = InnProduct(N, x, x, *sd);
  xy = InnProduct(N, x, y, *sd);
  xz = InnProduct(N, x, z, *sd);
  yy = InnProduct(N, y, y, *sd);
  yz = InnProduct(N, y, z, *sd);
  zz = InnProduct(N, z, z, *sd);
  // printf("cfunc: sx= %f, sy = %f, zz =%f\n", sx,sy,zz);

  (*varargout).c[0] = -sx + xz;
  (*varargout).c[1] = sy - yz;
  (*varargout).c[2] = 0.5 * (xx - yy);
  (*varargout).c[3] = -xy;

  (*varargout).v[0] = sx;
  (*varargout).v[1] = sy;
  (*varargout).v[2] = sz;
  (*varargout).v[3] = xx;
  (*varargout).v[4] = xy;
  (*varargout).v[5] = xz;
  (*varargout).v[6] = yy;
  (*varargout).v[7] = yz;
  (*varargout).v[8] = zz;

  free(x);
  free(y);
  free(z);
}
// EOF cfunc

void cfunc_raaptr(unsigned int N, double alpha, double delta, double alphaP,
                  double deltaP, double theta, double Amp, double omega,
                  double iota, double thetaN, double phi0, double *Phi,
                  double *s, double *sd, struct cfunc_OUTPUT *varargout) {
  unsigned int i;
  double alphatilde;
  double a, b, c, d, e, f;
  double Pp, Pc, Fp, Fc;
  double A, psi;
  double *x, *y, *z;
  double sx, sy, sz, xx, xy, xz, yy, yz, zz;

  // double res = 9.876;
  //(*varargout).fitVal = res;

  // printf("cfunc: alpha = %f, alphaP = %f, Phi[5] = %f, sd = %f\n", alpha,
  // alphaP, *(Phi+4), *sd);

  // for (i = 0; i < 6; i++) {
  //   printf("cfunc: s[%d] = %f\n", i, *(s+i));
  // }

  x = (double *)malloc(N * sizeof(double));
  y = (double *)malloc(N * sizeof(double));
  z = (double *)malloc(N * sizeof(double));

  // varargout->c = (double *)malloc(4 * sizeof(double));
  // varargout->v = (double *)malloc(9 * sizeof(double));

  // printf("cfunc: varargout.fitVal = %f\n", (*varargout).fitVal);
  alphatilde = alpha - alphaP;
  a = cos(deltaP);
  b = sin(deltaP);
  c = cos(alphatilde);
  d = sin(alphatilde);
  e = cos(delta);
  f = sin(delta);

  Pp = -pow(a, 2.0) * (1.0 - 2.0 * pow(c, 2.0) + pow(c, 2.0) * pow(e, 2.0)) +
       pow(b, 2.0) * pow(e, 2.0) -
       0.5 * sin(2.0 * deltaP) * c * sin(2.0 * delta);

  Pc = 2.0 * a * d * (a * c * f - b * e);
  // printf("cfunc: Pp = %f\n", Pp);
  Fp = Pp / (1.0 - cos(theta));
  Fc = Pc / (1.0 - cos(theta));
  // printf("cfunc: Fp = %f\n", Fp);
  // printf("cfunc: Amp = %f\n", Amp);
  A = 2.0 * Amp *
      sqrt(pow(1.0 + pow(cos(iota), 2.0), 2.0) *
               pow(Fp * cos(2.0 * thetaN) - Fc * sin(2.0 * thetaN), 2.0) +
           4.0 * pow(cos(iota), 2.0) *
               pow((Fp * sin(2.0 * thetaN) + Fc * cos(2.0 * thetaN)), 2.0));
  // printf("cfunc: A = %f\n", A);
  //  tmp=-2*cos(iota)/(1+cos(iota)^2)*(Fp*sin(2*thetaN)+Fc*cos(2*thetaN))/(Fp*cos(2*thetaN)-Fc*sin(2*thetaN));
  //  solve psi, atan or atan2 ?
  //  psi=atan(tmp);

  // psi = rt_atan2d_snf( -2.0 * cos(iota) * (Fp * sin(2.0 * thetaN) + Fc *
  // cos(2.0 * thetaN)),
  //       (1.0 + pow(cos(iota),2.0)) * (Fp * cos(2.0 * thetaN) - Fc * sin(2.0 *
  //       thetaN)) );
  //
  psi = atan2(-2.0 * cos(iota) *
                  (Fp * sin(2.0 * thetaN) + Fc * cos(2.0 * thetaN)),
              (1.0 + pow(cos(iota), 2.0)) *
                  (Fp * cos(2.0 * thetaN) - Fc * sin(2.0 * thetaN)));
  //("cfunc: psi = %f\n", psi);

  for (i = 0; i < N; i++) {
    x[i] = 0.5 * A * cos(psi + Phi[i]);
    y[i] = -0.5 * A * sin(psi + Phi[i]);
    z[i] = -0.5 * A * cos(2 * phi0 + psi + Phi[i]);
  }
  // printf("cfunc: x[0] = %f, x[1] = %f\n", *(x+0), *(x+1));

  // c are combination of inner weighted product of s with X,Y,Z
  // printf("cfunc: N= %d, s[5] = %f\n", N, *(s+5) );
  sx = InnProduct_raaptr(N, s, x, sd); //  scalar
  sy = InnProduct_raaptr(N, s, y, sd);
  sz = InnProduct_raaptr(N, s, z, sd);
  xx = InnProduct_raaptr(N, x, x, sd);
  xy = InnProduct_raaptr(N, x, y, sd);
  xz = InnProduct_raaptr(N, x, z, sd);
  yy = InnProduct_raaptr(N, y, y, sd);
  yz = InnProduct_raaptr(N, y, z, sd);
  zz = InnProduct_raaptr(N, z, z, sd);
  // printf("cfunc: sx= %f, sy = %f, zz =%f\n", sx,sy,zz);

  (*varargout).c[0] = -sx + xz;
  (*varargout).c[1] = sy - yz;
  (*varargout).c[2] = 0.5 * (xx - yy);
  (*varargout).c[3] = -xy;

  (*varargout).v[0] = sx;
  (*varargout).v[1] = sy;
  (*varargout).v[2] = sz;
  (*varargout).v[3] = xx;
  (*varargout).v[4] = xy;
  (*varargout).v[5] = xz;
  (*varargout).v[6] = yy;
  (*varargout).v[7] = yz;
  (*varargout).v[8] = zz;

  free(x);
  free(y);
  free(z);
}

/*! \brief Inner product function for MaxPhase codes. */
double InnProduct(unsigned int N, double *X, double *Y, double sd) {
  unsigned int i;
  double result;
  double c = 0.0;

  // printf("InnProduct: sd = %f\n", sd);
  // printf("InnProduct: X[0] = %f\n", *(X+0));

  for (i = 0; i < N; i++) {
    c += (X[i] * Y[i]);
  }
  // printf("InnProduct: c = %f\n", c);

  result = c / (sd * sd);
  // printf("InnProduct: result = %f\n", result);

  return result;
}
// EOF: InnProduct

/*! \brief Inner product function for RAAPTR codes. */
double InnProduct_raaptr(unsigned int N, double *X, double *Y, double *sd) {
  unsigned int i;
  double result;
  double c = 0.0;

  // printf("InnProduct: sd = %f\n", sd);
  // printf("InnProduct: X[0] = %f\n", *(X+0));

  for (i = 0; i < N; i++) {
    c += (X[i] * Y[i]) / (sd[i] * sd[i]);
  }
  // printf("InnProduct: c = %f\n", c);

  // result = c / (sd[i] * sd[i]);
  result = c;
  // printf("InnProduct: result = %f\n", result);

  return result;
}

/*! Allocate special parameter structure */
struct llr_pso_params *llrparam_alloc(unsigned int N, unsigned int Np) {
  struct llr_pso_params *llp =
      (struct llr_pso_params *)malloc(sizeof(struct llr_pso_params));
  llp->Np = Np;
  llp->N = N;
  llp->sd = (double *)malloc(Np * sizeof(double));
  llp->alphaP = (double *)malloc(Np * sizeof(double));
  llp->deltaP = (double *)malloc(Np * sizeof(double));
  llp->phiI = (double *)malloc(Np * sizeof(double));
  llp->yr = (double *)malloc(N * sizeof(double));

  double **s = (double **)malloc(Np * sizeof(double *));
  size_t lpc1;
  for (lpc1 = 0; lpc1 < Np; lpc1++) {
    s[lpc1] = (double *)malloc(N * sizeof(double));
  }

  llp->s = s;

  return llp;
}

/*! Deallocate special parameter structure specific to LLR_PSO fitness function.
 */
void llrparam_free(struct llr_pso_params *llp) {

  size_t lpc;
  size_t Np = llp->Np;

  for (lpc = 0; lpc < Np; lpc++) {
    free(llp->s[lpc]);
  }
  free(llp->s);
  free(llp->sd);
  free(llp->alphaP);
  free(llp->deltaP);
  free(llp->yr);
  free(llp->phiI);
  free(llp);
}
