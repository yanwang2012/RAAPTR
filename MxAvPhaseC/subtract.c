#include "subtract.h"
#include "maxphase.h"
#include "LLR_Mp_Av.h"
#include "hdf5.h"
#include "gslhdf5_io.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <sys/stat.h>
#include <stdio.h>
#include <string.h>
/*! \file
\brief Subtract estimated source out of simulated data.

\author Yiqian Qian.
*/

struct estSrcParams *file2Srcparam(char *outputFileName)
{
	herr_t status;
	hid_t SrcPar = H5Fopen(outputFileName, H5F_ACC_RDONLY, H5P_DEFAULT);
	if (SrcPar < 0)
	{
		printf("Error opening file %s \n", outputFileName);
	}

	gsl_vector *bestRealLoc = hdf52gslvector(SrcPar, "bestRealLoc");
	
	FILE * fPtr;
	fPtr = fopen("bestRealLoc.txt","w");
	gsl_vector_fprintf(fPtr,bestRealLoc,"%f");
	fclose(fPtr);
	
	size_t nDim = bestRealLoc->size;
	printf("nDim is: %zu \n", nDim);

	struct estSrcParams *srcp = srcp_alloc(nDim - 7);

	srcp->alpha = gsl_vector_get(bestRealLoc, 0);
	// printf("alpha is: %lf \n",srcp->alpha);
	srcp->delta = gsl_vector_get(bestRealLoc, 1);
	srcp->omega = gsl_vector_get(bestRealLoc, 2);
	srcp->phi0 = gsl_vector_get(bestRealLoc, 3);
	srcp->Amp = gsl_vector_get(bestRealLoc, 4);
	srcp->iota = gsl_vector_get(bestRealLoc, 5);
	srcp->thetaN = gsl_vector_get(bestRealLoc, 6);
	// printf("thetaN is: %lf \n",srcp->thetaN);

	size_t lpc1;
	for (lpc1 = 0; lpc1 < nDim - 7; lpc1++)
	{
		gsl_vector_set(srcp->psrPhase, lpc1, gsl_vector_get(bestRealLoc, lpc1 + 7));
	}
	
	FILE * Fpsr;
	Fpsr = fopen("psrPhase.txt","w");
	gsl_vector_fprintf(Fpsr,srcp->psrPhase,"%f");
	fclose(Fpsr);
	
	status = H5Fclose(SrcPar);
	if (status < 0)
	{
		fprintf(stdout, "Error closing file %s \n", outputFileName);
	}

	gsl_vector_free(bestRealLoc);

	return srcp;
}

gsl_matrix * timingResiduals(struct estSrcParams *srcp, struct llr_pso_params *splParams)
{
	//struct llr_pso_params *splParams = (struct llr_pso_params *)inParams->splParams;
	size_t Np;
	size_t N;
	/* estimated source parameters. */
	double alpha, delta, omega, phi0, Amp, iota, thetaN;

	gsl_vector *skyLocSrc = gsl_vector_calloc((size_t)3);
	gsl_vector *skyLocPsr = gsl_vector_calloc((size_t)3);
	/* pulsar parameters. */
	double *alphaP, *deltaP, *yr;
	double theta, res;


	alpha = srcp->alpha;
	delta = srcp->delta;
	omega = srcp->omega;
	phi0 = srcp->phi0;
	Amp = srcp->Amp;
	iota = srcp->iota;
	thetaN = srcp->thetaN;

	alphaP = splParams->alphaP;
	deltaP = splParams->deltaP;
	yr = splParams->yr;
	Np = (size_t)splParams->Np;
	N = (size_t)splParams->N;

	gsl_matrix * timResiduals = gsl_matrix_calloc(Np,N);
	gsl_vector *psrPhase = gsl_vector_calloc(Np);
	gsl_matrix * tmp = gsl_matrix_calloc(N,(size_t)1);

	size_t i,j,k;

	for (i == 0; i < Np; i++){
		gsl_vector_set(psrPhase, i, gsl_vector_get(srcp->psrPhase, i));
	}

	gsl_vector_set(skyLocSrc, 0, cos(delta) * cos(alpha));
	gsl_vector_set(skyLocSrc, 1, cos(delta) * sin(alpha));
	gsl_vector_set(skyLocSrc, 2, sin(delta));

	for (j == 0; j < Np; j++){
		gsl_vector_set(skyLocPsr, 0, cos(deltaP[j]) * cos(alphaP[j]));
		gsl_vector_set(skyLocPsr, 1, cos(deltaP[j]) * sin(alphaP[j]));
		gsl_vector_set(skyLocPsr, 2, sin(deltaP[j]));

		gsl_blas_ddot(skyLocSrc, skyLocPsr, &res);
		theta = acos(res);
		printf("timing_theta is: %e\n",theta);
		tmp = FullResiduals(srcp, alphaP[j], deltaP[j],gsl_vector_get(psrPhase,j),theta,yr);

		for (k = 0; k < N; k++){
			gsl_matrix_set(timResiduals,j,k,gsl_matrix_get(tmp,k,0));
		}
	}

	gsl_matrix_free(tmp);
	gsl_vector_free(psrPhase);
	gsl_vector_free(skyLocSrc);
	gsl_vector_free(skyLocPsr);
	
	return timResiduals;
}

gsl_matrix * FullResiduals(struct estSrcParams * srcp, double alphaP, double deltaP,double phiI, double theta, double *yr)
{


	size_t N = sizeof(yr)/sizeof(yr[0]);
	printf("Size of yr is: %zu\n", N);
	double alpha, delta, omega, phi0, Amp, iota, thetaN;
	alpha = srcp->alpha;
	delta = srcp->delta;
	omega = srcp->omega;
	phi0 = srcp->phi0;
	Amp = srcp->Amp;
	iota = srcp->iota;
	thetaN = srcp->thetaN;

	size_t i, j, k;
	gsl_matrix *C = gsl_matrix_calloc(8, 1);
	gsl_matrix *A = gsl_matrix_calloc(N, 8);
	gsl_matrix *r = gsl_matrix_calloc(N,1);
	double alphatilde;
	double a, b, c, d, e, f;
	double Pp, Pc, Fp, Fc, FpC, FpS, FcC, FcS;
	double CosIota, TwoThetaN, tmp1, tmp2, tmpC, tmpC2, tmpS, tmpS2;
	double *omegaT;
	omegaT = (double *)malloc(N * sizeof(double));

	alphatilde = alpha - alphaP;
	a = cos(deltaP);
	b = sin(deltaP);
	c = cos(alphatilde);
	d = sin(alphatilde);
	e = cos(delta);
	f = sin(delta);

	Pp = -pow(a, 2.0) * (1.0 - 2.0 * pow(c, 2.0) + pow(c, 2.0) * pow(e, 2.0)) +
		 pow(b, 2.0) * pow(e, 2.0) - 0.5 * sin(2.0 * deltaP) * c * sin(2.0 * delta);

	Pc = 2.0 * a * d * (a * c * f - b * e);
	//printf("cfunc: Pp = %f\n", Pp);
	Fp = Pp / (1.0 - cos(theta));
	Fc = Pc / (1.0 - cos(theta));

	CosIota = cos(iota);
	TwoThetaN = 2 * thetaN;
	tmp1 = Amp * (1 + pow(CosIota,2.0));
	tmp2 = Amp * 2 * CosIota;
	tmpC = cos(TwoThetaN);
	tmpS = sin(TwoThetaN);

	gsl_matrix_set(C,0,0, -tmp1 * tmpC);
	gsl_matrix_set(C,1,0, -tmp2 * tmpS);
	gsl_matrix_set(C,2,0, -gsl_matrix_get(C,0,0));
	gsl_matrix_set(C,3,0, gsl_matrix_get(C,1,0));
	gsl_matrix_set(C,4,0, tmp1 * tmpS);
	gsl_matrix_set(C,5,0, -tmp2 * tmpC);
	gsl_matrix_set(C,6,0, -gsl_matrix_get(C,4,0));
	gsl_matrix_set(C,7,0, gsl_matrix_get(C,5,0));

	tmpC2 = cos(2 * phi0) - cos(2 * phiI);
	tmpS2 = sin(2 * phi0) - sin(2 * phiI);
	FpC = Fp * tmpC2;
	FpS = Fp * tmpS2;
	FcC = Fc * tmpC2;
	FcS = Fc * tmpS2;

	for(i == 0; i < N; i++){
		omegaT[i] = omega * yr[i];
	}

	for(j == 0; j < 8; j++){
		for(k = 0; k < N; k++){
			gsl_matrix_set(A,k,j, FpC * cos(omegaT[k]));
		}
	}

	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, C, 0, r);

	FILE * fr;
	fr = fopen("r.txt","w");
	gsl_matrix_fprintf(fr,r,"%g");
	fclose(fr);

	gsl_matrix_free(A);
	gsl_matrix_free(C);

	return r;
}

void srcpara_free(struct estSrcParams *srcp)
{
	gsl_vector_free(srcp->psrPhase);
}

void printParam(struct estSrcParams *srcp)
{
	printf("alpha: %lf \n"
		   "delta: %lf \n"
		   "omega: %lf \n"
		   "phi0: %lf \n"
		   "Amp: %lf \n"
		   "iota: %lf \n"
		   "thetaN: %lf \n",
		   srcp->alpha, srcp->delta, srcp->omega, srcp->phi0,
		   srcp->Amp, srcp->iota, srcp->thetaN);
}

struct estSrcParams *srcp_alloc(size_t nDim)
{
	struct estSrcParams *srcp = (struct estSrcParams *)malloc(sizeof(struct estSrcParams));
	srcp->psrPhase = gsl_vector_calloc(nDim);
	return srcp;
}