#include "subtract.h"
//#include "maxphase.h"
//#include "LLR_Mp_Av.h"
//#include "ptapso.h"
#include "hdf5.h"
#include "gslhdf5_io.h"
#include <gsl/gsl_vector.h>
#include <sys/stat.h>
#include <stdio.h>
#include <string.h>
/*! \file
\brief Subtract estimated source out of simulated data.

\author Yiqian Qian.
*/

struct estSrcParams * file2Srcparam(char *outputFileName){
	herr_t status;
	hid_t SrcPar = H5Fopen(outputFileName, H5F_ACC_RDONLY, H5P_DEFAULT);
	if (SrcPar < 0)
	{
		printf("Error opening file %s \n", outputFileName);
	}

	gsl_vector * bestRealLoc = hdf52gslvector(SrcPar, "bestRealLoc");

	FILE * fPtr;
	fPtr = fopen("bestRealLoc.txt","w");
	gsl_vector_fprintf(fPtr,bestRealLoc,"%f");
	fclose(fPtr);

	size_t nDim = bestRealLoc->size;
	printf("nDim is: %zu \n",nDim);

	struct estSrcParams * srcp = srcp_alloc(nDim);

	srcp->alpha = gsl_vector_get(bestRealLoc,0);
	printf("alpha is: %lf \n",srcp->alpha);

	srcp->delta = gsl_vector_get(bestRealLoc,1);
	srcp->omega = gsl_vector_get(bestRealLoc,2);
	srcp->phi0 = gsl_vector_get(bestRealLoc,3);
	srcp->Amp = gsl_vector_get(bestRealLoc,4);
	srcp->iota = gsl_vector_get(bestRealLoc,5);
	srcp->thetaN = gsl_vector_get(bestRealLoc,6);
	printf("thetaN is: %lf \n",srcp->thetaN);

	FILE * f;
	f = fopen("psrPhase.txt","w");
	// srcp->psrPhase = gsl_vector_calloc(nDim-7); // initialize dimensionality of pulsar phase.
	printf("faaaak this");
	gsl_vector_fprintf(f,srcp->psrPhase,"%f");
	fclose(f);

	size_t lpc1;
	for(lpc1 = 0; lpc1 < nDim-7; lpc1++){
		gsl_vector_set(srcp->psrPhase,lpc1,gsl_vector_get(bestRealLoc,lpc1+7));
	}
	status = H5Fclose(SrcPar);
	if(status < 0){
		fprintf(stdout,"Error closing file %s \n",outputFileName);
	}

	gsl_vector_free(bestRealLoc);

	return srcp;
}

void srcpara_free(struct estSrcParams *srcp){
	gsl_vector_free(srcp->psrPhase);

}

void printParam(struct estSrcParams *srcp){
	printf("alpha: %lf \n"
			"delta: %lf \n"
			"omega: %lf \n"
			"phi0: %lf \n"
			"Amp: %lf \n"
			"iota: %lf \n"
			"thetaN: %lf \n",
			srcp->alpha, srcp->delta,srcp->omega,srcp->phi0,
			srcp->Amp,srcp->iota,srcp->thetaN);
}

struct estSrcParams * srcp_alloc(size_t nDim){
	struct estSrcParams * srcp = (struct estSrcParams *)malloc(sizeof(struct estSrcParams));
	srcp->psrPhase = gsl_vector_calloc(nDim-7);
}