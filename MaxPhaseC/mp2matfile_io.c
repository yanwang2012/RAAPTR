#include "maxphase.h"
#include "ptapso.h"
#include "mat.h"
#include "matrix.h"
#include <stdio.h>

/*! \file
\brief Functions for doing I/O with .mat files.

\author Soumya D. Mohanty
*/

/*! 
Read data from .mat file and fill out the fitness function 
parameter strucutre. See \ref fitFuncParams.
*/
struct fitFuncParams * file2ffparam(char *srchParamsFile){
	
	/* Read xmaxmin from srchParamsFile file */
	MATFile *srchPar = matOpen(srchParamsFile, "r");
	if (srchPar == NULL){
		printf("Error reading parameter file %s\n", srchParamsFile);
		abort();
	}
	
	mxArray *xmaxmin;
	xmaxmin = matGetVariable(srchPar,"xmaxmin");
	/* xmaxmin should be a 2 dimensional array */
	if (mxGetNumberOfDimensions(xmaxmin) !=2){
		printf("Parameter specification is not understood\n");
		mxDestroyArray(xmaxmin);
        abort();
	}
	/* Search Space dimensionality */
	size_t nDim = mxGetM(xmaxmin);
	
	/* transfer xmaxmin to fitness function parameter struct */
	double *xmaxminData = mxGetPr(xmaxmin);
	struct fitFuncParams *ffp = ffparam_alloc(nDim);
    size_t lpc1;
	for(lpc1 = 0; lpc1 < nDim; lpc1++){
		gsl_vector_set(ffp->rmin,lpc1,xmaxminData[lpc1+nDim]);
		gsl_vector_set(ffp->rangeVec,lpc1,xmaxminData[lpc1]-xmaxminData[lpc1+nDim]);
	}
	/* Close file */
	int stat;
	if ((stat =  matClose(srchPar)))
		printf("Error closing file %s %d\n", srchParamsFile, stat);
	
	mxDestroyArray(xmaxmin);
	
	return ffp;
}

