#include "maxphase.h"
#include <gsl/gsl_vector.h>
#include <stdio.h>
int main(){
    size_t nDim,lp;
	nDim = 5;
	gsl_vector *testVec = gsl_vector_calloc(nDim);
	struct fitFuncParams *ffp = ffparam_alloc(nDim);
	for (lp = 0; lp < nDim; lp++){
		gsl_vector_set(ffp->rmin,lp,-5);
		gsl_vector_set(ffp->rangeVec, lp, 10);
		gsl_vector_set(testVec,lp,lp*.1);
	}
	s2rvector(testVec,ffp->realCoord,ffp->rmin,ffp->rangeVec);
	for (lp = 0; lp < nDim; lp++){
		printf("%f ",gsl_vector_get(ffp->realCoord,lp));
	}
	printf("\n");
	
	ffparam_free(ffp);
}