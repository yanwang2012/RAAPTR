#include "maxphase.h"
#include <gsl/gsl_vector.h>
#include <stdio.h>
int main(){
	size_t lpc;
	double testVec[5] = {.1,.2,.6,.3,-.7};
	double min = -.2;
	double max = 0.58;
	gsl_vector *xVec = gsl_vector_alloc(5);
	for (lpc = 0; lpc < 5; lpc++)
		gsl_vector_set(xVec,lpc,testVec[lpc]);
	limitVecComponent(xVec,min,max);
	for (lpc = 0; lpc< 5; lpc++){
		printf("%f ", testVec[lpc]);
		printf("%f ", gsl_vector_get(xVec,lpc));
		printf(" min=%f max=%f\n", min, max);
	}
}