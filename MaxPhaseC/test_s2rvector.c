#include "maxphase.h"
#include <gsl/gsl_vector.h>
#include <stdio.h>
/*
Test function for s2rvector.c.

gcc -c test_s2rvector.c -I/usr/local/include
gcc test_s2rvector.o maxphaseutils.o -L/usr/local/lib -lgsl -lgslcblas -lm -o test_s2rvector.out
*/
int main(){

    /* Test matrix */
    int nrows = 5, ncols = 3;
    int lpr, lpc;
    double stdCoord[5][3]={{.1,.2,.5},
                          {1.0,.1,.2},
 					      {.1,.2,0.0},
 					      {2.0,.1,.2},
 					      {.1,.2,2.0}};
  	/* Range information to use for unstandardizing. */
  	double rmin[3] = {-5.0,-4.0,2.0};
  	double rmax[3] = {5.0,2.0,4.0};
						  
	gsl_vector *xVec = gsl_vector_alloc(ncols);
	gsl_vector *rcOut = gsl_vector_alloc(ncols);
	gsl_vector *rminVec = gsl_vector_alloc(ncols);
	gsl_vector *rngVec = gsl_vector_alloc(ncols);
	for(lpc = 0; lpc < ncols; lpc++){
		gsl_vector_set(rminVec, lpc, rmin[lpc]);
		gsl_vector_set(rngVec, lpc, rmax[lpc] - rmin[lpc]);
	}
	
	/* Load the pointer arrays */
	for (lpr=0; lpr < nrows; lpr++){
		for (lpc = 0; lpc < ncols; lpc++){
			gsl_vector_set(xVec, lpc, stdCoord[lpr][lpc]);			
		}
		s2rvector(xVec,rminVec,rngVec,rcOut);	
		for (lpc = 0; lpc < ncols; lpc++){
			printf("%f ", gsl_vector_get(xVec, lpc));			
		}
		printf("\n");
		for (lpc = 0; lpc < ncols; lpc++){
			printf("%f ",gsl_vector_get(rcOut, lpc));			
		}
		printf("\n");			
	}
	
	//Wrap up
	gsl_vector_free(xVec);
	gsl_vector_free(rcOut);
	gsl_vector_free(rminVec);
	gsl_vector_free(rngVec);
}