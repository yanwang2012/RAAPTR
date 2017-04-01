#include "ptapso_prep.h"
#include "maxphase.h"
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>

size_t nDim = 3;
double testLocs[3][3]={ {.1,.2,.3},
                 {1.5,.2,.3},
		         {0,.1,.2}
                };
double fitValVec[3];

void ptapso(double (*fitfunc)(const gsl_vector *, void *), void *fitfuncParams, 
            struct psoParamStruct psoParams, 
			struct returnData *psoResults){

	gsl_vector *xVec = gsl_vector_alloc(nDim);
	unsigned int lpr, lpc;
	unsigned int nrows = nDim;
	double x;
	
	for (lpr = 0; lpr < nrows; lpr++){
		for (lpc = 0; lpc < nDim; lpc++){
			gsl_vector_set(xVec,lpc,testLocs[lpr][lpc]);
			printf("%f, ",testLocs[lpr][lpc]);
		}
		printf("\n");
        fitValVec[lpr]=(*fitfunc)(xVec,fitfuncParams);
		printf("%f\n",fitValVec[lpr]);
	}

	
	gsl_rng *r = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(r,2571971);
	for (lpr = 0; lpr < nrows; lpr++){
		for (lpc = 0; lpc < nDim; lpc++){
			x = gsl_rng_uniform(r);
			gsl_vector_set(xVec,lpc,x);
			printf("%f, ",x);
		}
		printf("\n");
		fitValVec[lpr] = (*fitfunc)(xVec,fitfuncParams);
		printf("%f\n",fitValVec[lpr]);
	}

	
	/* Nelder Mead minimization */
	/* 1. Specify the function to minimize and the starting point.*/
	gsl_multimin_function func2minimz;
	func2minimz.n = nDim;
	func2minimz.f = fitfunc;
	func2minimz.params = fitfuncParams;
	gsl_vector *minimzStrtPt = gsl_vector_alloc(nDim);
	   /* Here we choose a random starting point */
	for(lpc = 0; lpc < nDim; lpc++){
		x = gsl_rng_uniform(r);
		gsl_vector_set(minimzStrtPt,lpc,x);
	}
	/* 2. Select and initialize the minimization algorithm.
	      Note the similarity of the interface with the randon number
	      generation.
	*/
	gsl_multimin_fminimizer *minimzrState = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, nDim);
	  /* Initial step size */
	gsl_vector *stpSz = gsl_vector_alloc(nDim);
	gsl_vector_set_all(stpSz,.1);
	  /* Tell the minimizer what the function to minimize is, the starting point
	     of the search and the initial step size. */
	gsl_multimin_fminimizer_set(minimzrState,&func2minimz,minimzStrtPt,stpSz);
	/* 3. Run a specified number of iterations (which is what we need in the PTAPSO case). */
	int status, conv_stat;
	double size, fitVal;
	for (lpc = 0; lpc < 20; lpc++){
		status = gsl_multimin_fminimizer_iterate(minimzrState);
		/* A non-zero value of status indicates failure */
		if (status)
			break;
		printf("%d, ",status);
		size = gsl_multimin_fminimizer_size(minimzrState);
		printf("%f, ",size);
		fitVal = gsl_multimin_fminimizer_minimum(minimzrState);
		printf("%f\n",fitVal);
		
		conv_stat = gsl_multimin_test_size (size, 1e-2);
		if (conv_stat == GSL_SUCCESS)
		{
			printf ("converged to minimum at (standardized coordinates)\n");
			for (lpr = 0; lpr < nDim; lpr++){
				printf("%f, ", gsl_vector_get(minimzrState->x,lpr));
			}
			printf("\n");
			/* This shows how to retrieve real coordinates. */
			printf ("converged to minimum at (real coordinates)\n");
			for (lpr = 0; lpr < nDim; lpr++){
				printf("%f, ", gsl_vector_get(((struct fitFuncParams *)fitfuncParams)->realCoord,lpr));
			}
			printf("\n-------\n");
			
		}
	}
	
	/* Free function minimizer state */
	gsl_multimin_fminimizer_free(minimzrState);
	/* Free random number generator state */
    gsl_rng_free(r);
	/* Deallocate vectors */
	gsl_vector_free(xVec);
	gsl_vector_free(minimzStrtPt);
}