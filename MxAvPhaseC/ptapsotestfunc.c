/*! \file ptapsotestfunc.c
\brief Example of a fitness function.
 */
#include <stdio.h>
#include "ptapsotestfunc.h"
#include "maxphase.h"
#include <gsl/gsl_math.h>
/*
fitFuncVal = ptapsotestfunc_gsl(xVec,P)
A benchmark test function for PTAPSO
that computes the Rastrigin fitness function for
each row of xVec.  The fitness value is returned in fitFuncVal.
xVec is standardized, that is 0<=xVec(i,j)<=1. 
The values used to convert xVec(i,j)
internally before computing fitness are given in P.rmin and
P.rangeVec: 
xVec(j) -> xVec(j)*rangevec(j)+rmin(j).
fitFuncVal = infty if the point xVec falls
outside the hypercube defined by 0<=xVec(j)<=1.
The real coordinates are returned in P.realCoord. 

Soumya D. Mohanty, Jan 2016
- Derived from ptapsotestfunc.c. Converts to gsl_vector inputs and 
  uses the interface needed by GSL multi-dimensional local minimization routines.
*/



double ptapsotestfunc(gsl_vector *xVec, void  *inParamsPointer){
	
	unsigned int validPt;
    unsigned int lpc;
	//! [Cast fit func params]
	struct fitFuncParams *inParams = (struct fitFuncParams *)inParamsPointer;
	//! [Cast fit func params]
	unsigned int ncols = inParams->nDim;
	double rangeVec, rmin, x;
	double fitFuncVal;
	//! [Cast special params]
	/* Shows how to retrieve special parameters (dummy ones are supplied for this particular
	fitness function).
	*/
	struct ptapsotestfunc_params *splParams = (struct ptapsotestfunc_params *)inParams->splParams;
		
	/* This fitness function knows what fields are given in the special parameters struct */
    int dummy = splParams->dummyParam;
	//printf("PTAPSOTESTFUNC called with special parameters %d\n",dummy);
	//! [Cast special params]
	
	s2rvector(xVec,inParams->rmin,inParams->rangeVec,inParams->realCoord);
	
	validPt = chkstdsrchrng(xVec);
	
	if (validPt){
		inParams->fitEvalFlag = 1;
		fitFuncVal = 0;
		for (lpc = 0; lpc < ncols; lpc++){
			x = gsl_vector_get(inParams->realCoord,lpc);
			fitFuncVal = fitFuncVal + 
			               	    gsl_pow_int(x,2) - 
						   		 10.0*cos(2*M_PI*x) +
						   		  10;
		}
    }
	else{
		fitFuncVal=GSL_POSINF;
		inParams->fitEvalFlag = 0;
	}
   return fitFuncVal;
}