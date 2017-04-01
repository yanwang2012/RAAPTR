/* test function for ptapso.c */

#include "ptapso_debug.h"
#include "maxphase.h"
/* The header file tells us what fitness function we are calling
and what the parameter structure for this function is.
*/
#include "ptapsotestfunc.h"
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>


int main(){
	/* Number of search dimensions */
	unsigned int nDim = 3, lpc;
	double rmin[3] = {-5,-5,-5};
	double rmax[3] = {5, 5, 5};
	double rangeVec[3];
	
    /* Set up fitness function parameter struct. A subset of the fields
	   is universal to all fitness function.
	 */
	struct fitFuncParams inParams;
	inParams.rmin = gsl_vector_alloc(nDim);
	inParams.rangeVec = gsl_vector_alloc(nDim);
	inParams.realCoord = gsl_vector_alloc(nDim);
	inParams.nDim = nDim;
	/* Set special parameters (which are dummy for this fitness function)*/
	struct ptapsotestfunc_params splParams;
	splParams.dummyParam = 78520;
	inParams.splParams = &splParams;
	
	/* Error handling off */
	gsl_error_handler_t *old_handler = gsl_set_error_handler_off ();
	
	/* Initialize random number generator */
	gsl_rng *rngGen = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(rngGen,2571971);
	
	for (lpc = 0; lpc < nDim; lpc++){
		rangeVec[lpc]=rmax[lpc]-rmin[lpc];
		gsl_vector_set(inParams.rmin,lpc,rmin[lpc]);
		gsl_vector_set(inParams.rangeVec,lpc,rangeVec[lpc]);		
	}
	/* Set up pointer to fitness function. Use the prototype
	declaration given in the header file for the fitness function. */
	double (*fitfunc)(const gsl_vector *, void *) = ptapsotestfunc;
	
	/* Set up storage for output from ptapso. */
	struct returnData psoResults;
	gsl_vector *bestLocation = gsl_vector_alloc(nDim);
	psoResults.bestLocation = bestLocation;
	
	/* Set up the pso parameter structure.*/
	struct psoParamStruct_debug psoParams;
	psoParams.popsize=4;
	psoParams.maxSteps= 5; 
	psoParams.c1=2;
	psoParams.c2=2;
	psoParams.max_velocity = 0.2;
	psoParams.dcLaw_a = 0.9;
	psoParams.dcLaw_b = 0.4;
	psoParams.dcLaw_c = psoParams.maxSteps;
	psoParams.dcLaw_d = 0.2;
	psoParams.locMinIter = 10;
	psoParams.locMinStpSz = 0.01;
	psoParams.rngGen = rngGen;
	psoParams.rngFile = fopen("matlabRng.txt","r");
	psoParams.debugDumpFile = fopen("test_ptapso_debug_Dump_matlabRng.txt","w");
	
	ptapso_debug(nDim, fitfunc, &inParams, psoParams, &psoResults);
	
	/* Test local minimization routines in gsl: replacement of fminsearch */
	/* Free allocated memory */
	gsl_vector_free(inParams.rmin);
	gsl_vector_free(inParams.realCoord);
	gsl_vector_free(inParams.rangeVec);
	gsl_vector_free(psoResults.bestLocation);
	gsl_rng_free(rngGen);
	
	fclose(psoParams.rngFile);
	fclose(psoParams.debugDumpFile);
}