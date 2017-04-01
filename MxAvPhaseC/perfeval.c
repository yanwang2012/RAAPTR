#include "maxphase.h"
#include "LLR_PSO.h"
#include "ptapso.h"
#include "perfeval.h"
#include "mat.h"
#include "matrix.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
/*! \file
\brief Run MaxPhase on a specified input data file.

The function \ref perfeval is the main one. 

\author Soumya D. Mohanty
*/
/*!
This function reads pulsar timing residual data stored in an input data file, 
searches for the maximum of the Log-Likelihood Ratio using Particle Swarm Optimization (PSO), 
and stores the results into an output file. PSO is run multiple times on the same 
input data using independent random positions and velocities.

The PSO function \ref ptapso can accept any fitness function. See \ref fitfunc_example to see the 
required interface for a fitness function. Here, the  Log-Likelihood Ratio fitness function \ref LLR_PSO
is sent to PSO.
*/
void perfeval(struct fitFuncParams *ffp, /*!< Parameters for the \ref LLR_PSO fitness function */
              char *inputFileName, /*!< Name of the file containing data to analyze*/
			  char *outputFileName /*!< Name of the file to store output results in*/){
	
	/* Number of independent PSO runs */
	size_t nRuns = 8;
	/* General purpose loop counters */
	size_t lpc1, lpc2;
	/*General purpose variables */
	int stat;

	
	/* Each run has an different random number seed.
	  {e,pi,
	   sqrt(2),euler's const gamma,
	   i^i, Feignbaum's number,
	   exp(pi), Champernowne const}
	*/
	unsigned long int rngSeeds[8] = {271828182,314159265,
	                                 141421356,577215664,
								     207879576,466920160,
								     231406926,123456789};
	
	double (*fitfunc)(gsl_vector *, void *) = LLR_PSO;
	
	/* Choose the random number generation method */
	gsl_rng *rngGen = gsl_rng_alloc(gsl_rng_taus);

	//Load data from specified .mat input file 
	MATFile *inputFilePr;
	inputFilePr =  matOpen(inputFileName, "r");
	mxArray *genHypothesis, *id;
	genHypothesis = matGetVariable(inputFilePr,"genHypothesis");
	id = matGetVariable(inputFilePr,"id");
	//Remaining data loads into special parameter structure
	struct llr_pso_params *llp;
	llp = loadfile2llrparam(inputFilePr);
	if((stat =  matClose(inputFilePr)))
		printf("Error closing file %s %d\n", inputFileName, stat);
	//Load special parameter struct into fitness function struct	
	ffp->splParams = llp;
		
	/* Configure PSO parameters*/
	struct psoParamStruct psoParams;
	psoParams.popsize= 40;
	psoParams.maxSteps= 2000; 
	psoParams.c1=2;
	psoParams.c2=2;
	psoParams.max_velocity = 0.2;
	psoParams.dcLaw_a = 0.9;
	psoParams.dcLaw_b = 0.4;
	psoParams.dcLaw_c = psoParams.maxSteps;
	psoParams.dcLaw_d = 0.2;
	psoParams.locMinIter = 0;
	psoParams.locMinStpSz = 0.01;
	psoParams.rngGen = rngGen;
	psoParams.debugDumpFile = NULL;
	
	
	/*----------------
	    RUN PSO
	-----------------*/
		//Number of search dimensions
	size_t nDim = ffp->nDim;
	struct returnData *psoResults[nRuns];
	for (lpc1 = 0; lpc1 < nRuns; lpc1++){
		psoResults[lpc1]= returnData_alloc(nDim);
	} 
    	//Number of extrinsic parameters = Number of pulsars
	size_t Np = llp->Np;
	    // Time to complete
	gsl_vector *wallClkTimes = gsl_vector_alloc(nRuns);
	    // PSO timer variables
	clock_t psoStartTime, psoStopTime;
		// Loop over runs
	for (lpc1 = 0; lpc1 < nRuns; lpc1++){
		/* Initialize random number generator */
		gsl_rng_set(rngGen,rngSeeds[lpc1]);
		
		/*Timed PSO run */
		psoStartTime = clock();
		ptapso(nDim, fitfunc, ffp, &psoParams, psoResults[lpc1]);
        psoStopTime = clock();
		
		gsl_vector_set(wallClkTimes,lpc1, (((double)(psoStopTime - psoStartTime))/CLOCKS_PER_SEC)/60.0);
				
	}
	
	/* Store results in output file */
	MATFile *outputFilePr;
	outputFilePr = matOpen(outputFileName,"w");
	
	/* copy input file name to mxArray for storing in output .mat file */	
	mxArray *inputFileNameMatStr;
	inputFileNameMatStr = mxCreateString(inputFileName);
	
	/*----------------
	   Store into .mat file
	-----------------*/
		// Input file name
	if (matPutVariable(outputFilePr,"inputFile",inputFileNameMatStr))
		printf("Error storing variable %s in file %s\n", "inputFile", outputFileName);
		// Metadata from input file
	if (matPutVariable(outputFilePr,"id",id))
		printf("Error storing variable %s in file %s\n", "id", outputFileName);
	if (matPutVariable(outputFilePr,"genHypothesis",genHypothesis))
		printf("Error storing variable %s in file %s\n", "genHypothesis", outputFileName);
		// Store all other variables
	perfevaloutput2matfile(nRuns, nDim, Np, wallClkTimes,
					       psoResults, ffp, fitfunc, 
					       outputFilePr);
			
	if ((stat =  matClose(outputFilePr)))
		printf("Error closing file %s %d\n", outputFileName, stat);

	/* 
	  Wrap up
	*/
	gsl_rng_free(rngGen);
	llrparam_free(llp);
	for (lpc1 = 0; lpc1 < nRuns; lpc1++){
		returnData_free(psoResults[lpc1]);
	} 
	mxDestroyArray(inputFileNameMatStr);
	mxDestroyArray(genHypothesis);
	mxDestroyArray(id);
	gsl_vector_free(wallClkTimes);
}


/*! Dump output from multiple pso runs in \ref perfeval to a .mat file */
void perfevaloutput2matfile(size_t nRuns, size_t nDim, size_t Np,  
                      gsl_vector *wallClckTimes,
					  struct returnData *psoResults[],
					  struct fitFuncParams *ffp,
                      double (*fitfunc)(gsl_vector *, void *), 
					  MATFile *outputFilePr){

	//Loop counters
	size_t lpc1, lpc2, lpc3;
	
	/* Create storage for results from multiple PSO runs. Since
	   these will be stored in the output .mat files, we use 
	   mxArray's for these variables.
	*/
	mxArray *nRuns4Matfile, *bestRun4MatFile;
	nRuns4Matfile = mxCreateDoubleScalar(nRuns);
	
	/* Variables to store results from all the runs */	
	mxArray *bestLocationVec, *bestLocRealC, *bestFitValVec;
	mxArray *nFuncEvalsVec, *nIterVec, *wallClkTimeVec;
	/* Variables to store results from the best run */
	mxArray *bestRunBestLocation4Mat, *bestRunRealC;

	
	/* Results from all the runs */
	bestLocationVec = mxCreateDoubleMatrix(nRuns,nDim,mxREAL);
	bestLocRealC  = mxCreateDoubleMatrix(nRuns,nDim+Np,mxREAL);
	bestFitValVec = mxCreateDoubleMatrix(nRuns,1,mxREAL);
	nFuncEvalsVec = mxCreateDoubleMatrix(nRuns,1,mxREAL);
	nIterVec = mxCreateDoubleMatrix(nRuns,1,mxREAL);
	wallClkTimeVec = mxCreateDoubleMatrix(nRuns,1,mxREAL);
	/* Results from the best run */
	bestRunBestLocation4Mat = mxCreateDoubleMatrix(1,nDim,mxREAL); //Unstandardized location
	bestRunRealC = mxCreateDoubleMatrix(1,nDim+Np,mxREAL); // Location in real coordinates
	
	
	double *bestLocationVecPr, *bestLocRealCPr, *bestFitValVecPr;
	double *nFuncEvalsVecPr, *nIterVecPr, *wallClkTimeVecPr;
	double *bestRunBestLocation4MatPr, *bestRunRealCPr;
	
	bestLocationVecPr = mxGetPr(bestLocationVec);
	bestLocRealCPr = mxGetPr(bestLocRealC);
	bestFitValVecPr = mxGetPr(bestFitValVec);
	nFuncEvalsVecPr = mxGetPr(nFuncEvalsVec);
	nIterVecPr = mxGetPr(nIterVec);
    wallClkTimeVecPr = mxGetPr(wallClkTimeVec);
	bestRunBestLocation4MatPr = mxGetPr(bestRunBestLocation4Mat);
	bestRunRealCPr = mxGetPr(bestRunRealC);
	
	//Get results from all the runs
	for(lpc1 = 0; lpc1 < nRuns; lpc1++){
		//Fitness value
		bestFitValVecPr[lpc1] = psoResults[lpc1]->bestFitVal;
		//standardized coordinates (intrinsic params only)
		for(lpc2 = 0; lpc2 < nDim; lpc2 ++){
			bestLocationVecPr[lpc1+lpc2*nRuns] = gsl_vector_get(psoResults[lpc1]->bestLocation,lpc2);
		}
		//Other quantities
		nFuncEvalsVecPr[lpc1] = psoResults[lpc1]->totalFuncEvals;
		nIterVecPr[lpc1] = psoResults[lpc1]->totalIterations;
		wallClkTimeVecPr[lpc1] = gsl_vector_get(wallClckTimes,lpc1);
	}
	//Get the unstandardized coordinates
	for(lpc1 = 0; lpc1 < nRuns; lpc1++){
		fitfunc(psoResults[lpc1]->bestLocation,ffp);
		for(lpc2 = 0; lpc2 < nDim; lpc2++){
			bestLocRealCPr[lpc1+lpc2*nRuns] = gsl_vector_get(ffp->realCoord,lpc2);
	    }
		//Append PhiI values (extrinsic parameters)
		for (lpc2 = nDim; lpc2 < nDim+Np; lpc2++){
			bestLocRealCPr[lpc1+lpc2*nRuns] = ((struct llr_pso_params *)ffp->splParams)->phiI[lpc2-nDim];
		}
	}
	
	//Find the best run
	size_t bestRun = 0;
	double minFitVal = INFINITY;
	for (lpc1 = 0; lpc1 < nRuns; lpc1++){
		if(bestFitValVecPr[lpc1] < minFitVal){
			minFitVal = bestFitValVecPr[lpc1];
			bestRun = lpc1;
		}
	}
	bestRun4MatFile = mxCreateDoubleScalar((double)(bestRun+1));
	
	//Get the standardized location for the best run
	for(lpc1 = 0; lpc1 < nDim; lpc1++){
		bestRunBestLocation4MatPr[lpc1] = bestLocationVecPr[bestRun+lpc1*nRuns];
	}
	
	//Get the unstandardized location for the best run (including phiI values)
	for(lpc1 = 0; lpc1 < nDim+Np; lpc1++){
		bestRunRealCPr[lpc1] = bestLocRealCPr[bestRun+lpc1*nRuns];
	}

	//Store in .mat file
	if (matPutVariable(outputFilePr,"nRuns",nRuns4Matfile))
		printf("Error storing variable %s\n", "nRuns");
	
	if (matPutVariable(outputFilePr,"bestLocVals",bestLocationVec))
		printf("Error storing variable %s\n", "bestLocVals");
	
	if (matPutVariable(outputFilePr,"bestLocRealCVals",bestLocRealC))
		printf("Error storing variable %s\n", "bestLocRealCVals");

	if (matPutVariable(outputFilePr,"fitnessVals",bestFitValVec))
		printf("Error storing variable %s\n", "fitnessVals");	
	
	if (matPutVariable(outputFilePr,"numFitEvals",nFuncEvalsVec))
		printf("Error storing variable %s\n", "numFitEvals");
	
	if (matPutVariable(outputFilePr,"numIter",nIterVec))
		printf("Error storing variable %s\n", "numIter");			
	
	if (matPutVariable(outputFilePr,"time2complete",wallClkTimeVec))
		printf("Error storing variable %s\n", "time2complete");	
	
	if (matPutVariable(outputFilePr,"bestRun",bestRun4MatFile))
		printf("Error storing variable %s\n", "bestRun");
	
	if (matPutVariable(outputFilePr,"bestLocation",bestRunBestLocation4Mat))
		printf("Error storing variable %s\n", "bestLocation");
	
	if (matPutVariable(outputFilePr,"bestRealLoc",bestRunRealC))
		printf("Error storing variable %s\n", "bestRealLoc");	

    //Wrap up
	mxDestroyArray(bestLocRealC);
	mxDestroyArray(bestRunRealC);
	mxDestroyArray(nRuns4Matfile);
	mxDestroyArray(bestRun4MatFile);
	mxDestroyArray(bestLocationVec);
	mxDestroyArray(bestFitValVec);
	mxDestroyArray(nFuncEvalsVec);
	mxDestroyArray(nIterVec);
	mxDestroyArray(wallClkTimeVec);
	mxDestroyArray(bestRunBestLocation4Mat);
}