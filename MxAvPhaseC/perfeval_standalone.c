#include "maxphase.h"
#include "LLR_PSO.h"
#include "ptapso.h"
#include "mat.h"
#include "matrix.h"
#include <gsl/gsl_vector.h>
#include <sys/stat.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

/* Usage: perfeval search_param_file input_data_dir */
int main(int argc, char *argv[]){
	
	/* General purpose loop counters */
	size_t lpc1, lpc2, lpc3;
	/*General purpose variables */
	int stat;
	double dummyFitVal, minFitVal;
	clock_t psoStartTime, psoStopTime;
	
	/* Number of independent PSO runs */
	size_t nRuns = 8;
	/* Index of best run */
	size_t bestRun;
	
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
	
	double (*fitfunc)(const gsl_vector *, void *) = LLR_PSO;
	
	/* Error handling off */
	gsl_error_handler_t *old_handler = gsl_set_error_handler_off ();
	
	/* Choose the random number generation method */
	gsl_rng *rngGen = gsl_rng_alloc(gsl_rng_taus);
	
	
	if (argc != 3){
		printf("Usage: %s parameter_file_path analysis_directory\n", argv[0]);
		gsl_rng_free(rngGen);
		return 1;
	}
	/* Full path to search parameter file */
	char *srchParamsFile = argv[1];
	/* Full path to directory containing input files */
	char *simDataDir = argv[2];
	
	/* Read xmaxmin from srchParamsFile file */
	MATFile *srchPar = matOpen(srchParamsFile, "r");
	if (srchPar == NULL){
		printf("Error reading parameter file %s\n", srchParamsFile);
		gsl_rng_free(rngGen);
		return 2;
	}
	mxArray *xmaxmin;
	xmaxmin = matGetVariable(srchPar,"xmaxmin");
	
	/* xmaxmin should be a 2 dimensional array */
	if (mxGetNumberOfDimensions(xmaxmin) !=2){
		printf("Parameter specification is not understood\n");
		gsl_rng_free(rngGen);
		mxDestroyArray(xmaxmin);
		return 3;
	}
	/* Search Space dimensionality */
	size_t nDim = mxGetM(xmaxmin);
	
	/* transfer xmaxmin to fitness function parameter struct */
	double *xmaxminData = mxGetPr(xmaxmin);
	struct fitFuncParams *ffp = ffparam_alloc(nDim);
	struct llr_pso_params llp; 
	for(lpc1 = 0; lpc1 < nDim; lpc1++){
		gsl_vector_set(ffp->rmin,lpc1,xmaxminData[lpc1+nDim]);
		gsl_vector_set(ffp->rangeVec,lpc1,xmaxminData[lpc1]-xmaxminData[lpc1+nDim]);
	}
	/* Close file */
	if ((stat =  matClose(srchPar)))
		printf("Error closing file %s %d\n", srchParamsFile, stat);
	
	/* Path to folder for storing outputs */
	char *outDataDir = (char *)malloc( (strlen(simDataDir)+
		                                strlen("results")+
									    1 + 1)*sizeof(char) );
	strcpy(outDataDir,simDataDir);
	strcat(outDataDir,"/");
	strcat(outDataDir,"results");
	if (mkdir(outDataDir,S_IRWXU)){
		printf("Error creating results directory %s\n", outDataDir);
		printf("Check that it does not already exist.\n");
		gsl_rng_free(rngGen);
		mxDestroyArray(xmaxmin);
		ffparam_free(ffp);
		free(outDataDir);
		return 4;
	}
	
	/* Create list of input data files */
	char **inputFileList;
	size_t nInputFiles;
	//Dummy variable for compatibility with new version of listfileswext
	size_t dummyFileNameLen;
	inputFileList = listfileswext("mat", simDataDir, &nInputFiles, &dummyFileNameLen);
	printf("Number of input files %zu\n",nInputFiles);
	
	/*Find longest filename */
	size_t fileNameLen;
	size_t maxFileNameLen = 0;
	for (lpc1 = 0; lpc1 < nInputFiles; lpc1++){
		fileNameLen = strlen(inputFileList[lpc1]);
		if ( fileNameLen > maxFileNameLen)
			maxFileNameLen = fileNameLen;
	}
	
	/* Sufficient storage for input and output file name strings: 
	directory name length +2 (for '/0' and '/' (filesep)).
	*/
	char *inputFileName, *outputFileName;
	MATFile *inputFilePr, *outputFilePr;
	outputFileName = (char *)malloc((strlen(outDataDir)+maxFileNameLen+2)*sizeof(char));
	inputFileName = (char *)malloc((strlen(simDataDir)+maxFileNameLen+2)*sizeof(char));
	
	/* Variables to hold input file content */
	mxArray *simParams, *timingResiduals, *yr, *genHypothesis, *id;
	mxArray *Np, *N, *sd, *alphaP, *deltaP;
	double **s;
	/* Create storage for results from multiple PSO runs. Since
	   these will be stored in the output .mat files, we use 
	   mxArray's for these variables.
	*/
	mxArray *inputFileNameMatStr;
	mxArray *nRuns4Matfile;
	mxArray *bestLocationVec, *bestLocRealC, *bestFitValVec;
	mxArray *nFuncEvalsVec, *nIterVec, *wallClkTimeVec;
	mxArray *bestFitIndx, *bestLocation4Mat, *realC;
	double *bestLocationVecPr, *bestLocRealCPr, *bestFitValVecPr;
	double *nFuncEvalsVecPr, *nIterVecPr, *wallClkTimeVecPr;
	double *bestFitIndxPr, *bestLocation4MatPr, *realCPr;
	nRuns4Matfile = mxCreateDoubleScalar(nRuns);
	bestLocationVec = mxCreateDoubleMatrix(nRuns,nDim,mxREAL);
	
	bestFitValVec = mxCreateDoubleMatrix(nRuns,1,mxREAL);
	nFuncEvalsVec = mxCreateDoubleMatrix(nRuns,1,mxREAL);
	nIterVec = mxCreateDoubleMatrix(nRuns,1,mxREAL);
	wallClkTimeVec = mxCreateDoubleMatrix(nRuns,1,mxREAL);
	bestFitIndx = mxCreateDoubleScalar(-1.0);
	bestLocation4Mat = mxCreateDoubleMatrix(1,nDim,mxREAL);
	
	
	bestLocationVecPr = mxGetPr(bestLocationVec);
	
	bestFitValVecPr = mxGetPr(bestFitValVec);
	nFuncEvalsVecPr = mxGetPr(nFuncEvalsVec);
	nIterVecPr = mxGetPr(nIterVec);
    wallClkTimeVecPr = mxGetPr(wallClkTimeVec);
	bestFitIndxPr = mxGetPr(bestFitIndx);
	bestLocation4MatPr = mxGetPr(bestLocation4Mat);
	
	
     /*nInputFiles = 4; printf("DEBUGGING ON!\n");*/

	/* Set up for ptapso. */
	struct returnData psoResults;
	gsl_vector *bestLocation = gsl_vector_alloc(nDim);
	psoResults.bestLocation = bestLocation;
	struct psoParamStruct psoParams;
	psoParams.popsize=40;
	psoParams.maxSteps= 2000; 
	psoParams.c1=2;
	psoParams.c2=2;
	psoParams.max_velocity = 0.2;
	psoParams.dcLaw_a = 0.9;
	psoParams.dcLaw_b = 0.4;
	psoParams.dcLaw_c = psoParams.maxSteps;
	psoParams.dcLaw_d = 0.2;
	psoParams.locMinIter = 200;
	psoParams.locMinStpSz = 0.01;
	psoParams.rngGen = rngGen;
	psoParams.debugDumpFile = NULL;

	/* Analyze each input file */
	for (lpc1 = 0; lpc1 < nInputFiles; lpc1++){
		
		/* Construct input file name along with path */
		strcpy(inputFileName,simDataDir);
		strcat(inputFileName,"/");
		strcat(inputFileName,inputFileList[lpc1]);
		/* copy file name to mxArray for storing in output .mat file */
		inputFileNameMatStr = mxCreateString(inputFileName);
		/* Construct output file name along with path */
		strcpy(outputFileName,outDataDir);
		strcat(outputFileName,"/");
		strcat(outputFileName,inputFileList[lpc1]);
		
		/* Load contents of inPutFileName into
		   fitness function parameter structures
		*/
		inputFilePr =  matOpen(inputFileName, "r");
		simParams = matGetVariable(inputFilePr,"simParams");
		yr = matGetVariable(inputFilePr,"yr");
		timingResiduals = matGetVariable(inputFilePr,"timingResiduals");
		genHypothesis = matGetVariable(inputFilePr,"genHypothesis");
		id = matGetVariable(inputFilePr,"id");
		/* Extract fields of simParams */
		Np = mxGetField(simParams, 0, "Np");
		N = mxGetField(simParams, 0, "N");
		sd = mxGetField(simParams, 0, "sd");
		alphaP = mxGetField(simParams, 0, "alphaP");
		deltaP = mxGetField(simParams, 0, "deltaP");
		if ((stat =  matClose(inputFilePr)))
			printf("Error closing file %s %d\n", inputFileName, stat);

		/* Load fitness function parameter structure */
		llp.Np = (unsigned int)mxGetScalar(Np);
		llp.N = (unsigned int)mxGetScalar(N);
		
		/* The variable 'kp' is not needed in the C code.
		   Some variables to be stored out in the .mat file require Np
		 */
		bestLocRealC  = mxCreateDoubleMatrix(nRuns,nDim+llp.Np,mxREAL);
		realC = mxCreateDoubleMatrix(1,nDim+llp.Np,mxREAL);
		bestLocRealCPr = mxGetPr(bestLocRealC);
		realCPr = mxGetPr(realC);
		
		/* Allocate storage for timing residuals field 's' in fitness
		function parameter struct. This is needed because
		timingResiduals is a 2D array while it is actually stored as 
		a 1 D array in the mxArray. The allocation must be done in the loop
		because Np is read from the input file.
		*/
		s = (double **)malloc(llp.Np*sizeof(double *));
		for (lpc2 = 0; lpc2 < llp.Np; lpc2++){
			s[lpc2] = (double *) malloc(llp.N*sizeof(double));
		}

		/* load timing residuals into fitness function param struct */
		double *trPr = mxGetPr(timingResiduals);
		unsigned int countElements = 0;
		for (lpc2 = 0; lpc2 < llp.N; lpc2++){
			for (lpc3 = 0; lpc3 < llp.Np; lpc3++){
				s[lpc3][lpc2] = trPr[countElements];
				countElements++;
			}
		}
		
		/* Allocate storage for phiI: special output returned by LLR_PSO */
		llp.phiI = (double *) malloc(llp.Np*sizeof(double));
		/* Transfering 1D arrays is easy through their pointers */
		llp.s = s;
		llp.sd = mxGetPr(sd);
		llp.alphaP = mxGetPr(alphaP);
		llp.deltaP = mxGetPr(deltaP);
		llp.yr = mxGetPr(yr);
		
				
		/* Finally pass the special fitness function struct into the 
		   standard fitness func struct */
		ffp->splParams = &llp;
		/* Independent PSO runs */
		printf("Input file %s\n",inputFileName);
		/* Track best run */
		bestRun = -1;
		minFitVal = INFINITY;
		for (lpc2 = 0; lpc2 < nRuns; lpc2++){
			
			/* Initialize random number generator */
			gsl_rng_set(rngGen,rngSeeds[lpc2]);
			
			/*Timed PSO run */
			psoStartTime = clock();
			ptapso(nDim, fitfunc, ffp, psoParams, &psoResults);
	        psoStopTime = clock();
			wallClkTimeVecPr[lpc2] = (((double)(psoStopTime - psoStartTime))/CLOCKS_PER_SEC)/60.0;
			
			bestFitValVecPr[lpc2] = psoResults.bestFitVal;
			nFuncEvalsVecPr[lpc2] = psoResults.totalFuncEvals;
			nIterVecPr[lpc2] = psoResults.totalIterations;
			
			for(lpc3 = 0; lpc3 < nDim; lpc3++){
				printf("%f ",gsl_vector_get(psoResults.bestLocation,lpc3));
				bestLocationVecPr[lpc2+lpc3*nRuns] = gsl_vector_get(psoResults.bestLocation,lpc3);
				dummyFitVal = fitfunc(psoResults.bestLocation,ffp);
				bestLocRealCPr[lpc2+lpc3*nRuns] = gsl_vector_get(ffp->realCoord,lpc3);
			}
			printf("\n");
			printf("--> ");
			/* Load phiI values after the nDim parameter values */
			for (lpc3 = nDim; lpc3 < nDim+llp.Np; lpc3++){
				printf("%f ",llp.phiI[lpc3-nDim]);
				bestLocRealCPr[lpc2+lpc3*nRuns]=llp.phiI[lpc3-nDim];
			}
			printf("\n");
			
						
			if(psoResults.bestFitVal < minFitVal){
				minFitVal = psoResults.bestFitVal;
				bestRun = lpc2+1;
				for(lpc3 = 0; lpc3 < nDim; lpc3++){
					bestLocation4MatPr[lpc3] = bestLocationVecPr[lpc2+lpc3*nRuns];
					/* FIXME May be more efficient to read off from the gsl_vector */
					realCPr[lpc3] = bestLocRealCPr[lpc2+lpc3*nRuns];
				}
				/* Store phiI values for best run */
				for (lpc3 = nDim; lpc3 < nDim+llp.Np; lpc3++){
					realCPr[lpc3]=llp.phiI[lpc3-nDim];
				}
				bestFitIndxPr[0] = (double) bestRun;				
			}
			

			
			/*
			printf("Run #%zu\n",lpc2+1);
			printf("Best Location found: \n");
			for (lpc3 = 0; lpc3 < nDim; lpc3++){
				printf("%.18f, ",gsl_vector_get(psoResults.bestLocation,lpc3));		
			}
			printf("\n");
			printf("Best Fitness Value: %e\n", psoResults.bestFitVal);
			printf("Phi Values:");
			for (lpc3 = 0; lpc3 < llp.Np; lpc3++){
				printf("%f ",llp.phiI[lpc3]);
			}
			printf("\n");
			*/
			
		}
		/* printf("***********************************************************\n");*/
		
        
		
		/* Store results in output file */
		outputFilePr = matOpen(outputFileName,"w");
		if (matPutVariable(outputFilePr,"inputFile",inputFileNameMatStr))
			printf("Error storing variable %s in file %s\n", "inputFile", outputFileName);
		
		if (matPutVariable(outputFilePr,"id",id))
			printf("Error storing variable %s in file %s\n", "id", outputFileName);
		
		if (matPutVariable(outputFilePr,"genHypothesis",genHypothesis))
			printf("Error storing variable %s in file %s\n", "genHypothesis", outputFileName);
		
		if (matPutVariable(outputFilePr,"nRuns",nRuns4Matfile))
			printf("Error storing variable %s in file %s\n", "nRuns", outputFileName);
		
		if (matPutVariable(outputFilePr,"bestLocVals",bestLocationVec))
			printf("Error storing variable %s in file %s\n", "bestLocVals", outputFileName);
		
		if (matPutVariable(outputFilePr,"bestLocRealCVals",bestLocRealC))
			printf("Error storing variable %s in file %s\n", "bestLocRealCVals", outputFileName);

		if (matPutVariable(outputFilePr,"fitnessVals",bestFitValVec))
			printf("Error storing variable %s in file %s\n", "fitnessVals", outputFileName);	
		
		if (matPutVariable(outputFilePr,"numFitEvals",nFuncEvalsVec))
			printf("Error storing variable %s in file %s\n", "numFitEvals", outputFileName);
		
		if (matPutVariable(outputFilePr,"numIter",nIterVec))
			printf("Error storing variable %s in file %s\n", "numIter", outputFileName);			
		
		if (matPutVariable(outputFilePr,"time2complete",wallClkTimeVec))
			printf("Error storing variable %s in file %s\n", "time2complete", outputFileName);	
		
		if (matPutVariable(outputFilePr,"bestRun",bestFitIndx))
			printf("Error storing variable %s in file %s\n", "bestRun", outputFileName);
		
		if (matPutVariable(outputFilePr,"bestLocation",bestLocation4Mat))
			printf("Error storing variable %s in file %s\n", "bestLocation", outputFileName);
		
		if (matPutVariable(outputFilePr,"bestRealLoc",realC))
			printf("Error storing variable %s in file %s\n", "bestRealLoc", outputFileName);		
				
		if ((stat =  matClose(outputFilePr)))
			printf("Error closing file %s %d\n", outputFileName, stat);

		/* Destroy dynamically allocated memory before reading
		   next file's contents.
		*/
		for (lpc2 = 0; lpc2 < llp.Np; lpc2++){
			free(s[lpc2]);
		}
		free(s);
		free(llp.phiI);
		mxDestroyArray(inputFileNameMatStr);
		mxDestroyArray(timingResiduals);
		mxDestroyArray(yr);
		mxDestroyArray(genHypothesis);
		mxDestroyArray(id);
		mxDestroyArray(bestLocRealC);
		mxDestroyArray(realC);
		/* Destroying an mxArray holding a structure also
		   destroys its fields. So no need to deallocate
		   Np, N, sd, deltaP, alphaP separately.
		*/
		mxDestroyArray(simParams);
	}
		
	/* ----------------------------
	Deallocate storage 
	-----------------------------*/
	free(outDataDir);
	gsl_rng_free(rngGen);
	gsl_vector_free(bestLocation);
	ffparam_free(ffp);
	mxDestroyArray(xmaxmin);
	/*mxDestroyArray(nRunsPr);*/
	for(lpc1 = 0; lpc1 < nInputFiles; lpc1++){
	    /* Free up dynamically allocated memory for file list*/
		free(inputFileList[lpc1]);
	}
	/* Free up dynamically allocated memory */
	free(inputFileList);
	free(inputFileName);
	free(outputFileName);
	/* Memory allocated to variables storing the output from 
	   ptapso
	*/
	mxDestroyArray(nRuns4Matfile);
	mxDestroyArray(bestLocationVec);
	mxDestroyArray(bestFitValVec);
	mxDestroyArray(nFuncEvalsVec);
	mxDestroyArray(nIterVec);
	mxDestroyArray(wallClkTimeVec);
	mxDestroyArray(bestFitIndx);
	mxDestroyArray(bestLocation4Mat);
	
	
	/* Everything executed successfully */
	return 0;
}