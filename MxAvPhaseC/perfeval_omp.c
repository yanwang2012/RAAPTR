#include "maxphase.h"
#include "LLR_Mp_Av.h"
#include "ptapso.h"
#include "perfeval_omp.h"
#include "gslhdf5_io.h"
#include "hdf5_hl.h"
#include "omp.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
/*! \file
\brief Run MaxPhase or AvPhase on a specified input data file.

\author Soumya D. Mohanty
*/
/*!
This function reads pulsar timing residual data stored in an input data file, 
searches for the maximum of the Log-Likelihood Ratio using Particle Swarm Optimization (PSO), 
and stores the results into an output file. PSO is run multiple times on the same 
input data using independent random positions and velocities. The independent runs are parallelized 
 using OpenMP (OMP).

The function ptapso() can apply PSO to any fitness function provided it has the interface
 specied in how_to_code_fitnessFunc.txt. In the present case,
 the  Log-Likelihood Ratio fitness function LLR_av() (or LLR_mp())
is sent to PSO.
*/
void perfeval_omp(struct fitFuncParams *ffp, /*!< Parameters for the fitness function */
              char *inputFileName, /*!< Name of the file containing data to analyze*/
			  char *outputFileName, /*!< Name of the file to store output results in*/
			  char *mp_av_select /*!< Select Max or AvPhase algorithm */){
	
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

	 /* Choose the algorithm to use. */
	 double (*fitfunc)(gsl_vector *, void *);
	 if (!strcmp(mp_av_select,"maxPhase")){
		 fitfunc = LLR_mp;
	 }
	 else if(!strcmp(mp_av_select,"avPhase")){
		 fitfunc = LLR_av;
	 }
	
	/* Choose the random number generation method.
	   Need independent random number generators 
	   inside an OMP for loop (?)*/
	gsl_rng *rngGen[8];
	for(lpc1 = 0; lpc1 < 8; lpc1++){
		rngGen[lpc1] = gsl_rng_alloc(gsl_rng_taus);
	}

	//Load data from specified .hdf5 input file 
	herr_t status;
	hid_t inFile = H5Fopen(inputFileName, H5F_ACC_RDONLY, H5P_DEFAULT);
	if (inFile < 0){
		fprintf(stdout,"Error opening file\n");
		abort();
	}
	char genHypothesis[10];
	status = H5LTread_dataset_string(inFile,"genHypothesis",genHypothesis);
	if (status < 0){
		fprintf(stdout,"Error reading genHypothesis\n");
		abort();
	}
	//fprintf(stdout,"genHypothesis: %s\n",genHypothesis);
	
	struct {   
		int snr_id;
        int loc_id;
        int omg_id;
        int rlz_id;} id;
    status = H5LTread_dataset_int(inFile,"snr_id",&id.snr_id);
	if (status < 0){
		fprintf(stdout,"Error reading snr_id\n");
		abort();
	}
    status = H5LTread_dataset_int(inFile,"loc_id",&id.loc_id);
	if (status < 0){
		fprintf(stdout,"Error reading loc_id\n");
		abort();
	}
    status = H5LTread_dataset_int(inFile,"omg_id",&id.omg_id);
	if (status < 0){
		fprintf(stdout,"Error reading omg_id\n");
		abort();
	}
    status = H5LTread_dataset_int(inFile,"rlz_id",&id.rlz_id);
	if (status < 0){
		fprintf(stdout,"Error reading rlz_id\n");
		abort();
	}

	//Remaining data loads into special parameter structure
	struct llr_pso_params *llp;
	llp = loadfile2llrparam(inFile);
	
	//Close file
	status = H5Fclose(inFile);
	if(status < 0){
		fprintf(stdout,"Error closing file %s \n", inputFileName);
	}

	//Load special parameter struct into fitness function struct	
	ffp->splParams = llp;
		
	
	
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
	// Loop over runs using omp. 
	/* We must use deep copies of 
	   the fitness parameter structure because otherwise it is 
       shared and different OMP workers will overwrite the realCoord
	   field. Similarly copies of other variables need to be used if there is a 
	   chance that sharing them will lead to intereference between the workers.
	   This is a not a problem in an MPI code since one must explicity pass
	   values of variables to the workers. 
	*/
	struct fitFuncParams *ffpCopy;
	struct psoParamStruct psoParams;
#pragma omp parallel for private(psoStartTime,psoStopTime,ffpCopy,psoParams)
	for (lpc1 = 0; lpc1 < nRuns; lpc1++){
		fprintf(stdout,"Running PSO run %zu on worker %d\n",lpc1,omp_get_thread_num());
		/* Initialize random number generator separately for each worker*/
		gsl_rng_set(rngGen[lpc1],rngSeeds[lpc1]);
	    /* Configure PSO parameters*/
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
		psoParams.debugDumpFile = NULL;
		psoParams.rngGen = rngGen[lpc1];
		
	    /*Clone fitness parameter struct */
		ffpCopy = ffparams_clone(ffp);
		
		/*Timed PSO run */
		psoStartTime = clock();
		ptapso(nDim, fitfunc, ffpCopy, &psoParams, psoResults[lpc1]);
        psoStopTime = clock();
	
		gsl_vector_set(wallClkTimes,lpc1, (((double)(psoStopTime - psoStartTime))/CLOCKS_PER_SEC)/60.0);
		
		//Delete local copy of fitness parameter struct
		llrparam_free(ffpCopy->splParams);
		ffparam_free(ffpCopy);		
	}
	
	/* Store results in output file */
	hsize_t dims;
	// MATFile *outputFilePr;
// 	outputFilePr = matOpen(outputFileName,"w");
	hid_t outFile = H5Fcreate(outputFileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if (outFile < 0){
		fprintf(stdout,"Error opening file %s \n", outputFileName);
		abort();
	}
	//Store input file name
	status =  H5LTmake_dataset_string(outFile, "inputFile", inputFileName);
	if (status < 0){
		fprintf(stdout, "Error writing variable inputFile in file %s\n", outputFileName);
	}	
	//Store hypothesis used for data generation
	status =  H5LTmake_dataset_string(outFile, "genHypothesis", genHypothesis);
	
	//Store out fields of 'id' struct
	/* 	struct {   
		int snr_id;
        int loc_id;
        int omg_id;
        int rlz_id;} id;
	*/
	dims = 1;
	status = H5LTmake_dataset_int(outFile,"snr_id",1,&dims,&id.snr_id);
	if (status < 0){
		fprintf(stdout, "Error writing variable snr_id in file %s\n", outputFileName);
	}	
	status = H5LTmake_dataset_int(outFile,"loc_id",1,&dims,&id.loc_id);
	if (status < 0){
		fprintf(stdout, "Error writing variable loc_id in file %s\n", outputFileName);
	}
	status = H5LTmake_dataset_int(outFile,"omg_id",1,&dims,&id.omg_id);
	if (status < 0){
		fprintf(stdout, "Error writing variable omg_id in file %s\n", outputFileName);
	}
	status = H5LTmake_dataset_int(outFile,"rlz_id",1,&dims,&id.rlz_id);
	if (status < 0){
		fprintf(stdout, "Error writing variable rlz_id in file %s\n", outputFileName);
	}
		
    // Store all other variables
	perfevalomp2hdf5file(nRuns, nDim, Np, wallClkTimes,
					       psoResults, ffp, fitfunc, mp_av_select,
					       outFile);
						   
    //Close output file
   	status = H5Fclose(outFile);
   	if(status < 0){
   		fprintf(stdout,"Error closing file %s\n", outputFileName);
   	}

	/* 
	  Wrap up
	*/
	for (lpc1 = 0; lpc1 < 8; lpc1++){
		gsl_rng_free(rngGen[lpc1]);		
	}
	llrparam_free(llp);
	for (lpc1 = 0; lpc1 < nRuns; lpc1++){
		returnData_free(psoResults[lpc1]);
	} 
	gsl_vector_free(wallClkTimes);
}


/*! Dump output from multiple pso runs in perfeval_omp() to a .mat file */
void perfevalomp2hdf5file(size_t nRuns, size_t nDim, size_t Np,  
                      gsl_vector *wallClckTimes,
					  struct returnData *psoResults[],
					  struct fitFuncParams *ffp,
                      double (*fitfunc)(gsl_vector *, void *), 
					  char *mp_av_select,
					  hid_t outFile){

	//Loop counters
	size_t lpc1, lpc2, lpc3;
	//Dummy variable
	double dummyFitVal, dummyAmp;
	
	/* Create storage for results from multiple PSO runs. 
	  (This is very inefficient because it is a translation
	   from the earlier code that used  mat file output.)
	*/
	gsl_matrix *bestLocationVecPr = gsl_matrix_alloc(nRuns,nDim);
	gsl_matrix *bestLocRealCPr = gsl_matrix_alloc(nRuns,nDim+Np);
	gsl_vector *bestFitValVecPr = gsl_vector_alloc(nRuns);
	gsl_vector *nFuncEvalsVecPr = gsl_vector_alloc(nRuns);
	gsl_vector *nIterVecPr = gsl_vector_alloc(nRuns);
	gsl_vector *wallClkTimeVecPr = gsl_vector_alloc(nRuns);
	gsl_vector *bestRunBestLocationPr = gsl_vector_alloc(nDim);
	gsl_vector *bestRunRealCPr = gsl_vector_alloc(nDim+Np);
	
	//Get results from all the runs
	for(lpc1 = 0; lpc1 < nRuns; lpc1++){
		//Fitness value
		gsl_vector_set(bestFitValVecPr, lpc1, psoResults[lpc1]->bestFitVal);
		//standardized coordinates (intrinsic params only)
		for(lpc2 = 0; lpc2 < nDim; lpc2 ++){
			gsl_matrix_set(bestLocationVecPr,lpc1,lpc2,gsl_vector_get(psoResults[lpc1]->bestLocation,lpc2));
		}
		//Other quantities
		gsl_vector_set(nFuncEvalsVecPr, lpc1, psoResults[lpc1]->totalFuncEvals);
		gsl_vector_set(nIterVecPr, lpc1, psoResults[lpc1]->totalIterations);
		gsl_vector_set(wallClkTimeVecPr, lpc1, gsl_vector_get(wallClckTimes,lpc1));
	}
	//Get the unstandardized coordinates
	for(lpc1 = 0; lpc1 < nRuns; lpc1++){
		fitfunc(psoResults[lpc1]->bestLocation,ffp);
		//For AvPhase, also get the pulsar phases using MaxPhase
		if(!strcmp(mp_av_select,"avPhase")){
			/* The previous call to fitfunc sets the real coordinate values,
			   allowing LogLikelihoodRatioMP5 to be called and pulsar phases to 
			   be estimated. However, the amplitude needs to be reconverted back to
			   log because its anti-log is taken inside LogLikelihoodRatioMP5 and
			   returned in the fitness function parameter structure.
			*/
			dummyAmp = gsl_vector_get(ffp->realCoord,4);
			dummyAmp = log10(dummyAmp);
			gsl_vector_set(ffp->realCoord,4,dummyAmp);
			dummyFitVal = LogLikelihoodRatioMP5(ffp); 
	    }
			 
		for(lpc2 = 0; lpc2 < nDim; lpc2++){
			gsl_matrix_set(bestLocRealCPr, lpc1, lpc2, gsl_vector_get(ffp->realCoord,lpc2));
	    }
		//Append PhiI values (extrinsic parameters)
		for (lpc2 = nDim; lpc2 < nDim+Np; lpc2++){
			gsl_matrix_set(bestLocRealCPr,lpc1, lpc2, ((struct llr_pso_params *)ffp->splParams)->phiI[lpc2-nDim]);
		}
	}
	
	//Find the best run
	size_t bestRun = 0;
	double minFitVal = INFINITY;
	for (lpc1 = 0; lpc1 < nRuns; lpc1++){
		if(gsl_vector_get(bestFitValVecPr,lpc1) < minFitVal){
			minFitVal = gsl_vector_get(bestFitValVecPr,lpc1);
			bestRun = lpc1;
		}
	}
	
	//Get the standardized location for the best run
	for(lpc1 = 0; lpc1 < nDim; lpc1++){
		gsl_vector_set(bestRunBestLocationPr,lpc1,gsl_matrix_get(bestLocationVecPr, bestRun, lpc1));
	}
	
	//Get the unstandardized location for the best run (including phiI values)
	for(lpc1 = 0; lpc1 < nDim+Np; lpc1++){
		gsl_vector_set(bestRunRealCPr, lpc1, gsl_matrix_get(bestLocRealCPr, bestRun, lpc1));
	}

	//Store in .hdf5 file
	dscalar2hdf5(outFile, "nRuns", nRuns);
	gslmatrix2hdf5(outFile,"bestLocVals",bestLocationVecPr);
	gslmatrix2hdf5(outFile, "bestLocRealCVals", bestLocRealCPr);
	gslvector2hdf5(outFile, "fitnessVals", bestFitValVecPr);
	gslvector2hdf5(outFile, "numFitEvals", nFuncEvalsVecPr);
	gslvector2hdf5(outFile, "numIter", nIterVecPr);
	gslvector2hdf5(outFile, "time2complete", wallClkTimeVecPr);
	dscalar2hdf5(outFile, "bestRun", bestRun+1);
	gslvector2hdf5(outFile, "bestLocation", bestRunBestLocationPr);	
	gslvector2hdf5(outFile, "bestRealLoc", bestRunRealCPr);

    //Wrap up
	gsl_matrix_free(bestLocationVecPr);
	gsl_matrix_free(bestLocRealCPr);
	gsl_vector_free(bestFitValVecPr);
	gsl_vector_free(nFuncEvalsVecPr);
	gsl_vector_free(nIterVecPr);
	gsl_vector_free(wallClkTimeVecPr);
	gsl_vector_free(bestRunBestLocationPr);
	gsl_vector_free(bestRunRealCPr);
}


/*! Load data from .hdf5 file into the special parameter structure
for LLR_PSO fitness function. Once used, the output should be 
destroyed using llrparam_free. */
struct  llr_pso_params * loadfile2llrparam(hid_t inFile){
	
	double **s;
	
	gsl_vector *yr_Pr = hdf52gslvector(inFile,"yr");
	
	gsl_matrix *trPr = hdf52gslmatrix(inFile,"timingResiduals");

	gsl_vector *sd_Pr = hdf52gslvector(inFile,"sd");
	
	gsl_vector *alphaP_Pr = hdf52gslvector(inFile,"alphaP");
		
	gsl_vector *deltaP_Pr = hdf52gslvector(inFile,"deltaP");
	
	/* Load fitness function parameter structure */
	size_t lpc1, lpc2;
	//struct llr_pso_params *llp = (struct llr_pso_params *)malloc(sizeof(struct llr_pso_params)); 
	size_t Np = (size_t)hdf52dscalar(inFile,"Np");
	size_t N = (size_t)hdf52dscalar(inFile,"N");
	struct llr_pso_params *llp = llrparam_alloc((unsigned int) N, (unsigned int) Np);
	// llp->Np = Np;
	// llp->N = N;
	// llp->sd = (double *)malloc(Np*sizeof(double));
	// llp->alphaP = (double *)malloc(Np*sizeof(double));
	// llp->deltaP = (double *)malloc(Np*sizeof(double));
	// llp->phiI = (double *)malloc(Np*sizeof(double));
	for (lpc1 = 0; lpc1 < Np; lpc1++){
		llp->sd[lpc1] = gsl_vector_get(sd_Pr,lpc1);	
		llp->alphaP[lpc1] = gsl_vector_get(alphaP_Pr,lpc1);
		llp->deltaP[lpc1] = gsl_vector_get(deltaP_Pr,lpc1);
	}
	// llp->yr = (double *) malloc(N*sizeof(double));
	for (lpc1 = 0; lpc1 < N; lpc1++){
		llp->yr[lpc1] = gsl_vector_get(yr_Pr,lpc1);
	}

	/* load timing residuals into fitness function param struct */
	// s = (double **)malloc(Np*sizeof(double *));
	//unsigned int countElements = 0;
	// for (lpc1 = 0; lpc1 < Np; lpc1++){
	// 	s[lpc1] = (double *) malloc(N*sizeof(double));
	// }
	for (lpc2 = 0; lpc2 < Np; lpc2++){
		for (lpc1 = 0; lpc1 < N; lpc1++){
			llp->s[lpc2][lpc1] = gsl_matrix_get(trPr, lpc2, lpc1);
		}
	}
	// llp->s = s;
	
	//Wrap up
	gsl_vector_free(yr_Pr); 
	gsl_matrix_free(trPr); 
	gsl_vector_free(sd_Pr);
	gsl_vector_free(alphaP_Pr);
	gsl_vector_free(deltaP_Pr);
	
	return llp;
}
