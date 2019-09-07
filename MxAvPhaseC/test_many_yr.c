#include "backcomp.h"
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
\brief Test harness for code that allows different observational epochs for different pulsars. 

\author Soumya D. Mohanty

# Usage
test_many_yr.out search_param_file input_data_file output_data_file mp_av_select
- search_param_file: Full path to a .hdf5 file containing the parameters for the run.
- input_data_file: Full path to the  the .hdf5 input file to analysis.
- output_data_file: Full path to the output file that will be created.
- mp_av_select: 'maxPhase' or 'avPhase' algorithm 

## Format of search_param_file
This is a .hdf5 file. It should contain a dataset called 'xmaxmin'.
This is a two column matrix. Each row contains the minimum and maximum values, in that order,
defining the search interval along a particular parameter for PSO. 

## Format of input data file
See the documentation for the simulation data generation codes in RAAPTR/GENSIMDATA.

*/
int main(int argc, char *argv[]){
	/* General purpose variables */
	size_t lpc1, lpc2, lpc3;
    int stat;
	
	if (argc != 5){
		fprintf(stdout,"Usage: %s parameter_file_path input_file_path output_file_path mp_av_select\n", argv[0]);
		return 1;
	}
	/* Full path to search parameter file */
	char *srchParamsFile = argv[1];
	/* Full path to input file */
	char *inputFileName = argv[2];	
	
	/* Path to output file */
	char *outputFileName = argv[3];
	
	/* Which algorithm to use */
	char *mp_av_select = argv[4];
	
	/* Error handling off */
	gsl_error_handler_t *old_handler = gsl_set_error_handler_off ();
	
	//Obtain fitness function parameter struct
	struct fitFuncParams *ffp;
	ffp = file2ffparam(srchParamsFile);

	/* Analyze input file*/
	fprintf(stdout, "Analyzing file %s \n",inputFileName);
	fprintf(stdout, "Output will be stored in %s\n",outputFileName);
	fprintf(stdout,"******************************************\n");	
		
    /* Number of independent PSO runs */
	size_t nRuns = 2;
	
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
	 else{
		 printf("Option %s is not recognized. Use maxPhase or avPhase.\n",mp_av_select);
		 return 1;
	 }
	
	/* Choose the random number generation method.
	   Need independent random number generators 
	   inside an OMP for loop (?)*/
	gsl_rng *rngGen[nRuns];
	for(lpc1 = 0; lpc1 < nRuns; lpc1++){
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
		abort();
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
	for (lpc1 = 0; lpc1 < nRuns; lpc1++){
		//fprintf(stdout,"Running PSO run %zu on worker %d\n",lpc1,omp_get_thread_num());
		/* Initialize random number generator separately for each worker*/
		gsl_rng_set(rngGen[lpc1],rngSeeds[lpc1]);
	    /* Configure PSO parameters*/
		psoParams.popsize= 4;
		psoParams.maxSteps= 20; 
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
        //printf("Starting PSO run %d\n",lpc1);
		psoStartTime = clock();
		ptapso(nDim, fitfunc, ffpCopy, &psoParams, psoResults[lpc1]);
        psoStopTime = clock();
	
		gsl_vector_set(wallClkTimes,lpc1, (((double)(psoStopTime - psoStartTime))/CLOCKS_PER_SEC)/60.0);
		
        //printf("Wrapping up PSO run %d\n",lpc1);
		//Delete local copy of fitness parameter struct
		llrparam_free(ffpCopy->splParams);
		ffparam_free(ffpCopy);		
        //printf("Completed pso run %d\n",lpc1);
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
	for (lpc1 = 0; lpc1 < nRuns; lpc1++){
		gsl_rng_free(rngGen[lpc1]);		
	}
	llrparam_free(llp);
	for (lpc1 = 0; lpc1 < nRuns; lpc1++){
		returnData_free(psoResults[lpc1]);
	} 
	gsl_vector_free(wallClkTimes);

	/* ----------------------------
	Deallocate storage 
	-----------------------------*/
	ffparam_free(ffp);
	
	/* Everything executed successfully */
	return 0;
}