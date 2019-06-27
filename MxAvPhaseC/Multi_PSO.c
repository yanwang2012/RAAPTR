#include "postprocess.h"
#include "hdf5_hl.h"
#include "gslhdf5_io.h"
#include "perfeval_omp.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
/*! \file
/brief Post processing program for estimated output files from perfeval_omp().

\author Yiqian Qian
 */
/*!
This function will call perfeval_omp() as many times as user want or meet the threshold set by user to subtract timing residual of 
estimated source recursively untill meets the user's demand.
 */
int main(int argc, char *argv[]){
    // will be in a for loop
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

    perfeval_omp(ffp, inputFileName, outputFileName, mp_av_select);

	/* ----------------------------
	Deallocate storage
	-----------------------------*/
	ffparam_free(ffp);

    // new part should be done
    /* This part should gather the information fo estimated source and calculate its corresponding timing residual
    and subtract it from the simulation data, then creat a new input data file pass to PSO(upper part). 
     */
   // char *whatisthis = OutFile[0];
    char *Filename = argv[1];
    printf("Output file is: %s \n", Filename);
   // fprintf(stdout,"What is adress 0 %s",whatisthis);
   struct llo_pso_params * llp;
   llp = Amp2Snr(Filename);
   printf("Year is %10.5f",llp->yr);

   
   /* Everyting excuted successfully */
    return 0;
}

struct llr_pso_params * Amp2Snr(char *inputFileName){

    // load the .hdf5 file
    herr_t status;
    hid_t inFile = H5Fopen(inputFileName, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (inFile < 0)
    {
        fprintf(stdout, "Error opening file\n");
        abort();
    }
    //Remaining data loads into special parameter structure
    size_t Np = (size_t)hdf52dscalar(inFile,"Np");
	size_t N = (size_t)hdf52dscalar(inFile,"N");
    struct llr_pso_params *llp = llrparam_alloc((unsigned int) N, (unsigned int) Np);
    llp = loadfile2llrparam(inFile);

    //Close file
    status = H5Fclose(inFile);
    if (status < 0)
    {
        fprintf(stdout, "Error closing file %s \n", inputFileName);
    }

    //Wrap up
	gsl_vector_free(yr_Pr);
	gsl_matrix_free(trPr);
	gsl_vector_free(sd_Pr);
	gsl_vector_free(alphaP_Pr);
	gsl_vector_free(deltaP_Pr);

    return llp;
}

struct fitFuncParams * file2ffparam(char *srchParamsFile){

	herr_t status;
	hid_t srchPar = H5Fopen(srchParamsFile, H5F_ACC_RDONLY, H5P_DEFAULT);
	if (srchPar < 0){
		fprintf(stdout,"Error opening file %s\n", srchParamsFile);
		abort();
	}

	/* Read xmaxmin from srchParamsFile file */
	gsl_matrix *xmaxmin = hdf52gslmatrix(srchPar,"xmaxmin");
	/* Search Space dimensionality */
	size_t nDim = xmaxmin->size1;

	/* transfer xmaxmin to fitness function parameter struct */
	struct fitFuncParams *ffp = ffparam_alloc(nDim);
    size_t lpc1;
	for(lpc1 = 0; lpc1 < nDim; lpc1++){
		gsl_vector_set(ffp->rmin,lpc1,gsl_matrix_get(xmaxmin,lpc1,1));
		gsl_vector_set(ffp->rangeVec,lpc1,gsl_matrix_get(xmaxmin,lpc1,0)-gsl_matrix_get(xmaxmin,lpc1,1));
		//fprintf(stdout,"%f %f\n",gsl_vector_get(ffp->rmin,lpc1), gsl_vector_get(ffp->rangeVec,lpc1));
	}
	/* Close file */
	status = H5Fclose(srchPar);
	if(status < 0){
		fprintf(stdout,"Error closing file %s \n", srchParamsFile);
	}

	gsl_matrix_free(xmaxmin);

	return ffp;
}