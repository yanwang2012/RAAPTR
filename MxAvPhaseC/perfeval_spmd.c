#include "maxphase.h"
#include "LLR_Mp_Av.h"
#include "ptapso.h"
#include "perfeval_omp.h"
#include "hdf5.h"
#include "gslhdf5_io.h"
#include <gsl/gsl_vector.h>
#include <sys/stat.h>
#include <stdio.h>
#include <string.h>

/*! \file 
\brief Run \ref perfeval_omp.c on an input data file. 

\author Soumya D. Mohanty

# Usage
perfeval_spmd.out search_param_file input_data_dir
- search_param_file: Full path to a .hdf5 file containing the parameters for the run.
- input_file: Full path to the  the .hdf5 input file to analysis.
- output_file: Full path to the output file that will be created.
- mp_av_select: 'maxPhase' or 'avPhase' algorithm 

## Format of search_param_file
This is a .hdf5 file. It should contain a dataset called 'xmaxmin'.
This is a two column matrix. Each row contains the minimum and maximum values, in that order,
defining the search interval along a particular parameter for PSO. 

## Format of input data file
See the documentation for the simulation data generation code.
*/
struct fitFuncParams * file2ffparam(char *);

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
		
    perfeval_omp(ffp, inputFileName, outputFileName, mp_av_select);

	/* ----------------------------
	Deallocate storage 
	-----------------------------*/
	ffparam_free(ffp);
	
	/* Everything executed successfully */
	return 0;
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