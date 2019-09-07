#include "maxphase.h"
#include "perfeval_omp.h"
#include <sys/stat.h>
#include <stdio.h>
#include <string.h>

/*! \file 
\brief Run perfeval_spmd() on an input data file. 

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


