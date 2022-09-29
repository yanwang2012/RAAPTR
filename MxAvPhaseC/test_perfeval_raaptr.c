#include "loadRAAPTR.h"
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

struct fitFuncParams *file2ffparam(char *); // decleration

int main(int argc, char *argv[])
{
	/* General purpose variables */
	size_t lpc1, lpc2, lpc3;
	size_t Np, N;
	double **tres;
	if (argc != 6)
	{
		fprintf(stdout, "Usage: %s parameter_file_path input_file_path output_file_path mp_av_select (maxPhase/avPhase) pulsar_catalog\n", argv[0]);
		return 1;
	}
	/* Full path to search parameter file */
	char *srchParamsFile = argv[1];
	/* Full path to input file */
	// printf("First argument is:%s\n",srchParamsFile);

	char *inputFileName = argv[2];
	/* Path to output file */
	char *outputFileName = argv[3];
	// size_t length = strlen(outputFileName);
	// printf("Length of outputFileName is %d\n",length);

	/* Which algorithm to use */
	char *mp_av_select = argv[4];

	/* Pulsar catalog */
	char *psrfile = argv[5];
	// printf("load pulsar file %s\n", psrfile);

	/* Error handling off */
	gsl_error_handler_t *old_handler = gsl_set_error_handler_off();

	// Obtain fitness function parameter struct
	struct fitFuncParams *ffp;
	ffp = file2ffparam(srchParamsFile);

	/* Analyze input file*/
	// fprintf(stdout,"inputFileNameis: %s \n", inputFileName);
	// fprintf(stdout,"outputFileName is: %s \n", outputFileName);
	fprintf(stdout, "Analyzing file %s \n", inputFileName);
	fprintf(stdout, "Output will be stored in %s\n", outputFileName);
	fprintf(stdout, "******************************************\n");
	/*Main function. Will be called several times according to user specified requirement.*/
	perfeval_omp_RAAPTR(ffp, inputFileName, outputFileName, mp_av_select, psrfile);

	/* ----------------------------
			Deallocate storage
		 -----------------------------*/
	ffparam_free(ffp);

	/* Everything executed successfully */
	fprintf(stdout, "All Done!\n");
	return 0;
}

struct fitFuncParams *file2ffparam(char *srchParamsFile)
{

	herr_t status;
	hid_t srchPar = H5Fopen(srchParamsFile, H5F_ACC_RDONLY, H5P_DEFAULT);
	if (srchPar < 0)
	{
		fprintf(stdout, "Error opening file %s\n", srchParamsFile);
		abort();
	}

	/* Read xmaxmin from srchParamsFile file */
	gsl_matrix *xmaxmin = hdf52gslmatrix(srchPar, "xmaxmin");
	/* Search Space dimensionality */
	size_t nDim = xmaxmin->size1;

	/* transfer xmaxmin to fitness function parameter struct */
	struct fitFuncParams *ffp = ffparam_alloc(nDim);
	size_t lpc1;
	for (lpc1 = 0; lpc1 < nDim; lpc1++)
	{
		gsl_vector_set(ffp->rmin, lpc1, gsl_matrix_get(xmaxmin, lpc1, 1));
		gsl_vector_set(ffp->rangeVec, lpc1, gsl_matrix_get(xmaxmin, lpc1, 0) - gsl_matrix_get(xmaxmin, lpc1, 1));
		// fprintf(stdout,"%f %f\n",gsl_vector_get(ffp->rmin,lpc1), gsl_vector_get(ffp->rangeVec,lpc1));
	}
	/* Close file */
	status = H5Fclose(srchPar);
	if (status < 0)
	{
		fprintf(stdout, "Error closing file %s \n", srchParamsFile);
	}

	gsl_matrix_free(xmaxmin);

	return ffp;
}