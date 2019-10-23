#include "subtract.h"
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
\brief Run perfeval_spmd() on an input data file for a specified 
requirement by user.

\author Soumya D. Mohanty & Yiqian Qian.

# Usage
Multi_PSO.exe search_param_file input_data_dir output_file mp_av_select number_of_iteration
- search_param_file: Full path to a .hdf5 file containing the parameters for the run.
- input_file: Full path to the  the .hdf5 input file to analysis.
- output_file: Full path to the output file that will be created.
- mp_av_select: 'maxPhase' or 'avPhase' algorithm
- number_of_iteration: An integer number describes how many times you want pso to run.
- threshold: lowest snr can approach.

## Format of search_param_file
This is a .hdf5 file. It should contain a dataset called 'xmaxmin'.
This is a two column matrix. Each row contains the minimum and maximum values, in that order,
defining the search interval along a particular parameter for PSO.

## Format of input data file
See the documentation for the simulation data generation code.
*/

struct fitFuncParams *file2ffparam(char *); //decleration
int main(int argc, char *argv[])
{
	/* General purpose variables */
	size_t lpc1, lpc2, lpc3;
	size_t Np, N;
	double **tres;
	if (argc != 6)
	{
		fprintf(stdout, "Usage: %s parameter_file_path input_file_path output_file_path mp_av_select number_of_iterations\n", argv[0]);
		return 1;
	}
	/* Full path to search parameter file */
	char *srchParamsFile = argv[1];
	/* Full path to input file */
	// printf("First argument is:%s\n",srchParamsFile);

	char *inputFileName = argv[2];
	/* Path to output file */
	char *outputFileName = argv[3];
	//size_t length = strlen(outputFileName);
	//printf("Length of outputFileName is %d\n",length);

	/* Which algorithm to use */
	char *mp_av_select = argv[4];

	/* Number of iterations */
	int num_ite = atoi(argv[5]); // transfer char to integer
	fprintf(stdout, "Number of iteration is: %d\n", num_ite);

	for (int ite = 0; ite < num_ite; ite++)
	{
		/* Multi PSO Process */

		/* Error handling off */
		gsl_error_handler_t *old_handler = gsl_set_error_handler_off();

		//Obtain fitness function parameter struct
		struct fitFuncParams *ffp;
		ffp = file2ffparam(srchParamsFile);

		/* Analyze input file*/
		fprintf(stdout, "Analyzing file %s \n", inputFileName);
		fprintf(stdout, "Output will be stored in %s\n", outputFileName);
		fprintf(stdout, "******************************************\n");
		/*Main function. Will be called several times according to user specified requirement.*/
		perfeval_omp(ffp, inputFileName, outputFileName, mp_av_select);

		/*------------------------------------------
	        Subtraction estimated timing residuals from source.
	        --------------------------------------------*/
		struct estSrcParams *srcp;
		struct llr_pso_params *llp;
		//Load data from specified .hdf5 input file
		srcp = file2Srcparam(outputFileName);
		fprintf(stdout, "Loading output file %s\n", outputFileName);

		herr_t status;
		hid_t inFile = H5Fopen(inputFileName, H5F_ACC_RDONLY, H5P_DEFAULT);
		if (inFile < 0)
		{
			fprintf(stdout, "Error opening file\n");
			abort();
		}

		llp = loadfile2llrparam(inFile);
		fprintf(stdout, "Loading input file %s\n", inputFileName);

		N = (size_t)llp->N;
		//printf("N: %zu\n",N);
		Np = (size_t)llp->Np;
		//printf("Np: %zu\n",Np);
		tres = llp->s;

		/*
            	FILE * fsrc;
            	fsrc = fopen("SrcRes.txt","w");
              for(int m = 0; m < Np; m++){
              	for(int n = 0; n < Np; n++){
              		fprintf(fsrc,"%e\t", tres[m][n]);
                 }
        		fprintf(fsrc,"\n");
               }
    			fclose(fsrc);
             */

		gsl_matrix *timResiduals = gsl_matrix_calloc(Np, N);
		gsl_matrix *estRes = gsl_matrix_calloc(Np, N);
		estRes = timingResiduals(srcp, llp);
		//printf("Dimension of timResiduals: %zu %zu\n", timResiduals->size1, timResiduals->size2);
		/*
    			FILE * fest;
    			fest = fopen("estRes.txt","w");
    			printMatrix(fest,timResiduals,Np,N);
    			fclose(fest);
			 */

		size_t i, j;
		for (i = 0; i < Np; i++)
		{
			for (j = 0; j < N; j++)
			{
				gsl_matrix_set(timResiduals, i, j, tres[i][j] - gsl_matrix_get(estRes, i, j));
			}
		}

		/* Put subtracted timing residuals into input file as the new input file. */
		/*double buffer[Np][N];

		size_t m, n;
		for (m = 0; m < Np; m++)
		{
			for (n = 0; n < N; n++)
			{
				buffer[m][n] = gsl_matrix_get(timResiduals, m, n);
			}
		}*/
		/*	hid_t dset_id = H5Dopen1(inFile, "timingResiduals"); // Open an existing dataset.
		status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
		// Close dataset.
		H5Dclose(dset_id);
		// close file.
		status = H5Fclose(inFile); 
		if (status < 0)
		{
			fprintf(stdout, "Error closing file: %s\n", inputFileName);
		}
	*/
		/* Create new input file.*/
		char purefilename[strlen(argv[2])];
		char newinputfile[strlen(purefilename) + strlen("_sub1.hdf5")];
		strncpy(purefilename, argv[2], strlen(argv[2]) - strlen(".hdf5"));
		purefilename[strlen(argv[2]) - strlen(".hdf5")] = '\0'; //null character manually added
		sprintf(newinputfile, "%s_sub%d.hdf5", purefilename, ite + 1);
		fprintf(stdout, "New input file is: %s\n", newinputfile);

		hid_t ninFile = H5Fcreate(newinputfile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		if (ninFile < 0)
		{
			printf("Error creating new input file %s \n", newinputfile);
		}
		/* ----Copy all the other parameters to new file.-------------*/
		H5Ocopy(inFile, "/Amp", ninFile, "/Amp", H5P_DEFAULT, H5P_DEFAULT);
		H5Ocopy(inFile, "/N", ninFile, "/N", H5P_DEFAULT, H5P_DEFAULT);
		H5Ocopy(inFile, "/Np", ninFile, "/Np", H5P_DEFAULT, H5P_DEFAULT);
		H5Ocopy(inFile, "/alpha", ninFile, "/alpha", H5P_DEFAULT, H5P_DEFAULT);
		H5Ocopy(inFile, "/alphaP", ninFile, "/alphaP", H5P_DEFAULT, H5P_DEFAULT);
		H5Ocopy(inFile, "/delta", ninFile, "/delta", H5P_DEFAULT, H5P_DEFAULT);
		H5Ocopy(inFile, "/deltaP", ninFile, "/deltaP", H5P_DEFAULT, H5P_DEFAULT);
		H5Ocopy(inFile, "/genHypothesis", ninFile, "/genHypothesis", H5P_DEFAULT, H5P_DEFAULT);
		H5Ocopy(inFile, "/iota", ninFile, "/iota", H5P_DEFAULT, H5P_DEFAULT);
		H5Ocopy(inFile, "/kp", ninFile, "/kp", H5P_DEFAULT, H5P_DEFAULT);
		H5Ocopy(inFile, "/loc_id", ninFile, "/loc_id", H5P_DEFAULT, H5P_DEFAULT);
		H5Ocopy(inFile, "/noise", ninFile, "/noise", H5P_DEFAULT, H5P_DEFAULT);
		H5Ocopy(inFile, "/omega", ninFile, "/omega", H5P_DEFAULT, H5P_DEFAULT);
		H5Ocopy(inFile, "/omg_id", ninFile, "/omg_id", H5P_DEFAULT, H5P_DEFAULT);
		H5Ocopy(inFile, "/perfect_fitness", ninFile, "/perfect_fitness", H5P_DEFAULT, H5P_DEFAULT);
		H5Ocopy(inFile, "/phi0", ninFile, "/phi0", H5P_DEFAULT, H5P_DEFAULT);
		H5Ocopy(inFile, "/rlz_id", ninFile, "/rlz_id", H5P_DEFAULT, H5P_DEFAULT);
		H5Ocopy(inFile, "/sd", ninFile, "/sd", H5P_DEFAULT, H5P_DEFAULT);
		H5Ocopy(inFile, "/snr", ninFile, "/snr", H5P_DEFAULT, H5P_DEFAULT);
		H5Ocopy(inFile, "/snr_chr", ninFile, "/snr_chr", H5P_DEFAULT, H5P_DEFAULT);
		H5Ocopy(inFile, "/snr_id", ninFile, "/snr_id", H5P_DEFAULT, H5P_DEFAULT);
		H5Ocopy(inFile, "/thetaN", ninFile, "/thetaN", H5P_DEFAULT, H5P_DEFAULT);
		//H5Ocopy(inFile, '/timingResiduals', ninFile, '/timingResiduals', H5P_DEFAULT, H5P_DEFAULT);
		H5Ocopy(inFile, "/timingResiduals_tmp", ninFile, "/timingResiduals_tmp", H5P_DEFAULT, H5P_DEFAULT);
		H5Ocopy(inFile, "/yr", ninFile, "/yr", H5P_DEFAULT, H5P_DEFAULT);
		/*-------------Hard coded------------------------------*/
		gslmatrix2hdf5(ninFile, "timingResiduals", timResiduals);
		status = H5Fclose(ninFile);
		if (status < 0)
		{
			printf("Error closing new input file %s \n", newinputfile);
		}
		status = H5Fclose(inFile);
		if (status < 0){
			printf("Error closing input file %s \n", inputFileName);
		}

		inputFileName = newinputfile;
		/*
   			 FILE * f;
   			 f = fopen("timingResiduals.txt", "w");
   			 printMatrix(f,timResiduals,Np,N);// print timing residual to file f. 
   			 fclose(f);
             */

		/* Creat new output file. */
		char purename[strlen(argv[3])];
		strncpy(purename, argv[3], strlen(argv[3]) - strlen(".hdf5")); // get the pure file name without extension.
		purename[strlen(argv[3]) - strlen(".hdf5")] = '\0';			   //null character manually added

		char newName[strlen(purename) + strlen("_0.hdf5")];

		sprintf(newName, "%s_%d.hdf5", purename, ite + 1);
		fprintf(stdout, "New outputFileName = %s\n", newName);
		fprintf(stdout, "******************************************\n");
		outputFileName = newName;

		/* ----------------------------
	        	Deallocate storage
	         -----------------------------*/

		srcpara_free(srcp);
		llrparam_free(llp);
		gsl_matrix_free(timResiduals);
		gsl_matrix_free(estRes);
		ffparam_free(ffp);
	}
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
		//fprintf(stdout,"%f %f\n",gsl_vector_get(ffp->rmin,lpc1), gsl_vector_get(ffp->rangeVec,lpc1));
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
