#include "maxphase.h"
#include "LLR_PSO.h"
#include "ptapso.h"
#include "mp2matfile_io.h"
#include "perfeval.h"
#include "mat.h"
#include "matrix.h"
#include <gsl/gsl_vector.h>
#include <sys/stat.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>

/*! \file 
\brief Run \ref perfeval on a set of input data files in parallel.

This is the main code that produces the simulation results in
the papers reporting analyses based on the MaxPhase algorithm (Wang, Mohanty, Jenet, ApJ, 2015).

\author Soumya D. Mohanty

# Usage
mpirun -n number_of_jobs perfeval_mpi.out search_param_file input_data_dir mpi_job_file
- search_param_file: Full path to a .mat file containing the parameters for the run.
- input_data_dir: Full path to the directory containing the input files to analysis.
- mpi_job_file: Full path to an ascii file containing information about how to split the 
input files among the workers.

## Format of search_param_file
This is a .mat file (Matlab file format). It should contain a variable called 'xmaxmin'.
This is a two column matrix. Each row contains the minimum and maximum values, in that order,
defining the search interval along a particular parameter for PSO.

## Format of mpi_job_file
This is a two column text file that lists the starting and ending file numbers to assign to 
an MPI worker. Use the Matlab script 'nfilesperworker.m' to divide a given total number
of files among a given number of workers as equitably as possible. 

## Format of input data files
See the documentation for the simulation data generation code.
*/
int main(int argc, char *argv[]){
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
	/*-----------------------------------------------------*/
	/* General purpose variables */
	size_t lpc1, lpc2, lpc3;
    int stat;
	
	if (argc != 4){
		printf("Usage: %s parameter_file_path analysis_directory mpi_job_specification_file_path \n", argv[0]);
		return 1;
	}
	/* Full path to search parameter file */
	char *srchParamsFile = argv[1];
	/* Full path to directory containing input files */
	char *simDataDir = argv[2];
	/* Full path to file containing start/stop file number pairs for each MPI worker */
	char *mpiJobFileName = argv[3];
	
	/* Read start and stop indices of files from MPI job specification file */
	FILE *mpiJbFile = fopen(mpiJobFileName,"r");
	size_t countlines = 0;
	unsigned int startFileIndx, stopFileIndx;
	if(mpiJbFile != NULL){
		while(fscanf(mpiJbFile,"%u %u",&startFileIndx,&stopFileIndx)==2){
			countlines++;
		}
    }
	rewind(mpiJbFile);
	gsl_vector_uint *startFile = gsl_vector_uint_alloc(countlines);
	gsl_vector_uint *stopFile = gsl_vector_uint_alloc(countlines);
	for(lpc1 = 0; lpc1 < countlines; lpc1++){
		fscanf(mpiJbFile,"%u %u", &startFileIndx, &stopFileIndx);
		//Job file uses indices starting from 1 to label files. Remove 1.
		gsl_vector_uint_set(startFile,lpc1,startFileIndx-1);
		gsl_vector_uint_set(stopFile,lpc1,stopFileIndx-1);
	}
	fclose(mpiJbFile);
	
	
	/* Path to folder for storing outputs */
	char *outDataDir = (char *)malloc( (strlen(simDataDir)+ 1 +
		                                strlen("results_")+
									    world_size + 1)*sizeof(char) );
	strcpy(outDataDir,simDataDir);
	strcat(outDataDir,"/");
	strcat(outDataDir,"results_");
	char outDataDirNum[world_size+1];
	sprintf(outDataDirNum,"%d",world_rank);
	strcat(outDataDir,outDataDirNum);
	if (mkdir(outDataDir,S_IRWXU)){
		printf("Error creating results directory %s\n", outDataDir);
		printf("Check that it does not already exist.\n");
        gsl_vector_uint_free(startFile);
		gsl_vector_uint_free(stopFile);
		free(outDataDir);
		MPI_Finalize();
		return 4;
	}
	
	/* Error handling off */
	gsl_error_handler_t *old_handler = gsl_set_error_handler_off ();
	
	//Obtain list of input file names
	/* Create list of input data files */
	char **inputFileList;
	size_t nInputFiles, maxFileNameLen;
	inputFileList = listfileswext("mat", simDataDir, &nInputFiles, &maxFileNameLen);
	printf("Analyzing %zu files across %d workers\n", nInputFiles, world_size);
	if (inputFileList == NULL){
		printf("Invalid input file list\n");
		free(outDataDir);
		return 1;
	}
	
	/* Sufficient storage for input and output file name strings: 
	directory name length +2 (for '/0' and '/' (filesep)).
	*/
	char *inputFileName, *outputFileName;
	outputFileName = (char *)malloc((strlen(outDataDir)+maxFileNameLen+2)*sizeof(char));
	inputFileName = (char *)malloc((strlen(simDataDir)+maxFileNameLen+2)*sizeof(char));	
	
	//Obtain fitness function parameter struct
	struct fitFuncParams *ffp;
	ffp = file2ffparam(srchParamsFile);

	/* Analyze input files in parallel */
	size_t startFileNum, endFileNum;
	startFileNum = (size_t)gsl_vector_uint_get(startFile,world_rank);
	endFileNum = (size_t)gsl_vector_uint_get(stopFile,world_rank);
	printf("Analyzing files %zu to %zu on worker %d\n",startFileNum+1,endFileNum+1,world_rank);
	printf("Output will be stored in directory %s\n",outDataDir);
	printf("******************************************\n");
	for (lpc1 = startFileNum; lpc1 <= endFileNum; lpc1++){
	    /* Construct input file name along with path */
		strcpy(inputFileName,simDataDir);
		strcat(inputFileName,"/");
		strcat(inputFileName,inputFileList[lpc1]);
		
		/* Construct output file name along with path */
		strcpy(outputFileName,outDataDir);
		strcat(outputFileName,"/");
		strcat(outputFileName,inputFileList[lpc1]);
		
		printf("Input file %s\n",inputFileName);
		
		perfeval(ffp, inputFileName, outputFileName);
	}

	/* ----------------------------
	Deallocate storage 
	-----------------------------*/
	ffparam_free(ffp);
	for(lpc1 = 0; lpc1 < nInputFiles; lpc1++){
	    /* Free up dynamically allocated memory for file list*/
		free(inputFileList[lpc1]);
	}
	/* Free up dynamically allocated memory */
	free(inputFileList);
	free(outputFileName);
	free(inputFileName);
	free(outDataDir);
	gsl_vector_uint_free(startFile);
	gsl_vector_uint_free(stopFile);
	/*----------------------------------------------------*/
    // Finalize the MPI environment.
    MPI_Finalize();
	
	/* Everything executed successfully */
	return 0;
}