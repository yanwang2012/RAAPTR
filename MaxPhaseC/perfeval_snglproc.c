#include "maxphase.h"
#include "mp2matfile_io.h"
#include "perfeval.h"
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>

int main(int argc, char *argv[]){
	
	if (argc != 3){
		printf("Usage: %s parameter_file_path analysis_directory\n", argv[0]);
		return 1;
	}
	
	/* Error handling off */
	gsl_error_handler_t *old_handler = gsl_set_error_handler_off ();
	
	/* Full path to search parameter file */
	char *srchParamsFile = argv[1];
	/* Full path to directory containing input files */
	char *simDataDir = argv[2];

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
		free(outDataDir);
		return 1;
	}
	
	//Obtain list of input file names
	/* Create list of input data files */
	char **inputFileList;
	size_t nInputFiles, maxFileNameLen;
	inputFileList = listfileswext("mat", simDataDir, &nInputFiles, &maxFileNameLen);
	if (inputFileList == NULL){
		printf("Invalid input file list\n");
		free(outDataDir);
		return 1;
	}
	printf("Number of input files %zu\n",nInputFiles);

	/* Sufficient storage for input and output file name strings: 
	directory name length +2 (for '/0' and '/' (filesep)).
	*/
	char *inputFileName, *outputFileName;
	outputFileName = (char *)malloc((strlen(outDataDir)+maxFileNameLen+2)*sizeof(char));
	inputFileName = (char *)malloc((strlen(simDataDir)+maxFileNameLen+2)*sizeof(char));	
	
	//Obtain fitness function parameter struct
	struct fitFuncParams *ffp;
	ffp = file2ffparam(srchParamsFile);
	
	//Process input files
	size_t lpc1;
	for (lpc1 = 0; lpc1 < nInputFiles; lpc1++){
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
	
	// Wrap up
	ffparam_free(ffp);
	for(lpc1 = 0; lpc1 < nInputFiles; lpc1++){
	    /* Free up dynamically allocated memory for file list*/
		free(inputFileList[lpc1]);
	}
	/* Free up dynamically allocated memory */
	free(inputFileList);
	free(outDataDir);
	free(outputFileName);
	free(inputFileName);
	
}