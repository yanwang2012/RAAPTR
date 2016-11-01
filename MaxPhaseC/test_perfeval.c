#include "perfeval.h"
#include "maxphase.h"
#include "mp2matfile_io.h"
#include <stdio.h>

/*Test suite for perfeval.c function
Compile using: 
make (To compile a bunch of functions)
gcc -c -g test_perfeval.c -I/usr/local/include 
	-I/Applications/MATLAB_R2016a.app/extern/include
gcc test_perfeval.o perfeval.o LLR_PSO.o  maxphaseutils.o ptapso.o mp2matfile_io.o  
    -L/usr/local/lib -lgsl -lgslcblas
	-L/Applications/MATLAB_R2016a.app/bin/maci64 -lmat -lmx 
	-Xlinker -rpath -Xlinker /Applications/MATLAB_R2016a.app/bin/maci64 
	-lm -o test_perfeval.out
*/
int main(int argc, char *argv[]){
	/*Supply searchparams file and input data file 
	  paths on the command line*/
	
	/* Full path to search parameter file */
	char *srchParamsFile = argv[1];
	/* Full path to directory containing input files */
	char *inputFileName = argv[2];
	
	char outputFileName[] = "test_perfeval_out.mat";
	
	/* Error handling off */
	gsl_error_handler_t *old_handler = gsl_set_error_handler_off ();
	
	//Obtain fitness function parameter struct
	struct fitFuncParams *ffp;
	ffp = file2ffparam(srchParamsFile);
	
	//Run a loop in order to detect memory leak
	size_t lpc;
	for(lpc = 0; lpc < 5; lpc++){
		printf("Loop # %zu\n",lpc);
		perfeval(ffp,inputFileName, outputFileName);
	}
	
	//Wrap up
	ffparam_free(ffp);

}
