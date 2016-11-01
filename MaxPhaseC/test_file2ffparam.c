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
int main(){
	/*Supply searchparams file and input data file 
	  paths on the command line*/
	
	/* Full path to search parameter file */
	char srchParamsFile[] = "searchParams_simDataX.mat";

	//Obtain fitness function parameter struct
	struct fitFuncParams *ffp;
	size_t lpc;
	for (lpc = 0; lpc < 500000; lpc++){
		ffp = file2ffparam(srchParamsFile);
		//Wrap up
		ffparam_free(ffp);
	}

}