#include "subtract.h"
//#include "maxphase.h"
//#include "LLR_Mp_Av.h"
//#include "ptapso.h"
#include "hdf5.h"
#include "gslhdf5_io.h"
#include <gsl/gsl_vector.h>
#include <sys/stat.h>
#include <stdio.h>
#include <string.h>

int main(int argc, char *argv[]){
	char *outputFileName = argv[0];
	struct estSrcParams * srcp;
	srcp = estSrcParams(outputFileName);
	printParam(srcp);
    srcpara_free(srcp);
	return 0;
}