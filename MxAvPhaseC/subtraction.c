#include "subtraction.h"
#include "maxphase.h"
#include "LLR_Mp_Av.h"
#include "ptapso.h"
#include "hdf5.h"
#include "gslhdf5_io.h"
#include <gsl/gsl_vector.h>
#include <sys/stat.h>
#include <stdio.h>
#include <string.h>
/*! \file
\brief Subtract estimated source out of simulated data.

\author Yiqian Qian.
*/
struct estSrcParams * file2Srcparam(char *outputFileName){
	herr_t status;
	hid_t SrcPar = H5Fopen(outputFileName, H5F_ACC_RDONLY, H5P_DEFAULT);
	if (SrcPar < 0)
	{
		printf("Error opening file %s\n", outputFileName);
	}
	struct estSrcParams *srcp; // need to be modified
	return srcp;
}
