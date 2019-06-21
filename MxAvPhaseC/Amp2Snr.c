#include "hdf5.h"
#include "gslhdf5_io.h"
#include "perfeval_omp.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
/*! \file
/brief Subtract timing residuals of estimated source from input data and creat a new input data file.

\author Yiqian Qian
 */
/*!
This function reads the estimated file and calculate the corresponding SNR and timing residual of the estimated source. 
Then subtract the timing residual from it and creat a new input data file.
 */
void Amp2Snr(char *inputFileName)
{

    // load the .hdf5 file
    herr_t status;
    hid_t inFile = H5Fopen(inputFileName, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (inFile < 0)
    {
        fprintf(stdout, "Error opening file\n");
        abort();
    }
    //Remaining data loads into special parameter structure
    struct llr_pso_params *llp;
    llp = loadfile2llrparam(inFile);

    //Close file
    status = H5Fclose(inFile);
    if (status < 0)
    {
        fprintf(stdout, "Error closing file %s \n", inputFileName);
    }
}