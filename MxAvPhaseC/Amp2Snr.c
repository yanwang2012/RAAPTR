#include "Amp2Snr.h"
#include "hdf5_hl.h"
#include "gslhdf5_io.h"
#include "perfeval_omp.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
/*! \file
/brief Subtract timing residual of the estimated source from output data and creat a new input data file.

\author Yiqian Qian
 */
/*!
This function calculates the estimated timing residual of the estimated source and 
Then subtract the timing residual from it to creat a new input data file.
 */
int main(int argc, char *OutFile[]){
   // char *whatisthis = OutFile[0];
    char *Filename = OutFile[1];
    printf("Output file is: %s \n", Filename);
   // fprintf(stdout,"What is adress 0 %s",whatisthis);
   Amp2Snr(Filename);
    return 0;
}

void Amp2Snr(char *inputFileName){

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