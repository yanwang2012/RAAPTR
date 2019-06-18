#include "maxphase.h"
#include "LLR_Mp_Av.h"
#include "ptapso.h"
#include "perfeval_omp.h"
#include "gslhdf5_io.h"
#include "hdf5_hl.h"
#include "omp.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

/*! \file
\brief Convert amplitude to SNR and return the timing residuals of estimated source.
\author Yiqian Qian
*/
/*!
This function reads the amplitude stored in estimated file and convert it back
to SNR, mean while returns the timing residuals of the estimated source.
*/
void Amp2Snr(char *inputFileName)
{

  // load data from .hdf5 file
  herr_t status;
  hid_t inFile = H5Fopen(inputFileName, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (inFile < 0)
  {
    fprintf(stdout, "Error openning file /n");
    abort();
  }

  struct llr_pso_params *llp;
  llp = loadfile2llrparam(inFile);

  // Close file
  status = H5Fclose(inFile);
  if (status < 0)
  {
    fprintf(stdout, "Error closing file %s /n", inputFileName);
  }
};