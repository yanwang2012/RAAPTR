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
\brief Convert amplitude to SNR and return the timing residuals of estimated source.
\author Yiqian Qian
*/
/*!
This function reads the amplitude stored in estimated file and convert it back
to SNR, mean while returns the timing residuals of the estimated source.
*/

struct srcParams {
  double Np;
  double N;
  double alpha;
  double delta;
  double omega;
  double phi0;
}srcparams;
