#include "hdf5_hl.h"
#include "gslhdf5_io.h"
#include "perfeval_omp.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
/*! \file
\brief Header file for Amp2Snr().
*/
struct estSrcParams{
    double alpha;
    double delta;
    double omega;
    double phi0;
    double Amp;
    double iota;
    double thetaN;
};

struct estSrcParams * file2Srcparam(char *);
struct fitFuncParams * file2ffparam(char *);