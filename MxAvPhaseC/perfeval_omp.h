#if !defined(PERFEVAL_OMP_HDR)
#define PERFEVAL_OMP_HDR

#include "maxphase.h"
#include "LLR_Mp_Av.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "ptapso.h"
#include "gslhdf5_io.h"
#include "hdf5.h"

/*! \file 
\brief Header file for perfeval_omp().
*/

void perfeval_omp(struct fitFuncParams *, char *, char *, char *);

/* Load special fitness function parameters from file*/
struct  llr_pso_params * loadfile2llrparam(hid_t);

/* Dump output from multiple pso runs in perfeval to .hdf5 file */
void perfevalomp2hdf5file(size_t, size_t, size_t, gsl_vector *,
					  struct returnData *[], struct fitFuncParams *,
                      double (*)(gsl_vector *, void *), hid_t );
					  

#endif
