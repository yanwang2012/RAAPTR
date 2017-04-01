#if !defined(PERFEVALHDR)
#define PERFEVALHDR

#include "maxphase.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "ptapso.h"
#include "mat.h"

/*! \file 
\brief Header file for \ref perfeval.c.
*/

void perfeval(struct fitFuncParams *, char *, char *);

/* Dump output from multiple pso runs in perfeval to .mat file */
void perfevaloutput2matfile(size_t, size_t, size_t, gsl_vector *,
					  struct returnData *[], struct fitFuncParams *,
                      double (*)(gsl_vector *, void *), MATFile *);

#endif