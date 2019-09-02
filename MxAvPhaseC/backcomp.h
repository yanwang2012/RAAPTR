#if !defined(BACKCOMP_HDR)
#define BACKCOMP_HDR

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <stddef.h>
#include <stdio.h>

gsl_matrix * raaptr_vec2mat(gsl_vector *, size_t);

#endif