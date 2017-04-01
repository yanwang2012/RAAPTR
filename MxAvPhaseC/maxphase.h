/*! \file maxphase.h
\brief Main header file for \ref maxphaseutils.c. 

\author Soumya D. Mohanty
*/
#if !defined(MAXPHASEHDR)
#define MAXPHASEHDR
#include <gsl/gsl_vector.h>

/*! 
A fitness function must use this
struct to ferry in parameter values needed to compute the fitness value, and
ferry out additional outputs, if any, from the fitness function besides the fitness value.
A fitness function can use the parameters supplied in this structure to translate
standardized input coordinates, each of which belongs to [0,1], into their real values.
*/
struct fitFuncParams{
size_t nDim;    /*!< Dimensionality of the fitness function. */
gsl_vector *rmin;     /*!< Minimum value of each coordinate. */
gsl_vector *rangeVec; /*!< Range of each coordinate. */
gsl_vector *realCoord;/*!< The unstandardized value of each coordinate is returned in this vector.*/
unsigned char fitEvalFlag; /*!< Set to 0 if fitness is infinity, else to 1.*/
/*! Pointer to a struct that carries additional parameters that are 
specific to a fitness function.
The header file of a fitness function can define this struct.
*/
void *splParams;
};

struct fitFuncParams * ffparam_alloc(size_t );
void ffparam_free(struct fitFuncParams *);

size_t chkstdsrchrng(const gsl_vector *);

void s2rvector(const gsl_vector *, const gsl_vector *, const gsl_vector *, gsl_vector *);

void limitVecComponent(gsl_vector *, double, double);

char ** listfileswext (const char *, const char *, size_t *, size_t *);

#endif