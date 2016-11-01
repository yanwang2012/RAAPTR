#if !defined(PTAPSODBGHDR)
#define PTAPSODBGHDR
#include "ptapso.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>

/* PSO parameter structure for debugging*/
struct psoParamStruct_debug {
	size_t popsize;
	size_t maxSteps;
	double c1;
	double c2;
	double max_velocity;
	double dcLaw_a;
	double dcLaw_b;
	double dcLaw_c;
	double dcLaw_d;
	/* 
	 Max number of iterations of 
	 local minimizer. Set to 0 to
	 switch off local minimization of
	 gbest.
	*/
	size_t locMinIter;
	/* Initial step size of the GSL
	   local minimizer (in standardized
	   coordinates) along each coordinate.
	*/
	double locMinStpSz;
	/* Pointer to GSL random number generator */
	gsl_rng *rngGen;
	/* Pointer to file if random numbers are to be
	   read from a file. If this is not 0, reading
	   random numbers from file will override 
	   generating random numbers on the fly. */
	FILE *rngFile;
	/* Pointer to file where to dump debug info */
	FILE *debugDumpFile;
};

/* ptapso function for debugging
*/
void ptapso_debug(size_t, /* Dimensionality of fitness function */
            double (*)(const gsl_vector *, void *), /* Pointer to fitness function */
		    void *, /* Fitness function parameter structure */
			struct psoParamStruct_debug, /* PSO parameters */
			struct returnData * /* Structure containing PSO output */
		   );
			
/* Function to initialize particle info */
void initPsoParticles_debug(struct particleInfo *, size_t , gsl_rng *, FILE *);

#endif