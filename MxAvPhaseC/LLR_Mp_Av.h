/*! \file LLR_Mp_Av.h
\brief Header file for \ref LLR_Mp_Av.c.
*/
#if !defined(LLR_MpAv_HDR)
#define LLR_MpAv_HDR

#include "maxphase.h"
#include <gsl/gsl_vector.h>

/*!
Special parameters struct to pass parameters specific to LLR_Mp_Av fitness function.
See how_to_code_fitnessFunc.txt for an example of a special parameter structure.
*/
struct llr_pso_params{
	/*! number of pulsars in the timing array */
	unsigned int Np;
	/*! number of observations */
	unsigned int N;
	/*! s[2]; (Np,N), GW signal */
	double **s;
	/*! (Np,1), standard deviation of noise for different pulsar */
	double *sd;
	/*! (Np,1), right ascension, in radian */
	double *alphaP;
	/*!  (Np,1), declination, in radian */
	double *deltaP;
	/*! (N,1), observation epoch, in year */
	double *yr;
	/*! Estimated phases (1, Np) */
	double *phiI;
};

/*! MISSING DOCUMENTATION */
struct cfunc_OUTPUT
{
//double fitVal;
double * c;
double * v;
};

/*! MISSING DOCUMENTATION */
struct lh_OUTPUT
{
double * phiI;
double * lhI;
};

// Creation of a new deep copy of the fitness function parameter struct
struct fitFuncParams * ffparams_clone(struct fitFuncParams *);

//-------AvPhase specific functions and structs--------
/* Structure for avPhase */
struct avPhase_param
{
	double p0; double p1; double p2; double p3; double p4; // p0-p4 is b
	double p5;  // p5 is norm
};

double f(double , void *);

double LLR_av(gsl_vector *, /*!< Data vector */
               void * /*!< Pointer to fitFuncParams struct */
			  );

//double LogLikelihoodRatioMP5(struct fitFuncParams *);
double AvPhaseLLR(struct fitFuncParams *);
//-----------------------------------------------------

double InnProduct(unsigned int N, double *, double *, double sd);


double LLR_mp(gsl_vector *, /*!< Data vector */
               void * /*!< Pointer to fitFuncParams struct */
			  );


double LogLikelihoodRatioMP5(struct fitFuncParams *);


void likelihood(double *, double *, struct lh_OUTPUT *);


void cfunc(unsigned int, double , double , double,  double, double,
	         double , double , double , double ,
				   double , double *, double *, double *, struct cfunc_OUTPUT *);
				   

struct llr_pso_params * llrparam_alloc(unsigned int, unsigned int);
				   
void llrparam_free(struct llr_pso_params *);


size_t  llrpsochkcoord(const size_t,  const gsl_vector *, const gsl_vector *, gsl_vector *,
                       gsl_vector *);
					   
double wraphalfcircangle(double);

#endif
