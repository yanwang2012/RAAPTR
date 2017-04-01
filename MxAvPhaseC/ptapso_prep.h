#if !defined(PTAPSOHDR)
#define PTAPSOHDR

#include <gsl/gsl_vector.h>

/* PSO parameter structure */
struct psoParamStruct {
	unsigned int popsize;
	unsigned int maxSteps;
	double c1;
	double c2;
	double max_initial_velocity;
	double max_velocity;
	double dcLaw_a;
	double dcLaw_b;
	double dcLaw_c;
	double dcLaw_d;
	/* 
	nmOptions = optimset('fminsearch');
	nmOptions.MaxFunEvals=200;
	nmOptions.Display = 'off';
	*/
};

/* Structure returned by ptapso.c */
struct returnData {
	unsigned long int totalSteps;
    unsigned long int totalFuncEvals;
	/* Number of search dimensions */
	unsigned int nDim;
    double *bestLocation;
    double bestSNR;
};

/* ptapso function 
First input argument is pointer to fitness function. Second input argument is
pointer to a structure containing parameters needed by the fitness function. 
The parameter structure should be defined in a header file associated with the
fitness function. 
*/
void ptapso(double (*)(const gsl_vector *, void *), 
					void *, struct psoParamStruct, struct returnData *);

#endif