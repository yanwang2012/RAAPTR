#include <stdio.h>
#include "ptapso.h"
#include <gsl/gsl_rng.h>

int main(){

	size_t lpc;
	struct particleInfo p[2];

	double rngNum;
	gsl_rng *rngGen = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(rngGen,2571971);

	size_t nDim = 3;

	for (lpc = 0; lpc < 2; lpc++){
		initPsoParticles(&p[lpc],nDim,rngGen);
		particleinfo_fwrite(stdout, &p[lpc]);
	}

	gsl_rng_free(rngGen);
	for (lpc = 0; lpc < 2; lpc++){
		particleinfo_free(&p[lpc]);
	}

}