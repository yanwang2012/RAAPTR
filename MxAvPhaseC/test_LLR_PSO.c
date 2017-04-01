/*! \file test_ptapsotestfunc.c
\brief Example of how to call a fitness function.

 If calling fitness function foo.c, include foo.h. Here, we
are calling ptapsotestfunc.c, so include ptapsotestfunc.h . */
#include <stdio.h>
#include "maxphase.h"
#include "LLR_PSO.h"
#include <gsl/gsl_math.h>

/*! Number of locations at which we will evaluate the fitness function.*/
const unsigned int nrows = 3;
/*! Dimensionality of the fitness function. */
const size_t nDim = 7;
/*!  Coordinates of the locations in standardized coordinates.
The first location is completely inside the unit hypercube.
The
second location is outside the unit hypercube (first coordinate is 1.5 > 1) and the third is right
on a hypercube face (first coordinate is = 0).  */
double testCoordPts[3][7] = { {.1,.2,.3,.4,.5,.6,.7},
                              {1.5,.2,.3,.4,.5,.6,.7},
							  {0,.1,.2,.4,.5,.6,.7}
                             };
/*! The minimum value of each coordinate in the real
    search space. */
double rmin[7] = {0,-M_PI_2, 0.001, 0, -8, 0, 0};
/*! The maximum value of each coordinate in the real search space. */
double rmax[7] = {2*M_PI, M_PI_2, 0.005, M_PI, -5, M_PI, M_PI};

unsigned int Np = 2, N = 6;
//double *s, *sd;
double sdVec[2] = {0.0001, 0.0002};
double aVec[2] = {M_PI/3.0, M_PI/6.0};
double dVec[2] = {M_PI/4.0, M_PI/4.0};
double sVec1[6] = {.001,.002,.003,.004,.005,.006};
double sVec2[6] = {0.007,.001,.002,.004,.005,.006};
double yrVec[6] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};

int main(){
	unsigned int lpr, lpc;

	double fitValVec[nrows];
	double realCoords[nrows][nDim];
	double rangeVec[nDim];

	gsl_vector *xVec = gsl_vector_alloc(nDim);

	//! [Initializing fitness param struct]
	/*
	  All fitness functions will use the fitFuncParams struct declared in
	  maxphase.h to ferry parameters needed by a fitness function.
	*/
	struct fitFuncParams inParams;
	inParams.nDim = nDim;
	inParams.rmin = gsl_vector_alloc(nDim);
	inParams.rangeVec = gsl_vector_alloc(nDim);
	inParams.realCoord = gsl_vector_alloc(nDim);

	struct llr_pso_params splParams;
  splParams.Np = Np;
  splParams.N = N;

  //inParams->s = (double **)malloc(3*sizeof(double*));
  splParams.sd = sdVec;
  printf("test_LLR_PSO: sd[0] = %f\n", splParams.sd[0]);
  printf("test_LLR_PSO: sd[1] = %f\n", splParams.sd[1]);

  splParams.alphaP = aVec;
  splParams.deltaP = dVec;

  splParams.s = (double **)malloc(2*sizeof(double *));
  splParams.s[0] = sVec1;
  splParams.s[1] = sVec2;
  printf("sVec[0][0] = %f\n", sVec1[0]);
  printf("inParams.s[0]+1 = %f\n", splParams.s[0][3]);
  printf("inParams.s[1]+5 = %f\n", splParams.s[1][5]);

  splParams.yr = yrVec;

  inParams.splParams = &splParams;
	//! [Initializing fitness param struct]

	/* Compute or supply the range for each coordinate */
	for (lpc = 0; lpc < nDim; lpc++){
		rangeVec[lpc]=rmax[lpc]-rmin[lpc];
		gsl_vector_set(inParams.rmin,lpc,rmin[lpc]);
		gsl_vector_set(inParams.rangeVec,lpc,rangeVec[lpc]);
	}


	/* Load test matrix, obtain fitness values and real coordinates*/
	for (lpr = 0; lpr < nrows; lpr++){
		for (lpc=0; lpc < nDim; lpc++){
			gsl_vector_set(xVec,lpc,testCoordPts[lpr][lpc]);
		}
		//fitValVec[lpr] = ptapsotestfunc(xVec,&inParams);
    fitValVec[lpr] = LLR_PSO(xVec,&inParams);
		for (lpc=0; lpc < nDim; lpc++){
			realCoords[lpr][lpc] = gsl_vector_get(inParams.realCoord,lpc);
		}
	}

	for (lpr = 0; lpr < nrows; lpr++){
		printf("Std. coord");
		for (lpc = 0; lpc < nDim; lpc++){
            printf("%f,", testCoordPts[lpr][lpc]);
		}
		printf("\n Real Coord:");
		for (lpc = 0; lpc < nDim; lpc++){
            printf("%f,", realCoords[lpr][lpc]);
		}
		printf("\n Fitness: %f \n",fitValVec[lpr]);
	}

	/* Free allocated memory */
	gsl_vector_free(xVec);
	gsl_vector_free(inParams.rmin);
	gsl_vector_free(inParams.realCoord);
	gsl_vector_free(inParams.rangeVec);

}
