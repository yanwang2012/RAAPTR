/*! \file test_ptapsotestfunc.c
\brief Example of how to call a fitness function.

 If calling fitness function foo.c, include foo.h. Here, we 
are calling ptapsotestfunc.c, so include ptapsotestfunc.h . */
#include <stdio.h>
#include "maxphase.h"
#include "ptapsotestfunc.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

/*! Number of locations at which we will evaluate the fitness function.*/
const unsigned int nrows = 3; 
/*! Dimensionality of the fitness function. */
const size_t nDim = 3;
/*!  Coordinates of the locations in standardized coordinates. 
The first location is completely inside the unit hypercube.
The
second location is outside the unit hypercube (first coordinate is 1.5 > 1) and the third is right
on a hypercube face (first coordinate is = 0).  */
double testCoordPts[3][3] = { {.1,.2,.3},
                              {1.5,.2,.3},
							  {0,.1,.2}
                             };	
/*! The minimum value of each coordinate in the real
    search space. */
double rmin[3] = {-5,-5,-5};
/*! The maximum value of each coordinate in the real search space. */
double rmax[3] = {5, 5, 5};	

						 
int main(){
	unsigned int lpr, lpc;

	double fitValVec[nrows];
	double realCoords[nrows][nDim];
	double rangeVec[nDim];
	unsigned char fitEvalFlag[nrows];
	
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
	/* Special parameters for a fitness function go into 
	   a special structure declared in the .h file for that
	   fitness function. Example: ptapsotestfunc.h since we 
	   are calling the ptapsotestfunc.c fitness function.
	*/
	struct ptapsotestfunc_params spclParams;
	spclParams.dummyParam = 10;
	inParams.splParams = &spclParams;
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
		fitValVec[lpr] = ptapsotestfunc(xVec,&inParams);
		fitEvalFlag[lpr] = inParams.fitEvalFlag;
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
		printf("\n Fitness: %f, ",fitValVec[lpr]);
		printf(" Function evaluated?: %d\n", (int)fitEvalFlag[lpr]);		
	}
	
	/* Evaluate fitness for the same points using vector views */
	printf("Testing Matrix Views\n");
	gsl_vector_view testCoordPtsRow;
	gsl_matrix *testCoordPtsMat = gsl_matrix_alloc(nrows, nDim);
	for (lpr = 0; lpr < nrows; lpr++){
		for (lpc = 0; lpc < nDim; lpc++){
			gsl_matrix_set(testCoordPtsMat,lpr,lpc,testCoordPts[lpr][lpc]);
			printf("%f, ",gsl_matrix_get(testCoordPtsMat,lpr,lpc));
		}
		printf("\n");
	}
	for (lpr = 0; lpr < nrows; lpr++){
		testCoordPtsRow = gsl_matrix_row(testCoordPtsMat,lpr);
		fitValVec[lpr] = ptapsotestfunc(&testCoordPtsRow.vector,&inParams);
		printf("Fitness: %f \n",fitValVec[lpr]);	
	}
	
	/* Free allocated memory */
	gsl_vector_free(xVec);
	gsl_vector_free(inParams.rmin);
	gsl_vector_free(inParams.realCoord);
	gsl_vector_free(inParams.rangeVec);
	gsl_matrix_free(testCoordPtsMat);
}