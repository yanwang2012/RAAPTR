// http://mpitutorial.com/tutorials/mpi-hello-world/
// Author: Wes Kendall
// Copyright 2011 www.mpitutorial.com
// This code is provided freely with the tutorials on mpitutorial.com. Feel
// free to modify it for your own use. Any distribution of the code must
// either provide a link to www.mpitutorial.com or keep this header intact.
//
// An intro MPI hello world program that uses MPI_Init, MPI_Comm_size,
// MPI_Comm_rank, MPI_Finalize, and MPI_Get_processor_name.
//
#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_vector.h>
#include "mat.h"
#include "matrix.h"

/* How to compile:
mpicc mpiHelloWorldMatFiles.c 
 -I/Applications/R2015bTAH/MATLAB_R2015b.app/extern/include 
 -L/Applications/R2015bTAH/MATLAB_R2015b.app/bin/maci64 
 -lmat -lmx 
 -Xlinker -rpath -Xlinker /Applications/R2015bTAH/MATLAB_R2015b.app/bin/maci64 
 -o mpihelloworldmatfiles

  How to run:
mpirun -n 4 ./mpihelloworldmatfiles testData/searchParams_simDataX.mat
*/

int main(int argc, char** argv) {
  // Initialize the MPI environment. The two arguments to MPI Init are not
  // currently used by MPI implementations, but are there in case future
  // implementations might need the arguments.
  MPI_Init(NULL, NULL);

  // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  // Get the name of the processor
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);
  
	/* Full path to search parameter file */
	char *srchParamsFile = argv[1];
	
	/* Read xmaxmin from srchParamsFile file */
	MATFile *srchPar = matOpen(srchParamsFile, "r");
	
	mxArray *xmaxmin, *dummy;
	dummy = mxCreateDoubleScalar(0.0);
	double *dummyPr = mxGetPr(dummy);
	
	xmaxmin = matGetVariable(srchPar,"xmaxmin");
	
	/* Search Space dimensionality */
	size_t nDim = mxGetM(xmaxmin);
	dummyPr[0] = (double)nDim;
	
    // Print off a hello world message
    printf("Hello world from processor %s, rank %d out of %d processors: nDim = %zu\n",
           processor_name, world_rank, world_size, nDim);
   /* Print the same message into a file */
   char outFilePrefix[] = "hworldMPIout_";
   char outFileNum[2];
   char *outFileName = (char *)malloc((strlen(outFilePrefix)+1)*sizeof(char));
   strcpy(outFileName,outFilePrefix);
   sprintf(outFileNum,"%d",world_rank);
   strcat(outFileName,outFileNum);
   printf("Out file name %s\n",outFileName);
   FILE *outFile = fopen(outFileName,"w");
   fprintf(outFile,"Hello world from processor %s, rank %d out of %d processors: nDim = %zu\n",
           processor_name, world_rank, world_size, nDim);
   fclose(outFile); 
   
   char *matOutFileName = (char *)malloc((strlen(outFileName)+4)*sizeof(char));
   strcpy(matOutFileName,outFileName);
   strcat(matOutFileName,".mat");
   MATFile *outMatFile = matOpen(matOutFileName,"w");
   if (matPutVariable(outMatFile,"test",dummy))
		printf("Error storing variable %s in file %s\n", "test", matOutFileName);
  
   if(matClose(outMatFile))
	   printf("Error closing file %s\n",matOutFileName);
	   
	if (matClose(srchPar))
		printf("Error closing file %s\n", srchParamsFile);
	
	//Dummy statements with GSL calls
	gsl_vector *gslDummyVec = gsl_vector_alloc(10);
	gsl_vector_free(gslDummyVec);
	
	mxDestroyArray(xmaxmin);
	free(outFileName);
	free(matOutFileName);

  // Finalize the MPI environment. No more MPI calls can be made after this
  MPI_Finalize();
}