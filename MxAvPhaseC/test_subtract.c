#include "subtract.h"
#include "perfeval_omp.h"
#include "hdf5_hl.h"
#include "gslhdf5_io.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <sys/stat.h>
#include <stdio.h>
#include <string.h>

int main(int argc, char *argv[]){
    char * inFileName = argv[1];
    char * outFileName = argv[2];

    size_t Np,N;
    double **tres; /*! timing residuals in input file */

    printf("Input file name is: %s\n", inFileName);
    printf("output file name is: %s\n", outFileName);

    struct llr_pso_params * llp;
    struct estSrcParams * srcp;

    //Load data from specified .hdf5 input file
	herr_t status;
	hid_t inFile = H5Fopen(inFileName, H5F_ACC_RDONLY, H5P_DEFAULT);
	if (inFile < 0){
		fprintf(stdout,"Error opening file\n");
		abort();
	}
	char genHypothesis[10];
	status = H5LTread_dataset_string(inFile,"genHypothesis",genHypothesis);
	if (status < 0){
		fprintf(stdout,"Error reading genHypothesis\n");
		abort();
    }

    llp = loadfile2llrparam(inFile);
    srcp = file2Srcparam(outFileName);

    N = (size_t)llp->N;
    //printf("N: %zu\n",N);
    Np = (size_t)llp->Np;
    //printf("Np: %zu\n",Np);
    tres = llp->s;

    gsl_matrix * timResiduals = gsl_matrix_calloc(Np,N);
    timResiduals = timingResiduals(srcp,llp);
    printf("Dimension of timResiduals: %zu %zu\n", timResiduals->size1, timResiduals->size2);

    /* Subtract estimated timing residual from source. */
    size_t i,j;
    for(i = 0; i < Np; i++){
        for(j = 0; j < N; j++){
            gsl_matrix_set(timResiduals,i,j, tres[i][j] - gsl_matrix_get(timResiduals,i,j));
        }
    }

    /* Put subtracted timing residuals into input file.*/
   /* herr_t status;
    hid_t inFile;
    inFile = H5Fopen(inFileName, H5F_ACC_RDONLY, H5P_DEFAULT);
    gslmatrix2hdf5(inFile,"timingResiduals",timResiduals);
    H5Fclose(inFile);*/

    
    FILE * f;
    f = fopen("timingResiduals.txt", "w");
    printMatrix(f,timResiduals,Np,N);// print timing residual to file f. 
    fclose(f);
    

    printParam(srcp);
    srcpara_free(srcp);
    gsl_matrix_free(timResiduals);

    /* success message */
    printf("All Done!");

    return 0;
}
