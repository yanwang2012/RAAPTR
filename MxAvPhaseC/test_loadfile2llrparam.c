#include "backcomp.h"
#include "maxphase.h"
#include "LLR_Mp_Av.h"
#include "perfeval_omp.h"
#include "gslhdf5_io.h"
#include "hdf5_hl.h"
#include "omp.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

/* Test that loadfile2llrparam can read matrix yr data from file correctly. 
   Also check backward compatibility with older files.
   
   gcc -c -I$(HDF5FLAGS) -I$(GSLFLAGS) perfeval_omp.c
   gcc -c -I$(HDF5FLAGS) -I$(GSLFLAGS) gslhdf5_io.c
   gcc -c -I$(GSLFLAGS) LLR_Mp_Av.c 
   gcc -c -I$(GSLFLAGS) backcomp.c 
   gcc $(OBJ) -L$(HDF5LIBS) -lhdf5 -lhdf5_hl -lz -L$(GSLLIBS) -lm -lgsl -lgslcblas -o test_loadfile2llrparam.out
   
   */
   int main(int argc, char *argv[]){
        //Input file 
        char *inputFileName = argv[1];

        //Load data from specified .hdf5 input file
        herr_t status;
        hid_t inFile = H5Fopen(inputFileName, H5F_ACC_RDONLY, H5P_DEFAULT);
        if (inFile < 0)
        {
            fprintf(stdout, "Error opening file %s\n",inputFileName);
            abort();
        }

        //Remaining data loads into special parameter structure
        struct llr_pso_params *llp;
        llp = loadfile2llrparam(inFile);

        //Display yr
        unsigned int N = llp->N;
        unsigned int Np = llp->Np;
        /* Display only a subset of values for each pulsar's yr. 
           Change dispBlckNum below to change block length. */
        size_t dispBlckNum = 5;
        size_t dispBlck[3*dispBlckNum];
        size_t i;
        for (i = 0; i < dispBlckNum; i++){
            dispBlck[i]=i;
        }
        for (i = dispBlckNum; i < 2*dispBlckNum; i++){
            dispBlck[i]=i-dispBlckNum+(N/2);
        }
        for (i = 2*dispBlckNum; i < 3*dispBlckNum; i++){
            dispBlck[i]=i-2*dispBlckNum+N-dispBlckNum;
        }
        printf("Printing these elements of yr\n");
        for (i = 0; i < 3*dispBlckNum; i++){
            printf("%zu,",dispBlck[i]); 
    
        }
        printf("\n");

        for (size_t lpc1 = 0; lpc1 < Np; lpc1++){
            printf("Pulsar %zu: ",lpc1);
            for (size_t lpc2 = 0; lpc2 < 3*dispBlckNum; lpc2++){
                printf("%f,",llp->yr[lpc1][dispBlck[lpc2]]);
            }
            printf("\n");
        }
        //Wrap up
        llrparam_free(llp);
   }