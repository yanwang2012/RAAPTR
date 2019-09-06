// Test harness for backcomp.c

#include "backcomp.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <stdio.h>

int main(int argc, char *argv[]){
    size_t lpc1, lpc2;
    //Test data
    double dataVec[] ={1,2,3,4,5};
    gsl_vector *data = gsl_vector_alloc(5);
    for(lpc1 = 0; lpc1 < 5; lpc1++){
        gsl_vector_set(data,lpc1,dataVec[lpc1]);
    }
    gsl_matrix *dataMat = raaptr_vec2mat(data, 5);
    for(lpc1 = 0; lpc1 < 5; lpc1++){
        for(lpc2 = 0; lpc2 < 5; lpc2++){
            printf("%f,",gsl_matrix_get(dataMat,lpc1,lpc2));
        }
        printf("\n");
    }
    gsl_vector_free(data);
    gsl_matrix_free(dataMat);
}