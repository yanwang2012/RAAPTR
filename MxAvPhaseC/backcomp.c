#include "backcomp.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <stdio.h>
#include <stddef.h>

/*! \file
\brief Functions for implementing backward compatibility. 

For example, in older input
data files, 'yr' is a vector because the observational epochs for all pulsars
were assumed to be identical. However, in newer code that handles variable epocs,
'yr' is a matrix. So, code is needed to seamlessly create a matrix from a vector.

\author Soumya D. Mohanty
*/

/*!
Creates and populates each row of a matrix with the same vector. 
Returns NULL if unsuccessful.
*/
gsl_matrix * raaptr_vec2mat(gsl_vector *vec, size_t nRows){
    if (vec == NULL){
        fprintf(stdout,'vec2mat: Received a null input\n');
        return NULL;
    }
    if (nRows < 1){
        fprintf(stdout,'vec2mat: Number of rows < 1\n')
        return NULL;
    }
    //Number of elements in vector = number of columns
    size_t nCols = vec->size;

    if (nCols < 1){
        fprintf(stdout, 'vec2mat: Number of columns < 1\n');
        return NULL;
    }
    
    gsl_matrix *mat = gsl_matrix_calloc(nRows,nCols);
    if (mat == NULL){
        fprintf(stdout, 'vec2mat: Could not allocate gsl_matrix\n');
        return NULL;
    }

    for (size_t lpc = 0; lpc < nRows; lpc++)
    {
        if(gsl_matrix_set_row(mat,lpc,vec)){
            fprintf(stdout,'vec2mat: Error copying vec to mat row %zd\n',lpc);
        };
    }

    return mat;
    

}

