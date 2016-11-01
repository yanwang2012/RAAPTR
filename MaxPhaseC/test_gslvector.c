/* test suite for gsl vector data type */
#include <stdio.h>
#include <gsl/gsl_vector.h>

int main(){
	/* create a gsl_vector */
	gsl_vector *vec = gsl_vector_calloc(5);
	printf("stride: %zu \n", vec->stride);
	/* formatted print out of vector elements */
	gsl_vector_fprintf(stdout,vec,"%f");
	/*----Always free allocated memory before exiting----*/
	gsl_vector_free(vec);
}