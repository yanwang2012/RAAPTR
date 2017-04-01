#if !defined(GSLHDF5HDR)
#define GSLHDF5HDR

#include "hdf5.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <stdio.h>

//Soumya D. Mohanty, 2016

//Input: hdf file, dataset name, gsl_vector
void gslvector2hdf5(hid_t, const char *, const gsl_vector *);

//Input: hdf file, dataset name, gsl_matrix
void gslmatrix2hdf5(hid_t, const char *, const gsl_matrix *);

//Put double scalar in HDF5 file
void dscalar2hdf5(hid_t, const char *, const double );

//Input: HDF5 File pointer, dataset name
//Output: gsl_vector *
gsl_vector * hdf52gslvector(hid_t, const char *);

//Read matrix from hdf5 file
gsl_matrix * hdf52gslmatrix(hid_t, const char *);

//Read double scalar from HDF5 file
double hdf52dscalar(hid_t, const char *);
	
#endif