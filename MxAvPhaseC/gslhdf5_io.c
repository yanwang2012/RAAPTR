#include "gslhdf5_io.h"
#include "hdf5_hl.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include <stdio.h>

//Soumya D. Mohanty, 2016

//Input: hdf file pointer, group name, dataset name, gsl_vector
void gslvector2hdf5(hid_t outFile, const char *dsetName, const gsl_vector *gslVec){
	//Size information
    size_t nSamples = gslVec->size;
	//Copy GSL vector data to ordinary array
	double outData[nSamples];
	size_t lpc;
	for(lpc = 0; lpc < nSamples; lpc++){
		outData[lpc] = gsl_vector_get(gslVec,lpc);
	}
	
	//Use HDF Lite function to write out data into file
	hsize_t dims = nSamples;
	herr_t status = H5LTmake_dataset_double (outFile, dsetName, 1, &dims, outData);
	if (status < 0){
		fprintf(stdout, "Error writing variable %s\n", dsetName);
	}
}

//Input: hdf file, dataset name, gsl_matrix
void gslmatrix2hdf5(hid_t outFile, const char *dsetName, const gsl_matrix *gslMat){
	//Size information
    size_t nRows = gslMat->size1;
    size_t nCols = gslMat->size2;
	printf("%zu %zu\n",nRows, nCols);
	
	//Copy GSL vector data to ordinary array
	double outData[nRows*nCols];
	size_t lpc1, lpc2;
	size_t countElements = 0;
	for(lpc1 = 0; lpc1 < nRows; lpc1++){
		for(lpc2 = 0; lpc2 < nCols; lpc2++){
			outData[countElements] = gsl_matrix_get(gslMat,lpc1,lpc2);
			countElements++;
		}
	}
	
	//Use HDF Lite function to write out data into file
	hsize_t dims[2];
	dims[0] = nRows;
	dims[1] = nCols;
	
	herr_t status = H5LTmake_dataset_double (outFile, dsetName, 2, dims, outData);
	if (status < 0){
		fprintf(stdout, "Error writing variable %s\n", dsetName);
	}	
}

//Put double scalar in HDF5 file
void dscalar2hdf5(hid_t outFile, const char *dsetName, const double data){
	hsize_t dims = 1;
	herr_t status = H5LTmake_dataset_double (outFile, dsetName, 1, &dims, &data);
	if (status < 0){
		fprintf(stdout, "Error writing variable %s\n", dsetName);
	}	
}

//Read double scalar from HDF5 file
double hdf52dscalar(hid_t inFile, const char *dsetName){
	herr_t status;
	//Check if the dataset exists
	if(!(status = H5LTfind_dataset(inFile, dsetName))){
		fprintf(stdout, "dataset %s does not exist\n", dsetName);
		return NAN;
	}
	//Get the dimensions of the dataset
	int nDims;
	status = H5LTget_dataset_ndims(inFile, dsetName, &nDims);
	if (status < 0){
		fprintf(stdout, "Error getting rank of dataset %s\n", dsetName);
		return NAN;
	}
	if (nDims > 1){
		//Get the number of elements along each dimension
		hsize_t *dims = (hsize_t *)malloc(nDims*sizeof(hsize_t));
		H5T_class_t classID;
		size_t typeSize;
		status = H5LTget_dataset_info(inFile, dsetName, dims, &classID, &typeSize);
		//TODO: Check the datatype also; assume for now that it is double
		if(status < 0){
			fprintf(stdout, "Could not get info on dataset %s\n", dsetName);
			free(dims);
			return NAN;
		}
		size_t dimProd=1;
		size_t lp;
		for (lp = 0; lp < nDims; lp++){
			dimProd = dimProd*dims[lp];
		}
		if (dimProd != 1){
			fprintf(stdout, "Dataset %s is not a scalar\n", dsetName);
			free(dims);
			return NAN;
	    }
		//Data is still a scalar
		free(dims);
	}
	
	//Read data		
	double data;
	status = H5LTread_dataset_double(inFile, dsetName, &data);
	
	return data;
}

//Input: HDF5 File pointer, dataset name
//Output: gsl_vector * (NULL for failure)
gsl_vector * hdf52gslvector(hid_t inFile, const char *dsetName){
	herr_t status;
	//Check if the dataset exists
	if(!(status = H5LTfind_dataset(inFile, dsetName))){
		fprintf(stdout, "dataset %s does not exist\n", dsetName);
		return NULL;
	}
	//Get the dimensions of the dataset
	int nDims;
	status = H5LTget_dataset_ndims(inFile, dsetName, &nDims);
	if (status < 0){
		fprintf(stdout, "Error getting rank of dataset %s\n", dsetName);
		return NULL;
	}
	if (nDims > 2){
		fprintf(stdout, "Dataset %s is not a vector\n", dsetName);
		return NULL;
	}
	//Get the number of elements along each dimension
	hsize_t *dims = (hsize_t *)malloc(nDims*sizeof(hsize_t));
	H5T_class_t classID;
	size_t typeSize;
	status = H5LTget_dataset_info(inFile, dsetName, dims, &classID, &typeSize);
	//TODO: Check the datatype also; assume for now that it is double
	if(status < 0){
		fprintf(stdout, "Could not get info on dataset %s\n", dsetName);
		free(dims);
		return NULL;
	}
	//Get the number of samples
	size_t nSamples;
	if (nDims == 1){
		nSamples = *dims;
	}
	else{
		if (dims[0] == 1 || dims[1] == 1) {
			//This is still a vector; 
			nSamples = dims[0];
			if (dims[1] > dims[0]) nSamples = dims[1];
		}
		else{
			fprintf(stdout, "Dataset %s does not contain a vector\n", dsetName);
			return NULL;
		}
	}

	
	// All checks passed
	double *dataVec = (double *)malloc(nSamples*sizeof(double));
	status = H5LTread_dataset_double(inFile, dsetName, dataVec);
	size_t lp;
	gsl_vector *gslVec = gsl_vector_alloc(nSamples);
	for (lp = 0; lp < nSamples; lp++){
		gsl_vector_set(gslVec,lp,dataVec[lp]);
	}
	
	//Wrap up
	free(dims);
	free(dataVec);
	
	return gslVec;
}

//Input: HDF5 File pointer, dataset name
//Output: gsl_vector * (NULL for failure)
gsl_matrix * hdf52gslmatrix(hid_t inFile, const char *dsetName){
	herr_t status;
	//Check if the dataset exists
	if(!(status = H5LTfind_dataset(inFile, dsetName))){
		fprintf(stdout, "dataset %s does not exist\n", dsetName);
		return NULL;
	}
	//Get the dimensions of the dataset
	int nDims;
	status = H5LTget_dataset_ndims(inFile, dsetName, &nDims);
	if (status < 0){
		fprintf(stdout, "Error getting rank of dataset %s\n", dsetName);
		return NULL;
	}
	if (nDims != 2){
		fprintf(stdout, "Dataset %s is not a 2D matrix\n", dsetName);
		return NULL;
	}
	//Get the number of elements along each dimension
	hsize_t *dims = (hsize_t *)malloc(nDims*sizeof(hsize_t));	
	H5T_class_t classID;
	size_t typeSize;
	status = H5LTget_dataset_info(inFile, dsetName, dims, &classID, &typeSize);
	//TODO: Check the datatype also; assume for now that it is double
	if(status < 0){
		fprintf(stdout, "Could not get info on dataset %s\n", dsetName);
		return NULL;
	}
	size_t nRows = (size_t)dims[0];
	size_t nCols = (size_t)dims[1];
	
	//Get the number of samples
	size_t nSamples = nRows * nCols;
	
	// All checks passed
	double *dataVec = (double *)malloc(nSamples*sizeof(double));
	status = H5LTread_dataset_double(inFile, dsetName, dataVec);
	
	//Fill gsl_matrix
	//HDF5 C API stores in row major order [r1, r2,...]; dims[0] is nRows
	size_t lp1, lp2;
	size_t countElements = 0;
	gsl_matrix *gslMat = gsl_matrix_alloc(nRows,nCols);
	for (lp1 = 0; lp1 < nRows ; lp1++){
		for (lp2 = 0; lp2 < nCols; lp2++){
			gsl_matrix_set(gslMat,lp1,lp2,dataVec[countElements]);
			countElements++;
		}
	}
	
	//Wrap up
	free(dims);
	free(dataVec);
	
	return gslMat;
}