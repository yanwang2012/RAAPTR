/* Code to test read/write of .mat files */
#include "mat.h"
#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <string.h>

int main(){

size_t lpc;

int stat;

const char file[] = "testMatFile.mat";

double data[] = {1.0, 2.0, 3.0, 4.0, 4.5, 3.5, 2.5, 1.5};

/* ------------
   We will test the translation from a gsl_vector to an mxArray 
   and then write out to a .mat file.
   ------------
*/

/* Create a gsl_vector */
gsl_vector *gsl_data = gsl_vector_calloc(8);

/* Load the gsl_vector */
for (lpc = 0; lpc < 8; lpc++)
  	gsl_vector_set(gsl_data,lpc,data[lpc]);

/* Create a temporary array */
double *tmpData = (double *)malloc(8*sizeof(double));

/* Copy gsl_vector data to the temporary array */
for (lpc = 0; lpc < 8; lpc++){
	tmpData[lpc] = gsl_vector_get(gsl_data,lpc);
}

/* An mxArray can only accept pointer to an array created using matlab
   memory allocation functions

double *gsl2mxDataPr = (double *)mxMalloc(8*sizeof(double));
*/


/* Load data from gsl_vector to memory pointed to by gsl2mxDataPr 
for (lpc = 0; lpc < 8; lpc++)
	gsl2mxDataPr[lpc] = gsl_vector_get(gsl_data,lpc);
*/

/* Create an mxArray */
mxArray *mxData = mxCreateDoubleMatrix(1, 8, mxREAL);
if (mxData == NULL){
	printf("Error creating mxArray\n");
	return(EXIT_FAILURE);
}

/* copy data from temporary array to mxArray */
memcpy(mxGetPr(mxData),tmpData,8*sizeof(double));

/* Make the mxArray point to the data pointed to by gsl2mxDataPr 
mxSetPr(mxData,gsl2mxDataPr);
*/

/* Open a mat file for writing data*/
MATFile *testFile = matOpen(file,"w");

if (testFile == NULL) {
    printf("Error creating file %s\n", file);
    return(EXIT_FAILURE);
  }

/* Save the mxArray to this file */
if (matPutVariable(testFile, "data", mxData))
	  printf("Error saving data to mat file\n");
  
/* Close the mat file */
if ((stat =  matClose(testFile)))
	printf("Error closing file %s %d\n", file, stat);

/* ------------
Now we will read the data back from the same .mat file and
load into a gsl_vector
   ------------
*/

/* Open the .mat file for reading data (we can use the 
same pointer to the file because we already closed the
associated file above)*/
testFile = matOpen(file,"r");

/* Print what variables are stored in the mat file */
int numVar;
char **varList;
if((varList=matGetDir(testFile,&numVar)) == NULL){
	printf("Error getting info about variables in %s.mat\n",file);
	if ((stat =  matClose(testFile)))
		printf("Error closing file %s %d\n", file, stat);
}
printf("Number of stored variables in %s.mat: %d\n",file,numVar);
printf("Names of stored variables \n");
for (lpc = 0; lpc < numVar; lpc++)
	printf("%s\n", varList[lpc]);

/* Load data into an mxArray */
mxArray *mxData2;
mxData2 = matGetVariable(testFile,"data");
/* Get info about the data */
mwSize nDims = mxGetNumberOfDimensions(mxData2);
printf("The variable data is a %d dimensional array\n",(int)nDims);
const mwSize *nRowsCols = mxGetDimensions(mxData2);
for (lpc = 0; lpc < nDims; lpc++){
	printf("Number of elements along dimension %zu is %d\n",lpc+1,(int)nRowsCols[lpc]);
}

/* Copy data from mxArray to ordinary dynamically allocated array (we
   do not assume that we know the size of the allocation here). However,
   we do assume that the data is double precision.
 */
size_t sizeData;
sizeData = mxGetNumberOfElements(mxData2);
printf("Number of elements in variable data %zu\n",sizeData);
printf("Printing content of mxArray loaded from mat file\n");
double *mx2Pr = mxGetPr(mxData2);
for (lpc = 0; lpc < sizeData; lpc++)
	printf("%f\n",mx2Pr[lpc]);

double *tmpData2 = (double *)malloc(sizeData*sizeof(double));

memcpy(tmpData2,mx2Pr,sizeData*sizeof(double));

/* Copy data from temporary array to a gsl_vector */
gsl_vector *gsl_data2 = gsl_vector_calloc(sizeData);

printf("Printing content of GSL vector \n");
for (lpc = 0; lpc < sizeData; lpc++){
	gsl_vector_set(gsl_data2,lpc,tmpData2[lpc]);
    printf("%f\n",gsl_vector_get(gsl_data2,lpc));
}

if ((stat =  matClose(testFile)))
	printf("Error closing file %s %d\n", file, stat);

/* Free dynamically allocated memory */
gsl_vector_free(gsl_data);
gsl_vector_free(gsl_data2);
free(tmpData);
free(tmpData2);
mxDestroyArray(mxData);
mxDestroyArray(mxData2);
mxFree(varList);
}

