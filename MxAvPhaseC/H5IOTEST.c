/*
 * @Author: Yiqian Qian
 * @Description: file content
 * @Date: 2022-08-31 15:56:01
 * @LastEditors: Yiqian Qian
 * @LastEditTime: 2022-09-19 22:41:40
 * @FilePath: /MxAvPhaseC/H5IOTEST.c
 */
#include "gslhdf5_io.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include <gsl/gsl_vector.h>
#include <stdio.h>
#include <string.h>

void main(int argc, char **argv)
{
  /*load data from .hdf5 file*/
  herr_t status;
  char *inputFile = argv[1];
  const char *psrnames[2] = {"B1855+09", "J1713+0747"};
  fprintf(stdout, "input filename is %s\n", inputFile);
  // Open HDF5 file
  hid_t inFile = H5Fopen(inputFile, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (inFile < 0)
  {
    fprintf(stdout, "Error openning file!\n");
    abort();
  }
  // declare yr array
  int numPsr = 2;
  double **yr_array = (double **)malloc(numPsr * sizeof(double *));
  // jagged array
  for (int i = 0; i < 2; i++)
  {
    fprintf(stdout, "Load pulsar %s\n", psrnames[i]);
    // Open group
    hid_t group = H5Gopen(inFile, psrnames[i], H5P_DEFAULT);
    if (group < 0)
    {
      fprintf(stdout, "Error openning group!\n");
      abort();
    }
    // Read dataset
    int N = hdf52dscalar(group, "N");
    fprintf(stdout, "Number of observations %d\n", N);

    gsl_vector *yr = hdf52gslvector(group, "yr");
    printf("Observations are:\n");
    gsl_vector_fprintf(stdout, yr, "%f");
    yr_array[i] = (double *)malloc(N * sizeof(double));
    // asign value
    for (int j = 0; j < N; j++)
    {
      yr_array[i][j] = gsl_vector_get(yr, j);
    }
    fprintf(stdout, "yr_array[%d][10] = %f\n", i, yr_array[i][10]);

    // free parameter
    gsl_vector_free(yr);
    // close group
    H5Gclose(group);
  }
  // close file
  H5Fclose(inFile);
  fprintf(stdout, "Pulsar %s yr is %.9f\n", psrnames[0], yr_array[0][10]);
  fprintf(stdout, "Pulsar %s yr is %.9f\n", psrnames[1], yr_array[1][10]);
  // free memory
  for (int i = 0; i < 2; i++)
    free(yr_array[i]);

  free(yr_array);
}