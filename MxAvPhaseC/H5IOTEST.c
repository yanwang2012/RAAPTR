/*
 * @Author: Yiqian Qian
 * @Description: file content
 * @Date: 2022-08-31 15:56:01
 * @LastEditors: Yiqian Qian
 * @LastEditTime: 2022-09-08 18:40:04
 * @FilePath: /MxAvPhaseC/H5IOTEST.c
 */
#include "gslhdf5_io.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include <gsl/gsl_vector.h>
#include <stdio.h>
#include <string.h>

void main(int argc, char **argv) {
  /*load data from .hdf5 file*/
  herr_t status;
  char * inputFile = argv[1];
  char * psrnames = argv[2];
  fprintf(stdout,"input filename is %s\n", inputFile);
  // Open HDF5 file
  hid_t inFile = H5Fopen(inputFile, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (inFile < 0) {
    fprintf(stdout, "Error openning file!\n");
    abort();
  }
  // Open group
  hid_t group = H5Gopen(inFile, psrnames, H5P_DEFAULT);
  if (group < 0) {
    fprintf(stdout, "Error openning group!\n");
    abort();
  }
  // Read dataset
  gsl_vector * yr = hdf52gslvector(group, "yr");
  int N = hdf52dscalar(group, "N");
  gsl_vector_fprintf(stdout, yr, "%f");
  fprintf(stdout, "Number of observations %d", N);

  // free memory
  gsl_vector_free(yr);
  // close group
  H5Gclose(group);
  // close file
  H5Fclose(inFile);
}
