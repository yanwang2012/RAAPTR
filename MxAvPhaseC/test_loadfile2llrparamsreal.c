/*
 * @Author: Yiqian Qian
 * @Description: test file for func loadfile2llrparam_real
 * @Date: 2022-09-19 19:05:41
 * @LastEditors: Yiqian Qian
 * @LastEditTime: 2022-09-21 15:23:10
 * @FilePath: /MxAvPhaseC/TEST_LoadFile2LLRParamsReal.c
 */
#include "gslhdf5_io.h"
#include "hdf5_hl.h"
#include "loadrealdata.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

void readpsrnames(char *, const char **);

void main(int argc, char *argv[])
{
  /*
  A function to load data from .hdf5 file for real data.
  */
  char *inputFile = argv[1];
  char *psrfile = argv[2];
  printf("pulsar file is %s\n", psrfile);

  const char *psrnames_dbug[] = {"B1855+09", "J1713+0747"};
  fprintf(stdout, "input filename is %s\n", inputFile);
  // printf("String length is %d\n", strlen(psrnames[0]));

  hid_t inFile = H5Fopen(inputFile, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (inFile < 0) {
    printf("Error: cannot open file %s", inputFile);
    abort();
  }
  // printf("inFile is %d\n", inFile);
  size_t Np = (size_t)hdf52dscalar(inFile, "Np");
  const char **psrnames = (const char **)malloc(Np * sizeof(const char *)); // jagged array to store pulsar names
  readpsrnames(psrfile, psrnames);
  // printf("Pulsar names are %s\n", psrnames[0]);

  struct real_data *llp;
  // printf("inFile is %d\n", inFile);
  printf("pulsarnames_dbug 0 %s, pulsarnames 0 %s", psrnames_dbug[0], psrnames[0]);
  llp = loadfile2llrparam_real(inFile, psrnames, Np);
  FILE *fptr = fopen("Output.txt", "w");
  fprintf(fptr, "Number of observations for psr %s is %f\n", psrnames[0],
          llp->N[0]);
  fprintf(fptr, "The location of pulsar %s is %f, %f\n", psrnames[0],
          llp->alphaP[0], llp->deltaP[0]);

  for (int psr = 0; psr < (int)Np; psr++)
  {
    // print yr for pulsar psr
    fprintf(fptr, "yr for pulsar %s is ", psrnames[psr]);
    for (int i = 0; i < llp->N[psr]; i++)
      fprintf(fptr, "%f\n", llp->yr[psr][i]);
    // print sd for pular psr
    fprintf(fptr, "sd for pulsar %s is ", psrnames[psr]);
    for (int i = 0; i < llp->N[psr]; i++)
      fprintf(fptr, "%f\n", llp->sd[psr][i]);
    // print timing residuals for pulsar psr
    fprintf(fptr, "timing residuals for pulsar %s is ", psrnames[psr]);
    for (int i = 0; i < llp->N[psr]; i++)
      fprintf(fptr, "%f\n", llp->s[psr][i]);
  }
  // close file
  H5Fclose(inFile);
  fclose(fptr);

  // free memory
  for (int i = 0; i < (int)Np; i++)
  {
    free(llp->yr[i]);
    free(llp->sd[i]);
    free(llp->s[i]);
  }

  free(llp->yr);
  free(llp->sd);
  free(llp->s);
  free(llp->alphaP);
  free(llp->deltaP);
  free(llp);
}

void readpsrnames(char *filename, const char **psrNames)
{
  FILE *fptr = fopen(filename, "r");
  FILE *fptr2 = fopen("PulsarNames.txt", "w");
  char *line = NULL;
  size_t len = 0;
  ssize_t read;
  int i = 0;
  while ((read = getline(&line, &len, fptr)) != -1)
  {
    psrNames[i] = (char *)malloc((strlen(line) + 1) * sizeof(char));
    strcpy(psrNames[i], line);
    printf("Read pulsar %s\n", psrNames[i]);
    fprintf(fptr2, "%s", psrNames[i]);
    i++;
  }
  fclose(fptr);
  fclose(fptr2);
}