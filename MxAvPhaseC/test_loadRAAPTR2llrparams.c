/*
 * @Author: Yiqian Qian
 * @Description: test file for func loadfile2llrparam_real
 * @Date: 2022-09-19 19:05:41
 * @LastEditors: Yiqian Qian
 * @LastEditTime: 2022-09-29 22:21:49
 * @FilePath: /MxAvPhaseC/test_loadRAAPTR2llrparams.c
 */
#include "gslhdf5_io.h"
#include "hdf5_hl.h"
#include "loadRAAPTR.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

void readpsrnames(const char *, char **, size_t);

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
  if (inFile < 0)
  {
    printf("Error: cannot open file %s", inputFile);
    abort();
  }
  // printf("inFile is %d\n", inFile);
  size_t Np = (size_t)hdf52dscalar(inFile, "Np");
  // size_t Np = 2;
  char **psrnames = (char **)malloc(
      Np * sizeof(char *)); // jagged array to store pulsar names
  readpsrnames(psrfile, psrnames, Np);
  // printf("Pulsar names are %s\n", psrnames[0]);

  struct RAAPTR_data *llp;
  // printf("inFile is %d\n", inFile);
  printf("pulsarnames_dbug 0 %s, pulsarnames 0 %s\n", psrnames_dbug[0],
         psrnames[0]);
  llp = loadRAAPTR2llrparams(inFile, psrnames, Np);
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
  double *ssd = llp->sd[1];
  printf("ssd[1] is %f\n sd[1][1] is %f", ssd[1], llp->sd[1][1]);
  // close file
  H5Fclose(inFile);
  fclose(fptr);

  // free memory
  for (int i = 0; i < (int)Np; i++)
  {
    free(llp->yr[i]);
    free(llp->sd[i]);
    free(llp->s[i]);
    // free(psrnames[i]);
  }

  free(llp->yr);
  free(llp->sd);
  free(llp->s);
  free(llp->alphaP);
  free(llp->deltaP);
  free(llp);
  // free(psrnames);
}

void readpsrnames(const char *filename, char **psrNames, size_t Np)
{
  FILE *fptr = fopen(filename, "r");
  FILE *fptr2 = fopen("PulsarNames.txt", "w");
  char *line = NULL;
  size_t len = 0;
  ssize_t read;
  int i = 0;

  while ((read = getline(&line, &len, fptr)) !=
         -1) // getline will store the newline character as part of the string.
  {
    psrNames[i] = (char *)malloc((strlen(line) + 1) * sizeof(char));
    // printf("The length of %s is %d\n", line, strlen(line));
    strcpy(psrNames[i], line);
    if (i != Np - 1)                        // the last line does not have a newline character
      psrNames[i][strlen(line) - 1] = '\0'; // change the last character to '\0'
    printf("Read pulsar %s\n", psrNames[i]);
    fprintf(fptr2, "%s loaded.\n", psrNames[i]);
    i++;
  }

  /*for (int i = 0; i < Np; i++ ){
     psrNames[i] = (char *)malloc(buffersize * sizeof(char));
     getline(&psrNames[i], &buffersize, fptr);
     fprintf(fptr2, "%s", psrNames[i]);
  }*/
  fclose(fptr);
  fclose(fptr2);
}
