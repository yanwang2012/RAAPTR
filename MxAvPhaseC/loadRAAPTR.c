/*
 * @Author: Yiqian Qian
 * @Description: Function to load real data.
 * @Date: 2022-09-20 13:12:27
 * @LastEditors: Yiqian Qian
 * @LastEditTime: 2022-10-01 13:11:22
 * @FilePath: /MxAvPhaseC/loadRAAPTR.c
 */
#include "loadRAAPTR.h"
#include "gslhdf5_io.h"
#include "hdf5_hl.h"
#include <gsl/gsl_vector.h>
#include <stdio.h>
#include <string.h>

struct RAAPTR_data *llrparam_alloc_RAAPTR(unsigned int Np) {
  struct RAAPTR_data *llp =
      (struct RAAPTR_data *)malloc(sizeof(struct RAAPTR_data));
  llp->Np = Np;
  llp->N = (double *)malloc(Np * sizeof(double));
  llp->sd = (double **)malloc(Np * sizeof(double *));
  llp->alphaP = (double *)malloc(Np * sizeof(double));
  llp->deltaP = (double *)malloc(Np * sizeof(double));
  llp->phiI = (double *)malloc(Np * sizeof(double));
  llp->yr = (double **)malloc(Np * sizeof(double *));
  llp->s = (double **)malloc(Np * sizeof(double *));

  return llp;
}

struct RAAPTR_data *loadRAAPTR2llrparams(hid_t inFile, const char **psrNames,
                                         size_t Np) {

  struct RAAPTR_data *llp = llrparam_alloc_RAAPTR((unsigned int)Np);
  printf("First pulsar in loadRAAPTR2llrparams is %s\n", psrNames[0]);
  for (int psr = 0; psr < Np; psr++) {
    hid_t inGroup = H5Gopen(inFile, psrNames[psr], H5P_DEFAULT);
    if (inGroup < 0) {
      printf("Error: cannot open group %s", psrNames[psr]);
      abort();
    }
    // read out the data
    double N = hdf52dscalar(inGroup, "N");
    gsl_vector *yr_Pr = hdf52gslvector(inGroup, "yr");
    gsl_vector *trPr = hdf52gslvector(inGroup, "timingResiduals");
    gsl_vector *sd_Pr = hdf52gslvector(inGroup, "sd");

    // Dynamically allocate memory
    llp->sd[psr] = (double *)malloc(N * sizeof(double));
    llp->s[psr] = (double *)malloc(N * sizeof(double));
    llp->yr[psr] = (double *)malloc(N * sizeof(double));

    // signing data
    llp->Np = Np;
    llp->N[psr] = N;
    llp->alphaP[psr] = hdf52dscalar(inGroup, "alphaP");
    llp->deltaP[psr] = hdf52dscalar(inGroup, "deltaP");
    for (int i = 0; i < N; i++) {
      llp->sd[psr][i] = gsl_vector_get(sd_Pr, i);
      llp->s[psr][i] = gsl_vector_get(trPr, i);
      llp->yr[psr][i] = gsl_vector_get(yr_Pr, i);
    }

    // Wrap up
    gsl_vector_free(yr_Pr);
    gsl_vector_free(trPr);
    gsl_vector_free(sd_Pr);

    // Close group
    herr_t status = H5Gclose(inGroup);
    if (status < 0) {
      printf("Error closing group %s)", psrNames[psr]);
    }
  }
  return llp;
}

void raaptr_free(struct RAAPTR_data * llp){
  size_t Np = llp->Np;
  
  for (int i = 0; i < (int)Np; i++) {
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
  free(llp->phiI);
  free(llp->N);
  free(llp);
}