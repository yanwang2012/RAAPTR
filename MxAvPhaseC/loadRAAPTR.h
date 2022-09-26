/*
 * @Author: Yiqian Qian
 * @Description: Header file for loadrealdata.c
 * @Date: 2022-09-20 14:35:43
 * @LastEditors: Yiqian Qian
 * @LastEditTime: 2022-09-21 14:42:05
 * @FilePath: /MxAvPhaseC/loadrealdata.h
 */
/*! \file loadrealdata.h
\brief Header file for \ref loadrealdata.c
*/
#if !defined(REALDATA)
#define REALDATA
#include "hdf5.h"
#include "hdf5_hl.h"

struct RAAPTR_data {
  /*! number of pulsars in the timing array */
  unsigned int Np;
  /*! number of observations, jagged array */
  double *N;
  /*! GW signal, jagged array */
  double **s;
  /*! standard deviation of noise for different pulsar, jagged array */
  double **sd;
  /*! (Np,1), right ascension, in radian */
  double *alphaP;
  /*!  (Np,1), declination, in radian */
  double *deltaP;
  /*! observation epoch, in year, jagged array */
  double **yr;
  /*! Estimated phases (1, Np) */
  double *phiI;
};

struct RAAPTR_data *loadRAAPTR2llrparams(hid_t, const char **, size_t);

struct RAAPTR_data *llrRAAPTR_alloc(unsigned int);
#endif
