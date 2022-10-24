/*
 * @Author: Yiqian Qian
 * @Description: header file contains subtraction functions
 * @Date: 2020-09-02 12:41:57
 * @LastEditors: Yiqian Qian
 * @LastEditTime: 2022-10-20 17:14:44
 * @FilePath: /MxAvPhaseC/subtract.h
 */
#include "gslhdf5_io.h"
#include "loadRAAPTR.h"
#include "perfeval_omp.h"
#include "maxphase.h"
#include <gsl/gsl_vector.h>
#include <stdio.h>
#include <string.h>
/*! \file
\brief Header file for function  subtraction().
*/
struct estSrcParams{
    double alpha;
    double delta;
    double omega;
    double phi0;
    double Amp;
    double iota;
    double thetaN;
    gsl_vector * psrPhase;
};

struct estSrcParams * file2Srcparam(char *);
struct estSrcParams * srcp_alloc(size_t);
void printParam(struct estSrcParams *);
void srcpara_free(struct estSrcParams *);
gsl_matrix * timingResiduals(struct estSrcParams *, struct llr_pso_params *);
double ** timingResiduals_raaptr(struct estSrcParams *, struct RAAPTR_data *);
gsl_matrix * FullResiduals(struct estSrcParams *, double, double, double, double, double *, size_t);
void printMatrix(FILE *, gsl_matrix *, size_t, size_t);