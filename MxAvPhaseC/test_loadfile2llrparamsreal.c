/*
 * @Author: Yiqian Qian
 * @Description: file content
 * @Date: 2022-09-19 19:05:41
 * @LastEditors: Yiqian Qian
 * @LastEditTime: 2022-09-20 00:36:31
 * @FilePath: /MxAvPhaseC/test_loadfile2llrparamsreal.c
 */
#include "LLR_Mp_Av.h"
#include "gslhdf5_io.h"
#include "hdf5_hl.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

void main(int argc, char *argv[])
{
    // load data from .hdf5 file
    //herr_t status;
    char *inputFile = argv[1];
    const char *psrnames[2] = {"B1855+09", "J1713+0747"};
    fprintf(stdout, "input filename is %s\n", inputFile);
    // printf("String length is %d\n", strlen(psrnames[0]));

    hid_t inFile = H5Fopen(inputFile, H5F_ACC_RDONLY, H5P_DEFAULT);

    struct llr_pso_params_real *llp = loadfile2llrparam_real(inFile, psrnames);
    printf("Number of observations is %f\n", llp->N);
    printf("yr of pulsar %s is %f\n", psrnames[0], llp->yr[0]);
    printf("first sd of psr %s is %f\n", psrnames[0], llp->sd[0][0]);
    printf("timing residuals for psr %s are %f\n", psrnames[1], llp->s[1]);

    // close file
    H5Fclose(inFile);
    
    // free memory
    for (int i = 0; i < 2; i++)
    {
        free(llp->yr[i]);
        free(llp->sd[i]);
        free(llp->s[i]);
    }

    free(llp->yr);
    free(llp->sd);
    free(llp->s);
    free(llp);
}