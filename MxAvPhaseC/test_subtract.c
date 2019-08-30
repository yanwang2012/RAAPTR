#include "subtract.h"
//#include "maxphase.h"
//#include "LLR_Mp_Av.h"
//#include "ptapso.h"
#include "hdf5.h"
#include "gslhdf5_io.h"
#include <gsl/gsl_vector.h>
#include <sys/stat.h>
#include <stdio.h>
#include <string.h>

int main(int argc, char *argv[]){
    char * outFileName = argv[1];
    printf("output file name: %s\n", outFileName);
    struct estSrcParams * srcp;
    srcp = file2Srcparam(outFileName);
    printParam(srcp);
    srcpara_free(srcp);
    
    /* success message */
    printf("All Done!");

    return 0;
}