# Breif script to compile test_subtract.out
gcc -c -I$TACC_HDF5_INC -I$TACC_GSL_INC gslhdf5_io.c
gcc -c -I$TACC_HDF5_INC -I$TACC_GSL_INC subtract.c
gcc -c -I$TACC_HDF5_INC -I$TACC_GSL_INC test_subtract.c
gcc -c -I$TACC_HDF5_INC -I$TACC_GSL_INC perfeval_omp.c
gcc -c -I$TACC_HDF5_INC -I$TACC_GSL_INC maxphaseutils.c
gcc -c -I$TACC_HDF5_INC -I$TACC_GSL_INC LLR_Mp_Av.c
gcc -c -I$TACC_HDF5_INC -I$TACC_GSL_INC ptapso.c
gcc -fopenmp perfeval_omp.o LLR_Mp_Av.o maxphaseutils.o ptapso.o gslhdf5_io.o subtract.o test_subtract.o -L$TACC_HDF5_LIB -lhdf5 -lhdf5_hl -lz -L$TACC_GSL_LIB -lm -lgsl -lgslcblas -o test_subtract.out
