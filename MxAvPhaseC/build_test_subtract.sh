gcc -O3 -c -I$TACC_HDF5_INC -I$TACC_GSL_INC gslhdf5_io.c
gcc -O3 -c -I$TACC_HDF5_INC -I$TACC_GSL_INC subtract.c
gcc -O3 gslhdf5_io.o subtract.o -L$TACC_HDF5_LIB -lhdf5 -lhdf5_hl -lz -L$TACC_GSL_LIB -lm -lgsl -lgslcblas -o subtract.out
