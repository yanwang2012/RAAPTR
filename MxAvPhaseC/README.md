This is the C-code that implements the MaxPhase and AvPhase algorithms. The codes are documented using Doxygen. The 'Doxyfile' is provided as part of the distribution, so run 'doxygen' to generate the documentation. The Doxyfile may be out of date.
A global overview of the codes and examples will be provided in a User manual that is under development.

The main function is 'perfeval_spmd', which is also the default name of the executable. 

To test the codes:
1. Create a directory called 'results' inside TESTPERFEVALOMP
2. Compile the code using 'makefile_omp_local' after editing paths to HDF5, GSL etc. ('makefile_omp' is used for compilation on TACC systems.) **Note:** You may need to edit file 'maxphaseutils.c' on some systems. Go to the lines (there are two such lines) containing the string FNM_FILE_NAME. Change FNM_FILE_NAME to FNM_PATHNAME.
3. Set the number of threads to use in OpenMP: export OMP_NUM_THREADS = [desired number of threads]
4. Run the command: ./perfeval_spmd.out searchParams_simDataSKA.hdf5 TESTPERFEVALOMP/noise1.hdf5 TESTPERFEVALOMP/results/noise1_RAAPTR_test_avp.hdf5 avPhase
