This is the C-code that implements the MaxPhase and AvPhase algorithms. The codes are documented using Doxygen. The 'Doxyfile' is provided as part of the distribution, so run 'doxygen' to generate the documentation. 
(A global overview of the codes and examples will be provided in a User manual that is under development.)
Note: Currently, the MPI based codes will not work. Use the OpenMP codes ('perfeval_spmd') instead. The Doxyfile may be out of date.

To test the codes:
1. Create a directory called 'results' inside TESTPERFEVALOMP
2. Compile the code using makefile_omp
3. Set the number of threads to use in OpenMP
4. Run the command: ./perfeval_spmd.out searchParams_simDataSKA.hdf5 TESTPERFEVALOMP/noise1.hdf5 TESTPERFEVALOMP/results/noise1_RAAPTR_test_avp.hdf5 avPhase
