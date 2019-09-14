## FIXME
•	There appear to be duplicate mpavinfile2hdf5.m, one in GENSIMDATA and one in MxAvPhase. The latter is the latest version. Investigate and remove the one in GENSIMDATA.

•	There appear to be duplicate mpavoutfile2mat.m, one in GENSIMDATA and one in MxAvPhase. Investigate and remove the one in GENSIMDATA.

•	Assign only as many gsl_rng * as the number of PSO runs (perfeval_omp.c). Right now, the number of gsl_rng * is hard coded to 8 (the number of hard coded seeds).

•	Single letter variables names in LLR_Mp_Av.c should be changed so that find can be used to trace the occurrence of these variables in the code. 

•	In LLR_Mp_Av.c: L138, func ‘f’ is the function handle to be integrated. 

•	Move out all the file I/O functions defined in perfeval_omp.c: Their presence in perfeval_omp.c prevents compilation of test functions that do not need OpenMP.

•	Handle the annoying issue related to FNM_FILE_NAME (on any machine besides LS5) vs. FNM_PATHNAME (on LS5) in maxphaseutils.c.

## TODO
•	Allow user specification of PSO parameters in perfeval_spmd and perfeval_omp.
•	Allow change of number of PSO runs in perfeval_spmd and perfeval_omp.
•	Change unsigned int to size_t for better portability.
•	Remove double **s in LLR_Mp_Av.c since it can be replaced with splParams->s[i] (see Aug 28, 2019 entry in "PTA SuperN Lab Notebook.docx").
•	Name of the search parameter file used should be stored in each output file.
