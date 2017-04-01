# makefile for MaxPhase project.  Feb. 22, 2016

OBJ = perfeval.o  LLR_PSO.o  maxphaseutils.o ptapso.o mp2matfile_io.o perfeval_snglproc.o
CC = gcc
GSLFLAGS :=  $(shell gsl-config --cflags)
GSLLIBS := $(shell gsl-config --libs)
MATLABFLAGS = /Applications/MATLAB_R2016a.app/extern/include
MATLABLIBS = /Applications/MATLAB_R2016a.app/bin/maci64

perfeval_snglproc.out: $(OBJ)
	$(CC) $(OBJ) -L$(MATLABLIBS) -lmat -lmx -Xlinker -rpath -Xlinker \
	$(MATLABLIBS) $(GSLLIBS) -lm -o perfeval_snglproc.out
	
perfeval_snglproc.o: perfeval_snglproc.c  
	$(CC) -g -c -I$(MATLABFLAGS) $(GSLFLAGS) perfeval_snglproc.c

mp2matfile_io.o: mp2matfile_io.c  
	$(CC) -g -c -I$(MATLABFLAGS) $(GSLFLAGS) mp2matfile_io.c
	
perfeval.o: perfeval.c  
	$(CC) -g -c -I$(MATLABFLAGS) $(GSLFLAGS) perfeval.c

LLR_PSO.o: LLR_PSO.c  
	$(CC) -g -c LLR_PSO.c -I$(MATLABFLAGS) $(GSLFLAGS) 

maxphaseutils.o: maxphaseutils.c  
	$(CC) -g -c maxphaseutils.c $(GSLFLAGS) 

ptapso.o: ptapso.c  
	$(CC) -g -c ptapso.c $(GSLFLAGS) 
		
.PHONY: clean
clean:
		-rm *.o
		-rm -r simDataX/results

# EOF
