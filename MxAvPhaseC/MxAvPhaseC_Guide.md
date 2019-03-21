# About
This is a program to find the local optimization of the data files by using PSO.

# Preparation

1. Use 'test_gensimdata.m' file to generate simulation data files and search parameter files.

2. Converting the data files and search parameter files from .mat file to HDF5 file using 'mpavinfile2hdf5.m'.

3. Compile the code use command

```
make -f makefile_omp
```

**Note:**
Edit file maxphaseutils.c: Go to the lines (there are two such lines) containing the string FNM_FILE_NAME. Change FNM_FILE_NAME to FNM_PATHNAME.

# Search parameters for PSO
By using **parameters.m** file we generate two parameter files, **Sim_Params_X.mat** file is for simulation while the other **searchParams_simDataSKA_X.mat** file contains the search range of PSO, users can change the value in **parameters.m** file.

# Run on LS5
First on LS5 load the module use

```
ml gcc/5.2.0
ml gsl/2.2.1
```

to compile the code.

For running a job use

```
ml gsl/2.2.1
ml launcher
export OMP_NUM_THREADS=8
```

## Create the LAUNCHER JOB FILE
The format of the line is:

```
<executable> <path-to-parameter-file>
<path-to-input-data-file> <path-to-output-folder>
<option: maxPhase or AvPhase>
```

Here is an example of a single line of the job file:

```
~/PULSARTIMING/MaxPhaseC/perfeval_spmd.out
/work/02580/soumya/lonestar/PULSARTIMING/
searchParams_simDataSKA.hdf5
/work/02580/soumya/lonestar/PULSARTIMING/
simDataSKA_snr0123_loc3_HDF5/noise25.hdf5
/work/02580/soumya/lonestar/PULSARTIMING/
simDataSKA_snr0123_loc3_HDF5/results_run2/noise25.hdf5
maxPhase
```

To submit a job to LS5 you should creat a slurm file, here is an example.
**LAUNCHER_SLURM_FILE_PTA.slurm**

```
#!/bin/bash

#

# Simple SLURM script for submitting multiple serial

# jobs (e.g. parametric studies) using a script wrapper

# to launch the jobs.

#

# To use, build the launcher executable and your

# serial application(s) and place them in your WORKDIR

# directory.  Then, edit the CONTROL_FILE to specify  

# each executable per process.

#-------------------------------------------------------

#-------------------------------------------------------

#  

#         <------ Setup Parameters ------>

#

#SBATCH -J PTA_MultiSrc_SKA_mxp_Apr12_2017
         (This is the name of the job)

#SBATCH -N 166 (This is the number of nodes requested.
                Each node on LS5 has 24 processors.)

#SBATCH -n 498 (This is the number of jobs to
  run = number of lines in the LAUNCHER job file.
  With 8 parallel PSO runs per data file,
  we can process 24/8 = 3 data files on each node.
  Hence the number of nodes is 498/3 = 166.)

#SBATCH -p normal (This is the type of queue requested.
  Do not change this.)

#SBATCH -o PTA_Apr12_2017_3.o%j (This is the output file
  where any screen output will be redirected.)

#SBATCH -e PTA_Apr12_2017_3.e%j (This is the file
  where error messages from LAUNCHER will be logged.
  Important for diagnostics if a job fails.)

#SBATCH -t 3:00:00 (Time requested for the job.
  Should be slightly in excess of the anticipated
  completion time of all the jobs. However,
  if any one node hangs, all nodes will wait
  until this time limit is reached and the
  SUs consumed will correspondingly be higher!
  So, it should not be too much higher than the
  anticipated completion time.)

#SBATCH --mail-user=soumya.mohanty@utrgv.edu
(Email address where alerts will be sent by sbatch.)

#SBATCH --mail-type=all

#          <------ Account String ----->

# <--- (Use this ONLY if you have MULTIPLE accounts) --->

##SBATCH -A

#------------------------------------------------------



export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins

export LAUNCHER_RMI=SLURM

export LAUNCHER_JOB_FILE=PULSARTIMING/...
        GWBsimDataSKA_MaxPhase_JOBFILE.txt
        (Name of the LAUNCHER job file.)



$LAUNCHER_DIR/paramrun
```

To submit the job, the command is:  
```
sbatch LAUNCHER_SLURM_FILE_PTA.slurm  
```
