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
#SBATCH -J TEST_Run
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p normal
#SBATCH -o PTA_Apr12_2017_3.o%j
#SBATCH -e PTA_Apr12_2017_3.e%j
#SBATCH -t 3:00:00
#SBATCH --mail-user=yiqian.qian01@utrgv.edu
#SBATCH --mail-type=all
#          <------ Account String ----->
# <--- (Use this ONLY if you have MULTIPLE accounts) --->
##SBATCH -A
#------------------------------------------------------

export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=PULSARTIMING/TEST_JOBFILE.txt

$LAUNCHER_DIR/paramrun
