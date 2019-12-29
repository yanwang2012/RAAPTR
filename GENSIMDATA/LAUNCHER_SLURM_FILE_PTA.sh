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
#SBATCH -J Mauritius_Uniform
#SBATCH -N 167
#SBATCH -n 501
#SBATCH -p normal
#SBATCH -o MauritiusUni_Apr22_2019.o%j
#SBATCH -e MauritiusUni_Apr22_2019.e%j
#SBATCH -t 3:00:00
#SBATCH --mail-user=yiqian.qian01@utrgv.edu
#SBATCH --mail-type=all
#          <------ Account String ----->
# <--- (Use this ONLY if you have MULTIPLE accounts) --->
##SBATCH -A
#----------:q--------------------------------------------
export LAUNCHER_WORKDIR=/work/05884/qyqstc/lonestar/PULSARTIMING/RAAPTR/MxAvPhaseC
export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=/work/05884/qyqstc/lonestar/PULSARTIMING/RAAPTR/Maricius_JOBFILE2.txt
$LAUNCHER_DIR/paramrun
