#!/bin/bash
#PBS -P v29
#PBS -l walltime=00:05:00
#PBS -l vmem=500MB
#PBS -l ngpus=2
#PBS -l ncpus=1
#PBS -l jobfs=1GB
#PBS -wd
#PBS -j oe
#PBS -l other=physmem

module load nvidia
module load python/2.7.3
module load cuda/4.2.9
module load cuda-cvp
module load gcc/4.4.4
module load boost/1.46.1
module load netcdf/4.2.1.1


export COMPUTE_PROFILE=1 
export COMPUTE_PROFILE_CSV=1
compute_profile=profile_compute_fluxes
export COMPUTE_PROFILE_LOG="$compute_profile.csv"
export COMPUTE_PROFILE_CONFIG=".cp_config"

compute_fluxes.py > result_compute_fluxes
