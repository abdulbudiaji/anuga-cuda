#!/bin/bash
#PBS -P v29
#PBS -j oe
##PBS -l ngpus=2
##PBS -l ncpus=1
#PBS -l walltime=00:20:00
#PBS -l vmem=1G
#PBS -l other=physmem
#PBS -wd

echo $cmdXe
echo $PBS_JOBID
echo $profile

outfile=`echo $cmdXe | tr ' ' '_' | tr -d ' /;<>&#()'`
jobid=`echo $PBS_JOBID | awk -F. '{print $1}'`
outfile=outputXe.$outfile.$jobid

if [ "$profile" = true ] ; then
	export COMPUTE_PROFILE=1 
	export COMPUTE_PROFILE_CSV=1
	compute_profile=cvp_output_0
	export COMPUTE_PROFILE_LOG="$compute_profile.csv"
	export COMPUTE_PROFILE_CONFIG=".cp_config"
fi

module load nvidia
module load python/2.7.3
module load cuda/4.2.9
module load gcc/4.4.4
module load boost/1.46.1
module load netcdf/4.2.1.1

$cmdXe > $outfile 2>&1

if [ "$profile" = true ] ; then
	cvp_output < $COMPUTE_PROFILE_LOG > $compute_profile.log
	cvp_summarize < $COMPUTE_PROFILE_LOG > $compute_profile.summ
	echo "CVP full output is in $compute_profile.log, summary in $compute_profile.summ" 
fi
