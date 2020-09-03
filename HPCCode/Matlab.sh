#!/bin/bash -l

# Batch script to run a multi-threaded Matlab job on Legion with the upgraded
# software stack under SGE.
#
# Based on openmp.sh by:
# Owain Kenway, Research Computing, 20/Sept/2010
# Updated for RHEL 7, Oct 2015
# Updated for R2016b Jan 2017

# 1. Force bash as the executing shell.
#$ -S /bin/bash

#PBS -l walltime=23:55:00
#PBS -l select=1:ncpus=1:mem=9000Mb
#PBS -J 0-10

# 5. Reserve one Matlab licence - this stops your job starting and failing when no
#    licences are available.
#$ -l matlab=1

# 6. The way Matlab threads work requires Matlab to not share nodes with other
# jobs.
#$ -ac exclusive

# 7. Set the name of the job.
#$ -N Matlab_job1

# 10. Your work *must* be done in $TMPDIR 
cd $TMPDIR

# 11. Copy main Matlab input file and any additional files to TMPDIR

cp ${infile}${PBS_ARRAY_INDEX}.m .
Matlab_infile=`basename ${infile}${PBS_ARRAY_INDEX}.m`
for file in `echo $addinfiles | tr ':' ' '` 
do
  cp $file .
done

# 12. Run Matlab job

module load matlab
echo ""
echo "Running matlab -nosplash -nodisplay < $Matlab_infile ..."
echo ""
matlab -nosplash -nodesktop -nodisplay < $Matlab_infile

# 13. Preferably, tar-up (archive) all output files onto the shared scratch area
cp *.mat /home/
# Make sure you have given enough time for the copy to complete!