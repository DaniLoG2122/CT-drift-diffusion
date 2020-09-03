#!/bin/bash

###############################################
# USAGE:

if [ $# -ne 1 ]; then
	    echo $0: USAGE: build_array.sh inputs_directory
	    exit 1
fi


echo "Input Files are located in: $1"

################################################
# CLEANING:

# Check if previous ranges.txt or jobscript files exists and remove they do
if [ -f ./ranges.txt ] && [[ -f ./jobscript* ]]; then
	echo "ranges.txt exists, removing"
	rm ./ranges.txt ./jobscript*
fi

################################################
# CALCULATE JOB SIZE:

# List number of files in the directory
FILES=$(ls -1 $1 | wc -l )
echo "FILES:$FILES"

# Cases per job, 700 cases at 2min each. 720min in 24 hour jobs.
CASES_PER_JOB=45

# Calculate number of jobs
ARRAY_JOBS=$(( $FILES / $CASES_PER_JOB))
REM="$(( $FILES % $CASES_PER_JOB))"

if [ "$REM" -gt 0 ]; then
	ARRAY_JOBS=$(( 1 + ARRAY_JOBS ))
	echo "Array size: $ARRAY_JOBS"
else
	ARRAY_JOBS=$ARRAY_JOBS
	echo "Array size: $ARRAY_JOBS"
fi

################################################
# CALCULATE RANGES:

for i in $(seq 1 $ARRAY_JOBS); do
	case "$i" in
		1)
			VAR1="JOB_${i}_START"
			declare $VAR1="1"
			echo ${!VAR1} > ranges.txt
			VAR1="JOB_${i}_END"
			declare $VAR1="$((1 * $CASES_PER_JOB))"
			sed -i '$s/$/ - '$(echo ${!VAR1})'/' ranges.txt
			;;
		"$ARRAY_JOBS")
			if [ "$REM" -gt 1 ]; then 
				VAR1="JOB_${i}_START"
				declare $VAR1="$(($(($FILES - $REM)) + 1))"
				echo ${!VAR1} >> ranges.txt
				VAR1="JOB_${i}_END"
				declare $VAR1="$FILES"
				sed -i '$s/$/ - '$(echo ${!VAR1})'/' ranges.txt
			else
				VAR1="JOB_${i}_START"
				declare $VAR1="$(($(( $((${i} - 1)) * $CASES_PER_JOB))+1))"
				echo ${!VAR1} >> ranges.txt
				VAR1="JOB_${i}_END"
				declare $VAR1="$((${i} * $CASES_PER_JOB))"
				sed -i '$s/$/ - '$(echo ${!VAR1})'/' ranges.txt
			fi
			;;

		*)
			VAR1="JOB_${i}_START"
			declare $VAR1="$(($(( $(( ${i} - 1)) * $CASES_PER_JOB))+1))"
			echo ${!VAR1} >> ranges.txt
			VAR1="JOB_${i}_END"
			declare $VAR1="$(( ${i} * $CASES_PER_JOB))"
			sed -i '$s/$/ - '$(echo ${!VAR1})'/' ranges.txt
			;;

	esac
done

################################################
# WRITE ARRAY JOB:

cat <<"EOF" > jobscript_1.sh
#!/bin/bash
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=1:mem=4gb
	
module load matlab/R2020a
export MATLAB_PREFDIR=$TMPDIR/prefs

cd $TMPDIR

cp $PBS_O_WORKDIR/*.csv $PBS_O_WORKDIR/*.m $PBS_O_WORKDIR/*.txt $TMPDIR

export RANGE=$(head -n $PBS_ARRAY_INDEX ranges.txt | tail -1)
export START=$(echo $RANGE | cut -f 1 -d " ")
export END=$(echo $RANGE | cut -f 3 -d " ")

for i in $(seq $START $END); do
	# run program;
	echo "Running file main$i.m"
	cp $PBS_O_WORKDIR/Mains/Int1/main$i.m ./
	matlab -nosplash -nodesktop -nodisplay < ./main$i.m;
	cp $TMPDIR/results$i.mat /rdsgpfs/general/user/dl719/home/Results/Int1


done

EOF

# Add array size to jobscript
sed -i '4i#PBS -J 1-'$(echo ${ARRAY_JOBS})'' jobscript_1.sh

