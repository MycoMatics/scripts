
# find line in files with extension job.sh and change into
find . -type f -name "*job.sh" -exec sed -i '/^#PBS -o /s/.*/#PBS -o stdout.$PBS_JOBNAME.$PBS_JOBID/' {} \;


In this updated command:

	#find . -type f -name "*job.sh" still searches for files with the "*job.sh" pattern.
	#-exec sed -i executes the sed command in-place for each matching file.
	#/^#PBS -o /s/.*/#PBS -o stdout.$PBS_JOBNAME.$PBS_JOBID/ is used to search for lines that start with 
	#PBS -o and replace the entire line with #PBS -o stdout.$PBS_JOBNAME.$PBS_JOBID.