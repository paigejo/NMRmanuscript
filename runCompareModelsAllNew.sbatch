#!/bin/bash
#SBATCH --job-name runCompareModelsNew      # Set a name for your job. This is especially useful if you have multiple jobs queued.
#SBATCH --partition medium     # Slurm partition to use
#SBATCH --ntasks 1          # Number of tasks to run. By default, one CPU core will be allocated per task
#SBATCH --time 2-00:00        # Wall time limit in D-HH:MM
#SBATCH --mem-per-cpu=20000     # Memory limit for each tasks (in MB)
#SBATCH -o runCompareModelsNew_%A_%a.out    # File to which STDOUT will be written
#SBATCH -e runCompareModelsNew_%A_%a.err    # File to which STDERR will be written
#SBATCH --mail-type=ALL       # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=paigejo@uw.edu # Email to which notifications will be sent
module load R
Rscript --verbose runCompareModelsAllNew.R ${SLURM_ARRAY_TASK_ID} > runCompareModelsAllNew_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.Rout