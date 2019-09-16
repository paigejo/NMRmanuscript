#!/bin/bash
#SBATCH --job-name runValidationEd      # Set a name for your job. This is especially useful if you have multiple jobs queued.
#SBATCH --partition largemem    # Slurm partition to use
#SBATCH --ntasks 4          # Number of tasks to run. By default, one CPU core will be allocated per task
#SBATCH --time 5-00:00        # Wall time limit in D-HH:MM
#SBATCH --mem-per-cpu=32000     # Memory limit for each tasks (in MB)
#SBATCH -o runValidationEd_%A_%a.out    # File to which STDOUT will be written
#SBATCH -e runValidationEd_%A_%a.err    # File to which STDERR will be written
#SBATCH --mail-type=ALL       # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=paigejo@uw.edu # Email to which notifications will be sent
module load R
Rscript --verbose runValidationEd.R ${SLURM_ARRAY_TASK_ID} > runValidationEd_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.Rout