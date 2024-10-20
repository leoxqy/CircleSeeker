#!/bin/bash
#SBATCH --job-name=conda_build          # Job name
#SBATCH --partition=64c512g             # Partition to use
#SBATCH --ntasks=1                      # Number of tasks (processes)
#SBATCH --cpus-per-task=48               # Number of CPU cores per task
#SBATCH --time=3:00:00                  # Time limit (hh:mm:ss)
#SBATCH --output=build2.log  # Log file path

conda activate CircleSeeker

# Execute the command
conda build conda-recipe -c conda-forge -c bioconda
