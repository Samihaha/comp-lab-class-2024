#!/bin/bash
#SBATCH --job-name=kalj_simulation
#SBATCH --output=kalj_output.log
#SBATCH --error=kalj_error.log
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
