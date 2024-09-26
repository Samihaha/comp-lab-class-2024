#!/bin/bash

#SBATCH --job-name=md_simulation        # Job name

#SBATCH --ntasks=1                      # Run a single task

#SBATCH --cpus-per-task=12               # Number of CPU cores per task

#SBATCH --mem=8GB                       # Memory allocation

#SBATCH --time=12:00:00                 # Time limit in hours, minutes, and seconds

#SBATCH --output=md_simulation.out      # Output file

#SBATCH --mail-user=sa8200@nyu.edu  # Email address for notifications

#SBATCH --mail-type=ALL                 # Email notifications (ALL, NONE, END, FAIL)



# Load the necessary module for GROMACS

module load gromacs/2020.4



# Run the MD simulation with GROMACS

gmx_mpi mdrun -deffnm md


