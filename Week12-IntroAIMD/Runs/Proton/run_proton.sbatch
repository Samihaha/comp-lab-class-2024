#!/bin/bash
#SBATCH --job-name=proton_aimd
#SBATCH --output=proton_aimd.out
#SBATCH --error=proton_aimd.err
#SBATCH --time=06:00:00   # Wall time (hh:mm:ss)
#SBATCH --ntasks=12       # Total tasks
#SBATCH --cpus-per-task=1 # CPUs per task
#SBATCH --mem=20GB        # Memory per node

# Load the CP2K module
module load cp2k/openmpi/intel/20201212

# Run the simulation
mpirun -np 12 cp2k.popt -i proton.inp -o proton_aimd_1ps_every0.5fs.log
