#!/bin/bash
#SBATCH --job-name=production_run          # Job name
#SBATCH --output=production_run_%j.out     # Output file, with job ID for easy tracking
#SBATCH --error=production_run_%j.err      # Error file, with job ID for easy tracking
#SBATCH --time=00:20:00                    # Set a shorter time limit (adjust if necessary)
#SBATCH --nodes=1                          # Run on a single node
#SBATCH --ntasks-per-node=4                # Number of tasks per node
#SBATCH --cpus-per-task=1                  # Number of CPUs per task



# Define temperature list for production runs
TEMPERATURES=(1.5 1.0 0.9 0.8 0.7 0.65 0.6 0.55 0.5 0.475)

# Production simulations loop
echo "Starting production simulations..."
for T in "${TEMPERATURES[@]}"; do
        echo "Starting production run at temperature $T..."
        mpirun lmp -var configfile ../Inputs/n360/kalj_n360_T${T}.lmp -var id 1 -in ../Inputs/production_3d_binary.lmp
        # Optionally create a file indicating the completion of this temperature
        touch "production_T${T}_complete.txt"
done

echo "Production simulations completed."
