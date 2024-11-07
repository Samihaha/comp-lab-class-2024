#!/bin/bash
#SBATCH --job-name=anneal_system      # Job name
#SBATCH --output=anneal_system.out     # Output file
#SBATCH --error=anneal_system.err      # Error file
#SBATCH --time=01:00:00                # Time limit (1 hour)
#SBATCH --nodes=1                      # Run on a single node
#SBATCH --ntasks-per-node=4            # Number of tasks per node
#SBATCH --cpus-per-task=1              # Number of CPUs per task


# Step 1: Create the system
echo "Creating the system..."
mpirun lmp -var configfile /home/sa8200/comp-lab-class/comp-lab-class-2024/Week9-LAMMPS-SupercooledLiquids/Inputs/n360/kalj_n360_create.lmp -var id 1 -in /home/sa8200/comp-lab-class/comp-lab-class-2024/Week9-LAMMPS-SupercooledLiquids/Inputs/create_3d_binary.lmp

# Define the list of temperatures for annealing
TEMPERATURES=(1.5 1.0 0.9 0.8 0.7 0.65 0.6 0.55 0.5 0.475)

# Step 2: Anneal the system by progressively equilibrating to lower temperatures
echo "Starting annealing process..."
for T in "${TEMPERATURES[@]}"; do
    echo "Equilibrating at temperature $T..."
    mpirun lmp -var configfile /home/sa8200/comp-lab-class/comp-lab-class-2024/Week9-LAMMPS-SupercooledLiquids/Inputs/n360/kalj_n360_T${T}.lmp -var id 1 -in /home/sa8200/comp-lab-class/comp-lab-class-2024/Week9-LAMMPS-SupercooledLiquids/Inputs/anneal_3d_binary.lmp
done

echo "Annealing process completed."
