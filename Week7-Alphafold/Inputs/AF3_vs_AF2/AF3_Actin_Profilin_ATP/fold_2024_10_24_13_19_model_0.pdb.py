# coding: utf-8
# Load CIF file
traj = md.load('/home/sa8200/comp-lab-class/comp-lab-class-2024/Week7-Alphafold/Inputs/AF3_vs_AF2/AF3_Actin_Profilin_ATP/*.cif')

# Save as PDB
traj.save('/home/sa8200/comp-lab-class/comp-lab-class-2024/Week7-Alphafold/Inputs/AF3_vs_AF2/AF3.pdb')
# Import MDTraj
import mdtraj as md

# Load CIF file
traj = md.load('/home/sa8200/comp-lab-class/comp-lab-class-2024/Week7-Alphafold/Inputs/AF3_vs_AF2/AF3_Actin_Profilin_ATP/*.cif')

# Save as PDB
traj.save('/home/sa8200/comp-lab-class/comp-lab-class-2024/Week7-Alphafold/Inputs/AF3_vs_AF2/AF3.pdb')
# Import MDTraj
import mdtraj as md

# Load CIF file
traj = md.load('/home/sa8200/comp-lab-class/comp-lab-class-2024/Week7-Alphafold/Inputs/AF3_vs_AF2/AF3_Actin_Profilin_ATP/fold_2024_10_24_13_19_model_4.cif')

# Save as PDB
traj.save('/home/sa8200/comp-lab-class/comp-lab-class-2024/Week7-Alphafold/Inputs/AF3_vs_AF2/AF3.pdb')
# Import MDTraj
import mdtraj as md

# Load CIF file
traj = md.load('/home/sa8200/comp-lab-class/comp-lab-class-2024/Week7-Alphafold/Inputs/AF3_vs_AF2/AF3_Actin_Profilin_ATP/fold_2024_10_24_13_19_model_4.cif')

# Save as PDB
traj.save('/home/sa8200/comp-lab-class/comp-lab-class-2024/Week7-Alphafold/Inputs/AF3_vs_AF2/AF3.pdb')
import mdtraj as md
import matplotlib.pyplot as plt
import numpy as np
# Load the AF2 and AF3 predictions
af2_traj = md.load('/home/sa8200/comp-lab-class/comp-lab-class-2024/Week7-Alphafold/Inputs/AF3_vs_AF2/AF2_Actin_Profilin/Actin_Profilin_923b1_0/Actin_Profilin_923b1_0_unrelaxed_rank_005_alphafold2_multimer_v3_model_2_seed_000.pdb'')
af3_traj = md.load('/home/sa8200/comp-lab-class/comp-lab-class-2024/Week7-Alphafold/Inputs/AF3_vs_AF2/AF3_Actin_Profilin_ATP/fold_2024_10_24_13_19_model_0.cif')
# Ensure both trajectories have the same number of atoms
#if af2_traj.n_atoms != af3_traj.n_atoms:
    #raise ValueError("The two trajectories do not have the same number of atoms.")
af2_atoms = af2_traj.topology.select('chainid 0')
af3_atoms = af3_traj.topology.select('chainid 0')
# Ensure both selections have the same number of atoms
if len(af2_atoms) != len(af3_atoms):
    raise ValueError("The selected atoms do not have the same number of atoms.")
# Align using the selected atoms
af2_traj.superpose(af3_traj, atom_indices=af2_atoms)
# Compute distances between corresponding atoms in the two models
distances = np.linalg.norm(af2_traj.xyz[:, af2_atoms, :] - af3_traj.xyz[:, af3_atoms, :], axis=2)
# Compute average distance for each residue
residue_distances = []
for residue in af2_traj.topology.residues:
    if residue.chain.index == 0:  # Only consider residues in chain A
        residue_indices = [atom.index for atom in residue.atoms if atom.index in af2_atoms]
        if residue_indices:
            residue_distance = np.mean(distances[:, np.isin(af2_atoms, residue_indices)])
            residue_distances.append(residue_distance)
# Plot the distances
plt.figure(figsize=(10, 6))
plt.bar(range(len(residue_distances)), residue_distances)
plt.xlabel('Residue Number')
plt.ylabel('Distance (nm)')
plt.title('Distance between AF2 and AF3 predictions')
plt.show()
import mdtraj as md
import matplotlib.pyplot as plt
import numpy as np
# Load the AF2 and AF3 predictions
af2_traj = md.load('/home/sa8200/comp-lab-class/comp-lab-class-2024/Week7-Alphafold/Inputs/AF3_vs_AF2/AF2_Actin_Profilin/Actin_Profilin_923b1_0/Actin_Profilin_923b1_0_unrelaxed_rank_005_alphafold2_multimer_v3_model_2_seed_000.pdb')
af3_traj = md.load('/home/sa8200/comp-lab-class/comp-lab-class-2024/Week7-Alphafold/Inputs/AF3_vs_AF2/AF3_Actin_Profilin_ATP/fold_2024_10_24_13_19_model_0.cif')
# Ensure both trajectories have the same number of atoms
#if af2_traj.n_atoms != af3_traj.n_atoms:
    #raise ValueError("The two trajectories do not have the same number of atoms.")
af2_atoms = af2_traj.topology.select('chainid 0')
af3_atoms = af3_traj.topology.select('chainid 0')
# Ensure both selections have the same number of atoms
if len(af2_atoms) != len(af3_atoms):
    raise ValueError("The selected atoms do not have the same number of atoms.")
# Align using the selected atoms
af2_traj.superpose(af3_traj, atom_indices=af2_atoms)
# Compute distances between corresponding atoms in the two models
distances = np.linalg.norm(af2_traj.xyz[:, af2_atoms, :] - af3_traj.xyz[:, af3_atoms, :], axis=2)
# Compute average distance for each residue
residue_distances = []
for residue in af2_traj.topology.residues:
    if residue.chain.index == 0:  # Only consider residues in chain A
        residue_indices = [atom.index for atom in residue.atoms if atom.index in af2_atoms]
        if residue_indices:
            residue_distance = np.mean(distances[:, np.isin(af2_atoms, residue_indices)])
            residue_distances.append(residue_distance)
# Plot the distances
plt.figure(figsize=(10, 6))
plt.bar(range(len(residue_distances)), residue_distances)
plt.xlabel('Residue Number')
plt.ylabel('Distance (nm)')
plt.title('Distance between AF2 and AF3 predictions')
plt.show()
import mdtraj as md
import matplotlib.pyplot as plt
import numpy as np
# Load the AF2 and AF3 predictions
af2_traj = md.load('/home/sa8200/comp-lab-class/comp-lab-class-2024/Week7-Alphafold/Inputs/AF3_vs_AF2/AF2_Actin_Profilin/Actin_Profilin_923b1_0/Actin_Profilin_923b1_0_unrelaxed_rank_005_alphafold2_multimer_v3_model_2_seed_000.pdb')
af3_traj = md.load('/home/sa8200/comp-lab-class/comp-lab-class-2024/Week7-Alphafold/Inputs/AF3_vs_AF2/AF3_Actin_Profilin_ATP/*.cif')
# Ensure both trajectories have the same number of atoms
#if af2_traj.n_atoms != af3_traj.n_atoms:
    #raise ValueError("The two trajectories do not have the same number of atoms.")
af2_atoms = af2_traj.topology.select('chainid 0')
af3_atoms = af3_traj.topology.select('chainid 0')
# Ensure both selections have the same number of atoms
if len(af2_atoms) != len(af3_atoms):
    raise ValueError("The selected atoms do not have the same number of atoms.")
# Align using the selected atoms
af2_traj.superpose(af3_traj, atom_indices=af2_atoms)
# Compute distances between corresponding atoms in the two models
distances = np.linalg.norm(af2_traj.xyz[:, af2_atoms, :] - af3_traj.xyz[:, af3_atoms, :], axis=2)
# Compute average distance for each residue
residue_distances = []
for residue in af2_traj.topology.residues:
    if residue.chain.index == 0:  # Only consider residues in chain A
        residue_indices = [atom.index for atom in residue.atoms if atom.index in af2_atoms]
        if residue_indices:
            residue_distance = np.mean(distances[:, np.isin(af2_atoms, residue_indices)])
            residue_distances.append(residue_distance)
# Plot the distances
plt.figure(figsize=(10, 6))
plt.bar(range(len(residue_distances)), residue_distances)
plt.xlabel('Residue Number')
plt.ylabel('Distance (nm)')
plt.title('Distance between AF2 and AF3 predictions')
plt.show()
# In PyMOL
load /home/sa8200/comp-lab-class/comp-lab-class-2024/Week7-Alphafold/Inputs/AF3_vs_AF2/AF3_Actin_Profilin_ATP/*.cif
save /home/sa8200/comp-lab-class/comp-lab-class-2024/Week7-Alphafold/Inputs/AF3_vs_AF2/AF3_Actin_Profilin_ATP/ATP.pdb
# In PyMOL
load /home/sa8200/comp-lab-class/comp-lab-class-2024/Week7-Alphafold/Inputs/AF3_vs_AF2/AF3_Actin_Profilin_ATP/fold_2024_10_24_13_19_model_0.cif

save /home/sa8200/comp-lab-class/comp-lab-class-2024/Week7-Alphafold/Inputs/AF3_vs_AF2/AF3_Actin_Profilin_ATP/ATP.pdb
get_ipython().run_line_magic('load', '/home/sa8200/comp-lab-class/comp-lab-class-2024/Week7-Alphafold/Inputs/AF3_vs_AF2/AF3_Actin_Profilin_ATP/fold_2024_10_24_13_19_model_0.cif')
