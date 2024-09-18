# Protein-Interface-Redesign-Pipeline
Automated multi-step analytical pipeline to analyse protein interfaces and measure binding affinities, through MD simulations and umbrella sampling.
Author: Ozan Aksakal

Overview
This project implements an automated multi-step pipeline designed to analyze the interface residues of a protein (Prot_A) in complex with two other proteins (Prot_B and Prot_C). The goal is to prepare Prot_A for redesign, specifically to abrogate its binding with Prot_B while maintaining or enhancing its binding affinity with Prot_C.

The pipeline includes:
- Protein mutagenesis
- Molecular Dynamics (MD) simulations
- Umbrella sampling to calculate binding affinities
- Polar contact analysis between the proteins

The pipeline is executed in a High-Performance Computing (HPC) environment using GROMACS, PyMOL, and Python scripts.

Example Use Case:
- Prot_A: SARS-CoV-2 spike protein
- Prot_B: Antibody P5A-3C8
- Prot_C: ACE2 receptor

Requirements
- HPC environment: Compatible with SLURM workload manager
- Software:
  - GROMACS 2021.5
  - PyMOL CLI (v2.5.0)
  - Python 3.7
  - Singularity/Apptainer 1.3.0
- Input: Two PDB files containing Prot_A with Prot_B and Prot_C. The pipeline will compare the interface residues between Prot_A and Prot_B, and Prot_A and Prot_C.

Directory Structure
- master_script.sh: Main pipeline wrapper script
- scripts/: Contains sub-scripts for mutagenesis, MD simulations, and analysis
- results/: Output directory for final results
- MDP/: GROMACS MDP configuration files for the simulations

Workflow
1. Polar Contact Analysis
- tsv_compare_polar_contacts.py: Identifies key residues involved in polar contacts between Prot_A and Prot_B/Prot_C using PyMOL.
- Results are stored in polar_contacts.tsv and common_residues.tsv for downstream analysis.

2. Mutagenesis and MD Simulations
- mutate.py: Performs mutagenesis on Prot_A to create variants.
- MD simulations: Uses GROMACS to run simulations of Prot_A with Prot_B and Prot_C.
- Box generation: The solvated box is dynamically adjusted to the z-dimension for umbrella sampling.
- Umbrella sampling: Performed along a defined reaction coordinate using WHAM to calculate binding affinities.

3. Umbrella Sampling
- Umbrella sampling setup: A dynamic set of configurations is generated from the MD trajectories to simulate the pulling of Prot_A from Prot_B or Prot_C.
- Analysis: WHAM is used to calculate the free energy profile, resulting in profile.xvg and hist.xvg files for further analysis.

Execution
To run the full pipeline:
1. Prepare Input Data: Ensure you have two PDB files of Prot_A complexed with Prot_B and Prot_C, respectively.
2. Run the Master Script:
   ./master_script.sh PDB1 CHAIN1 PDB2 CHAIN2 MUTATION
   
Example: ./master_script.sh 6M0J E 7Z0X R GLU

3. Review Results: The final results including the free energy profiles and histograms will be found in the results/ directory.

Notes
- This pipeline is computationally intensive, particularly the umbrella sampling, and is designed for HPC environments.
- The umbrella sampling simulations have been set to 1 ns for demonstration purposes. For more accurate results, it is recommended to run for longer durations.
