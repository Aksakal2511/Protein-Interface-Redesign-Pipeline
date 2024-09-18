#!/bin/bash
#SBATCH -o output.log
#SBATCH -e error.log
#SBATCH -J v7master
#SBATCH --partition=compute
#SBATCH --ntasks=40
#SBATCH --ntasks-per-node=40
#SBATCH --time=72:00:00 # Time limit
#SBATCH --exclusive
#SBATCH --account=scw1631

start="$(date +%s)"

module load pymol/2.5.0-cli
module load python/3.7.0
module load singularity/apptainer/1.3.0
module load mpi/intel/2020/4
module load compiler/gnu/9/2.0
module load compiler/intel/2020/4
module load lustre_getcwd_fix/LU-9735
module load gromacs/2021.5-single

# Base directory
export WORKING_DIR=$(pwd)
export SCRIPTS_DIR="${WORKING_DIR}/scripts"
export MDP_DIR="${SCRIPTS_DIR}/MDP"

# Script paths
export MUTATE_SCRIPT_PATH="${SCRIPTS_DIR}/mutate.py"
export TSC_COMPARE_SCRIPT="${SCRIPTS_DIR}/tsv_compare_polar_contacts.py"
export PBCATOM_SCRIPT="${SCRIPTS_DIR}/pbc_atoms.py"

# Usage: ./master_script.sh $1 $2 $3 $4 $5

# Input variables
PDB1=$1
CHAIN1=$2
PDB2=$3
CHAIN2=$4
MUTATION=$5

# Subscript 1: Find polar contacts
${SCRIPTS_DIR}/find_polar_contacts.sh ${PDB1} ${CHAIN1} ${PDB2} ${CHAIN2}

# Subscript 2: Mutate and run MD
${SCRIPTS_DIR}/mutate_and_md.sh ${MUTATION}

# Subscript 3: Generate plots
python3 ${SCRIPTS_DIR}/plot_pmf.py

echo "Master script completed."
