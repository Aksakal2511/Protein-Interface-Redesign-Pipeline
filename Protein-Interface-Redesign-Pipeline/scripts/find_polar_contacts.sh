#!/bin/bash

PDB1=$1
CHAIN1=$2
PDB2=$3
CHAIN2=$4

# Run tsv_compare_polar_contacts.py
pymol -cq ${TSC_COMPARE_SCRIPT} -- ${PDB1} ${CHAIN1} ${PDB2} ${CHAIN2}

echo "Polar contacts TSV generated"

# Clean the input file to remove carriage return characters
if [ -f ${WORKING_DIR}/common_residues.tsv ]; then
    sed -i 's/\r//' ${WORKING_DIR}/common_residues.tsv
else
    echo "Error: common_residues.tsv file not found"
    exit 1
fi


##GENERATE MDP FILES###
# Define paths
common_tsv="common_residues.tsv"
contacts_tsv="polar_contacts.tsv"
mdp_template="${MDP_DIR}/md_template.mdp"
npt_umbrella_template="${MDP_DIR}/npt_umbrella_template.mdp"
md_umbrella_template="${MDP_DIR}/md_umbrella_template.mdp"

# Function to update MDP file with chain information
update_mdp() {
    local pdb_id="$1"
    local chain1="$2"
    local chain2="$3"
    local template="$4"
    local output_file="$5"

    # Read the template MDP file and replace the placeholders
    sed -e "s/pull_group1_name\s*=.*/pull_group1_name        = Chain$chain1/" \
        -e "s/pull_group2_name\s*=.*/pull_group2_name        = Chain$chain2/" \
        $template > $output_file

    echo "Generated MDP file: $output_file" >&2
}

# Function to process a PDB ID and generate the MDP file
process_pdb() {
    local pdb_id=$1
    local primary_chain=$2

    # Count interactions for each chain
    declare -A chain_interactions
    while IFS=$'\t' read -r file chain residue_num residue_name; do
        base_file=${file%_Contacts}
        if [[ $base_file == $pdb_id ]]; then
            if [[ $chain != $primary_chain ]]; then
                chain_interactions[$chain]=$((chain_interactions[$chain] + 1))
            fi
        fi
    done < $contacts_tsv

    # Find the chain with the highest number of interactions
    max_interactions=0
    interacting_chain=""
    for chain in "${!chain_interactions[@]}"; do
        if (( chain_interactions[$chain] > max_interactions )); then
            max_interactions=${chain_interactions[$chain]}
            interacting_chain=$chain
        fi
    done

    # Debugging check for interacting chain
    if [[ -z $interacting_chain ]]; then
        echo "Error: Interacting chain not found for PDB ID $pdb_id"
        return
    fi

    # Update the MDP files with the primary chain and interacting chain
    update_mdp "$pdb_id" "$primary_chain" "$interacting_chain" "$mdp_template" "${MDP_DIR}/md_${pdb_id}.mdp"
    update_mdp "$pdb_id" "$primary_chain" "$interacting_chain" "$npt_umbrella_template" "${MDP_DIR}/npt_umbrella_${pdb_id}.mdp"
    update_mdp "$pdb_id" "$primary_chain" "$interacting_chain" "$md_umbrella_template" "${MDP_DIR}/md_umbrella_${pdb_id}.mdp"
}

# Process each PDB ID in the common residues TSV file
declare -A chain1_map
declare -A chain2_map

while IFS=$'\t' read -r pdb_id primary_chain residue; do
    if [[ -z "${chain1_map[$pdb_id]}" ]]; then
        # Process the PDB to get the chains
        interacting_chain=$(process_pdb "$pdb_id" "$primary_chain")
        if [[ $? -eq 0 ]]; then
            # Store the primary and interacting chains
            chain1_map[$pdb_id]="$primary_chain"
            chain2_map[$pdb_id]="$interacting_chain"
        fi
    fi
done < "$common_tsv"

mkdir -p ${WORKING_DIR}/results
mv *cif ${WORKING_DIR}/results
mv *.tsv ${WORKING_DIR}/results
