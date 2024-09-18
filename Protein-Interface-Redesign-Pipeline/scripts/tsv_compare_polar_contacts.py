import sys
import csv
from pymol import cmd

# Define the distance criteria for H-bonds
distance_criteria = 3.5

# Polar contacts function
def find_contacts_with_prot_A(pdb_id, prot_A_chain, output_prot_A, output_contacts):
    
    # Load the PDB structure
    cmd.fetch(pdb_id, async_=0)

    # Get all chain IDs
    chains = cmd.get_chains(pdb_id)

    # Remove duplicate Prot_A chain ID
    chains.remove(prot_A_chain)

    # Lists to store contacts and residues
    prot_A_residues = []
    contact_residues = []

    for chain in chains:
        # Intra-chain atoms selection in Prot_A within the distance criteria of atoms
        cmd.select(f'contacts_{prot_A_chain}_{chain}', f'chain {prot_A_chain} within {distance_criteria} of chain {chain}')
        cmd.select(f'contacts_{chain}_{prot_A_chain}', f'chain {chain} within {distance_criteria} of chain {prot_A_chain}')

        # Combine both selections to get all contact atoms
        cmd.select(f'contacts_{prot_A_chain}_all', f'contacts_{prot_A_chain}_{chain} or contacts_{chain}_{prot_A_chain}')

        # Select residues in Prot_A that are involved in the contacts
        cmd.select(f'{output_prot_A}_residues', f'byres (contacts_{prot_A_chain}_all and chain {prot_A_chain})')

        # Select residues in the other chain that are involved in the contacts
        cmd.select(f'{output_contacts}_residues', f'byres (contacts_{prot_A_chain}_all and chain {chain})')

        # Function to list residues in a table format
        def list_residues(selection):
            residues = cmd.get_model(selection).atom
            unique_residues = {(atom.chain, atom.resi, atom.resn) for atom in residues}
            return sorted(unique_residues, key=lambda x: (x[0], x[1]))

        # List residues for Prot_A and the other chain
        prot_A_list = list_residues(f'{output_prot_A}_residues')
        contact_list = list_residues(f'{output_contacts}_residues')

        prot_A_residues.extend(prot_A_list)
        contact_residues.extend(contact_list)

    return prot_A_residues, contact_residues

# Input variables
PDB1 = sys.argv[1]
CHAIN1 = sys.argv[2]
PDB2 = sys.argv[3]
CHAIN2 = sys.argv[4]

# Process PDB1
cmd.reinitialize()
chain1_list, contacts1_list = find_contacts_with_prot_A(PDB1, CHAIN1, f'chain_{PDB1}', f'contacts_{PDB1}')

# Process PDB2
cmd.reinitialize()
chain2_list, contacts2_list = find_contacts_with_prot_A(PDB2, CHAIN2, f'chain_{PDB2}', f'contacts_{PDB2}')

# Save residues and contacts to a TSV file
with open('polar_contacts.tsv', 'w', newline='') as file:
    writer = csv.writer(file, delimiter='\t')
    writer.writerow(['Structure', 'Chain', 'Residue Number', 'Residue Name'])

    for chain, resi, resn in chain1_list:
        writer.writerow([PDB1, chain, resi, resn])
    for chain, resi, resn in contacts1_list:
        writer.writerow([f"{PDB1}_Contacts", chain, resi, resn])
    for chain, resi, resn in chain2_list:
        writer.writerow([PDB2, chain, resi, resn])
    for chain, resi, resn in contacts2_list:
        writer.writerow([f"{PDB2}_Contacts", chain, resi, resn])

# Compare residues within Prot_A in both structures
residues_1 = {f"{resn} {resi}" for chain, resi, resn in chain1_list}
residues_2 = {f"{resn} {resi}" for chain, resi, resn in chain2_list}

common_residues = residues_1.intersection(residues_2)
unique_1 = residues_1.difference(residues_2)
unique_2 = residues_2.difference(residues_1)

# Print only the common residues in Prot_A
print(f"\nCommon residues in Prot_A (Chain {CHAIN1} in {PDB1} and Chain {CHAIN2} in {PDB2}):")
if not common_residues:
    print("None found.")
else:
    for residue in common_residues:
        print(residue)

# Save common residues to a TSV file
with open('common_residues.tsv', 'w', newline='') as file:
    writer = csv.writer(file, delimiter='\t')
    for residue in common_residues:
        resn, resi = residue.split()
        writer.writerow([PDB1, CHAIN1, resi])
        writer.writerow([PDB2, CHAIN2, resi])

cmd.quit()

