#!/bin/bash

MUTATION=$1

# Process common_residues.tsv file
while IFS=$'\t' read -r pdb chain residue_number
do
    for x in ${MUTATION}
    do
        # Create and enter the main directory if it does not exist
        mkdir -p "MD/${pdb}/${x}${residue_number}"
        cd "MD/${pdb}/${x}${residue_number}"

        # Execute the pymol command
        echo "Running pymol -qc ${MUTATE_SCRIPT_PATH} -d 'mutate(\"${pdb}\", \"${chain}\", \"${residue_number}\", \"${x}\", 1)'"
        if ! pymol -qc "${MUTATE_SCRIPT_PATH}" -d "mutate(\"${pdb}\", \"${chain}\", \"${residue_number}\", \"${x}\", 1)"; then
            echo "Error: pymol command failed for ${pdb} ${chain}/${residue_number}/ ${x}"
            exit 1
        fi

        echo "Mutation Complete ${pdb}/${x}${residue_number}"   
        sleep 2
        
        # Define chains and extract atom info for indexing [chain_atom_ranges.txt]
        grep '^ATOM' *.pdb > clean.pdb
        chains=$(grep "^ATOM" clean.pdb | awk '{print $5}' | grep -E '^[A-Za-z]+$' | sort | uniq)
        declare -A chain_atom_ranges
        for chain in $chains; do
            grep "^ATOM.* $chain " clean.pdb > chain${chain}.pdb
            rm chain[0-9]*.pdb
            start_atom=$(grep "^ATOM" chain${chain}.pdb | head -n 1 | awk '{print $2}')
            end_atom=$(grep "^ATOM" chain${chain}.pdb | tail -n 1 | awk '{print $2}')
            chain_atom_ranges[$chain]="${start_atom}-${end_atom}"
        done

        for chain in "${!chain_atom_ranges[@]}"; do
            echo "Chain $chain: ${chain_atom_ranges[$chain]}"
        done > chain_atom_ranges.txt


        # Generate topology
        { echo "8" ; echo "1"; } | gmx_mpi pdb2gmx -f clean.pdb -o processed.gro -ignh

        echo "Topology Generated ${pdb}/${x}${residue_number}"
        sleep 2

        # Automate box dimension and add solvate
        buffer=2.0
        box_size=$(tail -n 1 processed.gro)
        
        x_size=$(echo $box_size | awk '{print $1}')
        y_size=$(echo $box_size | awk '{print $2}')
        z_size=$(echo $box_size | awk '{print $3}')
        
        new_x=$(echo "$x_size + 2*$buffer" | bc)
        new_y=$(echo "$y_size + 2*$buffer" | bc)
        new_z=$(echo "2*$z_size + 5.0" | bc)
        
        x_center=$(echo "$new_x / 2" | bc -l)
        y_center=$(echo "$new_y / 2" | bc -l)
        z_center=$(echo "$new_z / 4" | bc -l)
        
        gmx_mpi editconf -f processed.gro -o newbox.gro -center $x_center $y_center $z_center -box $new_x $new_y $new_z
        gmx_mpi solvate -cp newbox.gro -cs spc216.gro -o solv.gro -p topol.top

        echo "Box Defined and Solvated ${pdb}/${x}${residue_number}"
        sleep 2

        # Add ions
        gmx_mpi grompp -f ${MDP_DIR}/ions.mdp -c solv.gro -p topol.top -o ions.tpr
        echo 13 | gmx_mpi genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral #15 w BV
        
        echo "Ions Added ${pdb}/${x}${residue_number}"
        sleep 2

        #Indexing
        atom_ranges_file="chain_atom_ranges.txt"
        index_file="index.ndx"
        temp_ndx_file="temp_index.ndx"
        rename_input_file="rename_input.txt"
        > $rename_input_file
        
        declare -A chain_map
        while read -r line; do
            chain=$(echo $line | awk '{print $2}' | tr -d ':')
            range=$(echo $line | awk '{print $3}')
            chain_map["$range"]="$chain"
        done < $atom_ranges_file

gmx_mpi make_ndx -f solv_ions.gro -o $temp_ndx_file <<EOF
$(for range in "${!chain_map[@]}"; do
    echo "a ${range//-/ - }"
done)
q
EOF

        {
            group_number=18   #Automate next ;22 w BV
            for range in "${!chain_map[@]}"; do
                chain="${chain_map[$range]}"
                echo "name $group_number chain$chain"
                echo ""
                group_number=$((group_number + 1))
            done
            echo "q"
        } > $rename_input_file
        echo "Contents of rename_input_file:"
        cat $rename_input_file
        
        gmx_mpi make_ndx -f solv_ions.gro -n $temp_ndx_file -o $index_file < $rename_input_file

       # Identifying PBC atoms
        pull_group1_name=$(grep 'pull_group1_name' ${MDP_DIR}/md_${pdb}.mdp | awk '{print $3}' | sed 's/Chain//')
        pull_group2_name=$(grep 'pull_group2_name' ${MDP_DIR}/md_${pdb}.mdp | awk '{print $3}' | sed 's/Chain//')
        
        echo "Pull group 1 name: $pull_group1_name"
        echo "Pull group 2 name: $pull_group2_name"

        pull_group1_range=$(grep "Chain $pull_group1_name:" chain_atom_ranges.txt | awk '{print $3}')
        pull_group2_range=$(grep "Chain $pull_group2_name:" chain_atom_ranges.txt | awk '{print $3}')

        echo "Atom range for group 1: $pull_group1_range"
        echo "Atom range for group 2: $pull_group2_range"

        calculate_com() {
            local atom_range=$1
            local PDB_file=$2
            local x_sum=0
            local y_sum=0
            local z_sum=0
            local count=0
            
            awk -v start_atom="${atom_range%-*}" -v end_atom="${atom_range#*-}" \
                '$2 >= start_atom && $2 <= end_atom && $1 == "ATOM" {
                    x_sum += $7;
                    y_sum += $8;
                    z_sum += $9;
                    count += 1
                }
                END {
                    print x_sum / count, y_sum / count, z_sum / count
                }' "$pdb_file"
            }
        
        find_closest_atom() {
            local com=($1)
            local atom_range=$2
            local pdb_file=$3

            awk -v start_atom="${atom_range%-*}" -v end_atom="${atom_range#*-}" \
                -v com_x="${com[0]}" -v com_y="${com[1]}" -v com_z="${com[2]}" \
                'BEGIN { min_dist = 1e10 }
                $2 >= start_atom && $2 <= end_atom && $1 == "ATOM" {
                    dist = sqrt(($7 - com_x)^2 + ($8 - com_y)^2 + ($9 - com_z)^2);
                    if (dist < min_dist) {
                        min_dist = dist;
                        closest_atom = $2
                    }
                }
                END {
                    print closest_atom
                }' "$pdb_file"
        }
        
        pdb_file="clean.pdb"
        
        com1=$(calculate_com "$pull_group1_range" "$pdb_file")
	echo "COM for group 1: $com1"
        com2=$(calculate_com "$pull_group2_range" "$pdb_file")
	echo "COM for group 2: $com2"

        pbcatom1=$(find_closest_atom "$com1" "$pull_group1_range" "$pdb_file")
 	echo "Closest atom to COM for group 1: $pbcatom1" 
        pbcatom2=$(find_closest_atom "$com2" "$pull_group2_range" "$pdb_file")
	echo "Closest atom to COM for group 2: $pbcatom2" 

        mdp_files=("${MDP_DIR}/md_${pdb}.mdp" "${MDP_DIR}/npt_umbrella_${pdb}.mdp" "${MDP_DIR}/md_umbrella_${pdb}.mdp")

        for mdp_file in "${mdp_files[@]}"; do
            echo "Updating $mdp_file with pbcatom values"
            echo "pull-group1-pbcatom   = $pbcatom1" >> $mdp_file
            echo "pull-group2-pbcatom   = $pbcatom2" >> $mdp_file
        done

        # Energy minimisation
        mkdir EM
        cd EM

        gmx_mpi grompp -f ${MDP_DIR}/min.mdp -c ../solv_ions.gro -p ../topol.top -o em.tpr -n ../index.ndx
        mpirun -np 40 gmx_mpi mdrun -v -deffnm em

        echo "Energy Minimised ${pdb}/${x}${residue_number}"
        sleep 2

        # Equilibration
         echo "Starting constant volume equilibration..."

         cd ../
         mkdir NPT
         cd NPT
        
         gmx_mpi grompp -f ${MDP_DIR}/npt.mdp -c ../EM/em.gro -r ../EM/em.gro -p ../topol.top -o npt.tpr -n ../index.ndx
         mpirun -np 40 gmx_mpi mdrun -v -deffnm npt

         echo "Equilibration Phase 2 Complete ${pdb}/${x}${residue_number}"
         sleep 2

        # Production MD - Change .MDP Accordingly to Length of Simulation
         echo "Starting production MD simulation..."

         cd ../
         mkdir Production_MD
         cd Production_MD

         gmx_mpi grompp -f ${MDP_DIR}/md_${pdb}.mdp -c ../NPT/npt.gro -t ../NPT/npt.cpt -p ../topol.top -n ../index.ndx -o md.tpr
         mpirun -np 40 gmx_mpi mdrun -v -deffnm md -pf pullf.xvg -px pullx.xvg

         echo "Simulation Complete ${pdb}/${x}${residue_number}"
        sleep 2
        
        # Umbrella sampling
        echo "Starting production Umbrella sampling..."
        
        mkdir umbrella
        cd umbrella
        
        echo 0 | gmx_mpi trjconv -s ../md.tpr -f ../md.xtc -o conf.gro -sep

        for (( i=0; i<501; i++ ))
        do
            gmx_mpi distance -s ../md.tpr -f conf${i}.gro -n ../../index.ndx -select "com of group 18 plus com of group 19" -oall dist${i}.xvg
        done

        touch summary_distances.dat
        for (( i=0; i<501; i++ ))
        do
            d=`tail -n 1 dist${i}.xvg | awk '{print $2}'`
            echo "${i} ${d}" >> summary_distances.dat
            rm dist${i}.xvg
        done
        
        #Find .gro structures with 0.2nm distance apart and remove all other .gro files
        input_file="summary_distances.dat"
        target_distance=0.2
        start_line=5
        temp_file=$(mktemp)
        awk -v target="$target_distance" -v start="$start_line" '
            NR >= start {
                if (prev == "" || ($2 - prev) >= target) {
                    print $0;
                    prev = $2;
                }
            }
        ' "$input_file" > "$temp_file"
        lines_to_keep=$(awk '{print $1}' "$temp_file")
        patterns=$(echo "$lines_to_keep" | awk '{printf "conf%d.gro\n", $1}')
        pattern_file=$(mktemp)
        echo "$patterns" > "$pattern_file"
        find . -type f -name 'conf*.gro' | grep -v -f "$pattern_file" | xargs rm
        rm "$temp_file" "$pattern_file"
        
        echo "Kept files:"
        echo "$patterns"

        #Simulate configurations (set to 10ns)
        gro_files=$(ls conf*.gro | sort -t'.' -k1.5n)
        counter=1

        for gro_file in $gro_files; do
            echo "Processing $gro_file..."

            # Generate the npt tpr file
            gmx_mpi grompp -f ${MDP_DIR}/npt_umbrella_${pdb}.mdp -c $gro_file -r $gro_file -p ../../topol.top -n ../../index.ndx -o npt${counter}.tpr
            mpirun -np 40 gmx_mpi mdrun -v -deffnm npt${counter}

            # Generate the md tpr file
            gmx_mpi grompp -f ${MDP_DIR}/md_umbrella_${pdb}.mdp -c npt${counter}.gro -t npt${counter}.cpt -p ../../topol.top -n ../../index.ndx -o md${counter}.tpr
            mpirun -np 40 gmx_mpi mdrun -v -deffnm md${counter}

            # Increment the counter
            counter=$((counter + 1))
        done

        echo "Sampling Complete ${pdb}/${x}${residue_number}"
        sleep 2
        
        # Analysis
        > tpr.dat
        > pullf.dat
        
        for file in md*.tpr; do
            echo "$file" >> tpr.dat
        done
        
        for file in md*_pullf.xvg; do
            echo "$file" >> pullf.dat
        done
        
        mkdir -p ${WORKING_DIR}/results
        
        gmx_mpi wham -it tpr.dat -if pullf.dat -o ${WORKING_DIR}/results/profile_${pdb}_${x}_${residue_number}.xvg -hist ${WORKING_DIR}/results/hist_${pdb}_${x}_${residue_number}.xvg -unit kCal
        
        # Go back to the main directory
        cd ${WORKING_DIR}
    done
done < ${WORKING_DIR}/results/common_residues.tsv

