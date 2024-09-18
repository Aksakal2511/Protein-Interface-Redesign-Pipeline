#!/bin/bash

# Define the files and directories to keep
keep_files=("rm.sh" "master_script.sh" "README.txt")
keep_dirs=("MDP" "scripts")

# Find and delete all other files and directories
for item in $(ls); do
    # Check if the current file or directory is in the keep list
    if [[ ! " ${keep_files[@]} " =~ " ${item} " ]] && [[ ! " ${keep_dirs[@]} " =~ " ${item} " ]]; then
        # Remove the file or directory
        rm -r "$item"
    fi
done

