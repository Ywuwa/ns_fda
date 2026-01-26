#!/usr/bin/env bash
# ===========================================================================
# Script for running 4 FDA-schemes under various conditions of vorticity 
# of the ABC flow. Creates residual graphs images (in a shared folder) and  
# text files containing numerical residual values (in a separate folder)
# ===========================================================================
set -e  # script will fail if any error occurs

# Sets of values
FDA_NUMBERS=(1 2 3 4)
VORTICITY_VALUES=(1 3 5 7)

# Create ABC_residual folder
mkdir -p ABC_residual

for i in "${FDA_NUMBERS[@]}"; do
    # Replace a number in config
    sed -i -E "/FDA_number/s/[0-9]+/${i}/" config

    for j in "${VORTICITY_VALUES[@]}"; do
        # Replace value in headers/mesh_n_model.hpp
        sed -i -E "/vorticity/s/[0-9]+/${j}/" headers/mesh_n_model.hpp

        # Start the program
        ./ns_fda_3D

        # Rename residual.txt
        mv residual.txt "residual_fda${i}_${j}.txt"
    done

    # After the inner loop is done
    python3 residual_plot.py

    # Rename residual_plot.png for each j
    for j in "${VORTICITY_VALUES[@]}"; do
        if [[ -f "residual_fda${i}_${j}.txt" ]]; then
            cp residual_plot.png "residual_plot_${i}_${j}.png"
        fi
    done

    # Move each residual_fda*.txt to ABC_residual/
    mv residual_fda${i}_*.txt ABC_residual/
done