#!/bin/bash

#SBATCH --job-name=Ndopt
#SBATCH --output=cluster_output/Ndopt%j.out
#SBATCH --error=cluster_output/Ndopt%j.err
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --mem=64GB

# Note: This file is specific to Benoit Pasquier's USC HPC cluster directory structure.
# You will likely need to ensure Julia is available and
# edit the paths to access Julia and the GNOM code

# Load the julia module
export PATH=~/Applications/julia-1.6.2/bin:$PATH
export LD_LIBRARY_PATH=~/Applications/julia-1.6.2/lib

# Cd to the root folder
cd /home/geovault-06/pasquier/Projects/GNOM

# Optimize it!
julia src/Nd_model/setup_and_optimization.jl