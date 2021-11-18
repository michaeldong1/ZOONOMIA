#!/bin/bash

#SBATCH -A project
#SBATCH -p core
#SBATCH -n 1
#SBATCH -J branchh241_%a
#SBATCH -e error_out/branchh241_%j_%a.err
#SBATCH -o error_out/branchh241_%j_%a.out
#SBATCH -t 00:10:00

# MAFspeciesLister version 1.0

### Description of the pipeline
# Script used to get the list and number of species per position, plus distance with tree_doctor
# Customized for Uppmax SLURM batch system.
# To run the pipeline : ./submit_MAF_speciesLister.sh

# Requires the following data:
# - MAF alignment files, compressed in gz or zst.
# - An output directory.
# - reference species name
# - a tree in Newick format.

# The programms used were already installed on the Uppmax system and were called as modules. 
# Requires following packages to be installed: 
# - biopython v1.73 : https://biopython.org/wiki/Download
# - Phast v1.5 : https://github.com/CshlSiepelLab/phast
# - gzip 1.5 : https://www.gnu.org/software/gzip

# Bash info:
# GNU bash, version 4.2.46(2)-release (x86_64-redhat-linux-gnu)
# Copyright (C) 2011 Free Software Foundation, Inc.
# License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>

# SLURM version:
# slurm 19.05.8

echo -n "Time started: "
date

set -ex

TmpID=${SLURM_JOB_ID:-$$}
datetime=$(date +%d_%b_%H_%M) 

# Create a text listing all files in there:
input_directory="files/input/MAF" # input directory path
list_file="listMAF_"$datetime".txt" # list of maf.gz files to process
ls -d $input_directory > $list_file # create a text file

# species
species_name="Homo_sapiens"

# tree
tree="input/TREE/Zoonomia_ChrX_lessGC40_241species_30Consensus.nh"

# output directory. Create if non-existant.
output_dir="files/output/"$species_name"_speciesList"
if [ ! -d "$output_dir" ]; then
 	mkdir -p "$output_dir" || exit
fi

sbatch --array=1-$number_files run_scripts/run_MAF_speciesLister.sh $list_file $output_dir $species_name $tree

echo -n "Time ended: "
date
