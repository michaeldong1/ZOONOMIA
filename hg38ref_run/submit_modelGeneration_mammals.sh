#!/bin/bash

#SBATCH -A projectName
#SBATCH -J RepeatModel_%j
#SBATCH -o error_out/RepeatModel_%j.out
#SBATCH -o error_out/RepeatModel_%j.err
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:10:00

# modelGeneration_mammals v1.0

### Description
# A pipeine developped specifically to build a substitution model based on repeats calculated on the most ancestral sequence of the 241 Mammals Alignment.
# Following recoomentations, we calculated repeats on the 2nd most ancestral sequence before reconverting coordinates to most ancestral consensus.
# Designed to run on the Uppmax	SLURM system. Can be simply ed
# Script customized to run everything in one run.
# if not avaliable as module on your system, the following packages/ apps are required to be installed :
# You will need to change the command lines in accordance.
# to run it: ./submit_modelGeneration_mammals.sh

# Required data
# - A tree in newick format, related to the HAL
# - A HAL alignment file
# - name of species anc ancestors to calcuate against.

# Required softwares: 
# - Repeatmasker : https://www.repeatmasker.org/RepeatMasker/
# - hal : https://github.com/ComparativeGenomicsToolkit/hal
# - BEDOPS v2.4.38 : https://bedops.readthedocs.io/en/latest
# - BEDTools v2.29.2 : https://bedtools.readthedocs.io/en/latest
# - Phast v1.5 : https://github.com/CshlSiepelLab/phast
# - mafTools from Dent Earl : https://github.com/dentearl/mafTools

# Bash info:
# GNU bash, version 4.2.46(2)-release (x86_64-redhat-linux-gnu)
# Copyright (C) 2011 Free Software Foundation, Inc.
# License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>

# SLURM version:
# slurm 19.05.8

echo -n "Time started: "
date

hal_file="HAL/241-mammalian-2020v2.hal" # pathway to HAL file
tree="input/TREE/Zoonomia_ChrX_lessGC40_241species_30Consensus.nh" # Newick format tree
species_rmsk="fullTreeAnc238" # Ancestor to calculate repeats on. Can be different from most ancestral one.
species_final="fullTreeAnc239" # most ancestal species. If identical to "species_rmsk", leave blank.
species_reference="Homo_sapiens"
output_dir="files/output/model_"$species_final 

sbatch run_scripts/run_modelGeneration_mammals.sh $hal_file $tree $species_rmsk $species_final $species_reference $output_dir

echo -n "Time ended: "
date
