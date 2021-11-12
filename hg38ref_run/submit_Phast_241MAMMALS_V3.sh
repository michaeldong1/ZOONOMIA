#!/bin/bash

#SBATCH -A <name_project>
#SBATCH -e error_out/phast_%j.err
#SBATCH -o error_out/phast_%j.out
#SBATCH -J PHAST
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 2:00:00

# Phast_ZOONOMIA version 3.0

### Description of the pipeline
# Runs both PhastCons and PhyloP on the Broad 250-way human-referenced MAF alignment file, primarily splitted into 10 Mb length chunks
# It is currently set to run the default PhastCons and PhyloP command (with same parameters from UCSC PhastCons and PhyloP 100-way alignment)
# If needed don't hesitate to custom the parameters.
# This script is a batch script ready to be submitted to SLURM.
# To submit job: ./submit_Phast_241MAMMALS_V3.sh

# Requires the following data:
# - substitution model files for autosome, chromosome X and chromosome Y.
# if you only have one model, leave both chrX and Y parameters at 0.
# - An output directory.
# - a folder with all the MAFs (can be gzipped or zst'd). 
# !!! The MAF files should have this name structure:
# <chromosome>.<start>.<length>.maf[.<zst/gz/none> if compressed]

# The programms used were already installed on the Uppmax system and were called as modules. 
# Requires following packages to be installed: 
# - BEDOPS v2.4.38 : https://bedops.readthedocs.io/en/latest
# - Phast v1.5 : https://github.com/CshlSiepelLab/phast
# - Python 2.7.15 : https://www.python.org/downloads/release/python-2715

# Bash info:
# GNU bash, version 4.2.46(2)-release (x86_64-redhat-linux-gnu)
# Copyright (C) 2011 Free Software Foundation, Inc.
# License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>

# SLURM version:
# slurm 19.05.8

# Print time of start in the programm standard output (the .out file)
echo -n "Time started: "
date
TmpID=${SLURM_JOB_ID:-$$}
datetime=$(date +%d_%b_%H_%M) 

# Reference species
species_name="Homo_sapiens"

# Get list of files.
input_MAF_directory="files/input/MAF"
list_file="files_"$datetime".txt" # list of maf.gz files to process
ls -d $input_MAF_directory > $list_file # create a text file

# Get number of files to process.
number_files=$(wc -l $list_file | awk '{print $1}')

# Output directories
output_dir="files/output"
out_Phast=$output_dir"/Phast_"$datetime
out_PhastCons=$out_Phast"/PHASTCONS"
out_PhyloP=$out_Phast"/PHYLOP"
if [ ! -d $out_PhastCons ]; then
 	mkdir -p $out_PhastCons
fi
if [ ! -d $out_PhyloP ]; then
 	mkdir -p $out_PhyloP
fi

# Get the chromosome number from the file name
model_General="input/MOD/Anc239_allARs_100kb_lessGC40_241species_30Consensus.mod"
if [ ! -f $model_General ]; then
 	echo "No main model file";
	exit
fi
Xmodel_yes=1 # leave to 1 if you want to use chromosome X specific model.
model_X="input/MOD/human_100kb_chrX_lessGC40_241species_30Consensus.mod"
if [ ! -f $model_X && $Xmodel_yes == 1 ]; then
 	echo "No X model despite specifying 1 to Xmodel_yes";
	exit
fi
Ymodel_yes=1 # leave to 1 if you want to use chromosome Y specific model.
model_Y="input/MOD/human_100kb_chrY_lessGC40_241species_30Consensus.mod"
if [ ! -f $model_Y && $Ymodel_yes == 1 ]; then
 	echo "No Y model despite specifying 1 to Ymodel_yes";
	exit
fi

# Submit job:
sbatch --array=1-$number_files run_scripts/run_Phast_241MAMMALS_V3.sh $list_file $out_PhastCons $out_PhyloP $model_General $Xmodel_yes $Ymodel_yes $model_X $model_Y

echo -n "Time ended: "
date
