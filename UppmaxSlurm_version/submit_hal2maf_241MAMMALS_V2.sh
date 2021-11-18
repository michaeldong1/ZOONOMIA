#!/bin/bash

#SBATCH -A projectName
#SBATCH -J hal2maf
#SBATCH -o submit_hal2maf_%j.out
#SBATCH -e submit_hal2maf_%j.err
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1:00:00

set -ex

echo -n "Time started: "
date

# hal2maf_ZOONOMIA version 2.0

### Description of the pipeline
# Script used to extract human-referenced MAF, on autosomes, chromosome X and chromosome Y on given regions.
# Customized for Uppmax SLURM batch system.
# To run the pipeline : ./submit_hal2maf_241MAMMALS_V2.sh

# Requires the following data:
# - a HAL file to extract MAF files from
# - a BED with one region to extract per line. For this project, I segemented the human genome 
#   by 10Mb regions, giving 322 regions in total. See "files/input/human_split_positions_10Mb_onlychr.bed" for an example.
# 	If no file is provided, interrupts script
# - An output directory.

# The programms used were already installed on the Uppmax system and were called as modules. 
# Requires following packages to be installed: 
# - hal v2.1 : https://github.com/ComparativeGenomicsToolkit/hal 
# - kent tools v378 : https://github.com/ucscGenomeBrowser/kent
# - mafTools v01 : https://github.com/dentearl/mafTools
# - gzip 1.5 : https://www.gnu.org/software/gzip

# Bash info:
# GNU bash, version 4.2.46(2)-release (x86_64-redhat-linux-gnu)
# Copyright (C) 2011 Free Software Foundation, Inc.
# License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>

# SLURM version:
# slurm 19.05.8

# Check whether you need the duplicates filtered out or not. remove_dups=1 means dupes will be removed.

# Note : we decided to extract each segement one by one with the hal2maf command instead of running 
# hal2maf2MP.py (part of the hal package), whcich would have been faster, due to storage limitations. 
# If you have enough space and memory, we would recommend you to run hal2maf2MP.py with the --splitBySequence parameter. 
# This will extract each chromosome, scaffold and contig in separate MAFs. The MAFs extracted will be the same.

# Before running, edit parameters to match your files.

TmpID=${SLURM_JOB_ID:-$$}
datetime=$(date +%d_%b_%H_%M) 

# Reference species
species_name="Homo_sapiens"

# HAL file
hal_file_input="files/input/HAL/241-mammalian-2020v2.hal"
if [ ! -f "$coordinates_file" ]; then
 	echo "No HAL file";
	exit
fi

# Coordinates file 
coordinates_file="files/input/human_split_positions_10Mb_onlychr.bed" # coordinates file of each segment.
if [ ! -f "$coordinates_file" ]; then
 	echo "No coordinates file";
	exit
fi

# Get number of files to process.
number_files=$(wc -l $coordinates_file | awk '{print $1}')

# Output directory. Create if it does not exist.
output_dir="MAF_"$datetime
if [ ! -d "$output_dir" ]; then
 	mkdir "$output_dir" || exit
fi

# filter out duplicates
remove_dups=1 # 1 if duplicates need to be removed. 0 if not. Requires mafDuplicateFilter
# Leaving it at 0 might shift the results on positions with a high number of paralogs across same species.
# While mafDuplicateFilter might not keep the best paralog for each alignment block, it seriously limits the
# drawbacks from species overalignment.

# Create error message directory if it does not exist yet.
error_dir="error_out/hal2maf"
if [ ! -d "$error_dir" ]; then
 	mkdir -p "$error_dir" || exit
fi

# run conversion / submit job
sbatch --array=1-$number_files run_scripts/run_hal2maf_241MAMMALS_V2.sh $species_name $hal_file_input $output_dir $coordinates_file $remove_dups

# End-message
echo -n "Time ended: "
date
