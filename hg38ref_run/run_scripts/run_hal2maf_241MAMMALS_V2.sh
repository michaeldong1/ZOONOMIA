#!/bin/bash

#SBATCH -A projectName
#SBATCH -J hal2maf
#SBATCH -o error_out/hal2maf/hal2maf_HUMAN10Mb_%j_%a.out
#SBATCH -e error_out/hal2maf/hal2maf_HUMAN10Mb_%j_%a.err
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 1-00:00:00

echo -n "Time started: "
date

set -ex

module load bioinfo-tools
module load hal
module load kent

### !!! run the "submit" script instead. !!!
# Script and parameters for converting hal to maf. Main script. 
# hal2maf_ZOONOMIA version 2
# Used to extract human-referenced MAF, on autosomes, chromosome X and chromosome Y.
# README info in the submit script.

TmpID=${SLURM_JOB_ID:-$$}
datetime=$(date +%d_%b_%H_%M) 

# Get parameters:
species_name=$1
hal_file_input=$2
output_dir=$3
coordinates_file=$4 # coordinates file of each segment.
remove_dups=$5

# Get coordinates
chr=$(head -$SLURM_ARRAY_TASK_ID $coordinates_file | tail -1 | awk -F '\t' '{print $1}')
start=$(head -$SLURM_ARRAY_TASK_ID $coordinates_file | tail -1 | awk -F '\t' '{print $2}')
length=$(head -$SLURM_ARRAY_TASK_ID $coordinates_file | tail -1 | awk -F '\t' '{print $3-$2}')

# Parameters for hal2maf
maf_file_output=$output_dir"/$chr.$start.$length.maf"
options="--refGenome $species_name --refSequence $chr --onlyOrthologs --noAncestors --start $start --length $length"

# run conversion
hal2maf $options $hal_file_input $maf_file_output

# if yes, remove the duplicates using mafDuplicateFilter
if [[ $remove_dups == 1 ]]; then
    mafTools/bin/mafDuplicateFilter -m $maf_file_output > $maf_file_output"_tmp"
	mv $maf_file_output"_tmp" $maf_file_output
fi

# Gunzip the file
gzip $maf_file_output

echo -n "Time ended: "
date
