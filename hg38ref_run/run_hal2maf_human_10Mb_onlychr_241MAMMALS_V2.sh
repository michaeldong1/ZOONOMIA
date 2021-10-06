#!/bin/bash

#SBATCH -A snic2099-21-9
#SBATCH -J hal2maf
#SBATCH -o error_out/hal2maf/hal2maf_HUMAN10Mb_%j_%a.out
#SBATCH -e error_out/hal2maf/hal2maf_HUMAN10Mb_%j_%a.err
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 3-00:00:00
#SBATCH --array=1-322%50

echo -n "Time started: "
date

set -ex

module load bioinfo-tools
module load hal
module load kent

# Script for converting hal to maf.

TmpID=${SLURM_JOB_ID:-$$}
datetime=$(date +%d_%b_%H_%M) 

# Get parameters:
coordinates_file="files/input/human_split_positions_10Mb_onlychr.bed"
chr=$(head -$SLURM_ARRAY_TASK_ID $coordinates_file | tail -1 | awk -F '\t' '{print $1}')
start=$(head -$SLURM_ARRAY_TASK_ID $coordinates_file | tail -1 | awk -F '\t' '{print $2}')
length=$(head -$SLURM_ARRAY_TASK_ID $coordinates_file | tail -1 | awk -F '\t' '{print $3-$2}')

# Parameters
hal_file_input="files/input/HAL/241-mammalian-2020v2.hal"
maf_file_output="files/input/MAF/$chr.$start.$length.maf"
options="--refGenome Homo_sapiens --refSequence $chr --onlyOrthologs --noAncestors --start $start --length $length"

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
