#!/bin/bash

#SBATCH -A snic2019-8-369
#SBATCH -J hal2maf
#SBATCH -o out/hal2maf_HUMAN10Mb_%j_%a.out
#SBATCH -e out/hal2maf_HUMAN10Mb_%j_%a.err
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 5-00:00:00
#SBATCH --array=83-322%100

echo -n "Time started: "
date

set -ex

module load bioinfo-tools
module load hal

# Script for converting hal to maf.
# Will modify it to generate chunks of 100kb

TmpID=${SLURM_JOB_ID:-$$}
datetime=$(date +%d_%b_%H_%M) 

# Get parameters:
coordinates_file="/proj/uppstore2017228/KLT.04.200M/200m_MD/scripts/250_MAMMALS_scripts/LISTS/HAL_chr_sizes/human_split_positions_10Mb_onlychr.txt"
chr=$(head -$SLURM_ARRAY_TASK_ID $coordinates_file | tail -1 | awk -F '\t' '{print $1}')
start=$(head -$SLURM_ARRAY_TASK_ID $coordinates_file | tail -1 | awk -F '\t' '{print $2}')
length=$(head -$SLURM_ARRAY_TASK_ID $coordinates_file | tail -1 | awk -F '\t' '{print $3}')

# Parameters
hal_file_input="/proj/uppstore2017228/KLT.04.200M/200m_MD/data/new_250_MAMMALS_v2_20201120/HAL/241-mammalian-2020v2.hal"
maf_file_output="/proj/uppstore2017228/KLT.04.200M/200m_MD/data/new_250_MAMMALS_v2_20201120/MAF/HUMAN/241MAMMALS/$chr.$start.$length.maf"
options="--refGenome Homo_sapiens --refSequence $chr --onlyOrthologs --noAncestors --start $start --length $length"

# run conversion
hal2maf $options $hal_file_input $maf_file_output

# Gunzip the file
# gzip $maf_file_output

#compress
/proj/uppstore2017228/KLT.04.200M/200m_MD/resources/zstd/zstd-dev/zstd --rm $maf_file_output

echo -n "Time ended: "
date
