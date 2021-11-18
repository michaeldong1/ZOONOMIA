#!/bin/bash

#SBATCH -A snic2021-5-28
#SBATCH -p core
#SBATCH -n 4
#SBATCH -J branch241_%a
#SBATCH -e error_out/speciesLister/branch241_%j_%a.err
#SBATCH -o error_out/speciesLister/branch241_%j_%a.out
#SBATCH -t 12:00:00

echo -n "Time started: "
date

set -ex

module load bioinfo-tools
module load biopython/1.73
module load phast

# Script for converting hal to maf.
# Will modify it to generate chunks of 1Mb

TmpID=${SLURM_JOB_ID:-$$}
datetime=$(date +%d_%b_%H_%M) 

# Get name of file nb $SLURM_ARRAY_TASK_ID
list_file=$1
FILE=$(head -$SLURM_ARRAY_TASK_ID $list_file | tail -1)

# Get only the base name
fnamebase=${FILE##*/}
fname=${fnamebase%.*}
fname2=${fname%.*}
mafFILE=${FILE%.*}

# uncompress target file, keep original
# need to install zstd if you use zstd-compressed files
dontdelete=0
if [[ "$fnamebase" =~ ^*.gz$ ]]; then
        gzip -d -c $FILE > $mafFILE
elif [[ "$fnamebase" =~ ^*.zst$ ]]; then
        zstd -d -k $FILE
elif [[ "$fnamebase" =~ ^*.maf$ ]]; then
		fnamebase=${FILE##*/}
		fname=${fnamebase}
		fname2=${fname%.*}
		mafFILE=$FILE
		dontdelete=1 # important
fi

output_dir=$2
output_file=$output_dir"/"$fname2".bed"
reference_species=$3
tree=$4

# Run script to get list of species + count for whole alignment file:
python MAF_allele_stats_get_names_v1.py $mafFILE $output_file $reference_species $tree
gzip $output_file

# Delete decompressed file and temp files
if [ $dontdelete == 0 ]; then
rm $mafFILE
fi

echo -n "Time ended: "
date
