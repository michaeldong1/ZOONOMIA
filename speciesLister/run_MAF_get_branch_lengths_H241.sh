#!/bin/bash

#SBATCH -A snic2021-5-28
#SBATCH -p core
#SBATCH -n 4
#SBATCH -e out/branchh241_%j_%a.err
#SBATCH -J branchh241_%a
#SBATCH -o out/branchh241_%j_%a.out
#SBATCH -t 1-00:00:00
#SBATCH --array=2-322%50

echo -n "Time started: "
date

set -ex

module load bioinfo-tools
module load biopython/1.73
module load BEDTools
module load BEDOPS
module load phast

# Script for converting hal to maf.
# Will modify it to generate chunks of 1Mb

TmpID=${SLURM_JOB_ID:-$$}
datetime=$(date +%d_%b_%H_%M) 

# Get name of file nb $SLURM_ARRAY_TASK_ID
list_file="/proj/uppstore2017228/KLT.04.200M/200m_MD/scripts/250_MAMMALS_scripts/LISTS/MAF/human_onlychr_v1_mdong_10Mb_241MAMMALS_V2_MTDF.txt"
FILE=$(head -$SLURM_ARRAY_TASK_ID $list_file | tail -1)

# Get only the base name
fnamebase=${FILE##*/}
fname=${fnamebase%.*}
fname2=${fname%.*}
mafFILE=${FILE%.*}

# uncompress target file, keep original
/proj/uppstore2017228/KLT.04.200M/200m_MD/resources/zstd/zstd-dev/zstd -d -k $FILE

output_file="/proj/uppstore2017228/KLT.04.200M/200m_MD/data/new_250_MAMMALS_v2_20201120/OUT/branch_length/human_onlychr_v2_mdong_10Mb_241MAMMALS_MTDF/"$fname2".bed"
reference_species="Homo_sapiens"
tree="/proj/uppstore2017228/KLT.04.200M/200m_MD/data/new_250_MAMMALS_v2_20201120/TREE/Zoonomia_ChrX_lessGC40_241species_30Consensus.nh"

# Run_script to get alleles for main species, second species and whole alignment + count for whole alignment:
python MAF_allele_stats_get_names_v1.py $mafFILE $output_file $reference_species $tree
gzip $output_file

# Delete decompressed file and temp files
rm $mafFILE

echo -n "Time ended: "
date
