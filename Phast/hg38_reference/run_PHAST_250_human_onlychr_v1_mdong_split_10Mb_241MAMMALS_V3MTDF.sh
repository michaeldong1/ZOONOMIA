#!/bin/bash

#SBATCH -A <name_project>
#SBATCH -e out/pH241MTDF_%j_%a.err
#SBATCH -J P_H241MTDF
#SBATCH -o out/pH241MTDF_%j_%a.out
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 12:00:00
#SBATCH --array=1-322

# This script is a batch script ready to be submitted to SLURM.
# Runs both PhastCons and PhyloP on the Broad 250-way human-referenced MAF alignment file, primarily splitted into 10 Mb length chunks
# It runs a series of jobs in the indicated directory.
# It is currently set to run the standard PhastCons command (with same parameters from UCSC PhastCons 100-way alignment)
# If needed don't hesitate to custom the parameters.
# Before running, don't forget to:
# - check/change the SBATCH project
# - check/change the SBATCH output and error message pathways
# - check/change the number of cores needed for each job (no need for too high if the alignment files were primarily splitted using maf_parse)
# - check/change the array range to match the number of files.
# - check/change the indicated input ($FILES), output directory, output name and options below. 
# Be conscious that this array job is quite heavy, depending on how many and how big your files are.

# Print time of start in the programm standard output (the .out file)
echo -n "Time started: "
date

# Call for modules already installed on Uppmax
module load bioinfo-tools
module load phast/1.5
module load python/2.7.15
module load BEDOPS

# Get name of file nb $SLURM_ARRAY_TASK_ID
list_file="/proj/uppstore2017228/KLT.04.200M/200m_MD/scripts/250_MAMMALS_scripts/LISTS/MAF/human_onlychr_v1_mdong_10Mb_241MAMMALS_V2_MTDF.txt"
FILE=$(head -$SLURM_ARRAY_TASK_ID $list_file | tail -1)

# Get only the base name (names of the files are in this format: <chromosome>.<start>.<length>.maf.<format: gz or zst>)
fnamebase=${FILE##*/}
fname=${fnamebase%.*}
fname2=${fname%.*}
mafFILE=${FILE%.*}

# uncompress target file, keep original. Adapt if the file is compressed in another format than gz
/proj/uppstore2017228/KLT.04.200M/200m_MD/resources/zstd/zstd-dev/zstd -d -k $FILE
if [[ "$fnamebase" =~ ^*.gz$ ]]; then
        gzip -d -c $FILE > $mafFILE
elif [[ "$fnamebase" =~ ^*.zst$ ]]; then
        zstd -d -k $FILE
fi

# Get the chromosome number from the file name
chromosome="$(cut -d'.' -f1 <<< $fnamebase)"
chromosome_nb=$(echo $chromosome | sed -r 's/chr//g')

# Model file pathway
model_conserved="/proj/uppstore2017228/KLT.04.200M/200m_MD/data/new_250_MAMMALS_v2_20201120/MOD/V3_100kb/Anc239_allARs_100kb_lessGC40_241species_30Consensus.mod"
if [[ "$chromosome_nb" =~ ^X.* ]]; then
	model_conserved="/proj/uppstore2017228/KLT.04.200M/200m_MD/data/new_250_MAMMALS_v2_20201120/MOD/V3_100kb/human_100kb_chrX_lessGC40_241species_30Consensus.mod"
elif [[ "$chromosome_nb" =~ ^Y.* ]]; then
	model_conserved="/proj/uppstore2017228/KLT.04.200M/200m_MD/data/new_250_MAMMALS_v2_20201120/MOD/V3_100kb/human_100kb_chrY_lessGC40_241species_30Consensus.mod"
fi

##### PhyloP

# Output directory. If it does not exist, creates it.
output_dir="/proj/uppstore2017228/KLT.04.200M/200m_MD/data/new_250_MAMMALS_v2_20201120/OUT/PHAST/human_onlychr_v3_mdong_10Mb_241MAMMALS_MTDF/PHYLOP"
if [ ! -d "$output_dir" ]; then
 	mkdir "$output_dir"
fi

# Output name format for list of conservation scores on each position
output_name=$fname2"_scoresPhyloP_250.wig"

# Set PhyloP options
options="-i MAF --method LRT --mode CONACC --wig-scores --chrom $chromosome"

# Run PhyloP
phyloP $options $model_conserved $mafFILE > $output_dir/$output_name
# Generate a BED version	
wig2bed < $output_dir/$output_name >> $output_dir/$output_name".bed"

##### PhastCons

# Output directory. If it does not exist, creates it.
output_dir="/proj/uppstore2017228/KLT.04.200M/200m_MD/data/new_250_MAMMALS_v2_20201120/OUT/PHAST/human_onlychr_v3_mdong_10Mb_241MAMMALS_MTDF/PHASTCONS"
if [ ! -d "$output_dir" ]; then
	mkdir "$output_dir"
fi

# Output name format for list of conservation scores on each position
output_name=$fname2"_scoresPhastCons_250.wig"

# Output name format for list of most conserved regions
output_name2="scoresPhastCons_"$fname2"_viterbi_mostconserved.bed"

# Set PhastCons options
options="--msa-format MAF --target-coverage 0.3 --expected-length 45 --rho 0.3 --viterbi $output_dir/$output_name2 --score --seqname $chromosome --idpref $chromosome"

# Run the command
phastCons $options $mafFILE $model_conserved > $output_dir/$output_name
# Generate a BED version
wig2bed < $output_dir/$output_name >> $output_dir/$output_name"_scores.bed"

##### END

# Delete uncompressed file
rm $mafFILE

echo -n "Time ended: "
date
