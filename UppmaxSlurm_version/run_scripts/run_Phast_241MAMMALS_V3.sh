#!/bin/bash

#SBATCH -A <name_project>
#SBATCH -e error_out/phast/phast_%j_%a.err
#SBATCH -o error_out/phast/phast_%j_%a.out
#SBATCH -J PHAST_%j_%a
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 12:00:00

### Description of the pipeline
# Runs both PhastCons and PhyloP on the Broad 250-way human-referenced MAF alignment file, primarily splitted into 10 Mb length chunks
# It is currently set to run the default PhastCons and PhyloP command (with same parameters from UCSC PhastCons and PhyloP 100-way alignment)
# If needed don't hesitate to custom the parameters.
# This script is a batch script ready to be submitted to SLURM.
# !!! Run the submit script instead. Readme in "submit_Phast_241MAMMALS_V3.sh" .


# Print time of start in the programm standard output (the .out file)
echo -n "Time started: "
date

# Call for modules already installed on Uppmax
module load bioinfo-tools
module load phast/1.5
module load python/2.7.15
module load BEDOPS

# Get name of file nb $SLURM_ARRAY_TASK_ID
list_file=$1 # list of maf files to process
FILE=$(head -$SLURM_ARRAY_TASK_ID $list_file | tail -1)

# Get only the base name (names of the files are in this format: <chromosome>.<start>.<length>.maf.<format: gz or zst>)
fnamebase=${FILE##*/}
fname=${fnamebase%.*}
fname2=${fname%.*}
mafFILE=${FILE%.*}

# uncompress target file, keep original. Adapt if the file is compressed in another format than gz
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

# Get the chromosome number from the file name
chromosome="$(cut -d'.' -f1 <<< $fnamebase)"
chromosome_nb=$(echo $chromosome | sed -r 's/chr//g')

# Model file pathway
model_conserved=$4

if [[ "$chromosome_nb" =~ ^X.* && $5 == 1 ]]; then
	model_conserved=$7
fi

if [[ "$chromosome_nb" =~ ^Y.* && $6 == 1 ]]; then
	model_conserved=$8
fi

##### PhyloP

# Output directory. If it does not exist, creates it.
output_dir=$3

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
output_dir=$2

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
if [ $dontdelete == 0 ]; then
rm $mafFILE
fi

echo -n "Time ended: "
date
