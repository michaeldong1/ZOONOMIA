# ZOONOMIA project scripts 

## Introduction

This repository regroups the pipelines used to generate the different datasets used for the generation and analysis of the 241-mammal Zoonomia dataset, for the upcoming publications: 
- "Evolutionary constraint and innovation across hundreds of placental mammals" UNDER REVIEW
- and "Evolutionary constraint and innovation across hundreds of placental mammals" UNDER REVIEW

They are set to run on the Uppmax SLURM batch system.

The programs used were already installed on the Uppmax system.

For the list of programs to install or inputs to provide, check the submit scripts for the README.

## hg38_run

Includes the pipelines used for processing the data for the zoonomia paper on human.

Please check the parameters in the submits before running. 

### hal2maf_241MAMMALS

script : submit_hal2maf_241MAMMALS_V2.sh , version : v2.0

Description : Convert the HAL unreferenced alignment to MAF-formatted species-referenced alignments, given a reference species and a BED file specifying which fractions to extract. This script was used to extract the human-referenced alignemnt for PhyloP and PhastCons conservation calculations. 

Depending on the parameters, you can specify if you want the duplicates filetered out or not.

Note: The MAF alignment file available at https://cglgenomics.ucsc.edu/data/cactus is the same alignment than to the one extracted via this pipeline. In our project we segmented the MAF into 10Mb alignment fractions.

### speciesLister

script : submit_MAF_speciesLister.sh

version : v1.0

Description : Extracts the list, number of species and genomic distance for each position of an alignment.

### modelGeneration_mammals

script : submit_modelGeneration_mammals.sh

version : v1.0

Description : Generate an ancestral sequence repeat-based model based on the alignment.

### Phast_MAMMALS

script : submit_Phast_241MAMMALS_V3.sh

version : v3.0

Description : Calculates PhastCons and PhyloP scores with standard parameters for each MAF file in a folder.

