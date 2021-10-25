#!/usr/bin/env python

# Import packages
from Bio import AlignIO
from Bio.AlignIO import MafIO
import os
import sys
import numpy
import subprocess

species = sys.argv[3]
tree_file = sys.argv[4]

# Initiate species list
species_list = ""
species_nb = 0

# Open file to write
file1 = open(sys.argv[2],"w")

# 1 - Create counting table

# For each alignment block in the MAF file:
for multiple_alignment in AlignIO.parse(sys.argv[1], "maf"):

	# For each line (species) in the alignment block:
	for seqrec in multiple_alignment:
		
		# If the species is the reference one, save the coordinates
		if seqrec.id.startswith(species):

			# Reset the species list
			species_list = ""
			species_nb = 0

			# Print chromosome of alignment block in output 
               		chr = seqrec.id.split(".")[1];
			file1.write("%s\t" % chr);

			# Print start coordinate of alignment block in output
			start = seqrec.annotations["start"];
			file1.write("%s\t" % start);
			
			# print end coordinate of alignment block in output
			end = seqrec.annotations["start"] + seqrec.annotations["size"];
			file1.write("%s\t" % end);

		# For each line in the block, extract name of each species
		name = seqrec.id.split(".")[0]
		#file1.write("%s," % name);

		# Save it into the species list
		species_list += name + ",";
		species_nb = species_nb + 1

	# Remove last comma, write into output file
	species_list = species_list[:-1];
	file1.write("%s\t" % species_list);
	file1.write("%s\t" % species_nb);
	
	# Run tree doctor to get total branch length with given list of species, and write the branch length in the output file
	#tree_file = "/proj/uppstore2017228/KLT.04.200M/200m_MD/data/new_250_MAMMALS_v2_20201120/TREE/Zoonomia_ChrX_lessGC40_241species_30Consensus.nh"
	bashCmd = ["tree_doctor","-P",species_list,"-b",tree_file]
	process = subprocess.Popen(bashCmd, stdout=subprocess.PIPE)
	output, error = process.communicate()
	file1.write("%s" % output.split(" ")[1]);
	
	# Skip line for next block
	# file1.write('\n');

# close output file
file1.close()
