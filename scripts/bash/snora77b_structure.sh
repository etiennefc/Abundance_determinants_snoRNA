#!/bin/bash


## Extract SNORA77B (ENSG00000264346) sequence and the sequence of its flanking intronic sequences
# The 1st argument given to this script is the sno_sequence input fasta file
# The 2nd argument given to this script is the flanking regions input file
# The 3rd argument given to this script is the output folder (for mfe and dot bracket)
# The structure (.ps file) will be redirected with all the other snoRNA structures

all_sno_seq=$1
all_flanking_seq=$2
output=$3
temp_output="temp_output.txt"

# Extract the sequence of SNORA77B (ENSG00000264346)
sno_seq=$(grep -A 1 "ENSG00000264346" $all_sno_seq | tail -n1)

# Extract SNORA77B flanking regions line and reverse the order of the characters
# in that line so that it is in the right order (not the order specific to RNAcofold)
flanking=$(grep -A 1 "ENSG00000264346" $all_flanking_seq | tail -n1)
rev_flanking=$(echo $flanking | rev)
left=$(echo ${rev_flanking: -15} | head -c 9)  # Select only intronic nucleotides that paired in the terminal stem created by terminal_stem.smk
right=$(echo ${rev_flanking:0:18} | tail -c +4) # Select only intronic nucleotides that paired in the terminal stem created by terminal_stem.smk

# Create a fasta of the snoRNA sequence plus its flanking regions
echo ">ENSG00000264346_terminal_stem" > temp_input
echo "${left}${sno_seq}${right}" >> temp_input

# Fold the structure and output its MFE in $output and move the
# structure figure (.ps file) to the right directory
RNAfold --infile=temp_input --outfile=$temp_output &&
mv $temp_output $output &&
mv *.ps data/structure/stability/ &&

rm temp_input
