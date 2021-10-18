#!/bin/bash


## Convert svg to png in the same directory 
# The first argument given to this script is the directory in which to find the svg file(s)

dir=$1

# Iterate through all svg files
files=$(find $dir -type f -name "*.svg")

for file in $files
do
	file_wo_svg=${file::-4}
	inkscape $file --export-png=${file_wo_svg}.png
done	
