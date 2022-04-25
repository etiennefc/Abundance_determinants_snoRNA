#!/bin/bash

link_gtex=$1
output=$2

# Download GTEx whole TPM dataset
wget -O gtex_temp_file $link_gtex &&
zcat gtex_temp_file | awk 'NR>2' > gtex_temp_file2


# gtex_var is predefined in the shell command (see rule download_gtex_host_data)
# It corresponds to all the columns we want to extract in the whole GTEx TPM dataset
# We select only these columns and concat them (paste) horizontally in the final df
for col in $gtex_var
do
	echo $col
	awk -v FS="\t" -v COL=$col 'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}} {print $(f[COL])}' gtex_temp_file2 > ${col}_GTEX_TEMPO_files
done

paste *_GTEX_TEMPO_files > $output &&
rm gtex_temp_file* && rm *_GTEX_TEMPO_files
