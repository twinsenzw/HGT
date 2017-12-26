#!/bin/bash

# run all species combinations for a sample, and clean up the folder after the run.

#==========prereq===========
source /fastscratch/zhouw/MinDB/scripts/env_var.sh
cd "$HGT_dir_remote"/"$1"

#==========run HGT identification===========
while read line1; do
	while read line2; do
		if [ "$line1" != "$line2" ]; then
			"$HGT_scripts_remote"/HGT_species_pair.sh "$line1" "$line2" "$1"
		fi
	done < present_species2.txt 
	sed -i "/\<$line1\>/d" present_species1.txt
	sed -i "/\<$line1\>/d" present_species2.txt 
done < present_species1.txt

#==========Clean up===========
rm *.faa
rm *.fna
rm full_table*
rm *.dS
rm *.uc
rm *.cutoff

