#!/bin/bash

#==========prereq===========
module load khmer/2.0
source /fastscratch/zhouw/MinDB/scripts/env_var.sh
cd "$Gen_recon_dir"; mkdir "$1"

#==========estract genomes===========
cd "$Pan_out_dir"/sam_pan
Patho_out=$(echo "$1.sam-sam-report.tsv")
cat "$Patho_out" | awk -F '\t' '{print $1}' | tail -n +3 > "$Gen_recon_dir"/"$1"/chr_list.txt

#==========estract reads of the most abundant 20 bacteria===========
cd "$Gen_recon_dir"/$1
cnt=0
rm bac_chr_list.txt
while read line; do
	if [ -s "$DB_dir/bacteria_$line" ]; then
		cat "$Pan_out_dir"/sam_pan/updated_"$1.sam" | grep -v ^@ | awk -v chr="$line" '{ if ($3==chr) print "@"$1"\n"$10"\n+\n"$11}' >> "read_$line.fq"
		((cnt++))
	fi
	if [ $cnt -gt 50 ]; then
		break
	fi
done < chr_list.txt

#==========clean up===========
cd "$Pan_out_dir"/sam_pan
rm updated_"$1.sam"


