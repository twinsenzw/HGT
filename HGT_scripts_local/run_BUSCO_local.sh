#!/bin/bash

#==========prereq===========
source /home/twinsenzw/scratch@jax/MinDB/scripts/env_var.sh
cd "$HGT_dir"; mkdir "$1"; cp "$HGT_dir"/COG_list.txt "$1"
cd "$BUSCO_dir"; mkdir HGT; cd HGT; rm -r run*; rm * # all genome files stored in the directories have to have extension x.fa; x being the species name/ID

#==========copy assembly to local===========
cd "$Gen_recon_dir_local"/"$1" 
for d in megahit_read*/; do
	taxid=$(echo $d | awk -F '_' '{print $3}' | awk -F '.' '{print $1}')
	Price_contig=$(ls -t "$d"/Price.contigs.*.fa | head -n 1)
	cp "$Price_contig" "$BUSCO_dir"/HGT/"$taxid".fa
done

#==========run BUSCO locally, then copy output to remote===========
cd "$BUSCO_dir"/HGT
for fa in *.fa; do
	base=$(echo "$fa" | awk -F '.' '{print $1}')
	"$prodigal_dir"/prodigal.linux -i "$fa" -d "$base"_dna.fna -a "$base"_protein.faa
	# skip if not enough information to train prodigal
	if [ ! -s "$base"_dna.fna ]; then rm "$fa"; rm "$base"_dna.fna; rm "$base"_protein.faa; continue; fi
	sed -i 's/*//g' "$base"_protein.faa
	python "$BUSCO_dir"/BUSCO_v1.22.py -in "$base"_protein.faa -o "$base" -l "$BUSCO_dir"/bacteria -m OGS
	incompness=$(cat ./run_"$base"/full_table_"$base" | grep Missing | wc -l)
	if [ $incompness -le 30 ]; then # only consider genomes that are >20% completeness
		mv "$base"_dna.fna "$HGT_dir"/"$1"/"$base"_dna.fna
		mv "$base"_protein.faa "$HGT_dir"/"$1"/"$base"_protein.faa
		mv ./run_"$base"/full_table_"$base" "$HGT_dir"/"$1"/full_table_"$base"
		#list present species that passes filter
		echo "$base" >> "$HGT_dir"/"$1"/present_species1.txt
		echo "$base" >> "$HGT_dir"/"$1"/present_species2.txt
	fi
done
