#!/bin/bash

#================prereq=================
source /fastscratch/zhouw/MinDB/scripts/env_var.sh
module load khmer/2.0


#================format fasta==================
cd "$HGT_dir_remote"/"$3"

#---------cleanup fasta header and truncate wrap----------
awk '{print $1}' "$1_dna.fna" | sed -e 's/\(^>.*$\)/#\1#/' | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' > "$1_dna_clean.fna"
awk '{print $1}' "$2_dna.fna" | sed -e 's/\(^>.*$\)/#\1#/' | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d'> "$2_dna_clean.fna"
awk '{print $1}' "$1_protein.faa" | sed -e 's/\(^>.*$\)/#\1#/' | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' | sed 's/\*//g' > "$1_protein_clean.faa"
awk '{print $1}' "$2_protein.faa" | sed -e 's/\(^>.*$\)/#\1#/' | tr -d "\r" | tr -d "\n" | sed -e 's/$/#/' | tr "#" "\n" | sed -e '/^$/d' | sed 's/\*//g' > "$2_protein_clean.faa"

#---------encode species information in the header----------
sed -i "/^>/ s/$/%$1/g" "$1"_dna_clean.fna
sed -i "/^>/ s/$/%$2/g" "$2"_dna_clean.fna
while read line; do
	header1=$(cat "full_table_$1" | grep -F "$line" | awk -F '\t' '{print $3}' | head -n1)
	header2=$(cat "full_table_$2" | grep -F "$line" | awk -F '\t' '{print $3}' | head -n1)
	if [[ ! -z $header1 && $header2 ]]; then
		cat "$1_protein_clean.faa" | grep -A1 -w "$header1" > "$1_$2_$line.protein.fasta"
		cat "$2_protein_clean.faa" | grep -A1 -w "$header2" >> "$1_$2_$line.protein.fasta"	
		cat "$1_dna_clean.fna" | grep -A1 -w "$header1" > "$1_$2_$line.dna.fasta"
		cat "$2_dna_clean.fna" | grep -A1 -w "$header2" >> "$1_$2_$line.dna.fasta"
	fi
done < COG_list.txt

#================fit null distribution==================

#---------check if enough loci to fit null distribution---------
pnum=$(ls "$1_$2_"*.dna.fasta | wc -l)
if [ $pnum -le 5 ]; then
	echo "Not enough information"
	exit
fi

#---------codon based alignment---------
"$R_bin"/Rscript "$HGT_scripts_remote"/cdn_align_null.r "$clustal"

#---------dS calculation---------
rm "$1_$2.dS"
cd "$paml_bin"
while read line; do
	echo "      seqfile = $HGT_dir_remote/$3/$1_$2_$line.dna.fasta.aln.fa * sequence data file name" > yn00.ctl
	echo "      outfile = yn           * main result file" >> yn00.ctl
	echo "      verbose = 0  * 1: detailed output (list sequences), 0: concise output" >> yn00.ctl
	echo "" >> yn00.ctl
        echo "        icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below" >> yn00.ctl
	echo "" >> yn00.ctl
	echo "    weighting = 0  * weighting pathways between codons (0/1)?" >> yn00.ctl
	echo "   commonf3x4 = 0  * use one set of codon freqs for all pairs (0/1)?" >> yn00.ctl
	echo "*       ndata = 1" >> yn00.ctl
	./yn00
	cat ./2YN.dS | awk '{print $2}' | grep -v '^$' >> "$HGT_dir_remote/$3/$1_$2.dS"
done < "$HGT_dir_remote"/"$3"/COG_list.txt

#---------dS cut-off calculation---------
cd "$HGT_dir_remote"/"$3"
"$R_bin"/Rscript "$HGT_scripts_remote"/cut_off.r "$1_$2.dS" "$1_$2.cutoff"
cutoff=$(cat "$1_$2.cutoff")


#================test for dS difference==================
rm *.dna.fasta; rm *.protein.fasta; rm *.aln.fa

#---------append species info to fasta header and cat---------
cat "$1_protein_clean.faa" | sed "/^>/ s/$/#$1/g" > "$1_$2.cat.faa"
cat "$2_protein_clean.faa" | sed "/^>/ s/$/#$2/g" >> "$1_$2.cat.faa"

#---------cluster, and generate all pairs of 1 vs 2 sequences from the same cluster---------
"$usearch_dir"/usearch -cluster_fast "$1_$2.cat.faa" -id 0.5 -uc "$1_$2.uc"
cat "$1_$2.uc" | awk -F '\t' '{if ( $1 == "C" ) print $0}' | awk -F '\t' '{if ( $3 >= 2 ) print $2}' > cluster.list
rm member.list
while read line; do
	cat "$1_$2.uc" | awk -v l=$line -F '\t' '{if ( $1 != "C" && $2 == l ) print $2,"\t",$9}' >> member.list
done < cluster.list
mkdir "$1_$2"
while read line; do # only keep members if both species are represented in cluster
	cat member.list | awk -v l=$line -F '\t' '{if ( $1 == l ) print $2}' | awk -v s=$1 -F '#' '{if ( $2 == s ) print $0}' > "$1_$2/$line.$1.list"
	cat member.list | awk -v l=$line -F '\t' '{if ( $1 == l ) print $2}' | awk -v s=$2 -F '#' '{if ( $2 == s ) print $0}' > "$1_$2/$line.$2.list"
done < cluster.list
cd "$1_$2"
find -empty | awk -F '/' '{print $2}' | awk -F '.' '{print $1}' | sort | uniq > empty.list
while read line; do
	rm $line.*.list
done < empty.list
rm empty.list

#---------write all pairwise sequences---------
ls *.list | awk -F '.' '{print $1}' | sort | uniq > fcluster.list
while read line; do 
cnt=0
	while read l1; do
		while read l2; do
			cl1=$(echo $l1 | cut -c 1- | awk -F '#' '{print $1}') 
			cl2=$(echo $l2 | cut -c 1- | awk -F '#' '{print $1}')
			cat ../"$1_protein_clean.faa" | grep -A1 -w "$cl1" > "$line.$cnt.faa"
			cat ../"$2_protein_clean.faa" | grep -A1 -w "$cl2" >> "$line.$cnt.faa"
			cat ../"$1_dna_clean.fna" | grep -A1 -w "$cl1" > "$line.$cnt.fna"
			cat ../"$2_dna_clean.fna" | grep -A1 -w "$cl2" >> "$line.$cnt.fna"
			((cnt++))
		done < "$line.$2.list"
	done < "$line.$1.list"
done < fcluster.list
rm *.list

#---------calculate synonymous distance for each pair---------
#cp "$HGT_scripts_remote"/cdn_align_test.r cdn_align_test.r
"$R_bin"/Rscript "$HGT_scripts_remote"/cdn_align_test.r "$clustal"
ls *.fna.aln.fa > alignment.list
rm "$HGT_dir_remote/$3/$1_$2/test_dS.list"
cd "$paml_bin"
while read line; do
	echo "      seqfile = $HGT_dir_remote/$3/$1_$2/$line * sequence data file name" > yn00.ctl
	echo "      outfile = yn           * main result file" >> yn00.ctl
	echo "      verbose = 0  * 1: detailed output (list sequences), 0: concise output" >> yn00.ctl
	echo "" >> yn00.ctl
        echo "        icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below" >> yn00.ctl
	echo "" >> yn00.ctl
	echo "    weighting = 0  * weighting pathways between codons (0/1)?" >> yn00.ctl
	echo "   commonf3x4 = 0  * use one set of codon freqs for all pairs (0/1)?" >> yn00.ctl
	echo "*       ndata = 1" >> yn00.ctl
	./yn00
	dS=$(cat ./2YN.dS)
	echo $dS | cut -c 3- >> "$HGT_dir_remote/$3/$1_$2/test_dS.list"
done < "$HGT_dir_remote/$3/$1_$2/alignment.list"

#---------output result and sequences---------
cd "$HGT_dir_remote"/"$3"; mkdir Final_output; cd Final_output
echo "$1~$2~$cutoff" >> "$HGT_dir_remote"/"$3"/Final_output/output.txt
while read line; do
	distance=$(echo "$line" | awk '{print $3}')
	if [ $(bc <<< "$distance < $cutoff") -eq 1 ]; then
		echo "$line" >> output.txt
		# the gene names
		gene1=$(echo "$line" | awk '{print $1}')
		gene2=$(echo "$line" | awk '{print $2}')
		# which species to find the gene?
		g1s=$(echo "$gene1" | awk -F '%' '{print $2}')
		g2s=$(echo "$gene2" | awk -F '%' '{print $2}')
		cat ../"$g1s"_dna_clean.fna | grep -A 1 "^>$gene1" >> "$g1s.$g2s.fa"
		cat ../"$g2s"_dna_clean.fna | grep -A 1 "^>$gene2" >> "$g1s.$g2s.fa"
	fi
done < "$HGT_dir_remote/$3/$1_$2/test_dS.list"
rm -r "$HGT_dir_remote"/"$3"/"$1_$2"
