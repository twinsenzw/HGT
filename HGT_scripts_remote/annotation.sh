#!/bin/bash

cd /home/zhouw/HGT/"$1"

blasting
for fa in *.fa; do
	/home/zhouw/blast+/bin/blastx -query "$fa" -db /home/zhouw/blast+/bin/uniref90 -out "$fa".out -outfmt "6 qseqid qlen sseqid slen length bitscore pident" -num_alignments 1
done

#mapping

cat *.out | sort | uniq > nonred.out
cp output.txt output_bk.txt
while read line; do

	query=$(echo "$line" | awk -F '\t' '{print $1}')
	subject=$(echo "$line" | awk -F '\t' '{print $3}')
	bit=$(echo "$line" | awk -F '\t' '{print $6}')

	go=$(cat "/home/zhouw/HGT/map_go_uniref90.txt" | grep "$subject" | awk '{print $1}' | tr '\n' ',' | sed 's/\(.*\),/\1/')
	WEGO=$(cat "/home/zhouw/HGT/map_go_uniref90.txt" | grep "$subject" | awk '{print $1}' | tr '\n' '\t')

	#======updating output file======
	sed -i "s/$query/$query|$go/" output.txt

	#======output for WEGO========
	echo -e "$query\t$WEGO" >> WEGO.out

done < nonred.out

