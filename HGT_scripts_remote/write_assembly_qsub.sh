#!/bin/bash

#==========prereq===========
source /fastscratch/zhouw/MinDB/scripts/env_var.sh
cd "$Gen_recon_dir"; mkdir qsub2

#==========write qsub jobs===========
cd qsub2
while read line; do
	echo "#PBS -N recon_$line" > $line.qsub
	echo "#PBS -l walltime=150:00:00,mem=128gb,nodes=1:ppn=1" >> $line.qsub
	echo "module load khmer/2.0" >> $line.qsub
	echo "cd $Gen_recon_dir" >> $line.qsub
	echo "$HGT_scripts_remote/patho_to_reads.sh $line" >> $line.qsub
	echo "cd $line" >> $line.qsub
	echo "for fq in *.fq; do" >> $line.qsub
	echo "$megahit_dir"'/megahit -o ./megahit_$fq --out-prefix megahit -r $fq' >> $line.qsub
	echo 'cd megahit_$fq' >> $line.qsub
	echo "fq1=$(ls /data/ohlab/mWGS_rawdata/unzipped/$line*mate1)" >> $line.qsub
	echo "fq2=$(ls /data/ohlab/mWGS_rawdata/unzipped/$line*mate2)" >> $line.qsub
	echo "fq3=$(ls /data/ohlab/mWGS_rawdata/unzipped/$line*singletons)" >> $line.qsub
	echo 'ln -s $fq1 mate1.fq' >> $line.qsub
	echo 'ln -s $fq2 mate2.fq' >> $line.qsub
	echo 'ln -s $fq3 single.fq' >> $line.qsub
	echo "$Price_dir/PriceTI -fpp mate1.fq mate2.fq 500 95 -icf megahit.contigs.fa 1 1 5 -nc 10 -dbmax 72 -mol 30 -tol 20 -mpi 90 -target 90 2 1 1 -o Price.contigs.fa" >> $line.qsub
	echo "cd .." >> $line.qsub
	echo "done" >> $line.qsub
done < $HGT_scripts_remote/Met_list_2.txt
