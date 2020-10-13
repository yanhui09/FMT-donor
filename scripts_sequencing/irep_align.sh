#!/bin/bash
# $1 directory for .fa fasta files of mags/bins
# $2 forward reads
# $3 reverse reads
# $4 tail-cut/suffix for R1 to get a clean sample using basename; e.g. sample_1.fastq > sample basename sample_1.fastq _1.fastq
# $5 threads
# $6 out
source activate drep
# get the sample name
NAME=$(basename $2 $4)
# concatenate the mags fasta
mkdir -p $6
if [ ! -f "$6/allGenomes.fasta" ]; then
	cat $1/*.fa > $6/allGenomes.fasta
# bt2 index
	mkdir $6/bt2
	bowtie2-build $6/allGenomes.fasta $6/bt2/allGenomes.fasta
fi	
# bt2 alignment
bowtie2 -p 8 -x $6/bt2/allGenomes.fasta -1 $2 -2 $3 --reorder | shrinksam  > $6/allGenomes-vs-$NAME.sam

# run irep on the sam sorted by read
#iRep -f $1/*.fa -s $6/allGenomes-vs-$NAME.sam -t 8 -o $6/$NAME


source deactivate
exit

