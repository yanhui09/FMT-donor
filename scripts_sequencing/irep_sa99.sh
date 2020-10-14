#!/bin/bash
# enter the working directory
#cd $WORKING_DIR
# dereplicate the refined MAGs at 99% similarity using drep
conda activate derep
parallel -j 1 'dRep dereplicate -p 8 -comp 50 -con 5 {}_drep_sa99 -g {}/*.fa -nc 0.25 -sa 0.99' ::: mags
# change name
mv mags_drep_sa99/dereplicated_genomes mags_sa99
parallel -j 1 --plus 'mv {} {/sinASSsinBIN/}' ::: mags_sa99/*.fa
parallel -j 1 --plus 'mv {} {/_bin./bin_}' ::: mags_sa99/*.fa

# alignment
parallel -j 1 --plus './irep_align.sh mags_sa99 {} {/_1.fastq/_2.fastq} _1.fastq 8 irep_sa99' ::: ~/sharedFatboy/FMT_metagenome/FMT_temp/FMT/CLEAN_READS/*_1.fastq 
# run irep on the sam sorted by read
# memory block, run separately
#iRep -f reassembled_bins_drep_cov10_sa99/*.fa -s irep_sa99/allGenomes-vs-*.sam -t 8 -o irep_sa99/FMT
mkdir -p irep_sa99/irep_out
parallel -j 1 'iRep -f mags_sa99/*.fa -s {} -t 8 -o irep_sa99/irep_out/{/.}' ::: irep_sa99/allGenomes-vs-*.sam
conda deactivate

exit
