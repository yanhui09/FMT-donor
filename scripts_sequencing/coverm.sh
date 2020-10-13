#!/bin/bash
# enter the workding directory
#cd $WORK_DIR
# using coverm to calculate genome-wide coverage and breadth of coverage
mkdir coverm_genomewide
# only analyze 39 (old), 40 (new) two donors
cat mags/*.fa > allGenomes.fasta
# build up contig-bin tab
cd mags
parallel 'cat {} | grep ">" | sed "s/>//;s/$/\t{/}/" > {/.}.tab' ::: *.fa
cd -
cat mags/*.tab > allGenomes_ctb.tab
# manual alignment
mkdir bt2
bowtie2-build ./allGenomes.fasta bt2/allGenomes.fasta
parallel -j 1 --plus 'bowtie2 -p 8 -x bt2/allGenomes.fasta -1 {} -2 {/_1/_2} > allGenomes-vs-{/}.sam' ::: fastqs/*_1.fastq.gz

#
parallel -j 1 --plus 'coverm genome --bam-files {} --genome-definition allgenome.b2c -t 6 --min-covered-fraction 0 -m mean covered_fraction covered_bases length count --output-format sparse > coverm_genomewide/{/.sorted.bam/}.tsv' ::: *.sorted.bam

# make tmp file out into present file; or a larger space to avoid no buffer
parallel --tmpdir $PWD -j 1 --plus 'coverm genome -1 {} -2 {/_1/_2} -d mags -x fa -t 6 --min-covered-fraction 0 -m mean covered_fraction covered_bases length count --output-format sparse > {/.}.tsv' ::: fastq/*_1.fastq.gz

# use coverm to calculate all stuff
# fastq_path for alignment
ls ~/sharedFatboy/FMT_metagenome/FMT_temp/FMT/CLEAN_READS/*_1.fastq | sed 's/_1\.fastq//' > fastq.path
mkdir coverm_genomewide
parallel --tmpdir $PWD -j 1 --plus 'coverm genome -1 {}_1.fastq -2 {}_2.fastq -d mags -x fa -t 6 --min-covered-fraction 0 -m mean covered_fraction covered_bases length count --output-format sparse > coverm_genomewide/{/.}.tsv' :::: fastq.path


