#!/bin/bash
# enter working dircotry
#cd $WORK_DIR
#######################################################################
#infer gene catalogue
# infer gene for co-assembly and single-assembly contigs, combine all contigs to a single file: final_contig.fasta
# rename each sample output
parallel -j 4 --plus 'mv {} {/_1/} ' ::: ASSEMBLY_INDIVIDUAL/*

# prepend sample names for respective contigs
parallel -j 4 --plus 'sed "s/>/>{/}_/g" {}/final_assembly.fasta > {}/final_assembly_adj.fasta ' ::: ASSEMBLY_INDIVIDUAL/*

# concatenate all the co-assembly and single-assemble contigs  together
cat ASSEMBLY_INDIVIDUAL/*/final_assembly_adj.fasta CO-ASSEMBLY/final_assembly.fasta > final_contig.fasta


# all assembly
gmhmmp -a -d -f G -m ../mgm/MetaGeneMark_v1.mod final_contig.fasta -A protein.fasta -D nucleotide.fasta
# cd-hit to cluster the FMT gene catalogue, discarded
#cd-hit-est -i MIX_nucleotide.fasta -o gene_catalogue.fasta -c 0.95 -aS 0.9 -G 0 -g 1 -d 0 -T 6 -M 0
# use linclust to facilitate clustering speed
mmseqs createdb nucleotide.fasta combined
mkdir tmp
mmseqs linclust combined combined_cluster tmp -c 0.9 --min-seq-id 0.95  --threads 8
mmseqs result2repseq combined combined_cluster combined_clu_rep
mmseqs result2flat combined combined combined_clu_rep combined.linclust.fna --use-fasta-header

# extact ORF > 100 bp with seqkit
seqkit seq -m 100 combined.linclust.fna -w 0 > combined.linclust_f.fna

################################################################
#count the gene abundance matrix
source activate metawrap-env
# bwa index
bwa index combined.linclust.fna
parallel -j 1 --link --plus 'bwa mem -t 8 -M combined.linclust_f.fna {1} {2} | samtools view -Sb - > {1/.}.bam ' ::: /home/yanhui/sharedFatboy/FMT_metagenome/FMT_temp/FMT/CLEAN_READS/17119*_1.fastq ::: /home/yanhui/sharedFatboy/FMT_metagenome/FMT_temp/FMT/CLEAN_READS/17119*_2.fastq
# rm _1 witin file name
parallel -j 1 --plus 'mv {} {/_1/}' ::: *.bam
#  filter pairs mapping to the same gene using the read_count_bam.pl script below and then only takes reads that are mapped with a mapping quality better than 30 (-q30 below). After that we sort the bam-files
parallel -j 8 'samtools view -h {} | read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o {.}.sort.bam -' ::: *.bam
# index sample-bams
parallel -j 8 'samtools index {}' ::: *.sort.bam
# counting using samtools idxstats. It will output three columns: gene-name, gene-length, no. mapped_reads, no. unmapped_reads.
# header
parallel -j 1 'printf "%s\t" {/.} >> header' ::: *[0-9].bam
# gene_names
# select first file to get the gene names
samtools idxstats $(ls *sort.bam | head -n 1) | grep -v "*" | cut -f1 > gene_names
# parallel to get all files' count
parallel -j 8 --plus 'samtools idxstats {} | grep -v "*" | cut -f3 > counts_{/\.sort\.bam/}' ::: *.sort.bam
# create the count matrix
paste gene_names counts_* | cat header - > count_matrix.tab

#############gene annotation
# translate to faa (amino acid sequence) to use ghostkoala for KEGG annotation
seqkit translate combined.linclust_f.fna -j 8 --trim -w 0 > combined.linclust_f.faa
# split to ~300M subfiles for batches submission
seqkit split -j 8 -s 1300000 combined.linclust_f.faa -w 0 -O faa.split


# eggNOG annotation
# download EGGNOG, proteomes, virus.faa, cog annotation and taxonomy
#nohup wget -c http://eggnog5.embl.de/download/eggnog_5.0/e5.proteomes.faa &
#nohup wget -c http://eggnog5.embl.de/download/eggnog_5.0/e5.viruses.faa &
#nohup wget -c http://eggnog5.embl.de/download/eggnog_5.0/e5.taxid_info.tsv &
#nohup wget -c http://eggnog5.embl.de/download/eggnog_5.0/e5.og_annotations.tsv &
# download relationship between faa id and NOG group
#nohup wget -c http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/2/2_members.tsv.gz
#gzip -d 2_members.tsv.gz (bacteria) can't fix the relationship, direct use eggnog_mapper (eggNOG5.0)
# splite the files to chunks with maximum 1M seqs
#mkdir eggnog5.0
#split -l 2000000 -a 3 -d combined.linclust_f.fna eggnog5.0/input_file.chunk_
# generate all the commands that should be distributed in the cluster
# fna sequence using --translate
#for f in eggnog5.0/*.chunk_*; do
#emapper.py -m diamond --translate --no_annot --no_file_comments --cpu 8 -i $f -o $f; 
#done
#cd eggnog5.0
#cat *.chunk_*.emapper.seed_orthologs > input_file.emapper.seed_orthologs
#othology and functional annotation
#emapper.py --annotate_hits_table input_file.emapper.seed_orthologs --no_file_comments -o output_file --cpu 8


# Cazy database can be moved to: http://bcb.unl.edu/dbCAN2/
# use the integrated python-workflow for automatic annotation: https://github.com/linnabrown/run_dbcan
# use docker for a easy installation: seemingly based on v7
# rebuild with conda
#source activate run_dbcan
# not work on meta mode used diamond for direct blast
#run_dbcan.py combined.linclust_f.fna meta --dia_cpu 8 --hmm_cpu 8 --hotpep_cpu 8 --tf_cpu 8 --stp_cpu 8 --db_dir /home/yanhui/Data/depot_tool/run_dbcan/db --out_dir output_dbcan
#source deactivate
#mkdir CAZy
#diamond blastx --query combined.linclust_f.fna --db /home/yanhui/Data/depot_tool/run_dbcan/db/CAZy --out CAZy/CAZy.tab --outfmt 6 --sensitive --max-target-seqs 1 --evalue 1e-5  --block-size 20.0 --tmpdir CAZy/temp --index-chunks 1
# ARDB database moved to: https://card.mcmaster.ca/
#load pre-build diamond database path

#diamond blastx --query card/combined.linclust_f.fna --db /home/yanhui/sharedFatboy/database_META/card/card-data_20191023/card_20191023_homolog --out card/card.tab --outfmt 6 --sensitive --max-target-seqs 1 --evalue 1e-5  --block-size 20.0 --tmpdir card/temp --index-chunks 1

