#!/bin/bash
# activate the humann2 env
source activate humann2
# run strainphlan
mkdir strainphlan2_out_count
parallel -eta -j 6 'metaphlan2.py {} strainphlan2_out_count/{/.}_profile.txt --nproc 1 -t rel_ab_w_read_stats \
        --bowtie2out strainphlan2_out_count/{/.}_bowtie2.txt --samout strainphlan2_out_count/{/.}.sam.bz2 --input_type fastq \
        --bowtie2db /home/yanhui/Documents/Database/db_metaphlan2/db_v20/mpa_v20_m200 --mpa_pkl /home/yanhui/Documents/Database/db_metaphlan2/db_v20/mpa_v20_m200.pkl' \
        ::: cat_reads/*fastq # directory for clean fastqs

source deactivate
source activate graphlan
# determine the consensus of unique marker genes for each species in each sample
mkdir consensus_marker
# specify the bcftools and samtools which Nicola Segate provided
parallel -eta -j 5 'sample2markers.py --nproc 1 \
	--samtools_exe /home/yanhui/anaconda3/envs/graphlan/extrabin/samtools \
	--bcftools_exe /home/yanhui/anaconda3/envs/graphlan/extrabin/bcftools \
	--ifn_samples {} --input_type sam --output_dir consensus_marker/' \
	::: strainphlan2_out_count/*.sam.bz2
# extract the markers specially for species with interst

extract_markers.py --mpa_pkl /home/yanhui/Documents/Database/db_metaphlan2/db_v20/mpa_v20_m200.pkl \
	--ifn_markers /home/yanhui/Documents/Database/db_metaphlan2/db_v20/mpa_v20_m200.fna \
	--clade s__Lactobacillus_reuteri --ofn_markers lactobacillus_reuteri/s__Lactobacillus_reuteri.fasta
# print clade
# strainphlan, adjust accordingly
strainphlan.py --ifn_samples consensus_marker/*.markers \
	--mpa_pkl /home/yanhui/Documents/Database/db_metaphlan2/db_v20/mpa_v20_m200.pkl \
	--ifn_markers lactobacillus_reuteri/s__Lactobacillus_reuteri.fasta --ifn_ref_genomes lactobacillus_reuteri/GCF*fna.bz2 \
	--output_dir lactobacillus_reuteri --clades s__Lactobacillus_reuteri --nprocs_main 2


# construct tree with ggtree
strainphlan_ggtree.R lactobacillus_reuteri/RAxML_bestTree.s__Lactobacillus_reuteri.tree mapping.tsv lactobacillus_reuteri/s__Lactobacillus_reuteri.fasta lactobacillus_reuteri/strainphlan_gtree_1.png lactobacillus_reuteri/strainphlan_gtree_2.png


#############Lactobacillus crispatus
extract_markers.py --mpa_pkl /home/yanhui/Documents/Database/db_metaphlan2/db_v20/mpa_v20_m200.pkl \
        --ifn_markers /home/yanhui/Documents/Database/db_metaphlan2/db_v20/mpa_v20_m200.fna \
        --clade s__Lactobacillus_crispatus --ofn_markers lactobacillus_crispatus/s__Lactobacillus_crispatus.fasta
# print clade
# strainphlan, adjust accordingly
strainphlan.py --ifn_samples consensus_marker/*.markers \
        --mpa_pkl /home/yanhui/Documents/Database/db_metaphlan2/db_v20/mpa_v20_m200.pkl \
        --ifn_markers s__Lactobacillus_crispatus.fasta --output_dir lactobacillus_crispatus --clades s__Lactobacillus_crispatus --nprocs_main 2

# construct tree with ggtree
strainphlan_ggtree.R lactobacillus_crispatus/RAxML_bestTree.s__Lactobacillus_crispatus.tree mapping.tsv lactobacillus_crispatus/s__Lactobacillus_crispatus.fasta lactobacillus_crispatus/strainphlan_gtree_1.png lactobacillus_crispatus/strainphlan_gtree_2.png


source deactivate


exit
