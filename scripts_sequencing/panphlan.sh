#!/bin/bash
#　use panphlan_database.sh built lcr lre database, 2019.11.21 refseq
source activate graphlan

#################################
#2019.11.21 Refseq lactobacillus reuteri
# lactobacillus reuteri
parallel -j 1 --plus 'cat {/_1.fastq/}*.fastq |panphlan_map.py -c lre19 --i_bowtie2_indexes /home/yanhui/sharedFatboy/database_META/panphlan_new/database_lre19 --fastx fastq -o {/.}.csv --verbose -p 6' ::: ~/sharedFatboy/FMT_metagenome/FMT_temp/FMT/CLEAN_READS/17119*_1.fastq
# adjust the sample name
parallel -j 1 --plus 'mv {} {/_1_lreuteri19/}' ::: ./*bz2
# clean_up
mkdir lre19
mv *bz2 lre19

# merge absence presence matrix # sensitive mode to include more samples
panphlan_profile.py -c lre19 --i_bowtie2_indexes /home/yanhui/sharedFatboy/database_META/panphlan_new/database_lre19 -i ./ -o gene_presence_absence.tsv --add_strains --o_cov data_genefamily_coverage.tsv --min_coverage 1 --left_max 1.70 --right_min 0.30
# quick heatmap
panphlan_profile.py -c lre19 --i_bowtie2_indexes /home/yanhui/sharedFatboy/database_META/panphlan_new/database_lre19 -i ./ -o gene_presence_absence.tsv  --min_coverage 1 --left_max 1.70 --right_min 0.30

./format2hclust2.R gene_presence_absence.tsv ../../mapping.tsv gene_presence_absence_group.tsv
hclust2.py -i gene_presence_absence_group.tsv -o heatmap.png --metadata_rows 1 --legend_file legend.png --image_size 18 --cell_aspect_ratio 0.01 --dpi 300 --f_dist_f jaccard --s_dist_f jaccard --no_slabels --no_flabels
#################################
#2019.11.21 Refseq lactobacillus crispatus
# lactobacillus crispatus
parallel -j 1 --plus 'cat {/_1.fastq/}*.fastq |panphlan_map.py -c lcr19 --i_bowtie2_indexes /home/yanhui/sharedFatboy/database_META/panphlan_new/database_lcr19 --fastx fastq -o {/.}.csv --verbose -p 6' ::: ~/sharedFatboy/FMT_metagenome/FMT_temp/FMT/CLEAN_READS/17119*_1.fastq
# adjust the sample name
parallel -j 1 --plus 'mv {} {/_1/}' ::: ./*bz2
# clean_up
mkdir lcr19
mv *bz2 lcr19

# merge absence presence matrix # sensitive mode to include more samples
panphlan_profile.py -c lcr19 --i_bowtie2_indexes /home/yanhui/sharedFatboy/database_META/panphlan_new/database_lcr19 -i ./ -o gene_presence_absence.tsv --add_strains --o_cov data_genefamily_coverage.tsv --min_coverage 1 --left_max 1.70 --right_min 0.30
# quick heatmap
panphlan_profile.py -c lcr19 --i_bowtie2_indexes /home/yanhui/sharedFatboy/database_META/panphlan_new/database_lcr19 -i ./ -o gene_presence_absence.tsv --min_coverage 1 --left_max 1.70 --right_min 0.30

./format2hclust2.R gene_presence_absence.tsv ../../mapping.tsv gene_presence_absence_group.tsv
hclust2.py -i gene_presence_absence_group.tsv -o heatmap.png --metadata_rows 1 --legend_file legend.png --image_size 18 --cell_aspect_ratio 0.01 --dpi 300 --f_dist_f jaccard --s_dist_f jaccard --no_slabels --no_flabels
source deactivate
