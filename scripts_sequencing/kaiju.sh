#!/bin/sh
# enter working directory
#cd $WORKING_DIR
mkdir Kaiju
# Run kaiju
parallel -j 1 --link 'kaiju -t /home/projects/cu_10168/apps/yanhui/kaiju/kaiju_nr_euk_2019-06-25/nodes.dmp -f /home/projects/cu_10168/apps/yanhui/kaiju/kaiju_nr_euk_2019-06-25/nr_euk/kaiju_db_nr_euk.fmi -i {1} -j {2} -z 28 -o Kaiju/kaiju_{1/.}.out' ::: CLEAN_READS/*_1.fastq ::: CLEAN_READS/*_2.fastq

# clean file format
parallel -j 1 --plus 'mv {} {/_1\.out/.out}' ::: Kaiju/kaiju*out

# summarize kaiju result on phylum, genus and species level
node_path=/home/projects/cu_10168/apps/yanhui/kaiju/kaiju_nr_euk_2019-06-25/nodes.dmp
name_path=/home/projects/cu_10168/apps/yanhui/kaiju/kaiju_nr_euk_2019-06-25/names.dmp
# species level 
kaiju2table -t $node_path -n $name_path -r species -e -l superkingdom,phylum,class,order,family,genus,species -o kaiju_summary_nr_species.tsv *out
# genus level
kaiju2table -t $node_path -n $name_path -r genus -e -l superkingdom,phylum,class,order,family,genus -o kaiju_summary_nr_genus.tsv *out
# phylum level
kaiju2table -t $node_path -n $name_path -r species -e -l superkingdom,phylum -o kaiju_summary_nr_phylum.tsv *out


exit
