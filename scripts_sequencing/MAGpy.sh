#!/bin/bash

#run MAGpy
snakemake --use-conda -s MAGpy -j 6    
#draw the phylogenetic tree
perl scripts/produce_tree_diamond.pl diamond_bin_report_plus_manual_source-manual.tsv tree/MAGpy/MAGpy.tree.nwk
