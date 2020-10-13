#!/bin/sh
# activate metawrap environment
conda activate metawrap-env
# ender the working directory
# cd $WOEKING_DIR
# start from clean-reads
# 2.Co-Assembly
PROJECTPRE='17119'
cat CLEAN_READS/${PROJECTPRE}*_1.fastq > CLEAN_READS/ALL_READS_1.fastq
cat CLEAN_READS/${PROJECTPRE}*_2.fastq > CLEAN_READS/ALL_READS_2.fastq
wait

metawrap assembly -1 CLEAN_READS/ALL_READS_1.fastq -2 CLEAN_READS/ALL_READS_2.fastq -m 990 -t 28 --megahit -o CO_ASSEMBLY
wait

# 3. binning with three algorithms
metawrap binning -o INITIAL_BINNING_COassembly -t 28 -l 2000 -m 120 -a CO_ASSEMBLY/final_assembly.fasta --metabat2 --maxbin2 --concoct CLEAN_READS/${PROJECTPRE}*fastq
wait

# 4. consolidate bin sets with the bin_refinement module
# 4.1 consolide the co-binning with metabat2, maxbin2, concoct
metawrap bin_refinement -o BIN_REFINEMENT_COassembly -t 28 -A INITIAL_BINNING_COassembly/metabat2_bins/ -B INITIAL_BINNING_COassembly/maxbin2_bins/ -C INITIAL_BINNING_COassembly/concoct_bins/ -c 50 -x 5 
# 4.2 consolide the co-binning with vamb
metawrap bin_refinement -o BIN_REFINEMENT_COassembly2_v2 -t 28 -A BIN_REFINEMENT_COassembly/mmc_50_5_bins -B CO_ASSEMBLY/vamb_bins -c 50 -x 5

conda deactivate

exit

