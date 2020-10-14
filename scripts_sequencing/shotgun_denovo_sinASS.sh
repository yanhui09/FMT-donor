#!/bin/sh
# activate metawrap conda environment
conda activate metawrap-env
# Enter the working directoy where stored your sequencing data
#cd $WORKING_DIR
# 1.Read-QC
# parallel unzip, all fastq.gz under a directory of RAW_READS 
parallel -j 28 --will-cite 'gzip -d {}' ::: RAW_READS/*fastq.gz

# parallel for read_qc
mkdir -p READ_QC
for F in RAW_READS/*_1.fastq; do
        R=${F%_*}_2.fastq
        BASE=${F##*/}
        SAMPLE=${BASE%_*}
        metawrap read_qc -1 $F -2 $R -t 1 -x pig -o READ_QC/$SAMPLE &
done

parallel -j 32 --plus --link --will-cite '
#metawrap read_qc -1 {1} -2 {2} -t 1 -x pig -o {1/_1.fastq/}' \
::: RAW_READS/*_1.fastq ::: RAW_READS/*_2.fastq
wait

# move the read_qc result
mkdir -p READ_QC
mv RAW_READS/*/ READ_QC

# get final clean reads
mkdir -p CLEAN_READS
for i in READ_QC/*; do 
	b=${i#*/}
	mv ${i}/final_pure_reads_1.fastq CLEAN_READS/${b}_1.fastq
	mv ${i}/final_pure_reads_2.fastq CLEAN_READS/${b}_2.fastq
done
wait

# 2.Individually Assembly
PROJECTPRE='17119'
mkdir ASSEMBLY_INDIVIDUAL
parallel -j 4 --link --will-cite '
metawrap assembly -1 {1} -2 {2} -m 100 -t 6 --metaspades -o ASSEMBLY_INDIVIDUAL/{1/.} ' \
::: CLEAN_READS/*_1.fastq ::: CLEAN_READS/*_2.fastq 
wait 
# change the file name
parallel -j 24 --will-cite '
mv {} {/_1/}' ::: ASSEMBLY_INDIVIDUAL/*_1
wait

# 3. binning with three algorithms with minimum contig length 2000 in order to keep consistent with vamb
# prepend sample names for respective contigs
parallel -j 4 --plus 'sed "s/>/>{/}_/g" {}/final_assembly.fasta > {}/final_assembly_adj.fasta ' ::: ASSEMBLY_INDIVIDUAL/*

# remove the contigs with length below 2000bp
# using the python script inside metawrap (python = 2.7)
parallel -j 8 './rm_short_contigs.py 2000 {}/final_assembly_adj.fasta > {}/long_scaffolds.fasta' ::: ASSEMBLY_INDIVIDUAL/*

# concatenate all the long scaffold together
cat ASSEMBLY_INDIVIDUAL/*/long_scaffolds.fasta > long_scaffolds.fasta


####### create multiple-level files for fasta
parallel -j 2 'mkdir -p INITIAL_BINNING/{//}' ::: ASSEMBLY_INDIVIDUAL/${PROJECTPRE}*/long_scaffolds.fasta

parallel -j 6 'metawrap binning -o INITIAL_BINNING/{//} -t 28 -l 2000 -m 120 -a {} --metabat2 --maxbin2 --concoct CLEAN_READS/${PROJECTPRE}*fastq' ::: ASSEMBLY_INDIVIDUAL/${PROJECTPRE}*/long_scaffolds.fasta

# 4. consolidate bin sets with the bin_refinement module
# Bin with metabat2, maxbin2, concoct

# first derep each binning and then do the refinement
# renaming all the bins
# metabat2
mkdir -p indiviBINS_drep/metabat2_bins
parallel -j 2 'mkdir -p indiviBINS_drep/metabat2_bins/{/}' ::: INITIAL_BINNING/ASSEMBLY_INDIVIDUAL/${PROJECTPRE}*
wait
parallel -j 2 'cp {}/metabat2_bins/*fa indiviBINS_drep/metabat2_bins/{/}' ::: INITIAL_BINNING/ASSEMBLY_INDIVIDUAL/${PROJECTPRE}*
wait
parallel -j 2 'mv {} {//}_{/}' ::: indiviBINS_drep/metabat2_bins/${PROJECTPRE}*/*fa
rmdir indiviBINS_drep/metabat2_bins/${PROJECTPRE}*[0-9]
rm indiviBINS_drep/metabat2_bins/*unbinned*fa
# maxbin2
mkdir -p indiviBINS_drep/maxbin2_bins
parallel -j 2 'mkdir -p indiviBINS_drep/maxbin2_bins/{/}' ::: INITIAL_BINNING/ASSEMBLY_INDIVIDUAL/${PROJECTPRE}*
wait
parallel -j 2 'cp {}/maxbin2_bins/*fa indiviBINS_drep/maxbin2_bins/{/}' ::: INITIAL_BINNING/ASSEMBLY_INDIVIDUAL/${PROJECTPRE}*
wait
parallel -j 2 'mv {} {//}_{/}' ::: indiviBINS_drep/maxbin2_bins/${PROJECTPRE}*/*fa
rmdir indiviBINS_drep/maxbin2_bins/${PROJECTPRE}*[0-9]
rm indiviBINS_drep/maxbin2_bins/*unbinned*fa
# concoct
mkdir -p indiviBINS_drep/concoct_bins
parallel -j 2 'mkdir -p indiviBINS_drep/concoct_bins/{/}' ::: INITIAL_BINNING_underep/ASSEMBLY_INDIVIDUAL/${PROJECTPRE}*
wait
parallel -j 2 'cp {}/concoct_bins/*fa indiviBINS_drep/concoct_bins/{/}' ::: INITIAL_BINNING_underep/ASSEMBLY_INDIVIDUAL/${PROJECTPRE}*
wait
parallel -j 2 'mv {} {//}_{/}' ::: indiviBINS_drep/concoct_bins/${PROJECTPRE}*/*fa
rmdir indiviBINS_drep/concoct_bins/${PROJECTPRE}*[0-9]
rm indiviBINS_drep/concoct_bins/*unbinned*fa

# run drep on all the metabat2, maxbin2 and concoct separately
parallel -j 1 'dRep dereplicate -p 28 -comp 50 -con 5 {}_drep -g {}/*.fa ' ::: indiviBINS_drep/metabat2_bins
parallel -j 1 'dRep dereplicate -p 28 -comp 50 -con 5 {}_drep -g {}/*.fa ' ::: indiviBINS_drep/maxbin2_bins
parallel -j 1 'dRep dereplicate -p 28 -comp 50 -con 5 {}_drep -g {}/*.fa ' ::: indiviBINS_drep/concoct_bins 
# run bin_refinement module on the dreped 3 bins
metawrap bin_refinement -o BIN_REFINEMENT_indiviASSBIN_Drep -t 28 -A indiviBINS_drep/metabat2_bins_drep/dereplicated_metabat2 -B indiviBINS_drep/maxbin2_bins_drep/dereplicated_maxbin2 -C indiviBINS_drep/concoct_bins_drep/dereplicated_concoct -c 50 -x 5
# run bin_refinement on consilated bins and vamb
metawrap bin_refinement -o BIN_REFINEMENT_indiviASSBIN_Drep2_vambv2 -t 28 -A BIN_REFINEMENT_indiviASSBIN_Drep/mmc_50_5_bins -B BIN_REFINEMENT_indiviASSBIN_Drep/vamb_bins -c 50 -x 5

# 5. quantify refined MAGsonly use single-assembly
mv BIN_REFINEMENT_indiviASSBIN_Drep2_vambv2/metawrap_50_5_bins mags 

metawrap quant_bins -b mags -a long_scaffolds.fasta -o QUANT_BINS -t 28 CLEAN_READS/${PROJECTPRE}*fastq
wait

# 6. functionally annotate bins ; Annotated by KEGG online server Ghostkoala
metaWRAP annotate_bins -o FUNCT_ANNOT -t 28 -b mags

conda deactivate
exit
