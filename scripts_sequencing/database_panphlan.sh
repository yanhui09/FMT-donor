#!/bin/bash

# download all the assembly from the refseq
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
# extract all genome download-links from file assembly_summary.txt related to a species
read -p "Please input the selected species name:" A
read -p "Please define the database name:" B
# include the complete genome, contig and scaffolds
grep -E "$A" assembly_summary.txt | cut -f 20 > ftpdirpaths
# create download scripts
awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget "ftpdir,file}' ftpdirpaths > download_fna_files.sh

#awk 'BEGIN{FS=OFS="/";filesuffix="genomic.gff.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget "ftpdir,file}' ftpdirpaths > download_gff_files.sh

head -3 download_*
# run download
source download_fna_files.sh
#source download_gff_files.sh
# uncompress .gz files
gzip -d *.gz
# rename fna files; keep only GCF-number(remove all behind first dot)
for FILE in *.fna; do
 mv ${FILE} ${FILE%%.*}.fna
done;

#ls
# move the gff and fna to respective groups
#mkdir ncbi_download_fna ncbi_download_gff
mkdir ncbi_download_fna
mv *.fna ncbi_download_fna
#mv *.gff ncbi_download_gff
# use prokka to reannotate all genomes
source activate graphlan # conda env with prokka
# prokka annotation
for i in ncbi_download_fna/*.fna; do
	FILE=$(basename $i)
	prokka --quiet --cpu 6 --outdir prokka_gff --prefix ${FILE%%.*} $i --force
done;

#mkdir prokka_gff
#mv prokka/*.gff prokka_gff

panphlan_pangenome_generation.py -c $B --i_gff prokka_gff/ -o database_$B
# adjust the directory name
mv prokka_gff prokka_gff_$B
mv ncbi_download_fna ncbi_download_fna_$B
mv ffn_from_gff ffn_from_gff_$B
mv fna_from_gff fna_from_gff_$B

