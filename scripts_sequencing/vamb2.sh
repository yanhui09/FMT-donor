# use long scaffolds to bin > 2000 bps

##########computerome ///faster alignment with more cores
##!! ADD SAMPLE source for the fasta's header to keep unique.
# using the adjusted long_scaffolds.fasta with prepending sample-name, keep the header of each contig name unique
# calculate the coverage
# build bwa index
bwa index long_scaffolds.fasta
# map back the reads
parallel -j 10 --link --plus 'bwa mem -t 2 long_scaffolds.fasta {1} {2} | samtools sort -o {1/_1.fastq/}.bam -' ::: CLEAN_READS/*_1.fastq ::: CLEAN_READS/*_2.fastq
# collect the .bam files
mkdir vamb/bams
mv CLEAN_READS/*bam vamb/bams
################################### Yan's sever with CUDA/GPU
# concatenated fasta (unique names); And bam files
vamb --outdir output_v2 --fasta long_scaffolds.fasta --bamfiles vamb/bams/*bam --cuda -m 2000
# vamb2 ouput to bin-set fastas, only keep genomic bins, bin size > 200, 000bp (revise accordingly)
python3 vamb_fasta_above200000.py

exit
