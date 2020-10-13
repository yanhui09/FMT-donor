#lre19 table_s1
# match the prokka annotation and the kegg
anno_tab <- read.table("data/panphlan/panphlan_lre19_annotations.tsv", sep = "\t", header = T, row.names = 1,quote = "")
anno_tab$geneID <- rownames(anno_tab)

kegg_tab <- read.table("data/panphlan/prokka2kegg_lre19_merged.out", sep = "\t", header = F, fill = T, quote = "")
# match the two table
colnames(kegg_tab) <- c("GFF_locus_tag", "koID")
anno_tab_c <- merge(x=anno_tab,y=kegg_tab, by = "GFF_locus_tag", all.x = TRUE)
#import the gene abundance table
abun_tab <- read.table("data/panphlan/gene_presence_absence_lre19.tsv", sep = "\t", header = T, row.names = 1, quote = "", check.names = F)

mapping <- read.table("data/panphlan/mapping.tsv",sep = "\t", header = T, row.names = 1, quote = "")
#ref detailed info (only select "complete genomes")
ref <- read.table("data/panphlan/lre19_ref_cat.txt", sep = "\t", header = F, row.names = 1, quote = "")
ref_cg <- subset(ref, V2=="Complete Genome")
ref_cg_lab <- sapply(rownames(ref_cg), function(x) paste0("REF_",x))
ref_cg_lab <- stringr::str_replace(ref_cg_lab,"\\.[0-9]","")
#only include the data and the ref

sample_infer <- rownames(mapping)[rownames(mapping) %in% colnames(abun_tab)]
abun_tab_f <- abun_tab[,c(sample_infer,ref_cg_lab)]
abun_tab_f_out <- abun_tab_f[rowSums(abun_tab_f)!=0,]
# write out the data frame
write.table(abun_tab_f_out,"data/panphlan/gene_presence_absence_lre19_cg.tsv", sep = "\t", col.names = NA)
  
# fisher exact test on the abun_tab_f on FMT1 and FMT2

library(tidyr)
abun_tab_f$geneID <- rownames(abun_tab_f)
abun_tab_f_long <- reshape2::melt(abun_tab_f,id.vars="geneID",variable.name="SampleID", value.name="Presence")
head(abun_tab_f_long)
# only take the FMT1 and FMT2 samples
index <- match(abun_tab_f_long$SampleID, rownames(mapping))
abun_tab_f_long$Group <- mapping$Group[index]
abun_tab_f_long$Presence <- factor(abun_tab_f_long$Presence,labels = c("absent","present"))

abun_tab_f_long_f <- subset(abun_tab_f_long, Group %in% c("FMT1","FMT2")) 
abun_tab_f_long_f$Group <- factor(abun_tab_f_long_f$Group)
 
#############################fisher test
fisher_p <- function(i){
abun_tab_f_long_ff <- subset(abun_tab_f_long_f,geneID==i)
test <- fisher.test(xtabs(~Presence+Group, abun_tab_f_long_ff))
return(test$p.value)
}

p_val <- fisher_p(unique(abun_tab_f_long_f$geneID)[1])
for (i in unique(abun_tab_f_long_f$geneID)[-1]) {
  p_val <- c(p_val, fisher_p(i))
}

p_tab <- data.frame(cbind(geneID=unique(abun_tab_f_long_f$geneID),p=p_val))

#add two donors
donor_name <- rownames(mapping)[mapping$Group %in% c("DONOR1", "DONOR2")]
index <- match(p_tab$geneID, abun_tab_f$geneID)
p_tab <- cbind.data.frame(p_tab, abun_tab_f[index,donor_name])
index <- match(colnames(p_tab)[c(3,4)],rownames(mapping))
colnames(p_tab) <- c("geneID","p",as.character(mapping$Group[index]))
# add mean relative abundance for FMT1 and FMT2
index <- match(p_tab$geneID, abun_tab_f$geneID)
FMT1_name <- rownames(mapping)[mapping$Group == "FMT1"]
FMT1 <- rowMeans(abun_tab_f[index,colnames(abun_tab_f) %in% FMT1_name])
FMT2_name <- rownames(mapping)[mapping$Group == "FMT2"]
FMT2 <- rowMeans(abun_tab_f[index,colnames(abun_tab_f) %in% FMT2_name])
p_tab$FMT1 <- FMT1
p_tab$FMT2 <- FMT2
#integrated with annotation tsv
p_tab <- cbind.data.frame(p_tab, anno_tab_c[index,])
write.table(p_tab, "data/panphlan/fisher_p_lre19.tsv", sep = "\t", col.names = NA)
# take significant hits (fisher p < 0.05)
p_tab <- as.data.frame(p_tab)
p_tab$p <- as.numeric(levels(p_tab$p))[p_tab$p]
p_tab_f <- subset(p_tab, p < 0.05)
write.table(p_tab_f, "data/panphlan/fisher_p_lre19_sig.txt", sep = "\t", col.names = NA)
# inlude the different genes and KEGG annotation to table
# load in the kegg data
#lre19
lre19 <- read.table("data/panphlan/fisher_p_lre19_sig.txt", sep = "\t", header = T, row.names = 1)
#use table download from keggrest to add the annotation and module pathway, information
ko_anno <- read.table("data/panphlan/ko_anno.tsv", sep = "\t", header = F, check.names = F, quote = "")
#ko_anno$V1 <- stringr::str_replace(ko_anno$V1,"ko:","")
# module
module_anno <- read.table("data/panphlan/module_anno.tsv", sep = "\t", header = F, check.names = F, quote = "")
#module_anno$V1 <- stringr::str_replace(module_anno$V1, "md:","")
# ko to module
ko2module <- read.table("data/panphlan/ko2module.tsv", sep = "\t", check.names = F, quote = "")
# pathway
pathway_anno <- read.table("data/panphlan/pathway_anno.tsv", sep = "\t", header = F, check.names = F, quote = "")
# ko to pathway
ko2pathway <- read.table("data/panphlan/ko2pathway.tsv", sep = "\t", header = F, check.names = F,quote = "")
# add ko, module, pathway annotation
# ko annotation
lre19$koID_m <- paste0("ko:",lre19$koID)
index <- match(lre19$koID_m, ko_anno$V1)
index
lre19$ko_anno <- ko_anno$V2[index]
# module
index <- match(lre19$koID_m, ko2module$V2)
lre19$module_ID <- ko2module$V1[index]
# module annotation
index <- match(lre19$module_ID,module_anno$V1)
lre19$module_anno <- module_anno$V2[index]
# pathway
index <- match(lre19$koID_m, ko2pathway$V2)
lre19$pathway_ID <- ko2pathway$V1[index]
# pathway annotation
index <- match(lre19$pathway_ID, pathway_anno$V1)
lre19$pathway_anno <- pathway_anno$V2[index]
lre19_out <- subset(lre19, select = c(geneID, p, DONOR1, DONOR2, FMT1, FMT2, Product, EC_number, koID, ko_anno, module_ID, module_anno, pathway_ID, pathway_anno))
write.table(lre19_out, "table/table_s1.tsv", sep="\t", row.names = F, quote = F)


#########################################################
#lcr19 table_s2
# match the prokka annotation and the kegg
anno_tab <- read.table("data/panphlan/panphlan_lcr19_annotations.tsv", sep = "\t", header = T, row.names = 1,quote = "")
anno_tab$geneID <- rownames(anno_tab)

kegg_tab <- read.table("data/panphlan/prokka2kegg_lcr19_merged.out", sep = "\t", header = F, fill = T, quote = "")
# match the two table

colnames(kegg_tab) <- c("GFF_locus_tag", "koID")
anno_tab_c <- merge(x=anno_tab,y=kegg_tab, by = "GFF_locus_tag", all.x = TRUE)
#import the gene abundance table
abun_tab <- read.table("data/panphlan/gene_presence_absence_lcr19.tsv", sep = "\t", header = T, row.names = 1, quote = "", check.names = F)
tail(abun_tab)
mapping <- read.table("data/panphlan/mapping.tsv",sep = "\t", header = T, row.names = 1, quote = "")
#ref detailed info (only select "complete genomes")
ref <- read.table("data/panphlan/lcr19_ref_cat.txt", sep = "\t", header = F, row.names = 1, quote = "")
ref_cg <- subset(ref, V2=="Complete Genome")
ref_cg_lab <- sapply(rownames(ref_cg), function(x) paste0("REF_",x))
ref_cg_lab <- stringr::str_replace(ref_cg_lab,"\\.[0-9]","")
#only include the data and the ref
sample_infer <- rownames(mapping)[rownames(mapping) %in% colnames(abun_tab)]
abun_tab_f <- abun_tab[,c(sample_infer,ref_cg_lab)]
abun_tab_f_out <- abun_tab_f[rowSums(abun_tab_f)!=0,]
# write out the data frame
write.table(abun_tab_f_out,"data/panphlan/gene_presence_absence_lcr19_cg.tsv", sep = "\t", col.names = NA)

# fisher exact test on the abun_tab_f on FMT1 and FMT2

library(tidyr)
abun_tab_f$geneID <- rownames(abun_tab_f)
abun_tab_f_long <- reshape2::melt(abun_tab_f,id.vars="geneID",variable.name="SampleID", value.name="Presence")
head(abun_tab_f_long)
# only take the FMT1 and FMT2 samples
index <- match(abun_tab_f_long$SampleID, rownames(mapping))
abun_tab_f_long$Group <- mapping$Group[index]
abun_tab_f_long$Presence <- factor(abun_tab_f_long$Presence,labels = c("absent","present"))

abun_tab_f_long_f <- subset(abun_tab_f_long, Group %in% c("FMT1","FMT2")) 
abun_tab_f_long_f$Group <- factor(abun_tab_f_long_f$Group)

#############################fisher test
fisher_p <- function(i){
  abun_tab_f_long_ff <- subset(abun_tab_f_long_f,geneID==i)
  test <- fisher.test(xtabs(~Presence+Group, abun_tab_f_long_ff))
  return(test$p.value)
}

p_val <- fisher_p(unique(abun_tab_f_long_f$geneID)[1])
for (i in unique(abun_tab_f_long_f$geneID)[-1]) {
  p_val <- c(p_val, fisher_p(i))
}

p_tab <- data.frame(cbind(geneID=unique(abun_tab_f_long_f$geneID),p=p_val))

#add two donors
donor_name <- rownames(mapping)[mapping$Group %in% c("DONOR1", "DONOR2")]
index <- match(p_tab$geneID, abun_tab_f$geneID)
p_tab <- cbind.data.frame(p_tab, abun_tab_f[index,donor_name])
index <- match(colnames(p_tab)[c(3,4)],rownames(mapping))
colnames(p_tab) <- c("geneID","p",as.character(mapping$Group[index]))
# add mean relative abundance for FMT1 and FMT2
index <- match(p_tab$geneID, abun_tab_f$geneID)
FMT1_name <- rownames(mapping)[mapping$Group == "FMT1"]
FMT1 <- rowMeans(abun_tab_f[index,colnames(abun_tab_f) %in% FMT1_name])
FMT2_name <- rownames(mapping)[mapping$Group == "FMT2"]
FMT2 <- rowMeans(abun_tab_f[index,colnames(abun_tab_f) %in% FMT2_name])
p_tab$FMT1 <- FMT1
p_tab$FMT2 <- FMT2
#integrated with annotation tsv
p_tab <- cbind.data.frame(p_tab, anno_tab_c[index,])
write.table(p_tab, "data/panphlan/fisher_p_lcr19.tsv", sep = "\t", col.names = NA)

##only take the different enriched gene families
# take significant hits (fisher p < 0.05)
p_tab <- as.data.frame(p_tab)
p_tab$p <- as.numeric(levels(p_tab$p))[p_tab$p]
p_tab_f <- subset(p_tab, p < 0.05)
write.table(p_tab_f, "data/panphlan/fisher_p_lcr19_sig.txt", sep = "\t", col.names = NA)

# inlude the different genes and KEGG annotation to table
# load in the kegg data
#lcr19
lcr19 <- read.table("data/panphlan/fisher_p_lcr19_sig.txt", sep = "\t", header = T, row.names = 1)
#use table download from keggrest to add the annotation and module pathway, information
ko_anno <- read.table("data/panphlan/ko_anno.tsv", sep = "\t", header = F, check.names = F, quote = "")
# module
module_anno <- read.table("data/panphlan/module_anno.tsv", sep = "\t", header = F, check.names = F, quote = "")
# ko to module
ko2module <- read.table("data/panphlan/ko2module.tsv", sep = "\t", check.names = F, quote = "")
# pathway
pathway_anno <- read.table("data/panphlan/pathway_anno.tsv", sep = "\t", header = F, check.names = F, quote = "")
# ko to pathway
ko2pathway <- read.table("data/panphlan/ko2pathway.tsv", sep = "\t", header = F, check.names = F,quote = "")
# add ko, module, pathway annotation
# ko annotation
lcr19$koID_m <- paste0("ko:",lcr19$koID)
index <- match(lcr19$koID_m, ko_anno$V1)
index
lcr19$ko_anno <- ko_anno$V2[index]
# module
index <- match(lcr19$koID_m, ko2module$V2)
lcr19$module_ID <- ko2module$V1[index]
# module annotation
index <- match(lcr19$module_ID,module_anno$V1)
lcr19$module_anno <- module_anno$V2[index]
# pathway
index <- match(lcr19$koID_m, ko2pathway$V2)
lcr19$pathway_ID <- ko2pathway$V1[index]
# pathway annotation
index <- match(lcr19$pathway_ID, pathway_anno$V1)
lcr19$pathway_anno <- pathway_anno$V2[index]
lcr19_out <- subset(lcr19, select = c(geneID, p, DONOR1, DONOR2, FMT1, FMT2, Product, EC_number, koID, ko_anno, module_ID, module_anno, pathway_ID, pathway_anno))
write.table(lcr19_out, "table/table_s2.tsv", sep="\t", row.names = F, quote = F)
