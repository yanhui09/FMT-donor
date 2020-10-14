library(phyloseq)
library(vegan)

#pairwise adonois
pairwise_adonis <- function(dis_mat, design, groups){
  sub_design <- subset(design, Group %in% groups)
  sub_design$Group <- as.factor(sub_design$Group)
  dis_mat_f <- dis_mat[rownames(sub_design),rownames(sub_design)]
  adonis_sub <- adonis(dis_mat_f~Group, data=sub_design, permutations = 999)
  return(adonis_sub$aov.tab$`Pr(>F)`[1])
}
#pairwise adonois R2
pairwise_adonis_R2 <- function(dis_mat, design, groups){
  sub_design <- subset(design, Group %in% groups)
  sub_design$Group <- as.factor(sub_design$Group)
  dis_mat_f <- dis_mat[rownames(sub_design),rownames(sub_design)]
  adonis_sub <- adonis(dis_mat_f~Group, data=sub_design, permutations = 999)
  return(adonis_sub$aov.tab$R2[1])
}

# wrapped Adnois_paiwise function
adonis_pairwise <- function(ps,trans,dis) {
ps.sample <- readRDS(ps)
if (trans=="hellinger") {
  #Hellinger transformation
  ps.norm <- transform_sample_counts(ps.sample, function(x) sqrt(x/sum(x)*100))
} else if (trans=="rarefaction") {
  ps.norm = rarefy_even_depth(ps.sample, rngseed=1, sample.size=11000, replace=F)
} else {
  ps.norm = ps.sample  
}
# calculate distance
if (dis=="jaccard") {
  disWuf <- distance(ps.norm, method = dis, binary = T)
  } else {
  disWuf <- distance(ps.norm, method = dis)
  }

sub_design <- data.frame(sample_data(ps.sample))
disWuf <- as.matrix(disWuf)
disWuf <- disWuf[rownames(sub_design),rownames(sub_design)]
p12 <- pairwise_adonis(disWuf, sub_design, c("FMT1","FMT2"))
p1C <- pairwise_adonis(disWuf, sub_design, c("FMT1","CON"))
p2C <- pairwise_adonis(disWuf, sub_design, c("FMT2","CON"))

p_pairwise <- c(p12,p1C, p2C)
p_pairwise_adj <- p.adjust(p_pairwise,method = "BH")
p_pairwise_adj <- data.frame(padj=p_pairwise_adj)
rownames(p_pairwise_adj) <- c("FMT1 vs FMT2", "FMT1 vs CON", "FMT2 vs CON")
p_pairwise_adj

r12 <- pairwise_adonis_R2(disWuf, sub_design, c("FMT1","FMT2"))
r1C <- pairwise_adonis_R2(disWuf, sub_design, c("FMT1","CON"))
r2C <- pairwise_adonis_R2(disWuf, sub_design, c("FMT2","CON"))
r_pairwise <- c(r12,r1C, r2C)
r_pairwise <- data.frame(r2=r_pairwise)
rownames(r_pairwise) <- c("FMT1 vs FMT2", "FMT1 vs CON", "FMT2 vs CON")
r_pairwise

tab_pairwise <- cbind.data.frame(r_pairwise,p_pairwise, p_pairwise_adj)
return(tab_pairwise)
}

############################################################################
set.seed(1234)
# 16s, rarefied at 11000 counts,weighted and uweighted unifrac
tab_16s_wunifrac <- adonis_pairwise(ps="data/phyloseq_16s.rds",trans = "rarefaction", dis = "wunifrac") 
tab_16s_unifrac <- adonis_pairwise(ps="data/phyloseq_16s.rds",trans = "rarefaction", dis = "unifrac") 

# shotgun kaiju species, hellinger normalization, bray curtis, binary jaccard
tab_species_bray <- adonis_pairwise(ps="data/ps_kaiju_nr_species.rds",trans = "hellinger", dis = "bray") 
tab_species_jaccard <- adonis_pairwise(ps="data/ps_kaiju_nr_species.rds",trans = "hellinger", dis = "jaccard") 

# shotgun kaiju genus hellinger normalization, bray curtis, binary jaccard
tab_genus_bray <- adonis_pairwise(ps="data/ps_kaiju_nr_genus.rds",trans = "hellinger", dis = "bray") 
tab_genus_jaccard <- adonis_pairwise(ps="data/ps_kaiju_nr_genus.rds",trans = "hellinger", dis = "jaccard") 

# build a table
text4 <- function(x,y) paste("R2=",signif(x, digits = 2),", p.adj=", signif(y, digits = 2))
tab_16s_wunifrac$tab <- mapply(text4, tab_16s_wunifrac$r2, tab_16s_wunifrac$padj)
tab_16s_unifrac$tab <- mapply(text4, tab_16s_unifrac$r2, tab_16s_unifrac$padj)
tab_species_bray$tab <- mapply(text4, tab_species_bray$r2, tab_species_bray$padj)
tab_species_jaccard$tab <- mapply(text4, tab_species_jaccard$r2, tab_species_jaccard$padj)
tab_genus_bray$tab <- mapply(text4, tab_genus_bray$r2, tab_genus_bray$padj)
tab_genus_jaccard$tab <- mapply(text4, tab_genus_jaccard$r2, tab_genus_jaccard$padj)

tab_out <- cbind.data.frame(`Weighed Unifac (16S rRNA)`= tab_16s_wunifrac$tab,
                            `Unweighed Unifac (16S rRNA)`= tab_16s_unifrac$tab,
                            `Bray Curtis (Shotgun species)`= tab_species_bray$tab,
                            `Binary Jaccard (Shotgun species)`= tab_species_jaccard$tab,
                            `Bray Curtis (Shotgun genus)`= tab_genus_bray$tab,
                            `Binary Jaccard (Shotgun genus)`= tab_genus_jaccard$tab)
rownames(tab_out) <- rownames(tab_16s_wunifrac)

write.table(tab_out,file = "table/table_1.tsv", sep = "\t", col.names = NA, quote = FALSE)
