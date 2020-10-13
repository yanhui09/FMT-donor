#Figure.4. heatmap for deseq2 comparison
########################################################
library(phyloseq)
library(DESeq2)
#kaiju deseq2
ps.sample <- readRDS("data/ps_kaiju_nr_species.rds")

#quick analysis with phyloseq
#remove donor samples
ps.sample <- prune_samples(sample_data(ps.sample)$Group!="DONOR1"&sample_data(ps.sample)$Group!="DONOR2", ps.sample)
ps.sample
# transform to despe2
library(stringr)
### loop for all the grouping compared with "Media-GM controls"
Groups <- levels(sample_data(ps.sample)$Group)
Groups
Groups <- Groups[Groups!="CON"]
path_table <- "data/desep2_kaiju/"
dir.create(path_table)
for (group in Groups){
  ps <- prune_samples(sample_data(ps.sample)$Group %in% c(group, "CON"), ps.sample)
  ps
  # remove all error taxa
  ps.ds <- phyloseq_to_deseq2(ps, ~Group)
  # solve rows without a zero, deseq need to calculate the geometric zero, 
  cts <- counts(ps.ds)
  geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
  dds <- estimateSizeFactors(ps.ds, geoMeans=geoMeans)
  ps.ds <-  DESeq(dds, test="Wald", fitType="parametric")
  # result
  res = results(ps.ds, cooksCutoff = FALSE)
  #alpha = 0.0001
  sigtab = res
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
  head(sigtab)
  
  library("ggplot2")
  theme_set(theme_bw())
  scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
  }
  # Phylum order
  x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
  x = sort(x, TRUE)
  sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
  # Genus order
  x = tapply(sigtab$log2FoldChange, sigtab$Species, function(x) max(x))
  x = sort(x, TRUE)
  sigtab$Species = factor(as.character(sigtab$Species), levels=names(x))
  write.table(data.frame(sigtab), paste0(path_table,str_replace(group," ",""),"_CON.tsv"), sep="\t", col.names = NA)
}
##############################################
#FMT1 and FMT2
Groups <- levels(sample_data(ps.sample)$Group)
Groups
Groups <- "FMT1"
path_table <- "data/desep2_kaiju/"
dir.create(path_table)
for (group in Groups){
  ps <- prune_samples(sample_data(ps.sample)$Group %in% c(group, "FMT2"), ps.sample)
  ps
  # remove all error taxa
  ps.ds <- phyloseq_to_deseq2(ps, ~Group)
  # solve rows without a zero, deseq need to calculate the geometric zero, 
  cts <- counts(ps.ds)
  geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
  dds <- estimateSizeFactors(ps.ds, geoMeans=geoMeans)
  ps.ds <-  DESeq(dds, test="Wald", fitType="parametric")
  # result
  res = results(ps.ds, cooksCutoff = FALSE)
  sigtab = res
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
  head(sigtab)
  
  library("ggplot2")
  theme_set(theme_bw())
  scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
  }
  # Phylum order
  x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
  x = sort(x, TRUE)
  sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
  # Genus order
  x = tapply(sigtab$log2FoldChange, sigtab$Species, function(x) max(x))
  x = sort(x, TRUE)
  sigtab$Species = factor(as.character(sigtab$Species), levels=names(x))
  write.table(data.frame(sigtab), paste0(path_table,str_replace(group," ",""),"_FMT2.tsv"), sep="\t", col.names = NA)
}

#################################################################################
#draw heatmap 
#Include deseq2 p.adj FMT2 vs CON, FMT2 vs FMT1, taxonomy for rows
#Include Litter Information and Group NEC Information for colors
#load two list of deseq2
tab_FMT2CON_all <- read.table("table/desep2_kaiju/FMT2_CON.tsv", sep = "\t",header = T, row.names = 1)
tab_FMT2FMT1_all <- read.table("table/desep2_kaiju/FMT1_FMT2.tsv", sep = "\t", header = T, row.names = 1)
tab_FMT1CON_all <- read.table("table/desep2_kaiju/FMT1_CON.tsv", sep = "\t",header = T, row.names = 1)

tab_FMT2CON <- subset(tab_FMT2CON_all, padj < 0.05)
tab_FMT2FMT1 <-subset(tab_FMT2FMT1_all, padj < 0.05)
tab_FMT1CON <- subset(tab_FMT1CON_all, padj < 0.05)

ASV <- unique(c(rownames(tab_FMT2CON), rownames(tab_FMT2FMT1)))

ps.sample <- readRDS("data/ps_kaiju_nr_species.rds")
ps.sample.rel <- transform_sample_counts(ps.sample, function(x) x/sum(x)*100)
ps.sample.rel.sig <- prune_taxa(rownames(otu_table(ps.sample.rel)) %in% ASV,ps.sample.rel)

#select the rel-abun > 0.1%

# at least 1% relative abundance appearance in 10% samples
mat <- as.matrix(otu_table(ps.sample.rel.sig))
species2keep <- rownames(mat)[rowSums(mat>=1)/length(colnames(mat))> 0.1]
ps.sample.rel.sig <- prune_taxa(species2keep,ps.sample.rel.sig)

otu_abun_select <- data.frame(otu_table(ps.sample.rel.sig), check.names = F)
#import relavant metadata
metadata <- read.table("data/Metadata.tsv", sep = "\t", header = T, row.names = 1, check.names = F)
tax.clean <- data.frame(tax_table(ps.sample.rel.sig))
# use complex heatmap for visualization, multiple row annotation
# extract the relative abundance matrix for thoes responding signals
library(ComplexHeatmap)
library(dplyr)
library(circlize)
# create a variable to define the subgroup, CON, FMT1, FMT2, DONOR1, DONOR2 ...
# order the matrix by the subgroup
metadata$Group <- factor(metadata$Group, levels = c("CON","FMT1","FMT2","DONOR1","DONOR2"))
meta_order <- metadata[order(metadata$Group),]
# re_order the col
mat <- otu_abun_select
mat <- as.matrix(mat[,rownames(meta_order)])

base_mean = rowMeans(mat)
mat_scaled = t(scale(t(mat)))
# calculate heatmap annotation
tax_heatmap <- tax.clean[rownames(mat_scaled),]
tax_heatmap$FMT2vsCON <- sapply(rownames(tax_heatmap), function(x) ifelse(x %in% rownames(tab_FMT2CON),"*","ns"))
tax_heatmap$FMT2vsFMT1 <- sapply(rownames(tax_heatmap), function(x) ifelse(x %in% rownames(tab_FMT2FMT1),"*","ns"))
tax_heatmap$FMT1vsCON <- sapply(rownames(tax_heatmap), function(x) ifelse(x %in% rownames(tab_FMT1CON),"*","ns"))

index <- match(rownames(tax_heatmap), rownames(tab_FMT2CON_all))
index
tax_heatmap$FMT2vsCON_p <- tab_FMT2CON_all$padj[index]


index <- match(rownames(tax_heatmap), rownames(tab_FMT2FMT1_all))
index
tax_heatmap$FMT2vsFMT1_p <- tab_FMT2FMT1_all$padj[index]

index <- match(rownames(tax_heatmap), rownames(tab_FMT1CON_all))
index
tax_heatmap$FMT1vsCON_p <- tab_FMT1CON_all$padj[index]

tax_heatmap <- tax_heatmap[order(tax_heatmap$FMT2vsCON_p,tax_heatmap$FMT2vsFMT1_p, tax_heatmap$FMT1vsCON_p),]

max(c(-log10(tax_heatmap$FMT2vsCON_p), -log10(tax_heatmap$FMT2vsFMT1_p),-log10(tax_heatmap$FMT1vsCON_p)))
mean(c(-log10(tax_heatmap$FMT2vsCON_p), -log10(tax_heatmap$FMT2vsFMT1_p), -log10(tax_heatmap$FMT1vsCON_p)))
pvalue_col_fun = colorRamp2(c(0, 6, 12), c("lightseagreen", "white", "red"))

my_palette <- c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey", "darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")
# adjust tax_heatmap genus and species, pasteurella
tax_heatmap$Genus <- as.character(tax_heatmap$Genus)
tax_heatmap$Species <- as.character(tax_heatmap$Species)
tax_heatmap$Genus[tax_heatmap$Species=='[Pasteurella] aerogenes'] <- 'Pasteurella'
tax_heatmap$Species[tax_heatmap$Species=='[Pasteurella] aerogenes'] <- 'Pasteurella aerogenes'

library(RColorBrewer)
genus <- unique(as.character(tax_heatmap$Genus))
genus_col <- colorRampPalette(my_palette)(length(genus))
names(genus_col) <- genus


ha_row <- HeatmapAnnotation(FMT2vsCON=anno_simple(-log10(tax_heatmap$FMT2vsCON_p),col = pvalue_col_fun,pch = na_if(tax_heatmap$FMT2vsCON,"ns")),
                            FMT2vsFMT1=anno_simple(-log10(tax_heatmap$FMT2vsFMT1_p),col = pvalue_col_fun,pch = na_if(tax_heatmap$FMT2vsFMT1,"ns")),
                            FMT1vsCON=anno_simple(-log10(tax_heatmap$FMT1vsCON_p),col = pvalue_col_fun,pch = na_if(tax_heatmap$FMT1vsCON,"ns")),
                            Genus=anno_simple(tax_heatmap$Genus, col = genus_col),
                            which = "row")

ha_row_txt <- rowAnnotation(labels = anno_text(tax_heatmap$Species, which = "row"))
meta_order$Litter <- factor(meta_order$Litter, levels = c("1","2","DONOR1","DONOR2"), labels = c("SOW1","SOW2","DONOR1","DONOR2"))

max(meta_order$NEC_severity)
mean(meta_order$NEC_severity)
min(meta_order$NEC_severity)

meta_order$NEC <- factor(meta_order$NEC, levels = c(0,1),labels = c("No","Yes"))

ha_col = HeatmapAnnotation(NEC=anno_simple(meta_order$NEC, col = c("Yes"="tomato","No"="lightseagreen")),
                           FMT=anno_simple(meta_order$Group,col = c("CON"="orange", "FMT1"="maroon1", "FMT2"="aquamarine3", "DONOR1"="red3", "DONOR2"="darkgreen")))

Hist <- Heatmap(mat_scaled[rownames(tax_heatmap),] , cluster_columns = FALSE, cluster_rows = TRUE,
                name="Z-score", col=colorRamp2(c(-2, 0, 2), c("dodgerblue4", "white","deeppink3")),
                top_annotation = ha_col, left_annotation = ha_row, right_annotation = ha_row_txt,
                show_row_names = FALSE, show_column_names = FALSE,
                heatmap_legend_param = list(legend_direction = "horizontal"))
Hist

# define the two legend
lgd_group = Legend(title = "FMT", legend_gp = gpar(fill = c("orange", "maroon1", "aquamarine3", "red3", "darkgreen")),
                 labels = c("CON", "FMT1", "FMT2", "DONOR1", "DONOR2"))

lgd_nec = Legend(title = "NEC", legend_gp = gpar(fill = c("tomato", "lightseagreen")),
                   labels = c("Yes", "No"))

lgd_genus = Legend(title = "Genus", legend_gp = gpar(fill = genus_col),labels = genus, ncol = 3)

lgd_sig = Legend(title= " ", pch = "*", type = "points", labels = "p.adj < 0.05")
lgd_pvalue = Legend(title = "Deseq2 p.adj value", col = pvalue_col_fun, at = c(0, 3, 6, 9, 12),
                    labels = c("1", "1e-3", "1e-6", "1e-9", "1e-12"),direction = "horizontal")


png("figure/figure_4.png",height = 8, width = 10, res = 300, units = "in")
draw(Hist, heatmap_legend_list=list(lgd_group,
                                    lgd_nec, lgd_genus, lgd_pvalue,lgd_sig), ht_gap =,
     heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

