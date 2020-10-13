# Figure.7B
################################################
# load data to calculate the matrix
ps <- readRDS("data/ps_denovo_mags_donorSource.rds")
ps.rel <- ps

metadata <- data.frame(sample_data(ps.rel),check.names = FALSE)

abun_matri <- data.frame(otu_table(ps.rel),check.names = FALSE)
abun_matri_t <- t(abun_matri)
library(reshape2)
rel_tab_l <- melt(abun_matri_t,variable.name = "mags")
colnames(rel_tab_l) <- c("sample","mags","abundance")
# match the Group and Litter, NEC, NEC_severity information
index <- match(rel_tab_l$sample, rownames(metadata))
table(is.na(index))
rel_tab_l$Group <- metadata$Group[index]


library(rstatix)

tab_FMT2CON_all <- rel_tab_l %>%
  filter(Group %in% c("FMT2", "CON") ) %>%
  group_by(mags) %>%
  kruskal_test(abundance~Group) %>%
  adjust_pvalue(method = "BH")%>%
  add_significance("p.adj") 

tab_FMT2FMT1_all <- rel_tab_l %>%
  filter(Group %in% c("FMT2", "FMT1") ) %>%
  group_by(mags) %>%
  wilcox_test(abundance ~ Group) %>%
  adjust_pvalue(method = "BH")%>%
  add_significance("p.adj") 

tab_FMT1CON_all <- rel_tab_l %>%
  filter(Group %in% c("FMT1", "CON") ) %>%
  group_by(mags) %>%
  wilcox_test(abundance ~ Group) %>%
  adjust_pvalue(method = "BH")%>%
  add_significance("p.adj") 
#Include wilcoxon p.adj FMT2 vs CON, FMT2 vs FMT1, taxonomy for rows
tab_FMT2CON <- subset(tab_FMT2CON_all, p.adj < 0.05)
tab_FMT2FMT1 <-subset(tab_FMT2FMT1_all, p.adj < 0.05)
tab_FMT1CON <- subset(tab_FMT1CON_all, p.adj < 0.05)
ASV <- unique(c(as.character(tab_FMT2CON$mags), as.character(tab_FMT2FMT1$mags),as.character(tab_FMT1CON$mags)))

ps.sample.rel <- readRDS("data/ps_denovo_mags_donorSource.rds")

#only take the donors belong mags
tax <- data.table::as.data.table(tax_table(ps.sample.rel),keep.rownames = TRUE)
index <- as.character(tax$source) != "Others"
tax_lac <- tax[index,]$rn 
ps.sample.rel.donor <- prune_taxa(tax_lac,ps.sample.rel)
ps.sample.rel.donor

#different by p different colonized in FMT1 and FMT2
ps.sample.rel.sig <- prune_taxa(as.character(tab_FMT2FMT1$mags),ps.sample.rel.donor)
ps.sample.rel.sig

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
tax_heatmap$FMT2vsCON <- sapply(rownames(tax_heatmap), function(x) ifelse(x %in% tab_FMT2CON$mags,"*","ns"))
tax_heatmap$FMT2vsFMT1 <- sapply(rownames(tax_heatmap), function(x) ifelse(x %in% tab_FMT2FMT1$mags,"*","ns"))
tax_heatmap$FMT1vsCON <- sapply(rownames(tax_heatmap), function(x) ifelse(x %in% tab_FMT1CON$mags,"*","ns"))

index <- match(rownames(tax_heatmap), tab_FMT2CON_all$mags)
index
tax_heatmap$FMT2vsCON_p <- tab_FMT2CON_all$p.adj[index]
tax_heatmap$FMT2vsCON_p[tax_heatmap$FMT2vsCON_p=="NaN"] <- 1

index <- match(rownames(tax_heatmap), tab_FMT2FMT1_all$mags)
index
tax_heatmap$FMT2vsFMT1_p <- tab_FMT2FMT1_all$p.adj[index]
tax_heatmap$FMT2vsFMT1_p[tax_heatmap$FMT2vsFMT1_p=="NaN"] <- 1

index <- match(rownames(tax_heatmap), tab_FMT1CON_all$mags)
index
tax_heatmap$FMT1vsCON_p <- tab_FMT1CON_all$p.adj[index]
tax_heatmap$FMT1vsCON_p[tax_heatmap$FMT1vsCON_p=="NaN"] <- 1

tax_heatmap <- tax_heatmap[order(tax_heatmap$FMT2vsCON_p,tax_heatmap$FMT2vsFMT1_p,tax_heatmap$FMT1vsCON_p),]

max(c(-log10(tax_heatmap$FMT2vsCON_p), -log10(tax_heatmap$FMT2vsFMT1_p),-log10(tax_heatmap$FMT1vsCON_p)))
mean(c(-log10(tax_heatmap$FMT2vsCON_p), -log10(tax_heatmap$FMT2vsFMT1_p),-log10(tax_heatmap$FMT1vsCON_p)))
pvalue_col_fun = colorRamp2(c(0, 1, 2), c("lightseagreen", "white", "red"))

my_palette <- c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey", "darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")
library(RColorBrewer)

tax_heatmap$species_adj <- stringr::str_replace(as.character(tax_heatmap$species), "\\ OX=[0-9]*\\ ","")
tax_heatmap$species_adj <- stringr::str_replace(as.character(tax_heatmap$species_adj), "\\ \\(strain\\ K12\\)","")
tax_heatmap$species_adj <- stringr::str_replace(as.character(tax_heatmap$species_adj), "\\ DSM\\ 5837\\ =\\ ATCC\\ 49540","")
tax_heatmap$species_adj <- stringr::str_replace(as.character(tax_heatmap$species_adj), "\\ CAG:511","")
tax_heatmap$species_adj <- stringr::str_replace(as.character(tax_heatmap$species_adj), "\\ CAG:511","")
tax_heatmap$species_adj <- stringr::str_replace(as.character(tax_heatmap$species_adj), "\\ \\(strain\\ ATCC\\ 53608\\)","")
tax_heatmap$species_adj <- stringr::str_replace(as.character(tax_heatmap$species_adj), "\\ DSM\\ 15176","")

species <- unique(as.character(tax_heatmap$species_adj))
species_col <- colorRampPalette(my_palette[1:12])(length(species))
names(species_col) <- species

magSource <- unique(as.character(tax_heatmap$source))
magSource_col <- colorRampPalette(my_palette[1:12])(length(magSource))
names(magSource_col) <- magSource

ha_row <- HeatmapAnnotation(FMT2vsCON=anno_simple(-log10(tax_heatmap$FMT2vsCON_p),col = pvalue_col_fun,pch = na_if(tax_heatmap$FMT2vsCON,"ns")),
                            FMT2vsFMT1=anno_simple(-log10(tax_heatmap$FMT2vsFMT1_p),col = pvalue_col_fun,pch = na_if(tax_heatmap$FMT2vsFMT1,"ns")),
                            FMT1vsCON=anno_simple(-log10(tax_heatmap$FMT1vsCON_p),col = pvalue_col_fun,pch = na_if(tax_heatmap$FMT1vsCON,"ns")),
                            Species=anno_simple(tax_heatmap$species_adj, col = species_col),
                            which = "row")

meta_order$Litter <- factor(meta_order$Litter, levels = c("1","2","DONOR1","DONOR2"), labels = c("SOW1","SOW2","DONOR1","DONOR2"))
max(meta_order$NEC_severity)
mean(meta_order$NEC_severity)
min(meta_order$NEC_severity)

ha_col = HeatmapAnnotation(FMT=anno_simple(meta_order$Group,col = c("CON"="orange", "FMT1"="maroon1", "FMT2"="aquamarine3", "DONOR1"="red3", "DONOR2"="darkgreen")))

Hist <- Heatmap(mat_scaled[rownames(tax_heatmap),] , cluster_columns = FALSE, cluster_rows = TRUE,
                name="Z-score", col=colorRamp2(c(-2, 0, 2), c("dodgerblue4", "white","deeppink3")),
                top_annotation = ha_col, left_annotation = ha_row,
                show_row_names = FALSE, show_column_names = FALSE,
                heatmap_legend_param = list(legend_direction = "horizontal"))
Hist

# define the two legend
lgd_group = Legend(title = "FMT", legend_gp = gpar(fill = c("orange", "maroon1", "aquamarine3", "red3", "darkgreen")),
                 labels = c("CON", "FMT1", "FMT2", "DONOR1", "DONOR2"))

lgd_genus = Legend(title = "Species", legend_gp = gpar(fill = species_col),labels = species, ncol = 2)

lgd_sig = Legend(title= " ", pch = "*", type = "points", labels = "p.adj < 0.05")
lgd_pvalue = Legend(title = "Wilcoxon p.adj value", col = pvalue_col_fun, at = c(0, 1, 2, 3),
                    labels = c("1", "0.1", "0.01", "0.001"),direction = "horizontal")

png("figure/figure_7B.png",height = 6, width = 12, units = "in", res = 100)
draw(Hist, heatmap_legend_list=list(lgd_group, lgd_genus, 
                                    lgd_pvalue,lgd_sig), 
     heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()


