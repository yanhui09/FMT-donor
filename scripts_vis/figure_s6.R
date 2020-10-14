#Figure.s6 functional capacities heatmap
################################################
# load data to calculate the matrix
ps <- readRDS("data/ps_gf_keggmodule.rds")
# rm unclassified
keep <- rownames(otu_table(ps))[rownames(otu_table(ps))!="Unclassified"]
ps <- prune_taxa(keep, ps)
ps
ps.rel <- ps

# use deseq2 to identify different modules # not work use wilcox
metadata <- data.frame(sample_data(ps.rel),check.names = FALSE)
abun_matri <- data.frame(otu_table(ps.rel),check.names = FALSE)

desep2east <- function(x,y){
metadata_f <- subset(metadata,Group %in% c(x,y))
metadata_f$Group <- factor(metadata_f$Group)
abun_matri_f <- abun_matri[,rownames(metadata_f)]
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = abun_matri_f, colData = metadata_f, design= ~ Group)
head(dds)
# normalize original dds
dds2 <- DESeq(dds)
# check resultname
resultsNames(dds2)
res <- results(dds2)
# get results
return(res)
}

tab_FMT2CON_all <- data.frame(desep2east("FMT2","CON"))
tab_FMT2FMT1_all <- data.frame(desep2east("FMT2","FMT1"))
tab_FMT1CON_all <- data.frame(desep2east("FMT1","CON"))
tab_FMT2CON <- subset(tab_FMT2CON_all, padj < 0.001)
tab_FMT2FMT1 <-subset(tab_FMT2FMT1_all, padj < 0.001)
tab_FMT1CON <- subset(tab_FMT1CON_all, padj < 0.001)

tab_FMT2CON_005 <- subset(tab_FMT2CON_all, padj < 0.001)
tab_FMT2FMT1_005 <-subset(tab_FMT2FMT1_all, padj < 0.001)
tab_FMT1CON_005 <- subset(tab_FMT1CON_all, padj < 0.001)
ASV <- unique(c(as.character(rownames(tab_FMT2CON)), as.character(rownames(tab_FMT2FMT1)),as.character(rownames(tab_FMT1CON))))
# only select mean abundance > 1


ps.sample <- readRDS("data/ps_gf_keggmodule.rds")

min(colSums(otu_table(ps.sample)))
#use relative abundance considering the varied depth
ps.sample.rel <- transform_sample_counts(ps.sample, function(x) x/sum(x)*100)
tax <- as.data.frame(tax_table(ps.sample.rel),keep.rownames = TRUE)

ps.sample.rel.sig <- prune_taxa(rownames(otu_table(ps.sample.rel)) %in% ASV,ps.sample.rel)
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
tax_heatmap <- tax.clean[rownames(mat_scaled),, drop=FALSE]
tax_heatmap$FMT2vsCON <- sapply(rownames(tax_heatmap), function(x) ifelse(x %in% rownames(tab_FMT2CON_005),"*","ns"))
tax_heatmap$FMT2vsFMT1 <- sapply(rownames(tax_heatmap), function(x) ifelse(x %in% rownames(tab_FMT2FMT1_005),"*","ns"))
tax_heatmap$FMT1vsCON <- sapply(rownames(tax_heatmap), function(x) ifelse(x %in% rownames(tab_FMT1CON_005),"*","ns"))

index <- match(rownames(tax_heatmap), rownames(tab_FMT2CON_all))
index
tax_heatmap$FMT2vsCON_p <- tab_FMT2CON_all$padj[index]
tax_heatmap$FMT2vsCON_p[is.na(tax_heatmap$FMT2vsCON_p)] <- 1

index <- match(rownames(tax_heatmap), rownames(tab_FMT2FMT1_all))
index
tax_heatmap$FMT2vsFMT1_p <- tab_FMT2FMT1_all$padj[index]
tax_heatmap$FMT2vsFMT1_p[is.na(tax_heatmap$FMT2vsFMT1_p)] <- 1

index <- match(rownames(tax_heatmap), rownames(tab_FMT1CON_all))
index
tax_heatmap$FMT1vsCON_p <- tab_FMT1CON_all$padj[index]
tax_heatmap$FMT1vsCON_p[is.na(tax_heatmap$FMT1vsCON_p)] <- 1

tax_heatmap <- tax_heatmap[order(tax_heatmap$FMT2vsCON_p,tax_heatmap$FMT2vsFMT1_p,tax_heatmap$FMT1vsCON_p),]

max(c(-log10(tax_heatmap$FMT2vsCON_p), -log10(tax_heatmap$FMT2vsFMT1_p),-log10(tax_heatmap$FMT1vsCON_p)))
mean(c(-log10(tax_heatmap$FMT2vsCON_p), -log10(tax_heatmap$FMT2vsFMT1_p),-log10(tax_heatmap$FMT1vsCON_p)))
pvalue_col_fun = colorRamp2(c(0, 6, 12), c("lightseagreen", "white", "red"))

my_palette <- c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey", "darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")
library(RColorBrewer)

ha_row <- HeatmapAnnotation(FMT2vsCON=anno_simple(-log10(tax_heatmap$FMT2vsCON_p),col = pvalue_col_fun,pch = na_if(tax_heatmap$FMT2vsCON,"ns")),
                            FMT2vsFMT1=anno_simple(-log10(tax_heatmap$FMT2vsFMT1_p),col = pvalue_col_fun,pch = na_if(tax_heatmap$FMT2vsFMT1,"ns")),
                            FMT1vsCON=anno_simple(-log10(tax_heatmap$FMT1vsCON_p),col = pvalue_col_fun,pch = na_if(tax_heatmap$FMT1vsCON,"ns")),
                            which = "row")

ha_row_txt <- rowAnnotation(labels = anno_text(tax_heatmap$module, which = "row"))
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
lgd_sig = Legend(title= " ", pch = "*", type = "points", labels = "p.adj < 0.05")
lgd_pvalue = Legend(title = "Deseq2 p.adj value", col = pvalue_col_fun, at = c(0, 3, 6, 9, 12),
                    labels = c("1", "1e-3", "1e-6", "1e-9", "1e-12"),direction = "horizontal")

png("figure/figure_s6.png",height = 8, width = 12, units = "in", res = 100)
draw(Hist, heatmap_legend_list=list(lgd_group, lgd_nec, lgd_pvalue,lgd_sig),
     heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()


