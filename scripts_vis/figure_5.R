# 16s and Kaiju whole genome sequencing
library(ggplot2)
library(phyloseq)
library(extrafont)
library(RColorBrewer)
library(dplyr)
loadfonts()
mytheme<- theme(plot.title = element_text(hjust=0.5, family = "Arial", size=12),
                legend.position ="right",
                legend.text = element_text(family = "Arial", size = 10),
                legend.background = element_blank(),
                strip.background = element_blank(),
                strip.placement = "outside",
                strip.text = element_text(family = "Arial", size = 10),
                axis.text = element_text(family = 'Arial', size = 10, color="black"),
                panel.grid = element_blank())
my_palette <- c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey", "darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(5)

# 16s
#load data
ps <- readRDS("data/phyloseq_16s.rds")
ps

# Figure.5A 16s
metadata <- read.table("data/Metadata_16s.tsv", sep = "\t", header = T, row.names = 1)
# # generated to the Lactobacillus genus abundance
ps.rel <- transform_sample_counts(ps, function(x) x/sum(x)*100)
# only keep Lactobacillus
tax <- data.table::as.data.table(tax_table(ps.rel),keep.rownames = TRUE)
tax_lac <- tax[Genus=="Lactobacilli",]$rn
ps.rel_lac <- prune_taxa(tax_lac,ps.rel)

psdat.gen_lac <- tax_glom(ps.rel_lac, taxrank = "Genus")
#
ps.melt <- psmelt(psdat.gen_lac)
ps.melt_f <- subset(ps.melt, Litter!="Donor")

library(ggpubr)
library(rstatix)

stat.test <- ps.melt_f %>%
  group_by(Genus) %>%
  wilcox_test(Abundance~Group) %>%
  adjust_pvalue(method = "BH")%>%
  add_significance(p.col = "p") %>%
  filter(p.adj<0.05) %>%
  mutate(y.position=seq(60,70,length.out = 3))
stat.test

ps.melt_f$NEC <- factor(ps.melt_f$NEC, levels = c(0,1), labels = c("No", "Yes"))
p_16s_lac_box <-ggplot(data = ps.melt_f, mapping = aes(x=Group, y=Abundance)) +
  stat_boxplot(geom ='errorbar', linetype=1, width=0.5) + 
  geom_boxplot(outlier.shape = NA, aes(fill=Group)) +
  geom_point(position = position_jitter(w=0.1),size=2, aes(color=NEC)) +
  scale_color_manual(values = c("black", "gray")) +
  scale_fill_manual(values = cols)+
  labs(x="", y="Lactobacilli (%)", title = "16S rRNA",color="NEC") +
  stat_pvalue_manual(stat.test,label = "p.adj.signif", tip.length = 0, size = 6)+ 
  theme_classic() +
  mytheme

p_16s_lac_box
#######################################################
# Figure.5A shotgun
metadata <- read.table("data/Metadata.tsv", sep = "\t", header = T, row.names = 1)
ps <- readRDS("data/ps_kaiju_nr_genus.rds")
# generated to the Lactobacillus genus abundance
ps.rel <- transform_sample_counts(ps, function(x) x/sum(x)*100)

# only keep Lactobacillus
tax <- data.table::as.data.table(tax_table(ps.rel),keep.rownames = TRUE)
tax_lac <- tax[Genus=="Lactobacilli",]$rn
ps.rel_lac <- prune_taxa(tax_lac,ps.rel)

psdat.gen_lac <- tax_glom(ps.rel_lac, taxrank = "Genus")
ps.melt <- psmelt(psdat.gen_lac)
ps.melt_f <- subset(ps.melt, Litter!="DONOR1"&Litter!="DONOR2")
# match nec severity
index <- match(ps.melt_f$FMT_ID, metadata$FMT_ID)
index
ps.melt_f$NEC_severity <- metadata$NEC_severity[index]

library(rstatix)
ps.melt_f$NEC <- factor(ps.melt_f$NEC, levels = c(0,1), labels = c("No","Yes"))

stat.test <- ps.melt_f %>%
  group_by(Genus) %>%
  wilcox_test(Abundance~Group) %>%
  adjust_pvalue(method = "BH")%>%
  add_significance(p.col = "p") %>%
  filter(p.adj<0.05) %>%
  mutate(y.position=seq(65,75,length.out = 3))
stat.test

p_kaiju_lac_box <-ggplot(data = ps.melt_f, mapping = aes(x=Group, y=Abundance)) +
  stat_boxplot(geom ='errorbar', linetype=1, width=0.5) + 
  geom_boxplot(outlier.shape = NA, aes(fill=Group)) +
  geom_point(position = position_jitter(w=0.1),size=2, aes(color=NEC)) +
  scale_color_manual(values = c("black", "gray")) +
  scale_fill_manual(values = cols)+
  labs(x="", y="Lactobacilli (%)", title = "Shotgun",color="NEC") +
  stat_pvalue_manual(stat.test,label = "p.adj.signif", tip.length = 0, size = 6)+ 
  theme_classic() +
  mytheme
p_kaiju_lac_box
#####################################
# Figure 5B
#load data
ps <- readRDS("data/ps_kaiju_nr_species.rds")
ps
# screen out the 3 group by the relative abundance
ps.rel <- transform_sample_counts(ps, function(x) x/sum(x)*100)
# only keep Lactobacillus
# filter first then cluster --- faster, data.table faster
tax <- data.table::as.data.table(tax_table(ps.rel),keep.rownames = TRUE)
tax_lac <- tax[Genus=="Lactobacilli",]$rn 

ps.rel_lac <- prune_taxa(tax_lac,ps.rel)

psdat.gen_lac <- tax_glom(ps.rel_lac, taxrank = "Species")

ps.melt <- psmelt(psdat.gen_lac)
# change to charater for easy adjust level
ps.melt$Species <- as.character(ps.melt$Species)
ps.melt <- ps.melt %>%
  group_by(Group,Species) %>%
  mutate(mean=mean(Abundance))

# select group mean > 1 speciese
keep <- unique(ps.melt$Species[ps.melt$mean > 1])
ps.melt$Species[!(ps.melt$Species %in% keep)] <- "< 1% mean abund."
#to get the same rows together
ps.melt_sum <- ps.melt %>%
  group_by(Sample,Group,Species) %>%
  summarise(Abundance=sum(Abundance))
summary(ps.melt)
#join all qualitative palettes
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
set.seed(1)
my_palette_large <- colorRampPalette(col_vector)(length(unique(ps.melt_sum$Species)))

# Scrren out mean relative abundance for each group
ps.melt_sum$Group <- factor(ps.melt_sum$Group, levels = c("CON","FMT1","FMT2","DONOR1","DONOR2"))
p_species_kaiju_lac <- ggplot(ps.melt_sum, aes(x = Group, y = Abundance, color = Species, fill = Species)) + 
  geom_bar(stat = "summary", position = "stack",aes(fill=Species), fun="mean") + 
  labs(x="", y="Relative abundance (%)", title = "Shotgun") +
  scale_fill_manual(values = my_palette_large) +
  scale_color_manual(values = my_palette_large) +
  theme_classic() + 
  mytheme + 
  theme(legend.position = "right") +
  guides(fill=guide_legend(ncol = 1))
p_species_kaiju_lac

########################################
# Figure.5C lre lcr tree
library(ggtree)
# build a function to draw tree file
plot_tree <- function(tree, meta, tit){
# read newick tree
e.sir.tre <- read.tree(tree)
# create ggtree object
e.sir.gg <- ggtree( e.sir.tre )
# read in metadata file
e.sir.meta <- read.delim(meta, header = T, sep = "\t" )
# add metadata to dendrogram plot
e.sir.gg <- e.sir.gg %<+% e.sir.meta

# strainphlan_tree_1.pdf
p <- e.sir.gg +
  geom_tippoint( size = 3, aes( color = Group ) ) +
  aes( branch.length = 'length' ) +
  labs(title = tit) +
  theme_tree2() + theme(legend.position="right") +
  mytheme
return(p)
}

tree_lre <- plot_tree(tree = "data/strainphlan/RAxML_bestTree.s__Lactobacillus_reuteri.tree",
                      meta = "data/strainphlan/mapping.tsv",
                      tit = "Limosilactobacillus reuteri")

tree_lcr <- plot_tree(tree = "data/strainphlan/RAxML_bestTree.s__Lactobacillus_crispatus.tree",
                      meta = "data/strainphlan/mapping.tsv",
                      tit = "Lactobacillus crispatus")

######################################################
# multi-panel figure
library(patchwork)

p_1 <- p_16s_lac_box + labs(tag = "A") + p_kaiju_lac_box + theme(plot.tag.position = "topleft") +
  plot_layout(guides = "collect")

p_1 <- p_1*theme(axis.text.x.bottom = element_text(angle = 90))

p_2 <- p_species_kaiju_lac*theme(axis.text.x.bottom = element_text(angle = 90)) + 
  labs(tag = "B") + theme(plot.tag.position = "topleft")

p_3 <-  tree_lre + labs(tag = "C") + tree_lcr + theme(plot.tag.position = "topleft") + 
  plot_layout(guides = "collect")

(p_1|p_2)/
  p_3

ggsave("figure/figure_5.png", height = 6, width = 10)

