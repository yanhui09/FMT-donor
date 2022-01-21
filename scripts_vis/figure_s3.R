library(ggplot2)
library(phyloseq)
library(extrafont)
library(vegan)
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

###########################################Kaiju phylum#####################################
ps.sample <- readRDS("data/ps_kaiju_nr_phylum.rds")
ps.sample.rel <- transform_sample_counts(ps.sample, function(x) x/sum(x)*100)

psdat.gen <- tax_glom(ps.sample.rel, taxrank = "Phylum")
ps.melt <- psmelt(psdat.gen)

# change to character for easy-adjusted level
ps.melt$Phylum <- as.character(ps.melt$Phylum)

ps.melt <- ps.melt %>%
  group_by(Group, Phylum) %>%
  mutate(mean=mean(Abundance))

# select group mean > 1
keep <- unique(ps.melt$Phylum[ps.melt$mean > 1])
ps.melt$Phylum[!(ps.melt$Phylum %in% keep)] <- "< 1% mean abund."
#to get the same rows together
ps.melt_sum <- ps.melt %>%
  group_by(Sample,Group,Phylum) %>%
  summarise(Abundance=sum(Abundance))


#join all qualitative palettes
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


# Scrren out mean relative abundance for each group
ps.melt_sum$Group <- factor(ps.melt_sum$Group, levels = c("CON","FMT1","FMT2","DONOR1","DONOR2"))
p_phylum_shotgun <- ggplot(ps.melt_sum, aes(x = Group, y = Abundance, color = Phylum, fill = Phylum)) + 
  geom_bar(stat = "summary", position = "stack",aes(fill=Phylum), fun="mean") + 
  labs(x="") +
  scale_fill_manual(values = col_vector) +
  scale_color_manual(values = col_vector) +
  theme_classic() + 
  mytheme + 
  theme(legend.position = "right") +
  guides(fill=guide_legend(ncol = 3))
p_phylum_shotgun

###################################16s genus#########################
ps.16s <- readRDS("data/phyloseq_16s.rds")
# 16s genera
ps.16s.rel <- transform_sample_counts(ps.16s, function(x) x/sum(x)*100)
psdat.gen <- tax_glom(ps.16s.rel, taxrank = "Genus")
ps.melt <- psmelt(psdat.gen)
# change to charater for easy adjust level
ps.melt$Genus <- as.character(ps.melt$Genus)
ps.melt <- ps.melt %>%
  group_by(Group, Genus) %>%
  mutate(mean=mean(Abundance))
# select group mean > 1 genera
keep <- unique(ps.melt$Genus[ps.melt$mean > 1])
ps.melt$Genus[!(ps.melt$Genus %in% keep)] <- "< 1% mean abund."
ps.melt_sum_16s <- ps.melt %>%
  group_by(Sample,Group,Genus) %>%
  summarise(Abundance=sum(Abundance))
#####################################Shotgun genus######################
ps.sample <- readRDS("data/ps_kaiju_nr_genus.rds")
ps.shotgun <- ps.sample
ps.shotgun.rel <- transform_sample_counts(ps.shotgun, function(x) x/sum(x)*100)
colSums(otu_table(ps.shotgun.rel))
psdat.gen <- tax_glom(ps.shotgun.rel, taxrank = "Genus")
ps.melt <- psmelt(psdat.gen)
#write.table(ps.melt,"data/kaiju_nr_all_melt_unrare_g.tsv", sep = "\t", col.names = NA)
#ps.melt <- read.table("data/kaiju_nr_all_melt_unrare_g.tsv", sep = "\t", header = T, row.names = 1)
# change to charater for easy adjust level
ps.melt$Genus <- as.character(ps.melt$Genus)

library(dplyr)
ps.melt <- ps.melt %>%
  group_by(Group, Genus) %>%
  mutate(mean=mean(Abundance))

# select group mean > 1 Genus
keep <- unique(ps.melt$Genus[ps.melt$mean > 1])
ps.melt$Genus[!(ps.melt$Genus %in% keep)] <- "< 1% mean abund."
#to get the same rows together
ps.melt_sum_shotgun <- ps.melt %>%
  group_by(Sample,Group,Genus) %>%
  summarise(Abundance=sum(Abundance))
#######################################find shared genus for coloring
genus_16s <- as.character(unique(ps.melt_sum_16s$Genus))
genus_shotgun <- as.character(unique(ps.melt_sum_shotgun$Genus))
genus_intersect <- intersect(genus_16s,genus_shotgun)
#reset levels
ps.melt_sum_16s$Genus <- factor(ps.melt_sum_16s$Genus, levels = c(genus_intersect,genus_16s[! genus_16s %in% genus_intersect]))
##############plot 16s genus #############################
# Scrren out mean relative abundance for each group
ps.melt_sum_16s$Group <- factor(ps.melt_sum_16s$Group, levels = c("CON","FMT1","FMT2","DONOR1","DONOR2"))
p_genus_16s <- ggplot(ps.melt_sum_16s, aes(x = Group, y = Abundance, color = Genus, fill = Genus)) + 
  geom_bar(stat = "summary", position = "stack",aes(fill=Genus), fun="mean") + 
  labs(x="") +
  scale_fill_manual(values = col_vector) +
  scale_color_manual(values = col_vector) +
  theme_classic() + 
  mytheme + 
  theme(legend.position = "right") +
  guides(fill=guide_legend(ncol = 3))
p_genus_16s  

###############plot shotgun genus########################################
ps.melt_sum_shotgun$Genus <- factor(ps.melt_sum_shotgun$Genus, levels = c(genus_intersect,genus_shotgun[! genus_shotgun %in% genus_intersect]))

# Scrren out mean relative abundance for each group
ps.melt_sum_shotgun$Group <- factor(ps.melt_sum_shotgun$Group, levels = c("CON","FMT1","FMT2","DONOR1","DONOR2"))
p_genus_shotgun <- ggplot(ps.melt_sum_shotgun, aes(x = Group, y = Abundance, color = Genus, fill = Genus)) + 
  geom_bar(stat = "summary", position = "stack",aes(fill=Genus), fun="mean") + 
  labs(x="") +
  scale_fill_manual(values = col_vector) +
  scale_color_manual(values = col_vector) +
  theme_classic() + 
  mytheme + 
  theme(legend.position = "right") +
  guides(fill=guide_legend(ncol = 3))
p_genus_shotgun

#########multi-panel figure using patchwork
library(patchwork)
p_merged_2 <- (p_phylum_shotgun / p_genus_16s / p_genus_shotgun) & theme(legend.justification = "left")
p_merged_2 + plot_annotation(tag_levels = 'a') 
ggsave("figure/figure_s3.png", height = 8, width = 11)

sessionInfo()