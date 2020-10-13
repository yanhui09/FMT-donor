library(ggplot2)
library(phyloseq)
library(extrafont)
library(vegan)
loadfonts()
mytheme<- theme(plot.title = element_text(hjust=0.5, family = "Arial", size=12),
                #legend.title = element_blank(),
                legend.position ="right",
                legend.text = element_text(family = "Arial", size = 10),
                legend.background = element_blank(),
                strip.background = element_blank(),
                strip.placement = "outside",
                strip.text = element_text(family = "Arial", size = 10),
                #axis.title = element_text(family = 'Arial', size = 8),
                axis.text = element_text(family = 'Arial', size = 10, color="black"),
                panel.grid = element_blank())
my_palette <- c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey", "darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")

ps.sample <- readRDS("data/ps_kaiju_nr_phylum.rds")

#

#quick analysis with phyloseq
# rarefy without replacement
 ps.rarefied <- ps.sample
# ps.rarefied.rel <- transform_sample_counts(ps.rarefied, function(x) x/sum(x)*100)
# colSums(otu_table(ps.rarefied.rel))
# psdat.gen <- tax_glom(ps.rarefied.rel, taxrank = "Phylum")
# ps.melt <- psmelt(psdat.gen)
# write.table(ps.melt,"data/kaiju_nr_all_melt_unrare_p.tsv", sep = "\t", col.names = NA)
ps.melt <- read.table("data/kaiju_nr_all_melt_unrare_p.tsv", sep = "\t", header = T, row.names = 1)
# change to charater for easy adjust level
ps.melt$Phylum <- as.character(ps.melt$Phylum)
head(ps.melt)
library(dplyr)
ps.melt <- ps.melt %>%
  group_by(Group, Phylum) %>%
  mutate(mean=mean(Abundance))

#ps.melt$Species[ps.melt$mean < 1] <- "< 1% mean abund."
# select group mean > 1 Genus
keep <- unique(ps.melt$Phylum[ps.melt$mean > 1])
ps.melt$Phylum[!(ps.melt$Phylum %in% keep)] <- "< 1% mean abund."
#to get the same rows together
ps.melt_sum <- ps.melt %>%
  group_by(Sample,Group,Phylum) %>%
  summarise(Abundance=sum(Abundance))


summary(ps.melt)
# move below 1% to others
library(RColorBrewer)
#join all qualitative palettes
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
set.seed(1)
my_palette_large <- sample(col_vector, length(unique(ps.melt_sum$Phylum)))

# Scrren out mean relative abundance for each group
ps.melt_sum$Group <- factor(ps.melt_sum$Group, levels = c("CON","FMT1","FMT2","DONOR1","DONOR2"))
p_phylum_bac_all <- ggplot(ps.melt_sum, aes(x = Group, y = Abundance, color = Phylum, fill = Phylum)) + 
  geom_bar(stat = "summary", position = "stack",aes(fill=Phylum), fun.y="mean") + 
  #facet_wrap(.~ Group, ncol = 5, scales = "free") + 
  labs(x="") +
  scale_fill_manual(values = my_palette_large) +
  scale_color_manual(values = my_palette_large) +
  theme_classic() + 
  mytheme + 
  theme(legend.position = "right") +
  guides(fill=guide_legend(ncol = 2))
p_phylum_bac_all
# ################### Domain
ps.rarefied.rel <- transform_sample_counts(ps.rarefied, function(x) x/sum(x)*100)
psdat.gen <- tax_glom(ps.rarefied.rel, taxrank = "Domain")
ps.melt <- psmelt(psdat.gen)
# change to charater for easy adjust level
ps.melt$Domain <- as.character(ps.melt$Domain)
library(dplyr)
ps.melt <- ps.melt %>%
  group_by(Domain) %>%
  mutate(mean=mean(Abundance))
#to get the same rows together
ps.melt_sum <- ps.melt %>%
  group_by(Sample,Group,Domain) %>%
  summarise(Abundance=sum(Abundance))
#ps.melt$Domain[ps.melt$mean < 1] <- "< 1% mean abund."
summary(ps.melt)
# move below 1% to others
library(RColorBrewer)
#join all qualitative palettes
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
set.seed(1)
my_palette_large <- sample(col_vector, length(unique(ps.melt$Domain)))
ps.melt$Group <- factor(ps.melt$Group, levels = c("CON","FMT1","FMT2","DONOR1","DONOR2"))
p_phylum_bac_all_domain <- ggplot(ps.melt, aes(x = Group, y = Abundance, color = Domain, fill = Domain)) +
  geom_bar(stat = "summary", position = "stack",aes(fill=Domain), fun.y="mean") +
  #facet_wrap(.~ Group, ncol = 5, scales = "free") +
  labs(x="") +
  scale_fill_manual(values = my_palette_large) +
  scale_color_manual(values = my_palette_large) +
  theme_classic() +
  mytheme +
  theme(legend.position = "right") +
  guides(fill=guide_legend(ncol = 1))
p_phylum_bac_all_domain
ggsave("figure/paper/nr_euk_domain_p.pdf", height= 4 ,width = 6)

# ###########################################################
# # Bacteria
# ## screen out bacteria and virus separately
# 
# ps.sample <- readRDS("data/ps_kaiju_nr.rds")
# min(colSums(otu_table(ps.sample)))
# # remove taxa
# tax_virus <- rownames(tax.clean[tax.clean$Domain == "Bacteria",])
# ps.virus <- prune_taxa(rownames(otu_table(ps.sample)) %in% tax_virus, ps.sample)
# ps.virus
# min(colSums(otu_table(ps.virus)))
# #quick analysis with phyloseq
# # rarefy without replacement
# ps.rarefied = rarefy_even_depth(ps.virus, rngseed=1, sample.size=1000000, replace=F)
# ps.rarefied
# # calculate alpha diversity 
# rich <- estimate_richness(ps.rarefied)
# rownames(rich) <- rownames(sample_data(ps.rarefied))
# tab <- subset(rich, select = c("Shannon"))
# tab
# index <- match(rownames(rich), rownames(sample_data(ps.rarefied)))
# tab$Group <-sample_data(ps.rarefied)$Group[index]
# tab$Litter <- sample_data(ps.rarefied)$Litter[index]
# tab$Sex <- sample_data(ps.rarefied)$Sex[index]
# tab$NEC <- sample_data(ps.rarefied)$NEC[index]
# #tab$Pattern <- sample_data(ps.rarefied)$Pattern[index]
# tab$Sample <- rownames(rich)
# 
# library(reshape2)
# tab_long <- melt(tab,id.vars = c("Sample","Group", "Litter", "Sex","NEC"),
#                  variable.name = "Parameter", value.name = "Abundance")
# tab_long
# tab_long$Litter <- factor(tab_long$Litter)
# tab_long$NEC <- factor(tab_long$NEC)
# library(tidyr)
# library(rstatix)
# library(ggpubr)
# stat.test <- tab_long %>%
#   filter(Litter!="DONOR2" & Litter!="DONOR1") %>%
#   group_by(Parameter) %>%
#   wilcox_test(Abundance~Group) %>%
#   adjust_pvalue(method = "BH")%>%
#   add_significance(p.col = "p") %>%
#   filter(p.adj < 0.05) %>%
#   mutate(y.position=c(seq(7,8,length.out = 2)))
# #stat.test
# #write.table(stat.test, file = "table/shannon_stat_bacteria.tsv", sep = "\t", col.names = NA)
# tab_long$NEC <- factor(tab_long$NEC)
# p_alpha_bac <- ggplot(tab_long, aes(x= Group, y= Abundance)) +
#   #geom_violin(position = "dodge",
#   #            alpha=1, outlier.size=0, size=0.7, width=0.5, color="gray", fill="gray") +
#   geom_crossbar(stat="summary", fun.y=mean, fun.ymax=mean, fun.ymin=mean, fatten=1, width=.5)+
#   geom_point(aes(color=NEC),position = position_jitter(w=0.05),size=2) +
#   scale_color_manual(values = c("black", "gray")) +
#   labs(x="", y="Shannon index", title="Kaiju Bacteria") +
#   stat_pvalue_manual(stat.test,label = "p.adj.signif", tip.length = 0, size = 4)+ 
#   #facet_wrap(.~Parameter, ncol = 1, scales = "free") +
#   theme_classic() +
#   mytheme 
# p_alpha_bac
# #dir.create("figure/diversity")
# #ggsave("figure/diversity/alpha_bacteria.pdf", p, height = 5, width = 6)
# ###### beta diversity
# # calculate distance
# disWuf <- distance(ps.rarefied, method = "bray")
# #ordWuf_MDS <- ordinate(ps.rarefied, method = "PCA", distance = disWuf)
# #plot_scree(ordWuf_MDS, "Scree Plot: Weighted UniFrac MDS")
# #ordWuf_NMDS <- ordinate(ps.rarefied, method = "NMDS", distance = disWuf)
# #plot_scree(ordWuf_NMDS, "Scree Plot: Weighted UniFrac NMDS")
# #plot_ordination(ps.rarefied, ordWuf_NMDS) + 
# #  geom_point(mapping = aes(color=Group)) +
# #  ggtitle("NMDS: Weighted Unifrac") + 
# #facet_wrap(.~ Design, nrow = 2, scales = "free") +
# #  theme_classic() +
# #  mytheme +
# #  scale_color_manual(values = my_palette)
# file_save="Bray Curtis"
# # pca analysis based on distance matrix
# sub_design <- data.frame(sample_data(ps.rarefied))
# disWuf <- as.matrix(disWuf)
# disWuf <- disWuf[rownames(sub_design),rownames(sub_design)]
# #disWuf <- as.dist(disWuf)
# pcoa = cmdscale(disWuf, k=3, eig=T)
# # k is dimension, 3 is recommended; eig is eigenvalues 
# points = as.data.frame(pcoa$points)
# # get coordinate string, format to dataframme 
# colnames(points) = c("x", "y", "z")
# eig = pcoa$eig
# points = cbind(points, sub_design[match(rownames(points), rownames(sub_design)),])
# #calculate adonis P
# sub_design_f <- subset(sub_design, Group!="DONOR1"&Group!="DONOR2")
# sub_design_f$Group <- factor(sub_design_f$Group)
# disWuf_f <- disWuf[rownames(sub_design_f),rownames(sub_design_f)]
# adonis_sub <- adonis(disWuf_f~Group, data=sub_design_f, permutations = 999)
# adonis_R2 <- adonis_sub$aov.tab$R2
# head(adonis_sub$aov.tab)
# adonis_p <- adonis_sub$aov.tab$`Pr(>F)`
# 
# #my_palette_small <- sample(my_palette,length(levels(points$Group)))
# # draw two dimensional pcoa analysis  ##change group
# points$NEC <- factor(points$NEC)
# p_bray_bacteria = ggplot(points, aes(x=x, y=y)) +
#   geom_point(alpha=.8, size=2,aes(color=Group)) +
#   # stat_ellipse(aes(fill=Group),geom="polygon",level=0.8, alpha=0.2, type="t") +
#   labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
#        y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
#        color="FMT")+
#   ggtitle(paste(file_save, " PCoA (Kaiju Bacteria)\n",
#                 "Group effect, R2=", format(adonis_R2[1],digits = 2),", P=",format(adonis_p[1], digits = 1),"\n")) +
#   #facet_wrap(.~ Design, nrow = 2, scales = "free") +
#   # scale_color_manual(values = my_palette_small)+
#   #guides(color=guide_legend(title = "AgeGroup"), shape=guide_legend(title = "Probiotics"), size=guide_legend(title = "NEC")) +
#   theme_classic()+
#   mytheme
# p_bray_bacteria
# # assess the effect size
# # utilized vegan to assess the effect size
# sub_design <- subset(sub_design, Litter!='DONOR1' & Litter!='DONOR2')
# tab_adonis <- function(disWuf,x) {
#   disWuf <- disWuf[rownames(sub_design),rownames(sub_design)]
#   adonis_sub <- adonis(as.formula(paste(quote(disWuf),"~", x)), data=sub_design[rownames(disWuf),], permutations = 999)
#   adonis_R2 <- adonis_sub$aov.tab$R2
#   adonis_sub$aov.tab
#   tab <- subset(adonis_sub$aov.tab, select = c("R2", "Pr(>F)"))
#   tab <- data.frame(tab)
#   tab$Covariate <- rownames(tab)
#   colnames(tab) <- c("R2", "P", "Covariate")
#   tab <- tab[1,]
#   return(tab)
# }
# sub_design$NEC <- factor(sub_design$NEC)
# tab_1 <- tab_adonis(disWuf, "Group") 
# tab_2 <- tab_adonis(disWuf, "Litter")
# tab_3 <- tab_adonis(disWuf, "Sex")
# tab_4 <- tab_adonis(disWuf, "NEC")
# tab <- rbind.data.frame(tab_1, tab_2, tab_3, tab_4)
# tab
# tab_f <- tab[tab$P < 0.05,]
# #tab_f <- tab_f[order(tab_f$R2,decreasing = T),]
# tab_f$Covariate <- factor(tab_f$Covariate, levels = tab_f$Covariate[order(tab_f$R2)])
# tab_f
# tab_f$Covariate <- factor(tab_f$Covariate, levels = c("Litter","Group"),labels=c("SOW","FMT"))
# p_effectsize_bray_bac<- tab_f %>%
#   ggplot(aes(y=R2*100, x=Covariate)) +
#   geom_bar(stat='identity', fill="darkgray") +
#   coord_flip() +
#   labs(title = "Clinical covarites associated with GM (Kaiju Bacteria)", y= "Effect size(R2, %)", x="") +
#   theme_classic()+
#   mytheme
# p_effectsize_bray_bac
# # ggsave(paste0("figure/diversity/PCoA", "_", file_save,"_kaiju_bacteria.pdf"), p, width = 6, height = 4)
# # ###############binary jaccard
# # # calculate distance
# # disWuf <- distance(ps.rarefied, method = "jaccard", binary=T)
# # #ordWuf_MDS <- ordinate(ps.rarefied, method = "PCA", distance = disWuf)
# # #plot_scree(ordWuf_MDS, "Scree Plot: Weighted UniFrac MDS")
# # #ordWuf_NMDS <- ordinate(ps.rarefied, method = "NMDS", distance = disWuf)
# # #plot_scree(ordWuf_NMDS, "Scree Plot: Weighted UniFrac NMDS")
# # #plot_ordination(ps.rarefied, ordWuf_NMDS) + 
# # #  geom_point(mapping = aes(color=Group)) +
# # #  ggtitle("NMDS: Weighted Unifrac") + 
# # #facet_wrap(.~ Design, nrow = 2, scales = "free") +
# # #  theme_classic() +
# # #  mytheme +
# # #  scale_color_manual(values = my_palette)
# # file_save="Binary Jaccard"
# # # pca analysis based on distance matrix
# # sub_design <- data.frame(sample_data(ps.rarefied))
# # disWuf <- as.matrix(disWuf)
# # disWuf <- disWuf[rownames(sub_design),rownames(sub_design)]
# # #disWuf <- as.dist(disWuf)
# # pcoa = cmdscale(disWuf, k=3, eig=T)
# # # k is dimension, 3 is recommended; eig is eigenvalues 
# # points = as.data.frame(pcoa$points)
# # # get coordinate string, format to dataframme 
# # colnames(points) = c("x", "y", "z")
# # eig = pcoa$eig
# # points = cbind(points, sub_design[match(rownames(points), rownames(sub_design)),])
# # #calculate adonis P
# # sub_design_f <- subset(sub_design, Group!="DONOR1"&Group!="DONOR2")
# # sub_design_f$Group <- factor(sub_design_f$Group)
# # disWuf_f <- disWuf[rownames(sub_design_f),rownames(sub_design_f)]
# # adonis_sub <- adonis(disWuf_f~Group, data=sub_design_f, permutations = 999)
# # adonis_R2 <- adonis_sub$aov.tab$R2
# # head(adonis_sub$aov.tab)
# # adonis_p <- adonis_sub$aov.tab$`Pr(>F)`
# # 
# # #my_palette_small <- sample(my_palette,length(levels(points$Group)))
# # # draw two dimensional pcoa analysis  ##change group
# # points$NEC <- factor(points$NEC)
# # p = ggplot(points, aes(x=x, y=y)) +
# #   geom_point(alpha=.8, size=2,aes(color=Group,shape=NEC)) +
# #   # stat_ellipse(aes(fill=Group),geom="polygon",level=0.8, alpha=0.2, type="t") +
# #   labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
# #        y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
# #        shape="NEC")+
# #   ggtitle(paste(file_save, " PCoA\n",
# #                 "Group effect, R2=", format(adonis_R2[1],digits = 2),", P=",format(adonis_p[1], digits = 1),"\n")) +
# #   #facet_wrap(.~ Design, nrow = 2, scales = "free") +
# #   # scale_color_manual(values = my_palette_small)+
# #   #guides(color=guide_legend(title = "AgeGroup"), shape=guide_legend(title = "Probiotics"), size=guide_legend(title = "NEC")) +
# #   theme_classic()+
# #   mytheme
# # p
# # ggsave(paste0("figure/diversity/PCoA", "_", file_save,"_bacteria.pdf"), p, width = 6, height = 4)
# #ps.rarefied.rel <- transform_sample_counts(ps.rarefied, function(x) x/sum(x)*100)
# #psdat.gen <- tax_glom(ps.rarefied.rel, taxrank = "Species")
# #ps.melt <- psmelt(psdat.gen)
# #write.table(ps.melt,"data/kaiju_nr_bacteria_melt.tsv", sep = "\t", col.names = NA)
# ps.melt <- read.table("data/kaiju_nr_bacteria_melt.tsv", sep = "\t", header = T, row.names = 1)
# # change to charater for easy adjust level
# ps.melt$Species <- as.character(ps.melt$Species)
# ps.melt <- ps.melt %>%
#   group_by(Species) %>%
#   mutate(mean=mean(Abundance))
# 
# ps.melt$Species[ps.melt$mean < 1] <- "< 1% mean abund."
# summary(ps.melt)
# # move below 1% to others
# library(RColorBrewer)
# #join all qualitative palettes
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# set.seed(1)
# my_palette_large <- sample(col_vector, length(unique(ps.melt$Species)))
# 
# # Scrren out mean relative abundance for each group
# p_species_bac <- ggplot(ps.melt, aes(x = Group, y = Abundance, color = Species, fill = Species)) + 
#   geom_bar(stat = "summary", position = "fill",aes(fill=Species), fun.y="mean") + 
#   #facet_wrap(.~ Group, ncol = 5, scales = "free") + 
#   labs(x="") +
#   scale_fill_manual(values = my_palette_large) +
#   scale_color_manual(values = my_palette_large) +
#   theme_classic() + 
#   mytheme + 
#   theme(legend.position = "right") +
#   guides(fill=guide_legend(ncol = 4))
# p_species_bac
# 
# ###################################################################
# # Virus
# ## screen out bacteria and virus separately
# #virus
# ps.sample <- readRDS("data/ps_kaiju_nr.rds")
# min(colSums(otu_table(ps.sample)))
# # remove taxa
# tax_virus <- rownames(tax.clean[tax.clean$Domain == "Viruses",])
# ps.virus <- prune_taxa(rownames(otu_table(ps.sample)) %in% tax_virus, ps.sample)
# ps.virus
# colSums(otu_table(ps.virus))
# min(colSums(otu_table(ps.virus)))
# #quick analysis with phyloseq
# # rarefy without replacement
# #ps.rarefied = rarefy_even_depth(ps.virus, rngseed=1, sample.size=900, replace=F)
# ps.rarefied <- ps.virus
# 
# # calculate alpha diversity 
# rich <- estimate_richness(ps.rarefied)
# rownames(rich) <- rownames(sample_data(ps.rarefied))
# tab <- subset(rich, select = c("Shannon"))
# tab
# index <- match(rownames(rich), rownames(sample_data(ps.rarefied)))
# tab$Group <-sample_data(ps.rarefied)$Group[index]
# tab$Litter <- sample_data(ps.rarefied)$Litter[index]
# tab$Sex <- sample_data(ps.rarefied)$Sex[index]
# tab$NEC <- sample_data(ps.rarefied)$NEC[index]
# #tab$Pattern <- sample_data(ps.rarefied)$Pattern[index]
# tab$Sample <- rownames(rich)
# 
# library(reshape2)
# tab_long <- melt(tab,id.vars = c("Sample","Group", "Litter", "Sex","NEC"),
#                  variable.name = "Parameter", value.name = "Abundance")
# tab_long
# tab_long$Litter <- factor(tab_long$Litter)
# 
# library(tidyr)
# library(rstatix)
# library(ggpubr)
# stat.test <- tab_long %>%
#   filter(Litter!='DONOR1'& Litter!='DONOR2') %>%
#   group_by(Parameter) %>%
#   wilcox_test(Abundance~Group) %>%
#   adjust_pvalue(method = "BH")%>%
#   add_significance(p.col = "p") %>%
#   filter(p.adj < 0.05) %>%
#   mutate(y.position=c(6.5))
# stat.test
# #write.table(stat.test, file = "table/shannon_stat_virus.tsv", sep = "\t", col.names = NA)
# tab_long$NEC <- factor(tab_long$NEC)
# p_alpha_vir <- ggplot(tab_long, aes(x= Group, y= Abundance)) +
#   #geom_violin(position = "dodge",
#   #            alpha=1, outlier.size=0, size=0.7, width=0.5, color="gray", fill="gray") +
#   geom_crossbar(stat="summary", fun.y=mean, fun.ymax=mean, fun.ymin=mean, fatten=1, width=.5)+
#   geom_point(aes(color=NEC),position = position_jitter(w=0.05),size=2) +
#   scale_color_manual(values = c("black", "gray")) +
#   labs(x="", y="Shannon Index", title="Kaiju Virus") +
#   stat_pvalue_manual(stat.test,label = "p.signif", tip.length = 0, size = 4)+ 
#   #facet_wrap(.~Parameter, ncol = 1, scales = "free") +
#   theme_classic() +
#   mytheme 
# p_alpha_vir
# 
# # dir.create("figure/diversity")
# # ggsave("figure/diversity/alpha_virus.pdf", p, height = 5, width = 6)
# # virus level
# #######################################################
# # screen out the 3 group by the relative abundance
# # genera
# tax_table(ps.rarefied)
# ps.rarefied.rel <- transform_sample_counts(ps.rarefied, function(x) x/sum(x)*100)
# psdat.gen <- tax_glom(ps.rarefied.rel, taxrank = "Genus")
# ps.melt <- psmelt(psdat.gen)
# # change to charater for easy adjust level
# ps.melt$Genus <- as.character(ps.melt$Genus)
# ps.melt <- ps.melt %>%
#   group_by(Genus) %>%
#   mutate(mean=mean(Abundance))
# ps.melt$Genus[ps.melt$mean < 1] <- "< 1% mean abund."
# summary(ps.melt)
# # move below 1% to others
# library(RColorBrewer)
# #join all qualitative palettes
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# set.seed(1)
# my_palette_large <- colorRampPalette(col_vector)(length(unique(ps.melt$Genus)))
# 
# # Scrren out mean relative abundance for each group
# p_species_vir <- ggplot(ps.melt, aes(x = Group, y = Abundance, color = Genus, fill = Genus)) + 
#   geom_bar(stat = "summary", position = "fill",aes(fill=Genus), fun.y="mean") + 
#   labs(x="") +
#   #facet_wrap(.~ Group, ncol = 5, scales = "free") + 
#   scale_fill_manual(values = my_palette_large) +
#   scale_color_manual(values = my_palette_large) +
#   theme_classic() + 
#   mytheme + 
#   theme(legend.position = "bottom") +
#   guides(fill=guide_legend(ncol = 7))
# p_species_vir 
# 
# ###### beta diversity
# # calculate distance
# disWuf <- distance(ps.rarefied, method = "bray")
# #ordWuf_MDS <- ordinate(ps.rarefied, method = "PCA", distance = disWuf)
# #plot_scree(ordWuf_MDS, "Scree Plot: Weighted UniFrac MDS")
# #ordWuf_NMDS <- ordinate(ps.rarefied, method = "NMDS", distance = disWuf)
# #plot_scree(ordWuf_NMDS, "Scree Plot: Weighted UniFrac NMDS")
# #plot_ordination(ps.rarefied, ordWuf_NMDS) + 
# #  geom_point(mapping = aes(color=Group)) +
# #  ggtitle("NMDS: Weighted Unifrac") + 
# #facet_wrap(.~ Design, nrow = 2, scales = "free") +
# #  theme_classic() +
# #  mytheme +
# #  scale_color_manual(values = my_palette)
# file_save="Bray Curtis"
# # pca analysis based on distance matrix
# sub_design <- data.frame(sample_data(ps.rarefied))
# disWuf <- as.matrix(disWuf)
# disWuf <- disWuf[rownames(sub_design),rownames(sub_design)]
# #disWuf <- as.dist(disWuf)
# pcoa = cmdscale(disWuf, k=3, eig=T)
# # k is dimension, 3 is recommended; eig is eigenvalues 
# points = as.data.frame(pcoa$points)
# # get coordinate string, format to dataframme 
# colnames(points) = c("x", "y", "z")
# eig = pcoa$eig
# points = cbind(points, sub_design[match(rownames(points), rownames(sub_design)),])
# #calculate adonis P
# sub_design_f <- subset(sub_design, Group!="DONOR1"&Group!="DONOR2")
# sub_design_f$Group <- factor(sub_design_f$Group)
# disWuf_f <- disWuf[rownames(sub_design_f),rownames(sub_design_f)]
# adonis_sub <- adonis(disWuf_f~Group, data=sub_design_f, permutations = 999)
# adonis_R2 <- adonis_sub$aov.tab$R2
# head(adonis_sub$aov.tab)
# adonis_p <- adonis_sub$aov.tab$`Pr(>F)`
# 
# #my_palette_small <- sample(my_palette,length(levels(points$Group)))
# # draw two dimensional pcoa analysis  ##change group
# points$NEC <- factor(points$NEC)
# p_bray_vir = ggplot(points, aes(x=x, y=y)) +
#   geom_point(alpha=.8, size=2,aes(color=Group)) +
#   # stat_ellipse(aes(fill=Group),geom="polygon",level=0.8, alpha=0.2, type="t") +
#   labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
#        y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
#        color="FMT")+
#   ggtitle(paste(file_save, " PCoA (Kaiju Virus)\n",
#                 "Group effect, R2=", format(adonis_R2[1],digits = 2),", P=",format(adonis_p[1], digits = 1),"\n")) +
#   #facet_wrap(.~ Design, nrow = 2, scales = "free") +
#   # scale_color_manual(values = my_palette_small)+
#   #guides(color=guide_legend(title = "AgeGroup"), shape=guide_legend(title = "Probiotics"), size=guide_legend(title = "NEC")) +
#   theme_classic()+
#   mytheme
# p_bray_vir
# # ggsave(paste0("figure/diversity/PCoA", "_", file_save,"_kaiju_virus.pdf"), p, width = 6, height = 4)
# sub_design <- subset(sub_design, Litter!='DONOR1' & Litter!='DONOR2')
# tab_adonis <- function(disWuf,x) {
#   disWuf <- disWuf[rownames(sub_design),rownames(sub_design)]
#   adonis_sub <- adonis(as.formula(paste(quote(disWuf),"~", x)), data=sub_design[rownames(disWuf),], permutations = 999)
#   adonis_R2 <- adonis_sub$aov.tab$R2
#   adonis_sub$aov.tab
#   tab <- subset(adonis_sub$aov.tab, select = c("R2", "Pr(>F)"))
#   tab <- data.frame(tab)
#   tab$Covariate <- rownames(tab)
#   colnames(tab) <- c("R2", "P", "Covariate")
#   tab <- tab[1,]
#   return(tab)
# }
# sub_design$NEC <- factor(sub_design$NEC)
# tab_1 <- tab_adonis(disWuf, "Group")
# tab_2 <- tab_adonis(disWuf, "Litter")
# tab_3 <- tab_adonis(disWuf, "Sex")
# tab_4 <- tab_adonis(disWuf, "NEC")
# tab <- rbind.data.frame(tab_1, tab_2, tab_3, tab_4)
# tab
# tab_f <- tab[tab$P < 0.05,]
# #tab_f <- tab_f[order(tab_f$R2,decreasing = T),]
# tab_f$Covariate <- factor(tab_f$Covariate, levels = tab_f$Covariate[order(tab_f$R2)])
# tab_f$Covariate <- factor(tab_f$Covariate, levels = c("Litter","Group"),labels=c("SOW","FMT"))
# p_effectsize_bray_vir<- tab_f %>%
#   ggplot(aes(y=R2*100, x=Covariate)) +
#   geom_bar(stat='identity', fill="darkgray") +
#   coord_flip() +
#   labs(title = "Clinical covarites associated with GM (Kaiju Virus)", y= "Effect size(R2, %)", x="") +
#   theme_classic()+
#   mytheme
# p_effectsize_bray_vir
# 
# # ###############binary jaccard
# # # calculate distance
# # disWuf <- distance(ps.rarefied, method = "jaccard", binary=T)
# # #ordWuf_MDS <- ordinate(ps.rarefied, method = "PCA", distance = disWuf)
# # #plot_scree(ordWuf_MDS, "Scree Plot: Weighted UniFrac MDS")
# # #ordWuf_NMDS <- ordinate(ps.rarefied, method = "NMDS", distance = disWuf)
# # #plot_scree(ordWuf_NMDS, "Scree Plot: Weighted UniFrac NMDS")
# # #plot_ordination(ps.rarefied, ordWuf_NMDS) + 
# # #  geom_point(mapping = aes(color=Group)) +
# # #  ggtitle("NMDS: Weighted Unifrac") + 
# # #facet_wrap(.~ Design, nrow = 2, scales = "free") +
# # #  theme_classic() +
# # #  mytheme +
# # #  scale_color_manual(values = my_palette)
# # file_save="Binary Jaccard"
# # # pca analysis based on distance matrix
# # sub_design <- data.frame(sample_data(ps.rarefied))
# # disWuf <- as.matrix(disWuf)
# # disWuf <- disWuf[rownames(sub_design),rownames(sub_design)]
# # #disWuf <- as.dist(disWuf)
# # pcoa = cmdscale(disWuf, k=3, eig=T)
# # # k is dimension, 3 is recommended; eig is eigenvalues 
# # points = as.data.frame(pcoa$points)
# # # get coordinate string, format to dataframme 
# # colnames(points) = c("x", "y", "z")
# # eig = pcoa$eig
# # points = cbind(points, sub_design[match(rownames(points), rownames(sub_design)),])
# # #calculate adonis P
# # sub_design_f <- subset(sub_design, Group!="DONOR1"&Group!="DONOR2")
# # sub_design_f$Group <- factor(sub_design_f$Group)
# # disWuf_f <- disWuf[rownames(sub_design_f),rownames(sub_design_f)]
# # adonis_sub <- adonis(disWuf_f~Group, data=sub_design_f, permutations = 999)
# # adonis_R2 <- adonis_sub$aov.tab$R2
# # head(adonis_sub$aov.tab)
# # adonis_p <- adonis_sub$aov.tab$`Pr(>F)`
# # points$NEC <- factor(points$NEC)
# # #my_palette_small <- sample(my_palette,length(levels(points$Group)))
# # # draw two dimensional pcoa analysis  ##change group
# # p = ggplot(points, aes(x=x, y=y)) +
# #   geom_point(alpha=.8, size=2,aes(color=Group,shape=NEC)) +
# #   # stat_ellipse(aes(fill=Group),geom="polygon",level=0.8, alpha=0.2, type="t") +
# #   labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
# #        y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""))+
# #   ggtitle(paste(file_save, " PCoA\n",
# #                 "Group effect, R2=", format(adonis_R2[1],digits = 2),", P=",format(adonis_p[1], digits = 1),"\n")) +
# #   #facet_wrap(.~ Design, nrow = 2, scales = "free") +
# #   # scale_color_manual(values = my_palette_small)+
# #   #guides(color=guide_legend(title = "AgeGroup"), shape=guide_legend(title = "Probiotics"), size=guide_legend(title = "NEC")) +
# #   theme_classic()+
# #   mytheme
# # p
# # ggsave(paste0("figure/diversity/PCoA", "_", file_save,"_virus.pdf"), p, width = 6, height = 4)
# # screen out the 3 group by the relative abundance
# # genera
# ps.rarefied.rel <- transform_sample_counts(ps.rarefied, function(x) x/sum(x)*100)
# psdat.gen <- tax_glom(ps.rarefied.rel, taxrank = "Genus")
# ps.melt <- psmelt(psdat.gen)
# # change to charater for easy adjust level
# ps.melt$Genus <- as.character(ps.melt$Genus)
# ps.melt <- ps.melt %>%
#   group_by(Genus) %>%
#   mutate(mean=mean(Abundance))
# ps.melt$Genus[ps.melt$mean < 1] <- "< 1% mean abund."
# summary(ps.melt)
# # move below 1% to others
# library(RColorBrewer)
# #join all qualitative palettes
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# set.seed(1)
# my_palette_large <- colorRampPalette(col_vector)(length(unique(ps.melt$Genus)))
# 
# # Scrren out mean relative abundance for each group
# p_genus_vir <- ggplot(ps.melt, aes(x = Group, y = Abundance, color = Genus, fill = Genus)) + 
#   geom_bar(stat = "summary", position = "fill",aes(fill=Genus), fun.y="mean") + 
#   #facet_wrap(.~ Group, ncol = 5, scales = "free") + 
#   scale_fill_manual(values = my_palette_large) +
#   scale_color_manual(values = my_palette_large) +
#   theme_classic() + 
#   mytheme + 
#   theme(legend.position = "right") +
#   guides(fill=guide_legend(ncol = 3))
# p_genus_vir 
# 
# library(patchwork)
# q <- p_alpha_vir / p_bray_vir / p_genus_vir +
#   plot_annotation(tag_levels = 'A')
# q
# ggsave("figure/paper/supple_virus_community.pdf", height = 8, width = 10)
