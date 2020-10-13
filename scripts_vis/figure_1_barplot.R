library(phyloseq)
library(dplyr)
#select proteobacteria/enterobacterial, clostridium, lactobacillus, top 5 -10 bacteria
#collect 16s data
ps_16s <- readRDS("data/phyloseq_16s.rds")
ps_16s
# build a function to exact the > 5% mean abundance genus
abun_tab <- function(ps_16s){
ps.rel_16s <- transform_sample_counts(ps_16s, function(x) x/sum(x)*100)

psdat.gen_16s <- tax_glom(ps.rel_16s, taxrank = "Genus")
ps.melt_16s <- psmelt(psdat.gen_16s)
ps.melt_16s$Genus <- as.character(ps.melt_16s$Genus)
ps.melt_16s <- ps.melt_16s %>%
  group_by(Group,Genus) %>%
  mutate(mean=mean(Abundance))
# select group mean > 1 speciese
#keep <- unique(ps.melt_16s$Genus[ps.melt_16s$mean > 1])
#table(keep)
#ps.melt_16s$Genus[!(ps.melt_16s$Genus %in% keep)] <- "< 1% mean abund."

#group the bacteria according to our manual selection
#Bacteriodetes: Bacteroides
#Proterobacteria: Enterobacteriaceae
#Firmicutes: Lactobacillus, Streptococcus, Subdoligranulum, Enterococcus, Clostridium, Staphylococcus
keep_genus <- c("Bacteroides","Lactobacillus","Streptococcus","Subdoligranulum", "Enterococcus", "Clostridium", "Staphylococcus")
keep_family <- "Enterobacteriaceae"
ps.melt_16s$Genus[!(ps.melt_16s$Genus %in% keep_genus)] <- "Residual"
ps.melt_16s$Genus[(ps.melt_16s$Family %in% keep_family)] <- "Enterobacteriaceae"
#to get the same rows together
ps.melt_sum_16s <- ps.melt_16s %>%
  group_by(Sample, Group,Genus) %>%
  summarise(Abundance=sum(Abundance))
return(ps.melt_sum_16s)
}
ps.group_sum_16s <- abun_tab(ps_16s)

ps.group_sum_16s$Genus
#####################################################
#collect metagenomics data
ps_meta <- readRDS("data/ps_kaiju_nr_genus.rds")
ps.group_sum_meta <- abun_tab(ps_meta)
#ps.group_sum_meta <- ps.group_mean_meta

ps.group_sum_meta$platform <- "Shotgun"
ps.group_sum_16s$platform <- "16s"
ps.group_sum <- rbind.data.frame(ps.group_sum_16s, ps.group_sum_meta)
write.table(ps.group_sum,"data/figure_1_barplot.tsv", sep = "\t", col.names = NA)
#########################################
#draw with ggplot2
ps.melt.sum <- read.table("data/figure_1_barplot.tsv", sep = "\t", header = T, row.names = 1)
library(ggplot2)
library(extrafont)
loadfonts()
mytheme<- theme(plot.title = element_text(hjust=0.5, family = "Arial", size=14),
                #legend.title = element_blank(),
                #legend.position ="right",
                legend.text = element_text(family = "Arial", size = 12),
                legend.title = element_text(family = "Arial", size = 14),
                strip.background = element_blank(),
                strip.placement = "outside",
                strip.text = element_text(family = "Arial", size = 12, vjust = 1),
                axis.text = element_text(family = 'Arial', size = 12, color = "black"),
                axis.text.x.bottom = element_text(angle = - 45),
                axis.title = element_text(family = "Arial", size = 14),
                axis.ticks.x = element_blank(),
                panel.grid = element_blank())
library(RColorBrewer)
levels(ps.melt.sum$Genus)
# adjust the Genus label level so as to to reorder legend labels
ps.melt.sum$Genus <- factor(ps.melt.sum$Genus, levels = c("Residual", "Bacteroides", "Enterobacteriaceae", "Lactobacillus", "Streptococcus","Subdoligranulum", "Clostridium", "Staphylococcus", "Enterococcus"))
cols = rev(colorRampPalette(brewer.pal(10, "Paired"))(length(levels(ps.melt.sum$Genus))-1))
cols =c("lightgray",cols)

ps.melt.sum$Group <- factor(ps.melt.sum$Group, levels = c("CON","FMT1","FMT2","DONOR1","DONOR2"), labels = c("CON", "FMT1", "FMT2", "DONOR1", "DONOR2"))
ps.melt.sum$platform <- factor(ps.melt.sum$platform, levels = c("16s", "Shotgun"), labels = c("16S", "Shotgun"))
p <- ggplot(ps.melt.sum, aes(x = Group, y = Abundance, color = Genus, fill = Genus)) + 
  geom_bar(stat = "summary", position = "stack", fun=mean) + 
  labs(x="", y="Relative abundance (%)", title = "",color="", fill="") +
  facet_wrap(.~ platform) + 
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  theme_classic() + 
  mytheme
p

ggsave("figure/figure_1_barplot.png", height = 4.5, width = 8)

