# Figure.7A load data
table  <- read.table("data/irep_merged_sa99.tab", sep = "\t", header = T, row.names = 1, check.names = FALSE)
table[table=="n/a"] <- NA
rowSums(is.na(table))
table <- table[rowSums(is.na(table))!=40,]

mag_tax <- read.table("data/mag_taxaSource.tsv", sep = "\t", header = T, row.names = 1)

#import the taxonomy
library(tidyr)
library(stringr)
tax.clean <- subset(mag_tax, select = c(Superkingdom, phylum,class, order, order, family, genus.1))
mag_tax <- tax.clean[rownames(table),]

mapping <- read.table("data/Metadata.tsv", sep = "\t", header = T, row.names = 1)

### average bacteria replication
library(tidyverse)
table_long <- table %>% 
  mutate(mags=rownames(table)) %>%
  gather(sample, irep, -mags)

table_long <- na.omit(table_long)
table_long$irep <- as.numeric(table_long$irep)
index <- match(table_long$sample, rownames(mapping))
table_long$Group <- mapping$Group[index]
  
table_long_mean <- table_long %>%
  group_by(sample) %>%
  summarise(mean=mean(irep))
index <- match(table_long_mean$sample, rownames(mapping))
table_long_mean$Group <- as.character(mapping$Group[index])

table_long_mean$NEC <- mapping$NEC[index]
table_long_mean$NEC <- as.factor(table_long_mean$NEC)
library(rstatix)
library(ggpubr)
library(extrafont)
library(RColorBrewer)
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

# Lactobacillus irep
names <- rownames(mag_tax)[mag_tax$genus.1=="Lactobacillus"]
table_long_enterbac <- table_long[table_long$mags %in% names,]
table_long_mean <- table_long_enterbac %>%
  group_by(sample) %>%
  summarise(mean=mean(irep))
index <- match(table_long_mean$sample, rownames(mapping))
table_long_mean$Group <- as.character(mapping$Group[index])

table_long_mean$NEC <- mapping$NEC[index]
table_long_mean$NEC <- factor(table_long_mean$NEC, levels = c(0,1), labels = c("No", "Yes"))

stat.test <- table_long_mean %>%
  wilcox_test(mean~Group) %>%
  adjust_pvalue(method = "BH")%>%
  add_significance(p.col = "p")  %>%
  filter(p.adj<0.05) %>%
  mutate(y.position=seq(3,3,length.out = 1))
stat.test

# lactobacillus irep
table_long_mean$Group <- factor(table_long_mean$Group, levels = c("CON","FMT1","FMT2", "DONOR1","DONOR2"))

table_long_mean <- table_long_mean %>%
  filter(Group!="DONOR1" & Group!="DONOR2")


p_irep_lac_mean <-ggplot(data = table_long_mean, mapping = aes(x=Group, y=mean)) +
  stat_boxplot(geom ='errorbar', linetype=1, width=0.5) + 
  geom_boxplot(outlier.shape = NA, aes(fill=Group)) +
  geom_point(position = position_jitter(w=0.1),size=2, aes(color=NEC)) +
  labs(x="", y="iRep",color="NEC",title = "Lactobacilli") +
  scale_color_manual(values = c("black", "gray")) +
  scale_fill_manual(values = c(cols[1], cols[2], cols[3])) +
  stat_pvalue_manual(stat.test,label = "p.adj.signif", tip.length = 0, size = 6)+ 
  theme_classic() +
  mytheme
p_irep_lac_mean
ggsave("figure/figure_7A.png",p_irep_lac_mean, height = 4.5, width = 4)
  