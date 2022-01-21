library(ggplot2)
library(phyloseq)
library(extrafont)
library(vegan)
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
ps.sample <- readRDS("data/ps_kaiju_nr_species.rds")
###################################################
# Figure.3B Shannon shotgun metagenomics
# calculate alpha diversity
rich <- estimate_richness(ps.sample)
rownames(rich) <- rownames(sample_data(ps.sample))
tab <- subset(rich, select = c("Shannon"))
tab
index <- match(rownames(rich), rownames(sample_data(ps.sample)))
tab$Group <-sample_data(ps.sample)$Group[index]
tab$Litter <- sample_data(ps.sample)$Litter[index]
tab$Sex <- sample_data(ps.sample)$Sex[index]
tab$NEC <- sample_data(ps.sample)$NEC[index]
tab$Sample <- rownames(rich)
#
library(reshape2)
tab_long <- melt(tab,id.vars = c("Sample","Group", "Litter", "Sex","NEC"),
                 variable.name = "Parameter", value.name = "Abundance")

tab_long$Litter <- factor(tab_long$Litter)
tab_long$NEC <- factor(tab_long$NEC)
library(tidyr)
library(rstatix)
library(ggpubr)
stat.test <- tab_long %>%
  filter(Group!="DONOR1"&Group!="DONOR2") %>%
  group_by(Parameter) %>%
  wilcox_test(Abundance~Group) %>%
  adjust_pvalue(method = "BH")%>%
  add_significance(p.col = "p") %>%
  mutate(y.position=c(seq(4.5,5.5,length.out = 3 )))
stat.test

stat.test <- stat.test %>%
  filter(p.adj<0.05)

tab_long$Group <- factor(tab_long$Group, levels = c("CON","FMT1","FMT2","DONOR1","DONOR2"))
tab_long$NEC <- factor(tab_long$NEC, levels = c(0,1), labels = c("No","Yes"))
tab_long$Group <- factor(tab_long$Group,levels = c("CON","FMT1","FMT2", "DONOR1", "DONOR2"))
#color palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(5)

p_alpha_all <- ggplot(tab_long, aes(x= Group, y= Abundance)) +
  stat_boxplot(geom ='errorbar', linetype=1, width=0.5) + 
  geom_boxplot(outlier.shape = NA, aes(fill=Group)) +
  geom_point(position = position_jitter(w=0.1, h=0),size=2, aes(color=NEC)) +
  scale_color_manual(values = c("black", "gray")) +
  scale_fill_manual(values = cols) +
  labs(x="", y="Shannon index", title="Shotgun") +
  stat_pvalue_manual(stat.test,label = "p.adj.signif", tip.length = 0, size = 6)+
  theme_classic() +
  mytheme
p_alpha_all

#Figure.3D bray curtis PCoA Shotgun
#Hellinger transformation
ps.rarefied.dis <- transform_sample_counts(ps.sample, function(x) sqrt(x/sum(x)*100))
# calculate distance
disWuf <- distance(ps.rarefied.dis, method = "bray")
file_save="Bray Curtis"
# pca analysis based on distance matrix
sub_design <- data.frame(sample_data(ps.sample))
disWuf <- as.matrix(disWuf)
disWuf <- disWuf[rownames(sub_design),rownames(sub_design)]
#disWuf <- as.dist(disWuf)
pcoa = cmdscale(disWuf, k=3, eig=T)
# k is dimension, 3 is recommended; eig is eigenvalues
points = as.data.frame(pcoa$points)
# get coordinate string, format to dataframme
colnames(points) = c("x", "y", "z")
eig = pcoa$eig
points = cbind(points, sub_design[match(rownames(points), rownames(sub_design)),])
#calculate adonis P

sub_design_f <- subset(sub_design, Group!="DONOR1"&Group!="DONOR2")
sub_design_f$Group <- factor(sub_design_f$Group)
disWuf_f <- disWuf[rownames(sub_design_f),rownames(sub_design_f)]
adonis_sub <- adonis(disWuf_f~Group, data=sub_design_f, permutations = 999)
adonis_R2 <- adonis_sub$aov.tab$R2
head(adonis_sub$aov.tab)
adonis_p <- adonis_sub$aov.tab$`Pr(>F)`

# draw two dimensional pcoa analysis  ##change group
points$NEC <- factor(points$NEC)
points$Group <- factor(points$Group, levels = c("CON","FMT1","FMT2","DONOR1","DONOR2"))
p_bray_all = ggplot(points, aes(x=x, y=y)) +
  geom_point(alpha=.8, size=2,aes(color=Group)) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       color="FMT")+
  ggtitle(paste(file_save, " PCoA (Shotgun)\n",
                "Group effect, R2=", format(adonis_R2[1],digits = 2),", P=",format(adonis_p[1], digits = 1),"\n")) +
  theme_classic()+
  mytheme
p_bray_all
# use adonis to test effect size

sub_design <- subset(sub_design, Litter!='DONOR1' & Litter!='DONOR2')
tab_adonis <- function(disWuf,x) {
  disWuf <- disWuf[rownames(sub_design),rownames(sub_design)]
  adonis_sub <- adonis(as.formula(paste(quote(disWuf),"~", x)), data=sub_design[rownames(disWuf),], permutations = 999)
  adonis_R2 <- adonis_sub$aov.tab$R2
  adonis_sub$aov.tab
  tab <- subset(adonis_sub$aov.tab, select = c("R2", "Pr(>F)"))
  tab <- data.frame(tab)
  tab$Covariate <- rownames(tab)
  colnames(tab) <- c("R2", "P", "Covariate")
  tab <- tab[1,]
  return(tab)
}
sub_design$NEC <- factor(sub_design$NEC)
tab_1 <- tab_adonis(disWuf, "Group")
tab_2 <- tab_adonis(disWuf, "Litter")
tab_3 <- tab_adonis(disWuf, "Sex")
tab_4 <- tab_adonis(disWuf, "NEC")
tab <- rbind.data.frame(tab_1, tab_2, tab_3, tab_4)
tab
tab$Covariate <- c("FMT", "SOW", "Sex", "NEC")
tab_f <- tab[tab$P < 0.05,]

tab_f$Covariate <- factor(tab_f$Covariate, levels = tab_f$Covariate[order(tab_f$R2)])

p_effectsize_bray_all <- tab_f %>%
  ggplot(aes(y=R2*100, x=Covariate)) +
  geom_bar(stat='identity', fill="darkgray") +
  coord_flip() +
  labs(title = "Clinical covarites associated with GM (Shotgun)", y= "Effect size(R2, %)", x="") +
  theme_classic()+
  mytheme
p_effectsize_bray_all

#Figure.S5C and S5D binary jaccard PCoA shotgun
# calculate distance
disWuf <- distance(ps.rarefied.dis, method = "jaccard", binary=T)
file_save="Binary Jaccard"

# pca analysis based on distance matrix
sub_design <- data.frame(sample_data(ps.sample))
disWuf <- as.matrix(disWuf)
disWuf <- disWuf[rownames(sub_design),rownames(sub_design)]
#disWuf <- as.dist(disWuf)
pcoa = cmdscale(disWuf, k=3, eig=T)
# k is dimension, 3 is recommended; eig is eigenvalues
points = as.data.frame(pcoa$points)
# get coordinate string, format to dataframme
colnames(points) = c("x", "y", "z")
eig = pcoa$eig
points = cbind(points, sub_design[match(rownames(points), rownames(sub_design)),])

#calculate adonis P
sub_design_f <- subset(sub_design, Group!="DONOR1"&Group!="DONOR2")
sub_design_f$Group <- factor(sub_design_f$Group)
disWuf_f <- disWuf[rownames(sub_design_f),rownames(sub_design_f)]
adonis_sub <- adonis(disWuf_f~Group, data=sub_design_f, permutations = 999)
adonis_R2 <- adonis_sub$aov.tab$R2
head(adonis_sub$aov.tab)
adonis_p <- adonis_sub$aov.tab$`Pr(>F)`

points$NEC <- factor(points$NEC, levels = c(0,1),labels = c("No","Yes"))
points$Group <- factor(points$Group, levels = c("CON","FMT1","FMT2","DONOR1","DONOR2"))
p_jaccard = ggplot(points, aes(x=x, y=y)) +
  geom_point(alpha=.8, size=2,aes(color=Group, shape=NEC)) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       shape="NEC")+
  ggtitle(paste(file_save, " PCoA (Shotgun)\n",
                "Group effect, R2=", format(adonis_R2[1],digits = 2),", P=",format(adonis_p[1], digits = 1),"\n")) +
  guides(color=guide_legend(title = "FMT")) +
  theme_classic()+
  mytheme
p_jaccard
# effect size
sub_design <- subset(sub_design, Litter!='DONOR1' & Litter!='DONOR2')
tab_adonis <- function(disWuf,x) {
  disWuf <- disWuf[rownames(sub_design),rownames(sub_design)]
  adonis_sub <- adonis(as.formula(paste(quote(disWuf),"~", x)), data=sub_design[rownames(disWuf),], permutations = 999)
  adonis_R2 <- adonis_sub$aov.tab$R2
  adonis_sub$aov.tab
  tab <- subset(adonis_sub$aov.tab, select = c("R2", "Pr(>F)"))
  tab <- data.frame(tab)
  tab$Covariate <- rownames(tab)
  colnames(tab) <- c("R2", "P", "Covariate")
  tab <- tab[1,]
  return(tab)
}
sub_design$NEC <- factor(sub_design$NEC)
tab_1 <- tab_adonis(disWuf, "Group")
tab_2 <- tab_adonis(disWuf, "Litter")
tab_3 <- tab_adonis(disWuf, "Sex")
tab_4 <- tab_adonis(disWuf, "NEC")
tab <- rbind.data.frame(tab_1, tab_2, tab_3, tab_4)
tab
tab$Covariate <- c("FMT", "SOW", "Sex", "NEC")
tab_f <- tab[tab$P < 0.05,]
tab_f$Covariate <- factor(tab_f$Covariate, levels = tab_f$Covariate[order(tab_f$R2)])

p_effectsize_jaccard <- tab_f %>%
  ggplot(aes(y=R2*100, x=Covariate)) +
  geom_bar(stat='identity', fill="darkgray") +
  coord_flip() +
  labs(title = "Clinical covarites associated with GM (Shotgun)", y= "Effect size(R2, %)", x="") +
  theme_classic()+
  mytheme
p_effectsize_jaccard
