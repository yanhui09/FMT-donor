# Figure 3E. db-RDA shotgun db-RDA
library(ggplot2)
library(phyloseq)
library(extrafont)
getwd()
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

library(phyloseq)
library(vegan)
#species
species <- readRDS("data/ps_kaiju_nr_species.rds")
species <- prune_samples(sample_data(species)$Group!="DONOR1"&sample_data(species)$Group!="DONOR2", species)

species_tab <- as.data.frame(otu_table(species))
tax_tab <- tax_table(species)
tax_tab <- data.table::as.data.table(tax_tab, keep.rownames=TRUE)
index <- match(rownames(species_tab), tax_tab$rn)
rownames(species_tab) <- as.character(tax_tab$Species)[index]

species_tab <- species_tab[rowSums(species_tab)>0,]
species_tab <- t(species_tab)
# hellinger transformation
species_tab_h <- decostand(species_tab, method = 'hellinger')
mapping <- data.frame(sample_data(species))
mapping$SOW <- as.factor(mapping$Litter)
mapping$NEC <- factor(mapping$NEC, levels = c(0,1), labels = c("N","Y"))
mapping$Group  
library(vegan)
species_tab_h <- as.matrix(species_tab_h)
dis_bray <- vegdist(species_tab_h, method = 'bray')

env <- subset(mapping, select = c(Group,SOW, Sex, NEC, NEC_severity))
db_rda <- capscale(species_tab_h~., env, distance = 'bray', add = TRUE)
dbanov <- anova.cca(db_rda, by="term", permutations = 999)
dbanov
db_rda <- capscale(species_tab~Group*SOW, env, distance = 'bray', add = TRUE)
dbanov <- anova.cca(db_rda, by="term", permutations = 999)
db_rda$`Pr(>F)` <- p.adjust(db_rda$`Pr(>F)`, method = 'bonferroni')
dbanov

db_rda <- capscale(species_tab~Group+SOW, env, distance = 'bray', add = TRUE)
dbanov <- anova.cca(db_rda, by="term", permutations = 999)
dbanov

# extract dbrda using scaling 1
db_rda.scaling1 <- summary(db_rda, scaling = 1)
#RsquareAdj() extract R2 ?RsquareAdj() 
r2 <- RsquareAdj(db_rda)
db_rda_noadj <- r2$r.squared #original R2
db_rda_adj <- r2$adj.r.squared       #adjusted R2
db_rda_noadj
db_rda_adj

# recalculate using adjusted r2
db_rda_exp_adj <- db_rda_adj * db_rda$CCA$eig/sum(db_rda$CCA$eig)
db_rda_eig_adj <- db_rda_exp_adj * db_rda$tot.chi

db_rda_exp_adj
sum(db_rda_exp_adj)
db_rda_eig_adj
# permutation test
#globel tests, on all axis 999 ? anova.cca
db_rda_test <- anova.cca(db_rda, permutations = 999)
db_rda_test
#separate test on each axises
db_rda_test_axis <- anova.cca(db_rda, by = 'axis', permutations = 999)
# adjust P（BH）
db_rda_test_axis$`Pr(>F)` <- p.adjust(db_rda_test_axis$`Pr(>F)`, method = 'BH')
db_rda_test_axis

# co-linear
vif.cca(db_rda)
#variation partioning ?varpart
db_rda_vp <- varpart(dis_bray, env['Group'], env['SOW'], scale = FALSE, add = TRUE)
db_rda_vp


tab <- data.frame(var=c("FMT", "SOW", "Residuals", "Intersect"), Adj_R2=round(c(db_rda_vp$part$indfract[c(1,3,4), 3], 0), digits = 3))
head(tab)
tab$var <- factor(tab$var, levels = c("Residuals", "Intersect", "SOW", "FMT" ))
# make a barplot for visualization
p_bar <- ggplot(data = tab, aes(x= var, y= Adj_R2)) +
  geom_bar(stat = "identity", color="lightgray", fill="lightgray", position=position_dodge()) +
  labs(x = "", y = "", title = "Adjusted R-squared") +
  geom_text(aes(label=Adj_R2, y=Adj_R2 + 0.04), vjust=1, color="black",
            position = position_dodge(0.9), size=3.5)+
  coord_flip() + 
  theme_classic() +
  mytheme +
  theme(legend.position = "none")
p_bar

# Group explained merely
anova.cca(capscale(dis_bray~Group+Condition(SOW), env, add = TRUE), permutations = 999)
anova.cca(capscale(dis_bray~SOW+Condition(Group), env, add = TRUE), permutations = 999)
# no co-explained

# visualizatiom with ggplot2
#extract r2 and enviroment factors，firtst 2 axies，scaling1
db_rda.scaling1 <- summary(db_rda, scaling = 1)
db_rda.site <- data.frame(db_rda.scaling1$sites)[1:2]
db_rda.env <- data.frame(db_rda.scaling1$biplot)[1:2]
db_rda.spe <- data.frame(db_rda.scaling1$species)[1:2]
# only take the co-microbiome mean abundance > 1%

species_tab.rel <- transform_sample_counts(species, function(x) x/sum(x)*100)
species_tab.rel
species_tab.rel <- as.data.frame(otu_table(species_tab.rel))
tax_tab.rel <- tax_table(species)
tax_tab.rel <- data.table::as.data.table(tax_tab.rel, keep.rownames=TRUE)
index <- match(rownames(species_tab.rel), tax_tab.rel$rn)
rownames(species_tab.rel) <- as.character(tax_tab.rel$Species)[index]
species_core <- species_tab.rel[rowMeans(species_tab.rel)>1,]

db_rda.spe.core <- db_rda.spe[rownames(species_core),]


#add names and groups
rownames(db_rda.env)
db_rda.env$name <- c("FMT1","FMT2","SOW2")
rownames(db_rda.site)
db_rda.site$name <- rownames(db_rda.site)
index <- match(db_rda.site$name, rownames(mapping))
db_rda.site$group <- mapping$Group[index]
db_rda.site$SOW <- mapping$SOW[index]
db_rda.site$SOW <- factor(db_rda.site$SOW, levels = c(1,2), labels = c("SOW1", "SOW2"))

db_rda.spe.core$name <- rownames(db_rda.spe.core)
#calculate the adujsted r squared
exp_adj <- db_rda$CCA$eig/sum(db_rda$CCA$eig)
rda1_exp <- paste('db-RDA1:', round(exp_adj[1]*100, 2), '%')
rda2_exp <- paste('db-RDA2:', round(exp_adj[2]*100, 2), '%')

# get ggplto2 default colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(5)

p_dbrda <- ggplot(db_rda.site, aes(CAP1, CAP2)) +
  geom_point(aes(color = group#,shape=SOW
                 ),size=4) +
  scale_color_manual(values = c(cols[1], cols[2], cols[3])) +
  xlim(-1,1.5) +
  ylim(-0.8,1.2) +
  theme_classic() +
  mytheme +
  labs(x = rda1_exp, y = rda2_exp) +
  guides(color=guide_legend(title="FMT")) +
  geom_vline(xintercept = 0, color = 'gray', size = 0.5) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
  geom_text(data = db_rda.spe.core, aes(CAP1*0.4, CAP2*0.4, label=name), color='black', size=3)

p_dbrda




