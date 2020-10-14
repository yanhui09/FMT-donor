scaffold2mag <- read.table("data/allGenomes_ctb.tab", sep = "\t", header = FALSE)
gene2ko <- read.table("data/mag_kegg_annotated.txt", sep = "\t", header = FALSE)

# extract before last underscore
gene2ko$contig <- stringr::str_replace(gene2ko$V1,"_[^_]+$", "")
scaffold2mag$mag <- stringr::str_replace(scaffold2mag$V2, "\\.fa", "")

gene2ko_s <- subset(gene2ko, select = c(V2, contig))
colnames(gene2ko_s) <- c("ko", "contig")
scaffold2mag_s <- subset(scaffold2mag, select = c(V1, mag))
colnames(scaffold2mag_s) <- c("contig", "mag")
head(scaffold2mag_s)

ko2mag_long <- merge(gene2ko_s, scaffold2mag_s, by = "contig") 
ko2mag_long_s <- subset(ko2mag_long, select=c(ko,mag))

# ko to module
module_anno <- read.table("data/module_anno.tsv", sep = "\t", header = F, check.names = F, quote = "")
ko2module <- read.table("data/ko2module.tsv", sep = "\t", check.names = F, quote = "")

# match module annoation 
index <- match(ko2module$V1, module_anno$V1)
ko2module$module <- module_anno$V2[index]
ko2module$module <- mapply(function(x,y) paste0(x,":",y),stringr::str_replace(ko2module$V1,"md:",""),as.character(ko2module$module))
ko2module$koID <- stringr::str_replace(ko2module$V2,"ko:","")
#match annotation for the different kos

index <- match(ko2mag_long_s$ko, ko2module$koID)
ko2mag_long_s$path <- ko2module$module[index]
ko2mag_long_s_path <- subset(ko2mag_long_s, select = c(mag,path))
ko2mag_long_s_path <- na.omit(ko2mag_long_s_path)
# spread on pathway

# to wide format
library(reshape2)
ko2mag_wide <- dcast(ko2mag_long_s_path, path~ mag)
table(rowSums(ko2mag_wide[,-1])==0)
rownames(ko2mag_wide) <- ko2mag_wide$path
ko2mag_wide <- ko2mag_wide[,-1]

# load list of the .csv files
dataFiles <- lapply(Sys.glob("data/coverm_genomewide/*.tsv"), function(x) read.table(file=x,sep="\t", header = T))
name <- list.files("data/coverm_genomewide/", pattern = "*.tsv")
name <- stringr::str_replace(name, "\\.tsv", "")

library(dplyr)
coverage_list <- lapply(dataFiles, function(x) x%>% select(c("Genome", "Mean")))
coverage_list <- lapply(coverage_list, function(x) mutate_all(x, as.character))
breadth_list <- lapply(dataFiles, function(x) x%>% select(c("Genome", "Covered.Fraction")))
breadth_list <- lapply(breadth_list, function(x) mutate_all(x, as.character))
library(tidyr)
# list to dataframe
ls2df <- function(coverage_list,x){
coverage_df <- coverage_list %>%
  bind_rows(.id="sample") 
coverage_df <- coverage_df %>%
  spread(key = "sample", value = x)

rownames(coverage_df) <- coverage_df$Genome
coverage_df <- coverage_df[,as.character(c(1:40))]
colnames(coverage_df) <- name
#coverage_df[is.na(coverage_df)] <- 0
return(coverage_df)
}

coverage_df <- ls2df(coverage_list,"Mean")
breadth_df <- ls2df(breadth_list,"Covered.Fraction")
# mapping # group the mags according to the breadth, breath> 0.25 will be consider as detected
mapping <- read.table("data/Metadata.tsv", sep = "\t", header = T, row.names = 1)
donor1 <- rownames(mapping)[mapping$Group=="DONOR1"]
donor2 <- rownames(mapping)[mapping$Group=="DONOR2"]
breadth_donor1 <- breadth_df[,donor1, drop=FALSE]
breadth_donor1$`17119-05-40-243200090` <- as.numeric(breadth_donor1$`17119-05-40-243200090`)
mags_donor1 <- rownames(breadth_donor1)[breadth_donor1> 0.5]

breadth_donor2 <- breadth_df[,donor2, drop=FALSE]
breadth_donor2$`17119-05-39-243200084` <- as.numeric(breadth_donor2$`17119-05-39-243200084`)
mags_donor2 <- rownames(breadth_donor2)[breadth_donor2> 0.5]

shared <- intersect(mags_donor1,mags_donor2)
donor1_specific <- mags_donor1[!mags_donor1%in%shared]
donor2_specific <- mags_donor2[!mags_donor2%in%shared]

# preparing environment file
mags_group <- data.frame(mag=colnames(ko2mag_wide)) 
mags_group$source <- case_when(mags_group$mag %in% shared ~ 'Donors shared',
                               mags_group$mag %in% donor1_specific ~ 'Donor1 specific',
                               mags_group$mag %in% donor2_specific ~ 'Donor2 specific',
                               TRUE ~ 'Others')
rownames(mags_group) <- mags_group$mag
mags_group$source <- as.factor(mags_group$source)
summary(mags_group)

# metawrap taxonomy
tax <- read.table("data/mag_taxaSource.tsv", sep = "\t", header = T, row.names = 1)
tax.clean <- subset(tax, select = c("genus","species"))
#match mag_groups
index <- match(rownames(mags_group), rownames(tax.clean))
mags_group <- cbind.data.frame(mags_group, tax.clean[index,])
write.table(mags_group,"data/mag_donorSource.tsv", sep="\t", quote = FALSE)
# only focus the mags in the donors
mags_group <- mags_group[mags_group$source!="Others",]

######CCA 
library(vegan)
ko_tab <- t(ko2mag_wide)
colnames(ko_tab)
#filtering
ko_tab <- ko_tab[rownames(mags_group),]

#filtering colSums==0
ko_tab <- ko_tab[, colSums(ko_tab)>0]
colnames(ko_tab)

mags_group$source <- factor(mags_group$source, levels = c("Donors shared", "Donor1 specific", "Donor2 specific"))
db_rda <- capscale(ko_tab~source, mags_group,distance='bray', add = TRUE)
dbanov <- anova.cca(db_rda, by="term", permutations = 999)
dbanov
cca_wa <- db_rda$CCA$wa
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(5)
# match the source
index <- match(rownames(cca_wa), mags_group$mag)
cca_source <- mags_group$source[index]

# extract dbrda using scaling 1
db_rda.scaling1 <- summary(db_rda)
#RsquareAdj() extract R2 ?RsquareAdj()
r2 <- RsquareAdj(db_rda)
db_rda_noadj <- r2$r.squared #original R2
db_rda_adj <- r2$adj.r.squared       #adjusted R2
db_rda_noadj
db_rda_adj

# recalculate using adjusted r2
db_rda_exp_adj <- db_rda_adj * db_rda$CCA$eig/sum(db_rda$CCA$eig)
db_rda_eig_adj <- db_rda_exp_adj * db_rda$tot.chi

# permutation test
#globel tests, on all axis 999 ? anova.cca
db_rda_test <- anova.cca(db_rda, permutations = 999)
db_rda_test


# visualizatiom with ggplot2
#extract r2 and enviroment factors，firtst 2 axies，scaling1
db_rda.scaling1 <- summary(db_rda, scaling = 1)
db_rda.site <- data.frame(db_rda.scaling1$sites)[1:2]
db_rda.env <- data.frame(db_rda.scaling1$biplot)[1:2]
db_rda.spe <- data.frame(db_rda.scaling1$species)[1:2]
# use fisher test to find different kos between donor1 specific and donor2 specific
ko_tab_donor12 <- as.data.frame(t(ko_tab))
mags_group_f <- mags_group[mags_group$source!='Donors shared',]

ko_tab_donor12 <- ko_tab_donor12[,as.character(mags_group_f$mag)]
#using binary
ko_tab_donor12[ko_tab_donor12>1] <- 1 
ko_tab_donor12 <- as.data.frame(ko_tab_donor12)
#remove all 0s
ko_tab_donor12 <- ko_tab_donor12[rowSums(ko_tab_donor12)> 0,]  
#remove all 1s
ko_tab_donor12 <- ko_tab_donor12[!(rowMeans(ko_tab_donor12)==1),]
ko_tab_donor12$koID <- rownames(ko_tab_donor12)

ko_tab_donor12_long <- reshape2::melt(ko_tab_donor12,id.vars="koID",variable.name="MAGs", value.name="Presence")
mags_group_f$source <- as.character(mags_group_f$source)
index <- match(ko_tab_donor12_long$MAGs, mags_group_f$mag)
ko_tab_donor12_long$source <- mags_group_f$source[index]
#############################fisher test
fisher_p <- function(i){
  abun_tab_f_long_ff <- subset(ko_tab_donor12_long,koID==i)
  test <- fisher.test(xtabs(~Presence+source, abun_tab_f_long_ff))
  return(test$p.value)
}


p_val <- fisher_p(unique(ko_tab_donor12_long$koID)[1])
for (i in unique(ko_tab_donor12_long$koID)[-1]) {
  p_val <- c(p_val, fisher_p(i))
}
############################# show which has more kos
enrich <- function(i){
abun_tab_f_long_ff <- subset(ko_tab_donor12_long,koID==i)
tab<- xtabs(~Presence+source, abun_tab_f_long_ff)
Donor1 <- tab[rownames(tab)=="1",colnames(tab)=="Donor1 specific"]/colSums(tab)["Donor1 specific"]
Donor2 <- tab[rownames(tab)=="1",colnames(tab)=="Donor2 specific"]/colSums(tab)["Donor2 specific"]
out <- ifelse(as.numeric(Donor2)>as.numeric(Donor1),"Donor2","Donor1")
return(out)
}

Enrichment <- c()
for (i in unique(ko_tab_donor12_long$koID)) {
  Enrichment <- c(Enrichment,enrich(i))
}

#p_val <- p_val[-1]
p_tab <- data.frame(cbind(koID=unique(ko_tab_donor12_long$koID),p=p_val,enrichment=Enrichment))
p_tab$p <- as.numeric(levels(p_tab$p))[p_tab$p]
p_tab$p_adjusted <- p.adjust(p_tab$p, method = "BH")

p_tab_f <- p_tab[p_tab$p<0.05,]
########################################################
db_rda.spe.core <- db_rda.spe[as.character(p_tab_f$koID),]
rownames(db_rda.spe.core)
#add names and groups
rownames(db_rda.env)
db_rda.env$name <- c("Donor1 specific","Donor2 specific")
rownames(db_rda.site)
db_rda.site$name <- rownames(db_rda.site)
index <- match(db_rda.site$name, mags_group$mag)
db_rda.site$source <- mags_group$source[index]
db_rda.site$Genus <- mags_group$genus[index]
#only keep Lactobacillus
db_rda.site$Genus <- as.character(db_rda.site$Genus)
db_rda.site$Genus[db_rda.site$Genus!="Lactobacillus"] <- 'Others'
db_rda.site$Genus <- factor(db_rda.site$Genus, levels = c('Others','Lactobacillus'))

db_rda.spe.core$name <- rownames(db_rda.spe.core)
# only include M00079 M00078 CAG metabolism
selected <- c('M00078:Heparan sulfate degradation','M00079:Keratan sulfate degradation')
db_rda.spe.core <- subset(db_rda.spe.core, db_rda.spe.core$name %in% selected)

#calculate the adujsted r squared
exp_adj <-  db_rda$CCA$eig/sum(db_rda$CCA$eig)
rda1_exp <- paste('db-RDA1:', round(exp_adj[1]*100, 2), '%')
rda2_exp <- paste('db-RDA2:', round(exp_adj[2]*100, 2), '%')

#ggplot2
library(ggplot2)
library(extrafont)
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
# get ggplto2 default colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(5)

p_dbrda_mag <- ggplot(db_rda.site, aes(CAP1, CAP2)) +
  geom_point(aes(color = source, shape = Genus),size=2) +
  scale_color_manual(values = c(cols[1], cols[2], cols[3])) +
  theme_classic() +
  mytheme +
  labs(x = rda1_exp, y = rda2_exp) +
  guides(color=guide_legend(title="MAGs")) +
  geom_vline(xintercept = 0, color = 'gray', size = 0.5) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
  geom_segment(data = db_rda.env, aes(x = 0, y = 0, xend = CAP1, yend = CAP2), arrow = arrow(length = unit(0.2, 'cm')), size = 0.3, color = 'blue') +
  geom_text(data = db_rda.env, aes(CAP1 * 2, CAP2 * 1.2, label = name), color = 'blue', size = 3) +
  geom_text(data = db_rda.spe.core, aes(CAP1*0.1, CAP2*0.1, label=name), color='black', size=3) +
  theme(legend.position = "top",
        legend.box="vertical",
        legend.spacing.y = unit(0, 'cm')
        )

p_dbrda_mag
ggsave("figure/figure_7D.png", height = 5, width = 6)
#################################
#table_s3
module_out <- subset(p_tab, select = c(koID, p,p_adjusted, enrichment))
colnames(module_out) <- c("Module", "P", "Adjusted P","Enrichment")
module_out_f <- module_out[module_out$`Adjusted P`< 0.1,]
module_out_f <- module_out_f[order(module_out_f$`Adjusted P`),]
write.table(module_out_f,"table/table_s3.tsv", sep = "\t", row.names = F, quote = F)
#table_s4
#M00079 and M00078 level in Donor-specific MAGs and its phylogen
ko_tab_s <- ko_tab[,c('M00078:Heparan sulfate degradation','M00079:Keratan sulfate degradation')]
ko_tab_s <- as.data.frame(ko_tab_s)
ko_tab_s$mag <- rownames(ko_tab_s)
ko_tab_s_detailed <- merge(ko_tab_s, mags_group, by = "mag")
write.table(ko_tab_s_detailed, "table/table_s4.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
