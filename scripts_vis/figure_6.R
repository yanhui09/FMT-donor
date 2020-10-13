library(phyloseq)

# load data to calculate the matrix
ps <- readRDS("data/ps_kaiju_nr_species.rds")
ps.rel <- transform_sample_counts(ps, function(x) x/sum(x)*100)
ps.sample_f <- prune_samples(!(sample_data(ps.rel)$Group %in% c("DONOR1", "DONOR2")), ps.rel)
ps.sample_f

# remove species less 0.1% in all samples
mat <- as.matrix(otu_table(ps.sample_f))
head(mat)
#species2keep <- rownames(mat)[rowMeans(mat)>1]
# species appear in 10% samples with relative abundance > 1%
species2keep <- rownames(mat)[rowSums(mat >= 1)/length(colnames(mat))> 0.1]
#species2keep <- rownames(mat)

ps.sample_f <- prune_taxa(species2keep,ps.sample_f)
ps.sample_f

mat_tax <- data.frame(tax_table(ps.sample_f))
mat_tax$Genus <- as.character(mat_tax$Genus)
mat_tax$Species <- as.character(mat_tax$Species)
mat_tax$Genus[mat_tax$Species=='[Pasteurella] aerogenes'] <- 'Pasteurella'
mat_tax$Species[mat_tax$Species=='[Pasteurella] aerogenes'] <- 'Pasteurella aerogenes'

mat <- as.matrix(otu_table(ps.sample_f))

mat_scaled_orig = t(scale(t(mat)))
mapping <- data.frame(sample_data(ps.sample_f))
nec <- mapping[,"NEC_severity",drop=FALSE]

####
#include the diarhea scores, MPO colon, ABPAS colon
library(readxl)
list.data <- list()
for (i in sheets) {
  list.data[[i]] <- read_excel("data/phenotype.xlsx",sheet=i)
}

sheets
clinic_table <- as.data.frame(list.data$Clinic)
diarhea <- clinic_table[,c("ID","Diarrhea_severity")]
histo_table <- as.data.frame(list.data$Histology)
colnames(histo_table)
histo <- histo_table[,c("ID","ABPAS_colon_fraction", "CD3_colon_fraction", "MPO_colon_score", "FISH_SI_score")]
merge_tab <- merge(diarhea, histo, by="ID")
#match mapping
index <- merge_tab$ID %in% mapping$FMT_ID 
table(index)
merge_tab <- merge_tab[index,]
index <- match(merge_tab$ID, mapping$FMT_ID)
rownames(merge_tab) <- rownames(mapping)[index]
colnames(merge_tab) <- c("ID", "Diarrhea", "Goblet cell density", "CD3+ cell density", "MPO score", "FISH score")

### function to calculate respective scaled perason's correlation
cor_tab <- function(mat_scaled, nec){
#  mat_scaled <- mat_scaled_orig
mat_scaled <- mat_scaled[,rownames(nec)]
nec <- as.matrix(nec)
nec_scaled <- scale(nec)

cor_tab <- cor(t(mat_scaled), nec_scaled, method = "pearson")
cor_p_tab <- psych::corr.test(t(mat_scaled), nec_scaled, method = "pearson")$p
max(cor_tab)
min(cor_tab)
index <- match(rownames(cor_tab), rownames(mat_tax))
cor_tab_full <- cbind.data.frame(cor_tab, mat_tax[index,])
#includ p
index <- match(rownames(cor_tab_full), rownames(cor_p_tab))
cor_tab_full_p <- cbind(cor_tab_full, cor_p=cor_p_tab[index,])

return(cor_tab_full_p)
}


# nec
cor_tab_full_nec <- cor_tab(mat_scaled_orig,nec)
cor_tab_full_nec_f <- cor_tab_full_nec[abs(cor_tab_full_nec$NEC_severity)>0.2,]

# dirrhea
diarrhea <- merge_tab[,"Diarrhea",drop=FALSE]
diarrhea <- na.omit(diarrhea)
cor_tab_full_diarrhea <- cor_tab(mat_scaled_orig,diarrhea)
cor_tab_full_diarrhea_f <- cor_tab_full_diarrhea[abs(cor_tab_full_diarrhea$Diarrhea)>0.2,]

# ABPAS mucin
mucin <- merge_tab[,"Goblet cell density",drop=FALSE]
mucin <- na.omit(mucin)
cor_tab_full_mucin <- cor_tab(mat_scaled_orig, mucin)
cor_tab_full_mucin_f <- cor_tab_full_mucin[abs(cor_tab_full_mucin$`Goblet cell density`)>0.2,]

# CD3+ cell density
CD3 <- merge_tab[,"CD3+ cell density",drop=FALSE]
CD3 <- na.omit(CD3)
cor_tab_full_CD3 <- cor_tab(mat_scaled_orig, CD3)
cor_tab_full_CD3_f <- cor_tab_full_CD3[abs(cor_tab_full_CD3$`CD3+ cell density`)>0.2,]

# MPO
MPO <- merge_tab[,"MPO score",drop=FALSE]
MPO <- na.omit(MPO)
cor_tab_full_MPO <- cor_tab(mat_scaled_orig,MPO)
cor_tab_full_MPO_f <- cor_tab_full_MPO[abs(cor_tab_full_MPO$`MPO score`)>0.2,]

# FISH score
FISH <- merge_tab[,"FISH score",drop=FALSE]
FISH <- na.omit(FISH)
cor_tab_full_FISH <- cor_tab(mat_scaled_orig,FISH)
cor_tab_full_FISH_f <- cor_tab_full_FISH[abs(cor_tab_full_FISH$`FISH score`)>0.2,]

# create a network for the nec severity
library(igraph)
library(ggraph)

tidy_cors_nec <- data.frame(x="NEC severity", y=cor_tab_full_nec_f$Species, r=cor_tab_full_nec_f$NEC_severity)
tidy_cors_diarrhea <- data.frame(x="Diarrhea severity", y=cor_tab_full_diarrhea_f$Species, r=cor_tab_full_diarrhea_f$Diarrhea)
tidy_cors_mucin <- data.frame(x="Goblet cell density", y=cor_tab_full_mucin_f$Species, r=cor_tab_full_mucin_f$`Goblet cell density`)
tidy_cors_CD3 <- data.frame(x="CD3+ cell density", y=cor_tab_full_CD3_f$Species, r=cor_tab_full_CD3_f$`CD3+ cell density`)
tidy_cors_MPO <- data.frame(x="MPO score", y=cor_tab_full_MPO_f$Species, r=cor_tab_full_MPO_f$`MPO score`)
tidy_cors_FISH <- data.frame(x="FISH score", y=cor_tab_full_FISH_f$Species, r=cor_tab_full_FISH_f$`FISH score`)


tidy_cors <- rbind.data.frame(tidy_cors_nec, tidy_cors_diarrhea, tidy_cors_mucin, 
                              tidy_cors_CD3,tidy_cors_MPO, tidy_cors_FISH)
colnames(tidy_cors) <- c("x","y","Coefficient")

graph_cors <- tidy_cors %>%
  graph_from_data_frame(directed = FALSE)

p_NEC <- ggraph(graph_cors, layout = 'linear', circular=TRUE) +
  geom_edge_link(aes(edge_alpha = abs(Coefficient), edge_width = abs(Coefficient), color = Coefficient)) +
  guides(edge_alpha = "none", edge_width = "none") +
  scale_edge_colour_gradientn(limits = c(-1, 1), colors = c("dodgerblue2","firebrick2")) +
  geom_node_point(color = "white", size = 5) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_graph(base_family="sans") 
p_NEC

ggsave("figure/figure_6.png", height = 7, width=10)
