library(reshape2)
library(tidyr)
library(phyloseq)
library(stringr)
#construct phyloseq r project for analysis
# kaiju species level
####################################################################
#load data
abun_tab <- read.table("data/kaiju_summary_species_nr.tsv", sep = "\t", header = T, quote = "", check.names = F)
head(abun_tab[abun_tab$taxon_name=="cannot be assigned to a (non-viral) species",])
table(abun_tab$taxon_name=="cannot be assigned to a (non-viral) species")

# clean up the format
abun_tab$file <- stringr::str_replace(abun_tab$file,"\\.out", "")
head(abun_tab)

abun_tab$taxon_name <- stringr::str_replace(abun_tab$taxon_name, "cannot be assigned to a \\(non-viral\\) species","Unassigned at species level")
abun_tab$taxon_name <- stringr::str_replace(abun_tab$taxon_name, "unclassified", "Unclassified")
# revise lactobacillus species name
abun_tab$taxon_name <- stringr::str_replace(abun_tab$taxon_name, "Lactobacillus reuteri", "Limosilactobacillus reuteri")
abun_tab$taxon_name <- stringr::str_replace(abun_tab$taxon_name, "Lactobacillus crispatus", "Lactobacillus crispatus")
abun_tab$taxon_name <- stringr::str_replace(abun_tab$taxon_name, "Lactobacillus amylovorus", "Lactobacillus amylovorus")
abun_tab$taxon_name <- stringr::str_replace(abun_tab$taxon_name, "Lactobacillus helveticus", "Lactobacillus helveticus")
abun_tab$taxon_name <- stringr::str_replace(abun_tab$taxon_name, "Lactobacillus agilis", "Ligilactobacillus agilis")
abun_tab$taxon_name <- stringr::str_replace(abun_tab$taxon_name, "Lactobacillus vaginalis", "Limosilactobacillus vaginalis")
abun_tab$taxon_name <- stringr::str_replace(abun_tab$taxon_name, "Lactobacillus equicursoris", "Lactobacillus equicursoris")
abun_tab$taxon_name <- stringr::str_replace(abun_tab$taxon_name, "Lactobacillus salivarius", "Ligilactobacillus salivarius")
# use Lactobacilli for all previous saying for Lactobacillus genus
abun_tab$taxon_name <- stringr::str_replace(abun_tab$taxon_name, "Lactobacillus;", "Lactobacilli;")

# calculate the unclassified percentage
abun_tab_unclassified <- subset(abun_tab, abun_tab$taxon_name %in% c("Unclassified"))

# clean up unclassified reads
abun_tab <- subset(abun_tab, abun_tab$taxon_name!="Unclassified")
abun_tab <- subset(abun_tab,abun_tab$taxon_name!= "Unassigned at species level")

# clean up the
colnames(abun_tab) <- c("SampleID", "percent","reads", "taxon_id", "taxon_name")

# only select the SampleID, reads, taxon_name to create a long-format table
abun_tab_f <- abun_tab[,c("SampleID","reads","taxon_name")]
abun_tab_wide <- dcast(abun_tab_f, taxon_name ~ SampleID, value.var="reads",fun.aggregate = sum)
rownames(abun_tab_wide) <- abun_tab_wide$taxon_name

taxonomy_dt <- data.frame(abun_tab_wide[,1, drop=F])
head(taxonomy_dt)
tax <-
  taxonomy_dt %>% separate(taxon_name, c("Domain","Phylum", "Class", "Order", "Family", "Genus", "Species"), ";")
tax <- data.frame(tax, stringsAsFactors = FALSE)
tax.clean <- tax
tax.clean[tax.clean=="NA"] <- ""
tax.clean[is.na(tax.clean)] <- ""
for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,1] == "") {
    tax.clean[i, 1] <- "Unclassified"
  }}
for (i in 1:nrow(tax.clean)){
  for (j in 2:7) {
    if (tax.clean[i,j] == ""){
      tax.clean[i, j] <- tax.clean[i,j-1]
    }
  }
}

rownames(tax.clean) <- sapply(c(1:dim(tax.clean)[1]), function(x) paste0("S_",x))
rownames(abun_tab_wide) <- rownames(tax.clean)
abun_tab_wide <- abun_tab_wide[,-1]
colnames(abun_tab_wide)
colnames(abun_tab_wide) <- stringr::str_replace(colnames(abun_tab_wide),"kaiju_", "")
head(tax.clean)
table(tax.clean$Species %in% c("Unclassified"))
# metadata
metadata <- read.table(file = "data/Metadata.tsv", sep = "\t", header = T, row.names = 1)
metadata$SeqID <- rownames(metadata)
#transform for matrix for phyloseq construction

OTU = as.matrix(abun_tab_wide)
TAX = as.matrix(tax.clean)

OTU = otu_table(OTU, taxa_are_rows = TRUE)
TAX = tax_table(TAX)
sampledata <- sample_data(metadata)
# merge the data
ps_all <- phyloseq(OTU,TAX,sampledata)
ps_all
saveRDS(ps_all, file = "data/ps_kaiju_nr_species.rds")

#############################################################################
###################kaiju genus phyloseq
#load data
abun_tab <- read.table("data/kaiju_summary_genus_nr.tsv", sep = "\t", header = T, quote = "", check.names = F)
head(abun_tab[abun_tab$taxon_name=="cannot be assigned to a (non-viral) genus",])
table(abun_tab$taxon_name=="cannot be assigned to a (non-viral) genus")

# clean up the format
abun_tab$file <- stringr::str_replace(abun_tab$file,"\\.out", "")
head(abun_tab)

abun_tab$taxon_name <- stringr::str_replace(abun_tab$taxon_name, "cannot be assigned to a \\(non-viral\\) genus","Unassigned at genus level")
abun_tab$taxon_name <- stringr::str_replace(abun_tab$taxon_name, "unclassified", "Unclassified")
# use Lactobacilli for all previous saying for Lactobacillus genus
abun_tab$taxon_name <- stringr::str_replace(abun_tab$taxon_name, "Lactobacillus;", "Lactobacilli;")

# calculate the unclassified percentage
abun_tab_unclassified <- subset(abun_tab, abun_tab$taxon_name %in% c("Unclassified"))
abun_tab_unclassified
summary(abun_tab_unclassified)


# clean up unclassified reads
abun_tab <- subset(abun_tab, abun_tab$taxon_name!="Unclassified")
abun_tab <- subset(abun_tab,abun_tab$taxon_name!= "Unassigned at genus level")


# clean up the
colnames(abun_tab) <- c("SampleID", "percent","reads", "taxon_id", "taxon_name")
head(abun_tab)
# only select the SampleID, reads, taxon_name to create a long-format table
abun_tab_f <- abun_tab[,c("SampleID","reads","taxon_name")]

abun_tab_wide <- dcast(abun_tab_f, taxon_name ~ SampleID, value.var="reads",fun.aggregate = sum)
rownames(abun_tab_wide) <- abun_tab_wide$taxon_name

taxonomy_dt <- data.frame(abun_tab_wide[,1, drop=F])
head(taxonomy_dt)
tax <-
  taxonomy_dt %>% separate(taxon_name, c("Domain","Phylum", "Class", "Order", "Family", "Genus"), ";")
tax <- data.frame(tax, stringsAsFactors = FALSE)
tax.clean <- tax
tax.clean[tax.clean=="NA"] <- ""
tax.clean[is.na(tax.clean)] <- ""
for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,1] == "") {
    tax.clean[i, 1] <- "Unclassified"
  }}
for (i in 1:nrow(tax.clean)){
  for (j in 2:6) {
    if (tax.clean[i,j] == ""){
      tax.clean[i, j] <- tax.clean[i,j-1]
    }
  }
}

rownames(tax.clean) <- sapply(c(1:dim(tax.clean)[1]), function(x) paste0("S_",x))
rownames(abun_tab_wide) <- rownames(tax.clean)
abun_tab_wide <- abun_tab_wide[,-1]
colnames(abun_tab_wide)
colnames(abun_tab_wide) <- stringr::str_replace(colnames(abun_tab_wide),"kaiju_", "")
head(tax.clean)
table(tax.clean$Genus %in% c("Unclassified"))
# metadata

metadata <- read.table(file = "data/Metadata.tsv", sep = "\t", header = T, row.names = 1)
metadata$SeqID <- rownames(metadata)
#transform for matrix for phyloseq construction

OTU = as.matrix(abun_tab_wide)
TAX = as.matrix(tax.clean)

OTU = otu_table(OTU, taxa_are_rows = TRUE)
TAX = tax_table(TAX)
sampledata <- sample_data(metadata)
# merge the data
ps_all <- phyloseq(OTU,TAX,sampledata)
ps_all
saveRDS(ps_all, file = "data/ps_kaiju_nr_genus.rds")

#############################################################################
#################kaiju phylum phyloseq
#load data
abun_tab <- read.table("data/kaiju_summary_phylum_nr.tsv", sep = "\t", header = T, quote = "", check.names = F)
head(abun_tab[abun_tab$taxon_name=="cannot be assigned to a (non-viral) phylum",])
table(abun_tab$taxon_name=="cannot be assigned to a (non-viral) genus")

# clean up the format
abun_tab$file <- stringr::str_replace(abun_tab$file,"\\.out", "")
head(abun_tab)

abun_tab$taxon_name <- stringr::str_replace(abun_tab$taxon_name, "cannot be assigned to a \\(non-viral\\) phylum","Unassigned at phylum level")
abun_tab$taxon_name <- stringr::str_replace(abun_tab$taxon_name, "unclassified", "Unclassified")
# calculate the unclassified percentage
abun_tab_unclassified <- subset(abun_tab, abun_tab$taxon_name %in% c("Unclassified"))
abun_tab_unclassified
summary(abun_tab_unclassified)


# clean up unclassified reads
abun_tab <- subset(abun_tab, abun_tab$taxon_name!="Unclassified")
abun_tab <- subset(abun_tab,abun_tab$taxon_name!= "Unassigned at phylum level")


# clean up the
colnames(abun_tab) <- c("SampleID", "percent","reads", "taxon_id", "taxon_name")
head(abun_tab)
# only select the SampleID, reads, taxon_name to create a long-format table
abun_tab_f <- abun_tab[,c("SampleID","reads","taxon_name")]

abun_tab_wide <- dcast(abun_tab_f, taxon_name ~ SampleID, value.var="reads",fun.aggregate = sum)
rownames(abun_tab_wide) <- abun_tab_wide$taxon_name

taxonomy_dt <- data.frame(abun_tab_wide[,1, drop=F])
head(taxonomy_dt)
tax <-
  taxonomy_dt %>% separate(taxon_name, c("Domain","Phylum"), ";")
tax <- data.frame(tax, stringsAsFactors = FALSE)
tax.clean <- tax
tax.clean[tax.clean=="NA"] <- ""
tax.clean[is.na(tax.clean)] <- ""
for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,1] == "") {
    tax.clean[i, 1] <- "Unclassified"
  }}
for (i in 1:nrow(tax.clean)){
  for (j in 2:2) {
    if (tax.clean[i,j] == ""){
      tax.clean[i, j] <- tax.clean[i,j-1]
    }
  }
}

rownames(tax.clean) <- sapply(c(1:dim(tax.clean)[1]), function(x) paste0("S_",x))
rownames(abun_tab_wide) <- rownames(tax.clean)
abun_tab_wide <- abun_tab_wide[,-1]
colnames(abun_tab_wide)
colnames(abun_tab_wide) <- stringr::str_replace(colnames(abun_tab_wide),"kaiju_", "")
head(tax.clean)
table(tax.clean$Domain %in% c("Unclassified"))
# metadata
library(phyloseq)
metadata <- read.table(file = "data/Metadata.tsv", sep = "\t", header = T, row.names = 1)
metadata$SeqID <- rownames(metadata)

OTU = as.matrix(abun_tab_wide)
TAX = as.matrix(tax.clean)

OTU = otu_table(OTU, taxa_are_rows = TRUE)
TAX = tax_table(TAX)
sampledata <- sample_data(metadata)
# merge the data
ps_all <- phyloseq(OTU,TAX,sampledata)
ps_all
saveRDS(ps_all, file = "data/ps_kaiju_nr_phylum.rds")

###############################################################
# KEGG module for phyloseq
####load module abundance matrix
ko_matrix <- read.table("data/ko_matrix.tsv", sep = "\t", header = T,row.names = 1,check.names = F)
# module
module_anno <- read.table("data/module_anno.tsv", sep = "\t", header = F, check.names = F, quote = "")
# ko to module
ko2module <- read.table("data/ko2module.tsv", sep = "\t", check.names = F, quote = "")

# match module annoation 
index <- match(ko2module$V1, module_anno$V1)
ko2module$module <- module_anno$V2[index]
ko2module$module <- mapply(function(x,y) paste0(x,":",y),stringr::str_replace(ko2module$V1,"md:",""),as.character(ko2module$module))
index <- match(rownames(ko_matrix), stringr::str_replace(ko2module$V2,"ko:",""))
KEGG_l3 <- ko2module$module[index]
# replace the unannotated
KEGG_l3 <- c(as.character(KEGG_l3[1:length(KEGG_l3)-1]),"Unannotated")

# replace NA with unclassified
KEGG_l3 <- replace(KEGG_l3, is.na(KEGG_l3),"Unclassified")
ko_matrix$module <- KEGG_l3

library(plyr)
l3_matrix <- ddply(ko_matrix,"module",numcolwise(sum)) 
rownames(l3_matrix) <- l3_matrix$module
l3_matrix <- l3_matrix[,colnames(l3_matrix)!="module"]
abun_tab <-l3_matrix

tax.clean <- data.frame(module=rownames(abun_tab)) 
rownames(tax.clean) <- tax.clean$module

# build phyloseq project
metadata <- read.table("data/Metadata.tsv", sep = "\t", header = T, row.names = 1, check.names = F)
# merge the data
OTU <- otu_table(as.matrix(abun_tab), taxa_are_rows = T)
TAX <- tax_table(as.matrix(tax.clean))
sampledata <- sample_data(metadata)
ps <- phyloseq(OTU,TAX,sampledata)
ps
saveRDS(ps, file = "data/ps_gf_keggmodule.rds")

#############################################################################################################
# MAG-abundance phyloseq project
####load mags abundance matrix
abun_tab <- read.table("data/bin_abundance_table_sinASSsinBIN_c50t5.tab", sep = "\t", header = T, row.names = 1, check.names = F)
rownames(abun_tab) <- stringr::str_replace(rownames(abun_tab),"sinASSsinBIN_","")
rownames(abun_tab) <- stringr::str_replace(rownames(abun_tab),"bin\\.","bin_")
#load donor source and tax
tax_donor <- read.table("data/mag_taxaSource.tsv", sep = "\t", header = T, row.names = 1)
abun_tab <- t(t(abun_tab)/colSums(abun_tab))

# remove bins without reads mapping
abun_tab <- abun_tab[rownames(tax_donor),]
#####################
#import the taxonomy
tax.clean <- subset(tax_donor, select = c(source.1,genus.1, species))
colnames(tax.clean) <- c("source", "genus", "species")
# revise lactobacillus species name
tax.clean$species <- stringr::str_replace(tax.clean$species, "Lactobacillus reuteri", "Limosilactobacillus reuteri")
tax.clean$species <- stringr::str_replace(tax.clean$species, "Lactobacillus crispatus", "Lactobacillus crispatus")
tax.clean$species <- stringr::str_replace(tax.clean$species, "Lactobacillus amylovorus", "Lactobacillus amylovorus")
tax.clean$species <- stringr::str_replace(tax.clean$species, "Lactobacillus helveticus", "Lactobacillus helveticus")
tax.clean$species <- stringr::str_replace(tax.clean$species, "Lactobacillus agilis", "Ligilactobacillus agilis")
tax.clean$species <- stringr::str_replace(tax.clean$species, "Lactobacillus vaginalis", "Limosilactobacillus vaginalis")
tax.clean$species <- stringr::str_replace(tax.clean$species, "Lactobacillus equicursoris", "Lactobacillus equicursoris")
tax.clean$species <- stringr::str_replace(tax.clean$species, "Lactobacillus salivarius", "Ligilactobacillus salivarius")
# use Lactobacilli for all previous saying for Lactobacillus genus
tax.clean$genus <- stringr::str_replace(tax.clean$genus, "Lactobacillus;", "Lactobacilli;")

# build phyloseq project
#import relavant metadata
metadata <- read.table("data/Metadata.tsv", sep = "\t", header = T, row.names = 1, check.names = F)
# merge the data
OTU <- otu_table(as.matrix(abun_tab), taxa_are_rows = T)
TAX <- tax_table(as.matrix(tax.clean))
sampledata <- sample_data(metadata)

ps <- phyloseq(OTU,TAX,sampledata)
ps
saveRDS(ps, file = "data/ps_denovo_mags_donorSource.rds")

###############################################################################
# 16S rRNA amplicon sequencing phyloseq
# import data through phyloseq
otu <- read.table(file = "data/feature-table_16s.tsv", sep = "\t", header = T, row.names = 1)
taxonomy <- read.table(file = "data/taxonomy_16s.tsv", sep = "\t", header = T ,row.names = 1)
head(taxonomy)
taxonomy_f <- data.frame(taxonomy[rownames(otu),,drop=FALSE])

t_tab <- 
  taxonomy_f %>% separate(Taxon, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), ";")
#green gene format change
tax <- data.frame(t_tab)
tax.clean <- data.frame(row.names = row.names(tax),
                        Kingdom = str_replace(tax[,1], "k__",""),
                        Phylum = str_replace(tax[,2], "p__",""),
                        Class = str_replace(tax[,3], "c__",""),
                        Order = str_replace(tax[,4], "o__",""),
                        Family = str_replace(tax[,5], "f__",""),
                        Genus = str_replace(tax[,6], "g__",""),
                        Species = str_replace(tax[,7], "s__",""),
                        stringsAsFactors = FALSE)

tax.clean[is.na(tax.clean)] <- ""
tax.clean[tax.clean=="__"] <- ""

for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,7] != ""){
    tax.clean$Species[i] <- paste(tax.clean$Genus[i], tax.clean$Species[i], sep = " ")
  } else if (tax.clean[i,2] == ""){
    kingdom <- paste("Kingdom_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Phylum_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Class_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Order_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Family_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Genus",tax.clean$Genus[i], sep = "_")
  }
}


metadata <- read.table(file = "data/Metadata_16s.tsv", sep = "\t", header = T, row.names = 1)
mapping_new <- read.table(file = "data/Metadata.tsv", header = T, row.names = 1)
index <- match(metadata$FMT_ID, mapping_new$FMT_ID)
index
metadata$NEC <-  mapping_new$NEC[index]
metadata$NEC[is.na(metadata$NEC)] <- 0
metadata$SeqID <- rownames(metadata)

#transform for matrix for phyloseq construction
OTU = as.matrix(otu)
TAX = as.matrix(tax.clean)

OTU = otu_table(OTU, taxa_are_rows = TRUE)
TAX = tax_table(TAX)
sampledata <- sample_data(metadata)
TREE = read_tree("data/tree_16s.nwk")
# merge the data
ps_all <- phyloseq(OTU,TAX,sampledata,TREE)
ps_all
ps <- prune_taxa(rowSums(otu_table(ps_all))>0,ps_all)
ps
# change genus name lactobacillus to lactobacilli
tax_lab <- tax_table(ps)
tax_lab[,6] <- stringr::str_replace(tax_lab[,6], "Lactobacillus", "Lactobacilli")
tax_table(ps) <- tax_lab
saveRDS(ps, "data/phyloseq_16s.rds")
##################################################################################
sessionInfo()
