####################################################
# build mapping file for the vcontact2 network
votu <- read.table("data/votus/depth.txt", sep = "\t", header = TRUE, row.names = 1)
# donor2: 39 donor1: 40
colnames(votu)
votu.donor <- subset(votu, select = c("X17119.05.39.243200084.sort.bam", "X17119.05.40.243200090.sort.bam"))
library(tidyverse)
colnames(votu.donor) <- c("donor2", "donor1")
votu.donor <- votu.donor %>%
  mutate(source = case_when(donor2 > 0 & donor1 == 0 ~ "Donor 2 specific",
                               donor2 == 0 & donor1 > 0 ~ "Donor 1 specific",
                               donor2 > 0 & donor1 > 0 ~ "Donors shared",
                               TRUE ~ "Recipients specific"))
table(votu.donor$source)
# include the vcontact2 genome2genoem file
vct2_mapping <- read.table("data/votus/genome_by_genome_overview_tax_predictions.txt", sep = "\t", header = TRUE, row.names = 1)
vct2_mapping <- vct2_mapping %>%
  tibble::rownames_to_column() %>%
  select(rowname, Order, Family, Genus)

vct2_mapping.sample <- vct2_mapping[stringr::str_detect(vct2_mapping$rowname, "\\|\\|"),]
vct2_mapping.ref <- vct2_mapping[! stringr::str_detect(vct2_mapping$rowname, "\\|\\|"),]
# specify data source
vct2_mapping.sample$source <- "FMT study"
vct2_mapping.ref$source <- "vConTact2"
# create mapping
vct2_mapping.ref$mapping <- gsub("[[:punct:]].*", "", vct2_mapping.ref$rowname)
# merge the tax to the order level
unique(vct2_mapping.ref$mapping)

# # import greengene taxonomy to help find its order
# gg_tax <- read.table("data/votus/gg_taxU_l6.csv", sep = ";", quote = "")
# index <- match(vct2_mapping.ref$mapping, gg_tax$V6)
# vct2_mapping.ref$mapping_ <- gg_tax$V4[index]
# 
# # put archeae and actinobacteria
# archaea_l <- unique(gg_tax$V4[gg_tax$V1 == "Archaea"])
# actinobacteria_l <- unique(gg_tax$V4[gg_tax$V2 == "Actinobacteria"])
# vct2_mapping.ref$mapping_ <- sapply(vct2_mapping.ref$mapping_, function(x) ifelse(x %in% archaea_l, "Archaea", x))
# vct2_mapping.ref$mapping_ <- sapply(vct2_mapping.ref$mapping_, function(x) ifelse(x %in% actinobacteria_l, "Actinobacteria", x))
# # include the na
# vct2_mapping.ref$mapping_ <- mapply(function(x,y) ifelse(is.na(x), y, x), x=vct2_mapping.ref$mapping_, y=vct2_mapping.ref$mapping)
# 
# unique(vct2_mapping.ref$mapping_)

n=15
top_n <- names(sort(table(vct2_mapping.ref$mapping),decreasing = T))[1:n]
vct2_mapping.ref$mapping2 <- sapply(vct2_mapping.ref$mapping, function(x) ifelse(x %in% top_n, x, "Others"))
unique(vct2_mapping.ref$mapping2)

index <- match(vct2_mapping.sample$rowname, rownames(votu.donor))
vct2_mapping.sample$mapping <- votu.donor$source[index]
vct2_mapping.sample$mapping2 <- vct2_mapping.sample$mapping

vct2_mapping_combined <- rbind.data.frame(vct2_mapping.ref, vct2_mapping.sample)
# import color palette
#join all qualitative palettes
library(RColorBrewer)
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- unique(col_vector)
# build color palette
num_ <- length(unique(vct2_mapping_combined$mapping2))
#col_vector.c  <- sample(col_vector, num_)
col_vector.c <- colorRampPalette(col_vector[1:num_])(num_)
vct2_mapping_combined$color <- factor(vct2_mapping_combined$mapping2, labels = col_vector.c)

write.table(vct2_mapping_combined, "table/vct2_mapping.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
unique(vct2_mapping_combined$color)
unique(vct2_mapping_combined$mapping2)
# Add a legend
pdf("figure/vct2_legend.pdf",height = 8, width = 10)
#palette(unique(vct2_mapping_combined$color))
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center", 
       legend = c(unique(vct2_mapping_combined$mapping2)[c(18:19,17,20,2:16,1)]), 
       col = c(as.character(unique(vct2_mapping_combined$color))[c(18:19,17,20,2:16,1)]), 
       pch = c(rep(19,4),rep(17,16)), 
       bty = "n", 
       ncol =2,
       pt.cex = 2, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1))
dev.off()

######################################################3
# create a flow chart between the phage and 
crt <- read.table("data/votus/predict_CRT.tsv", sep = "\t")
crisprdetect <- read.table("data/votus/predict_CRISPRDetect.tsv", sep = "\t")
#spacepharer <- read.table("data/votus/predict_SpacePHARER.tsv", sep = "\t")
consensus <- merge(crt,crisprdetect, by = c("V1", "V2"))
# add information for mags and donors
# add coverage information for viral contigs
index <- match(consensus$V2, rownames(votu.donor))
consensus$vi_coverage <- votu.donor$source[index]
# add MAGs genus information
library(phyloseq)
ps <- readRDS("data/ps_denovo_mags_donorSource.rds")
tax <- data.frame(tax_table(ps),keep.rownames = TRUE)
head(tax)
tax$source <- sapply(tax$source, function(x) ifelse(x == "Others", "Recipients specific", x))

index <- match(consensus$V1, rownames(tax))
consensus$host_species <- tax$species[index]
# uncultured Blautia sp. OX=765821 
# uncultured Clostridium sp. OX=59620 
consensus$host_species[consensus$host_species == "uncultured Blautia sp. OX=765821 "] <- "Blautia sp. OX=765821"
consensus$host_species[consensus$host_species == "uncultured Clostridium sp. OX=59620 "] <- "Clostridium sp. OX=59620"
consensus$host_genus <- gsub("\ .*","", consensus$host_species)
consensus$host_genus <- gsub("\\[","", consensus$host_genus)
consensus$host_genus <- gsub("\\]","", consensus$host_genus)

others_list <- c("Clostridiales","Ruminococcaceae", "Ruminococcaceae", "Planctomycetaceae", "Erysipelotrichaceae", "Firmicutes","Rikenellaceae")
consensus$host <- sapply(consensus$host_genus, function(x) ifelse(x %in% others_list, "Others", x))
unique(consensus$host)
# add MAGs detection information
consensus$host_source <- tax$source[index]
consensus_f <- consensus %>%
  select(vi_coverage, host,host_source) %>%
  group_by(host, vi_coverage, host_source) %>%
  summarise(n=n())
colnames(consensus_f) <- c("Host", "Phage", "MAG", "n")

library(ggalluvial)
library(ggplot2)
library(extrafont)

loadfonts()
mytheme<- theme(plot.title = element_text(hjust=0.5, family = "Arial", size=12),
                #legend.position ="right",
                legend.text = element_text(family = "Arial", size = 10),
                legend.background = element_blank(),
                strip.background = element_blank(),
                strip.placement = "outside",
                strip.text = element_text(family = "Arial", size = 10),
                axis.title = element_text(family = 'Arial', size = 10),
                axis.text = element_text(family = 'Arial', size = 10, color="black"),
                panel.grid = element_blank())
# build new color-palette
library(RColorBrewer)
#join all qualitative palettes
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector = unique(col_vector)

my_palette <- col_vector[1:length(unique(consensus_f$Host))]


p <- ggplot(as.data.frame(consensus_f),
       aes(y = n, axis1 = Phage, axis2 = MAG)) +
  geom_alluvium(aes(fill = Host)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Phage", "MAG"), expand = c(.05, .05)) +
  scale_fill_manual(values = my_palette) +
  labs(y="") +
  guides(fill=guide_legend(ncol= 5)) + 
  theme_classic() + 
  mytheme +
  theme(axis.ticks.x.bottom = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank())
   
ggsave("figure/figure_viSpacer.png",p, height = 8, width = 7.5)
ggsave("figure/figure_viSpacer.pdf",p, height = 8, width = 7.5, device = cairo_pdf)
