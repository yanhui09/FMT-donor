crt <- read.table("data/votus/predict_CRT.tsv", sep = "\t")
crisprdetect <- read.table("data/votus/predict_CRISPRDetect.tsv", sep = "\t")
#spacepharer <- read.table("data/votus/predict_SpacePHARER.tsv", sep = "\t")

consensus <- merge(crt,crisprdetect, by = c("V1", "V2"))

# votus
votu_tax <- read.table("data/votus/tax_predict_table.txt", sep = "\t", header = TRUE, row.names = 1)
votu_tax <- votu_tax[rownames(votu),]
# load donor source
library(phyloseq)
ps <- readRDS("data/ps_denovo_mags_donorSource.rds")
#only take the donors belong mags
tax <- data.frame(tax_table(ps),keep.rownames = TRUE)
# add bin information
index <- match(consensus$V1,rownames(tax)) 
consensus <- cbind.data.frame(consensus,tax[index,])
# add viral information
index <- match(consensus$V2, rownames(votu_tax))
consensus <- cbind(consensus, votu_tax[index,])


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
                               donor2 > 0 & donor1 > 0 ~ "Shared by donors",
                               TRUE ~ "Not from donors"))
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
#vct2_mapping.ref$mapping <- vct2_mapping.ref$Order
vct2_mapping.ref$mapping <- gsub("[[:punct:]].*", "", vct2_mapping.ref$rowname)
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
col_vector.c <- colorRampPalette(col_vector[1:num_])(num_)
vct2_mapping_combined$color <- factor(vct2_mapping_combined$mapping2, labels = col_vector.c)

write.table(vct2_mapping_combined, "table/vct2_mapping.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
