#coASS coBinning
metabat2 <- read.table("data/BIN_REFINEMENT/coASS/metabat2_bins.stats", sep = "\t", header = T, row.names = 1)
maxbin2 <- read.table("data/BIN_REFINEMENT/coASS/maxbin2_bins.stats", sep = "\t", header = T, row.names = 1)
concoct <- read.table("data/BIN_REFINEMENT/coASS/concoct_bins.stats", sep = "\t", header = T, row.names = 1)
mmc_50_5 <- read.table("data/BIN_REFINEMENT/coASS/mmc_50_5_bins.stats", sep = "\t", header = T, row.names = 1)
vamb <- read.table("data/BIN_REFINEMENT/coASS/vamb_bins.stats", sep = "\t", header = T, row.names = 1)
metawrap_50_5 <- read.table("data/BIN_REFINEMENT/coASS/metawrap_50_5_bins.stats", sep = "\t", header = T, row.names = 1)

#function to adjust table 
metabat2$bin <- rownames(metabat2)
metabat2$type <- rep(deparse(substitute(metabat2)), dim(metabat2)[1])

maxbin2$bin <- rownames(maxbin2)
maxbin2$type <- rep(deparse(substitute(maxbin2)), dim(maxbin2)[1])

concoct$bin <- rownames(concoct)
concoct$type <- rep(deparse(substitute(concoct)), dim(concoct)[1])

mmc_50_5$bin <- rownames(mmc_50_5)
mmc_50_5$type <- rep(deparse(substitute(mmc_50_5)), dim(mmc_50_5)[1])

vamb$bin <- rownames(vamb)
vamb$type <- rep(deparse(substitute(vamb)), dim(vamb)[1])

metawrap_50_5$bin <- rownames(metawrap_50_5)
metawrap_50_5$type <- rep(deparse(substitute(metawrap_50_5)), dim(metawrap_50_5)[1])

bin_tab1 <- rbind.data.frame(metabat2,maxbin2,concoct,mmc_50_5,vamb,metawrap_50_5)
bin_tab1$strategy <- rep("coASScoBIN",dim(bin_tab1)[1])
  
#singleASS singleBinning
metabat2 <- read.table("data/BIN_REFINEMENT/singleASSBIN/metabat2_bins.stats", sep = "\t", header = T, row.names = 1)
maxbin2 <- read.table("data/BIN_REFINEMENT/singleASSBIN/maxbin2_bins.stats", sep = "\t", header = T, row.names = 1)
concoct <- read.table("data/BIN_REFINEMENT/singleASSBIN/concoct_bins.stats", sep = "\t", header = T, row.names = 1)
mmc_50_5 <- read.table("data/BIN_REFINEMENT/singleASSBIN/mmc_50_5_bins.stats", sep = "\t", header = T, row.names = 1)
vamb <- read.table("data/BIN_REFINEMENT/singleASSBIN/vamb_bins.stats", sep = "\t", header = T, row.names = 1)
metawrap_50_5 <- read.table("data/BIN_REFINEMENT/singleASSBIN/metawrap_50_5_bins.stats", sep = "\t", header = T, row.names = 1)

#function to adjust table 
metabat2$bin <- rownames(metabat2)
metabat2$type <- rep(deparse(substitute(metabat2)), dim(metabat2)[1])

maxbin2$bin <- rownames(maxbin2)
maxbin2$type <- rep(deparse(substitute(maxbin2)), dim(maxbin2)[1])

concoct$bin <- rownames(concoct)
concoct$type <- rep(deparse(substitute(concoct)), dim(concoct)[1])

mmc_50_5$bin <- rownames(mmc_50_5)
mmc_50_5$type <- rep(deparse(substitute(mmc_50_5)), dim(mmc_50_5)[1])

vamb$bin <- rownames(vamb)
vamb$type <- rep(deparse(substitute(vamb)), dim(vamb)[1])

metawrap_50_5$bin <- rownames(metawrap_50_5)
metawrap_50_5$type <- rep(deparse(substitute(metawrap_50_5)), dim(metawrap_50_5)[1])

bin_tab4 <- rbind.data.frame(metabat2,maxbin2,concoct,mmc_50_5,vamb,metawrap_50_5)
bin_tab4$strategy <- rep("sinASSsingleBIN",dim(bin_tab4)[1])


#########################################################################
#cobine all tha table
#bin_tab <- rbind.data.frame(bin_tab1,bin_tab2,bin_tab3,bin_tab4)
bin_tab <- rbind.data.frame(bin_tab1,bin_tab4)

#set completeness 50%, contamination 5% as threshold
library(dplyr)
bin_tab$quality <- rep("Waste(completeness<50% & contamination>10%)", dim(bin_tab)[1])
bin_tab <- bin_tab %>%
  filter(!grepl("unbin",bin)) %>% #filter out bin name with "unbin"
  mutate(quality=replace(quality,completeness<50&contamination<=10,"LQ-II(completeness<50% & contamination 5-10%))")) %>%
  mutate(quality=replace(quality,completeness<50&contamination<=5,"LQ-I(completeness<50% & contamination<=5%)")) %>%
  mutate(quality=replace(quality,completeness>=50&contamination<=10,"MQ-II(completeness>=50% & contamination 5-10%)")) %>%
  mutate(quality=replace(quality,completeness>=50&contamination<=5,"MQ-I(completeness>=50% & contamination<=5%)")) %>%
  mutate(quality=replace(quality,completeness>=90&contamination<=10,"HQ-II(completeness>=90% & contamination 5-10%)")) %>%
  mutate(quality=replace(quality,completeness>=90&contamination<=5,"HQ-I(completeness>=90% & contamination<=5%)"))

#visualization with ggplot2
library(extrafont)
library(ggplot2)
loadfonts()
mytheme<- theme(plot.title = element_text(hjust=0.5, family = "Arial", size=12),
                legend.title = element_blank(),
                legend.position ="right",
                legend.text = element_text(family = "Arial", size = 8),
                legend.background = element_blank(),
                strip.background = element_blank(),
                strip.placement = "outside",
                strip.text = element_text(family = "Arial", size = 8),
                #axis.title = element_text(family = 'Arial', size = 8),
                axis.text = element_text(family = 'Arial', size = 8, color="black"),
                panel.grid = element_blank())
#contig numbers
#adjust the order of assembly type and legend
bin_tab$type <- factor(bin_tab$type, levels = c("metabat2", "maxbin2","concoct","mmc_50_5","vamb","metawrap_50_5"))
bin_tab$quality <- factor(bin_tab$quality, levels = c("Waste(completeness<50% & contamination>10%)",
                                                      "LQ-II(completeness<50% & contamination 5-10%))",
                                                      "LQ-I(completeness<50% & contamination<=5%)",
                                                      "MQ-II(completeness>=50% & contamination 5-10%)",
                                                      "MQ-I(completeness>=50% & contamination<=5%)",
                                                      "HQ-II(completeness>=90% & contamination 5-10%)",
                                                      "HQ-I(completeness>=90% & contamination<=5%)"))
#rm bin_tab vamb_split
bin_tab <- subset(bin_tab, strategy!="sinASScoBIN")
bin_tab$strategy[bin_tab$strategy=="sinASScoBIN(vamb_split)"] <- "sinASScoBIN"
p <- ggplot(bin_tab) +
  geom_bar(stat='count', aes(x=type,fill=quality)) +
  facet_wrap(.~ strategy,ncol = 1,scales = 'free_x')+
  scale_fill_brewer(palette = "Spectral") +
  labs(x="", y="Number of bins") +
  theme_classic() +
  mytheme +
  theme(axis.ticks.x.bottom = element_blank(),
        axis.text.x.bottom = element_text(size = 10),
        strip.text = element_text(size = 10)
        )
p
ggsave("figure/figure_s2.png",p, width = 9, height = 6)
