library(ggradar)
library(readxl)
library(dplyr)
library(scales)

library(ggplot2)
library(tidyr)
library(extrafont)
library(ggpubr)
loadfonts()
mytheme<- theme(plot.title = element_text(hjust=0.5, family = "Arial", size=12),
                legend.position ="right",
                legend.text = element_text(family = "Arial", size = 8),
                legend.background = element_blank(),
                strip.background = element_blank(),
                strip.placement = "outside",
                strip.text = element_text(family = "Arial", size = 8),
                axis.text = element_text(family = 'Arial', size = 8, color="black"),
                panel.grid = element_blank())
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(5)
#import histological data
# mapping for sequencing
mapping <- read.table("data/Metadata.tsv", sep="\t", header = T, row.names = 1)
head(mapping)
# input for all the excels

sheets <- excel_sheets("data/phenotype.xlsx")
sheets

histo <- read_excel("data/phenotype.xlsx",sheet=6)
clinic <- read_excel("data/phenotype.xlsx",sheet=4)
histo <- subset(histo, ID %in% mapping$FMT_ID)
clinic <- subset(clinic, ID %in% mapping$FMT_ID)
colnames(clinic)
#NEC
#############################
nec_group <- clinic %>%
  select(Group, NEC, NEC_prox, NEC_mid, NEC_dist, NEC_colon, NEC_stomach) %>%
  filter(!is.na(NEC)) %>%
  group_by(Group) %>%
  summarise_if(is.numeric, mean, na.rm =TRUE ) %>%
  as.data.frame()
nec_group

colnames(nec_group) <- c("Group","NEC", "Prox", "Mid", "Dist", "Colon", "Stomach")
nec_group$Group <- c("CON",  "FMT1", "FMT2")

radar_NEC <- nec_group %>%
  select(-NEC) %>%
  mutate_each(rescale, -Group) %>%
  ggradar(plot.title = "NEC", legend.title = "Group",  font.radar = "Arial", group.colours = cols, #, grid.label.size = 5, axis.label.size = 4
          ) +
  theme(plot.title = element_text(hjust=0.5, family = "Arial", size = 14),
        legend.text = element_text(family = "Arial", size = 12),
        legend.title = element_text(family = "Arial", size = 14),
        strip.text = element_text(family = "Arial", size = 12, vjust = 1),
        axis.text = element_text(family = 'Arial', size = 12, color = "black"),
        plot.tag = element_text(family = "Arial", size = 14),
        )
radar_NEC

#Diarrhea
Diarrhea_group <- clinic %>%
  select(Group,NEC, Feces_d1,Feces_d2, Feces_d3, Feces_d4, Feces_d5) %>%
  filter(!is.na(NEC)) %>%
  group_by(Group) %>%
  summarise_if(is.numeric, mean, na.rm =TRUE ) %>%
  as.data.frame()
Diarrhea_group

colnames(Diarrhea_group) <- c("Group","NEC", "D1", "D2", "D3", "D4", "D5")
Diarrhea_group$Group <- c("CON", "FMT1", "FMT2")


radar_Diarrhea <- Diarrhea_group %>%
  select(-NEC) %>%
  mutate_each(rescale, -Group) %>%
  ggradar(plot.title = "Diarrhea",legend.title = "Group", font.radar = "Arial", group.colours = cols,#, grid.label.size = 5, axis.label.size = 4
          ) +
  theme(plot.title = element_text(hjust=0.5, family = "Arial", size = 14),
        legend.text = element_text(family = "Arial", size = 12),
        legend.title = element_text(family = "Arial", size = 14),
        strip.text = element_text(family = "Arial", size = 12, vjust = 1),
        axis.text = element_text(family = 'Arial', size = 12, color = "black"),
        plot.tag = element_text(family = "Arial", size = 14))

radar_Diarrhea

######################################################
# histology
index <- match(histo$ID, clinic$ID)
histo$NEC <- clinic$NEC[index]  
colnames(histo)
histo_group <- histo %>%
  select(Group, NEC,ABPAS_SI_fraction, ABPAS_colon_fraction, CD3_SI_fraction, CD3_colon_fraction, MPO_SI_score, MPO_colon_score,FISH_SI_score) %>%
  filter(!is.na(NEC)) %>%
  group_by(Group) %>%
  summarise_if(is.numeric, mean, na.rm =TRUE ) %>%
  as.data.frame()
histo_group

colnames(histo_group) <- c("Group","NEC", "Goblet cell SI", "Goblet cell Colon", "CD3+ cell SI", "CD3+ cell Colon", "MPO score SI", "MPO score Colon", "FISH SI")
histo_group$Group <- c("CON", "FMT1", "FMT2")

radar_histo <- histo_group %>%
  select(-NEC) %>%
  mutate_each(rescale, -Group) %>%
  ggradar(plot.title = "Histology", legend.title = "Group",  font.radar = "Arial", group.colours = cols#, grid.label.size = 5, axis.label.size = 4
          ) +
  theme(plot.title = element_text(hjust=0.5, family = "Arial", size = 14),
        legend.text = element_text(family = "Arial", size = 12),
        legend.title = element_text(family = "Arial", size = 14),
        strip.text = element_text(family = "Arial", size = 12, vjust = 1),
        axis.text = element_text(family = 'Arial', size = 12, color = "black"),
        plot.tag = element_text(family = "Arial", size = 14))

radar_histo


# all histological data
colnames(histo)
histo_boxplot <- histo %>%
  select(Group, NEC,ABPAS_SI_fraction, ABPAS_colon_fraction, CD3_SI_fraction, CD3_colon_fraction, MPO_SI_score, MPO_colon_score,FISH_SI_score) %>%
  filter(!is.na(NEC))
  
histo_boxplot$NEC <- factor(histo_boxplot$NEC, levels = c(0,1), labels = c("No","Yes"))
histo_boxplot$Group <- factor(histo_boxplot$Group, levels = c("CON","NEW","OLD"), labels = c("CON","FMT1","FMT2"))
colnames(histo_boxplot) <- c("Group","NEC", "Goblet cell SI", "Goblet cell Colon", "CD3+ cell SI", "CD3+ cell Colon", "MPO score SI", "MPO score Colon", "FISH SI")
#MPO score
histo_mpo <- histo_boxplot %>%
  select(Group, NEC, "MPO score SI", "MPO score Colon") %>%
  gather(Condition, Measurement, 3:4, factor_key = TRUE)
histo_mpo$Condition <- factor(histo_mpo$Condition, levels = c("MPO score SI", "MPO score Colon"), labels = c("SI","Colon"))

p_mpo <- ggplot(data = histo_mpo, aes(x=Condition,y=Measurement,fill=Group)) +
  stat_boxplot(geom ='errorbar', linetype=1, width=0.25, position = position_dodge(.85)) + 
  geom_boxplot(position=position_dodge(0.85),outlier.shape = NA) +
  geom_point(size=2,position = position_jitterdodge(dodge.width = 0.9), aes(color=NEC, group = Group)) +
  scale_color_manual(values = c("black", "gray")) +
  scale_fill_manual(values = cols) +
  labs(x="", y="MPO score", title="MPO score") +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5, family = "Arial", size=14),
        legend.text = element_text(family = "Arial", size = 12),
        legend.title = element_text(family = "Arial", size = 14),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(family = "Arial", size = 12, vjust = 1),
        axis.text = element_text(family = 'Arial', size = 12, color = "black"),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())

p_mpo

#CD3 density (fraction)
histo_cd3 <- histo_boxplot %>%
  select(Group, NEC, "CD3+ cell SI", "CD3+ cell Colon") %>%
  gather(Condition, Measurement, 3:4, factor_key = TRUE)
histo_cd3$Condition <- factor(histo_cd3$Condition, levels = c("CD3+ cell SI", "CD3+ cell Colon"), labels = c("SI","Colon"))

p_cd3 <- ggplot(data = histo_cd3, aes(x=Condition,y=Measurement,fill=Group)) +
  stat_boxplot(geom ='errorbar', linetype=1, width=0.25, position = position_dodge(.85)) + 
  geom_boxplot(position=position_dodge(0.85),outlier.shape = NA) +
  geom_point(size=2,position = position_jitterdodge(dodge.width = 0.9), aes(color=NEC, group = Group)) +
  scale_color_manual(values = c("black", "gray")) +
  scale_fill_manual(values = cols) +
  labs(x="", y="CD3+ cell density (%)", title="CD3+ cell density") +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5, family = "Arial", size=14),
        legend.text = element_text(family = "Arial", size = 12),
        legend.title = element_text(family = "Arial", size = 14),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(family = "Arial", size = 12, vjust = 1),
        axis.text = element_text(family = 'Arial', size = 12, color = "black"),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())

p_cd3

#ABPAS density (fraction)
histo_abpas <- histo_boxplot %>%
  select(Group, NEC, "Goblet cell SI", "Goblet cell Colon") %>%
  gather(Condition, Measurement, 3:4, factor_key = TRUE)
histo_abpas$Condition <- factor(histo_abpas$Condition, levels = c("Goblet cell SI", "Goblet cell Colon"), labels = c("SI","Colon"))

p_abpas <- ggplot(data = histo_abpas, aes(x=Condition,y=Measurement,fill=Group)) +
  stat_boxplot(geom ='errorbar', linetype=1, width=0.25, position = position_dodge(.85)) + 
  geom_boxplot(position=position_dodge(0.85),outlier.shape = NA) +
  geom_point(size=2,position = position_jitterdodge(dodge.width = 0.9), aes(color=NEC, group = Group)) +
  scale_color_manual(values = c("black", "gray")) +
  scale_fill_manual(values = cols) +
  labs(x="", y="Goblet cell density (%)", title="Goblet cell density") +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5, family = "Arial", size=14),
        legend.text = element_text(family = "Arial", size = 12),
        legend.title = element_text(family = "Arial", size = 14),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(family = "Arial", size = 12, vjust = 1),
        axis.text = element_text(family = 'Arial', size = 12, color = "black"),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())

p_abpas


#Fish colon is missing
histo_fish <- histo_boxplot %>%
  select(Group, NEC, "FISH SI") %>%
  gather(Condition, Measurement, 3, factor_key = TRUE)
histo_fish$Condition <- factor(histo_fish$Condition, levels = c("FISH SI"), labels = c("SI"))

p_fish <- ggplot(data = histo_fish, aes(x=Condition,y=Measurement,fill=Group)) +
  stat_boxplot(geom ='errorbar', linetype=1, width=0.25, position = position_dodge(.85)) + 
  geom_boxplot(position=position_dodge(0.85),outlier.shape = NA) +
  geom_point(size=2,position = position_jitterdodge(dodge.width = 0.9), aes(color=NEC, group = Group)) +
  scale_color_manual(values = c("black", "gray")) +
  scale_fill_manual(values = cols) +
  labs(x="", y="FISH score", title="FISH score") +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5, family = "Arial", size=14),
        legend.text = element_text(family = "Arial", size = 12),
        legend.title = element_text(family = "Arial", size = 14),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(family = "Arial", size = 12, vjust = 1),
        axis.text = element_text(family = 'Arial', size = 12, color = "black"),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())

p_fish

######################################################################################
#nec correlation and histology
NEC_tab <- clinic %>%
  select(ID, NEC, NEC_max, NEC_severity)


histo_tab <- histo %>%
  select(ID,  ABPAS_SI_fraction, ABPAS_colon_fraction, MPO_SI_score, MPO_colon_score,CD3_SI_fraction, CD3_colon_fraction,MPO_SI_score,MPO_colon_score,FISH_SI_score)

tab <- merge(NEC_tab, histo_tab, by="ID")
# add group information
index <- match(tab$ID, mapping$FMT_ID)
tab$Group <- as.character(mapping$Group[index])

tab <- tab %>%
  filter(ID %in% mapping$FMT_ID)
colnames(tab)

#transform to long-format
tab_long <- tab %>%
  gather(histo_category,x, 5:11, factor_key = TRUE) %>%
  gather(nec_category,y,3:4, factor_key = TRUE) %>%
  group_by(nec_category,histo_category) %>%
  mutate_each(rescale,-c(ID, NEC,Group))

#Rename facet labels
tab_long$NEC <- factor(tab_long$NEC, levels = c(0,1), labels = c("No","Yes"))
tab_long$Group <- factor(tab_long$Group, levels = c("CON","FMT1","FMT2"), labels = c("CON","FMT1","FMT2"))
tab_long$nec_category <- factor(tab_long$nec_category, levels = c("NEC_severity", "NEC_max"), labels = c("NEC severity", "NEC max"))
tab_long$histo_category <- factor(tab_long$histo_category, levels = c("ABPAS_SI_fraction","ABPAS_colon_fraction","CD3_SI_fraction", "CD3_colon_fraction", "MPO_SI_score", "MPO_colon_score", "FISH_SI_score"),
                                  labels = c("Goblet cell density SI", "Goblet cell density Colon", "CD3+ cell density SI", "CD3+ cell density Colon", "MPO score SI", "MPO score Colon", "FISH score SI"))
#exclude the goblet cell density
tab_long_f <- subset(tab_long, !histo_category %in% c("Goblet cell density SI","Goblet cell density Colon"))
library(ggpmisc)
myformula <- y ~ x
p_lm <- ggplot(tab_long_f, aes(y=y, x=x)) +
  geom_jitter(size=2, aes(color=Group, shape=NEC)) +
  #add linear regression
  stat_smooth(method = 'lm',formula = myformula, color="black") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "*`,`~")), size=rel(2.5), coef.digits = 2, rr.digits = 2, p.digits = 2,
               formula = myformula, parse = TRUE)  +
  facet_grid(nec_category ~ histo_category) +
  scale_color_manual(values = cols) +
  scale_y_continuous(breaks=c(0,0.5,1)) +
  scale_x_continuous(breaks=c(0,0.5,1)) +
  labs(x="Scaled histological data", y="Scaled NEC scores", title="") +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5, family = "Arial", size=14),
        legend.text = element_text(family = "Arial", size = 12),
        legend.title = element_text(family = "Arial", size = 14),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(family = "Arial", size = 12, vjust = 1),
        axis.text = element_text(family = 'Arial', size = 12, color = "black"),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())

p_lm  

library(patchwork)
q_1 <- (radar_NEC + radar_Diarrhea) + plot_layout(guides = "collect")
q_1

q_2 <- (p_abpas +  p_cd3) / (p_mpo + (p_fish + plot_spacer())) + plot_layout(guides = "collect")
q_2

q <- q_1 / q_2 / p_lm + plot_annotation(tag_levels = "A")
q

ggsave("figure/figure_s3.png", height = 12, width = 14)
