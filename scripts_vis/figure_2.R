# load package and data import
##################################################
library(ggplot2)
library(ggbeeswarm)
library(extrafont)
library(ggpubr)
library(readxl)
library(reshape2)
library(rstatix)
library(dplyr)
library(lme4)
library(lsmeans)
require(data.table)

getwd()
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

# import all the data sheet
# mapping for sequencing
mapping <- read.table("data/Metadata.tsv", sep="\t", header = T, row.names = 1)
# input for all the excels
sheets <- excel_sheets("data/phenotype.xlsx")

list.data <- list()
for (i in sheets) {
  list.data[[i]] <- read_excel("data/phenotype.xlsx",sheet=i)
}
############################
# Figure.2A, NEC incidence, chi-squre
data_tab <- data.frame(list.data$Clinic)
data_tab <- subset(data_tab, select = c(ID,Group, NEC))
data_tab_f <- data_tab[data_tab$ID %in% mapping$FMT_ID,]

data_tab_f$Group <- factor(data_tab_f$Group, levels = c("CON","NEW", "OLD"), labels = c("CON","FMT1","FMT2"))
data_tab_long <- melt(data_tab_f, id.vars = c("ID","Group"))
summary(data_tab_long$Group)

data_table_long<- data_tab_long %>%
  group_by(variable) %>%
  mutate(NEC_incidence = case_when(Group=="CON" ~ value/as.numeric(summary(data_tab_long$Group)[1])*100,
                                   Group=="FMT1" ~ value/as.numeric(summary(data_tab_long$Group)[2])*100,
                                   Group=="FMT2" ~ value/as.numeric(summary(data_tab_long$Group)[3])*100))

data_tab_long$value <-factor(data_tab_long$value)
data_tab_long$Group <- factor(data_tab_long$Group)

tab_chiqs <- data_tab_f %>%
  select(Group, NEC) %>%
  table() 
tab_chiqs<- t(tab_chiqs)
tab_chiqs
stat.test <- pairwise_fisher_test(tab_chiqs, p.adjust.method = "BH")%>%
             filter(p.adj<0.05) %>%
             mutate(y.position=c(43))
stat.test

p_NEC <- ggplot(data_table_long, aes(x= Group, y=NEC_incidence)) +
  geom_bar(stat = "identity", aes(fill=Group)) +
  labs(x="", y="NEC Incidence (%)", title=paste("")) +
  scale_fill_manual(values = cols) +
  stat_pvalue_manual(stat.test,label = "p.adj.signif", tip.length = 0,size = 6)+
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5, family = "Arial", size=12),
        legend.title = element_blank(),
        legend.position ="none",
        legend.text = element_text(family = "Arial", size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(family = "Arial", size = 10, vjust = 1),
        axis.text = element_text(family = 'Arial', size = 10, color = "black"),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())
p_NEC
#############################################################
#Figure.2B NEC severity lm, litter as random effect, #tukey_hsd
data_tab <- data.frame(list.data$Clinic)
data_tab <- subset(data_tab, select = c(ID,Group, NEC_max,Litter))
data_tab_f <- data_tab[data_tab$ID %in% mapping$FMT_ID,]
data_tab_f$Group <- factor(data_tab_f$Group, levels = c("CON","NEW", "OLD"), labels = c("CON","FMT1","FMT2"))
data_tab_long <- melt(data_tab_f, id.vars = c("ID","Group", "Litter"))

# lm function to list litter as random effect and adjusted by tukey_hsd
lm_fun <- function(data_tab_long) {
lmer.fit <- lmer(value ~ Group + (1|Litter), data = data_tab_long) 
sum_tukey <- lsmeans(lmer.fit, pairwise~Group, adjust="tukey")
tab <- as.data.frame(sum_tukey$contrasts)

groups <- setDT(tstrsplit(as.character(tab$contrast), "\ -\ ", fixed=TRUE))[]

stat.test <- tibble(group1=groups$V1, group2=groups$V2, p.adj=tab$p.value)
stat.test$p.adj.signif <- case_when(stat.test$p.adj >= 0.05 ~ 'na',
                                    stat.test$p.adj < 0.05 ~ "*",
                                    stat.test$p.adj < 0.01 ~ "**",
                                    stat.test$p.adj < 0.001 ~ "***")
return(stat.test)
}

stat.test <- lm_fun(data_tab_long)

stat.test <- stat.test %>%
    filter(p.adj<0.05) %>%
    mutate(y.position=c(seq(6.3, length.out = 1)))

stat.test

index <- match(data_tab_long$ID,as.character(mapping$FMT_ID))
data_tab_long$NEC <- mapping$NEC[index]
data_tab_long$NEC <- factor(data_tab_long$NEC, levels = c(0,1), labels = c("No","Yes"))


p_NECseverity <- ggplot(data_tab_long, aes(x= Group, y= value)) +
  geom_violin(aes(fill=Group)) +
  geom_beeswarm(cex = 2, size=2, aes(color=NEC)) +
  geom_hline(yintercept=3.5, linetype="dashed", color = "black") +
  scale_color_manual(values = c("black", "gray")) +
  scale_fill_manual(values = cols) +
  labs(x="", y="NEC max scores", title=paste("")) +
  stat_pvalue_manual(stat.test,label = "p.adj.signif", tip.length = 0,size = 6)+
  scale_y_continuous(breaks = seq(1,6,1)) +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5, family = "Arial", size=12),
        legend.text = element_text(family = "Arial", size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(family = "Arial", size = 10, vjust = 1),
        axis.text = element_text(family = 'Arial', size = 10, color = "black"),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())
p_NECseverity

#################
# Figure.2C lm_ratio
data_tab <- data.frame(list.data$LM_ratio)
data_tab <- subset(data_tab, select = c(ID,Group, LM_ratio, Litter))
data_tab_f <- data_tab[data_tab$ID %in% mapping$FMT_ID,]
data_tab_f <- na.omit(data_tab_f)
data_tab_f$Group <- factor(data_tab_f$Group, levels = c("CON","NEW", "OLD"), labels = c("CON","FMT1","FMT2"))
data_tab_long <- melt(data_tab_f, id.vars = c("ID","Group", "Litter"))

stat.test <- lm_fun(data_tab_long)

stat.test <- stat.test %>%
  filter(p.adj<0.05) %>%
  mutate(y.position=c(seq(58.5, length.out = 1)))
stat.test 

index <- match(data_tab_long$ID,as.character(mapping$FMT_ID))
data_tab_long$NEC <- mapping$NEC[index]
data_tab_long$NEC <- factor(data_tab_long$NEC, levels = c(0,1), labels = c("No","Yes"))


p_LM <- ggplot(data_tab_long, aes(x= Group, y= value)) +
  stat_boxplot(geom ='errorbar', linetype=1, width=0.5) + 
  geom_boxplot(outlier.shape = NA, aes(fill=Group)) +
  scale_color_manual(values = c("black", "gray")) +
  scale_fill_manual(values = cols) +
  labs(x="", y="LM ratio (%)", title=paste("")) +
  stat_pvalue_manual(stat.test,label = "p.adj.signif", tip.length = 0,size = 6)+
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5, family = "Arial", size=12),
        legend.title = element_blank(),
        legend.position ="none",
        legend.text = element_text(family = "Arial", size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(family = "Arial", size = 10, vjust = 1),
        axis.text = element_text(family = 'Arial', size = 10, color = "black"),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())
p_LM
# Figure.2D daily weight gain
data_tab <- data.frame(list.data$Body_dimensions)
data_tab <- subset(data_tab, select = c(ID,Group, Daily_weight_grain, Litter))
data_tab_f <- data_tab[data_tab$ID %in% mapping$FMT_ID,]
data_tab_f <- na.omit(data_tab_f)

data_tab_f$Group <- factor(data_tab_f$Group, levels = c("CON","NEW", "OLD"), labels = c("CON","FMT1","FMT2"))
data_tab_long <- melt(data_tab_f, id.vars = c("ID","Group", "Litter"))

stat.test <- lm_fun(data_tab_long)

stat.test <- stat.test %>%
  filter(p.adj<0.05) %>%
  mutate(y.position=c(seq(43, length.out = 1)),p.adj = format.pval(p.adj, digits = 2))
stat.test 
summary(data_tab_long)

index <- match(data_tab_long$ID,as.character(mapping$FMT_ID))
data_tab_long$NEC <- mapping$NEC[index]
data_tab_long$NEC <- factor(data_tab_long$NEC, levels = c(0,1), labels = c("No","Yes"))

p_bw <- ggplot(data_tab_long, aes(x= Group, y= value)) +
  stat_boxplot(geom ='errorbar', linetype=1, width=0.5) + 
  geom_boxplot(outlier.shape = NA, aes(fill=Group)) +
  geom_point(position = position_jitter(w=0.1, h=0),size=2, aes(color=NEC)) +
  scale_color_manual(values = c("black", "gray")) +
  scale_fill_manual(values = cols) +
  labs(x="", y="Growth rate (g/kg/d)", title=paste("")) +
  stat_pvalue_manual(stat.test,label = "p.adj.signif", tip.length = 0,size = 6)+
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5, family = "Arial", size=12),
        legend.title = element_blank(),
        legend.position ="none",
        legend.text = element_text(family = "Arial", size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(family = "Arial", size = 10, vjust = 1),
        axis.text = element_text(family = 'Arial', size = 10, color = "black"),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())
p_bw
###############################
#Figure.2E BBenzyme ##P_DPPIV
data_tab <- data.frame(list.data$BB_enzymes)
data_tab <- subset(data_tab, select = c(ID,Group, P_DPPIV, Litter))
data_tab_f <- data_tab[data_tab$ID %in% mapping$FMT_ID,]
data_tab_f <- na.omit(data_tab_f)

data_tab_f$Group <- factor(data_tab_f$Group, levels = c("CON","NEW", "OLD"), labels = c("CON","FMT1","FMT2"))
data_tab_long <- melt(data_tab_f, id.vars = c("ID","Group", "Litter"))

stat.test <- lm_fun(data_tab_long)

stat.test <- stat.test %>%
  filter(p.adj<0.06) %>%
  mutate(y.position=c(seq(3.8, length.out = 1)),p.adj = format.pval(p.adj, digits = 2))
stat.test 

summary(data_tab_long)
index <- match(data_tab_long$ID,as.character(mapping$FMT_ID))
data_tab_long$NEC <- mapping$NEC[index]
data_tab_long$NEC <- factor(data_tab_long$NEC, levels = c(0,1), labels = c("No","Yes"))

p_P_DPPIV <- ggplot(data_tab_long, aes(x= Group, y= value)) +
  stat_boxplot(geom ='errorbar', linetype=1, width=0.5) + 
  geom_boxplot(outlier.shape = NA, aes(fill=Group)) +
  geom_point(position = position_jitter(w=0.1, h=0),size=2, aes(color=NEC)) +
  scale_color_manual(values = c("black", "gray")) +
  scale_fill_manual(values = cols) +
  labs(x="", y="Proximal SI DPPIV (U/g)", title=paste("")) +
  stat_pvalue_manual(stat.test,label = "p.adj={p.adj}", tip.length = 0,size = 3)+
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5, family = "Arial", size=12),
        legend.title = element_blank(),
        legend.position ="none",
        legend.text = element_text(family = "Arial", size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(family = "Arial", size = 10, vjust = 1),
        axis.text = element_text(family = 'Arial', size = 10, color = "black"),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())
p_P_DPPIV

###############################
#Figure.2F BBenzyme ##M_maltase
data_tab <- data.frame(list.data$BB_enzymes)
data_tab <- subset(data_tab, select = c(ID,Group, M_maltase, Litter))
data_tab_f <- data_tab[data_tab$ID %in% mapping$FMT_ID,]
data_tab_f <- na.omit(data_tab_f)
data_tab_f$Group <- factor(data_tab_f$Group, levels = c("CON","NEW", "OLD"), labels = c("CON","FMT1","FMT2"))
data_tab_long <- melt(data_tab_f, id.vars = c("ID","Group", "Litter"))

stat.test <- lm_fun(data_tab_long)

stat.test <- stat.test %>%
  filter(p.adj<0.05) %>%
  mutate(y.position=c(seq(9, length.out = 1)))
stat.test 
summary(data_tab_long)

index <- match(data_tab_long$ID,as.character(mapping$FMT_ID))
data_tab_long$NEC <- mapping$NEC[index]
data_tab_long$NEC <- factor(data_tab_long$NEC, levels = c(0,1), labels = c("No","Yes"))
p_M_maltase <- ggplot(data_tab_long, aes(x= Group, y= value)) +
  stat_boxplot(geom ='errorbar', linetype=1, width=0.5) + 
  geom_boxplot(outlier.shape = NA, aes(fill=Group)) +
  geom_point(position = position_jitter(w=0.1, h=0),size=2, aes(color=NEC)) +
  scale_color_manual(values = c("black", "gray")) +
  scale_fill_manual(values = cols) +
  labs(x="", y="Median SI Maltase (U/g)", title=paste("")) +
  stat_pvalue_manual(stat.test,label = "p.adj.signif", tip.length = 0,size = 6)+
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5, family = "Arial", size=12), 
        legend.title = element_blank(),
        legend.position ="none",
        legend.text = element_text(family = "Arial", size = 10),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(family = "Arial", size = 10, vjust = 1),
        axis.text = element_text(family = 'Arial', size = 10, color = "black"),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())
p_M_maltase

# integrate all the sub-figure together wiith patchwork
library(patchwork)
p_merge <- p_NEC + p_NECseverity + p_LM  +
  p_bw  + p_P_DPPIV + p_M_maltase + plot_annotation(tag_levels = 'a') + plot_layout(guides = 'collect', ncol=3)
p_merge
ggsave("figure/figure_2.png",p_merge, height = 8, width = 13)
#ggsave("figure/figure_2.pdf",p_merge, height = 8, width = 13, device = cairo_pdf)
###########################sessionInfo to keep track
sessionInfo()