library(patchwork)
source("scripts_vis/phyloseq_figure_unoise.R")
source("scripts_vis/kaiju_beta_nr_euk_species.R")
# Figure.S5, based on qualitative dissimilarity metrics

p_merged <- (p_uifrac_16s | p_effectsize_unifrac)/ (p_jaccard| p_effectsize_jaccard)
p_merged + plot_annotation(tag_levels = 'a')
ggsave("figure/figure_s4.png", height = 6, width = 10)
