library(patchwork)
source("scripts_vis/phyloseq_figure_unoise.R")
source("scripts_vis/kaiju_beta_nr_euk_species.R")
source("scripts_vis/dbrda_kaiju.R")

# Figure. 3A
p_alpha_16s2 <- p_alpha_16s + labs(tag = "A") + theme(plot.tag.position = "topleft")
# Figure. 3B
p_alpha_all2 <- p_alpha_all + labs(tag = "B") + theme(plot.tag.position = "topleft")
# Figure.3C
p_wunifrac_16s2 <- p_wunifrac_16s + labs(tag = "C") + theme(plot.tag.position = "topleft")
# Figure.3D
p_bray_all2 <- p_bray_all + labs(tag = "D") + theme(plot.tag.position = "topleft")
# Figure.3E
p_dbrda2 <- p_dbrda+ labs(tag = "E") + theme(plot.tag.position = "topleft")
p_dbrda3 <- p_dbrda2 + p_bar  + plot_layout(design = layout_man)

#multi-panel figure
p_merged <- p_alpha_16s2  + p_alpha_all2 +
  p_wunifrac_16s2 + p_bray_all2 + plot_layout(ncol = 2)

p_merged / p_dbrda3

ggsave("figure/figure_3.png", height = 12, width=10)

sessionInfo()