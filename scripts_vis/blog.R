library(patchwork)
library(showtext)
font_add(family = "Arial",
 regular = "/usr/share/fonts/truetype/msttcorefonts/Arial.ttf")
font_add(family = "Arial Italic",
 regular = "/usr/share/fonts/truetype/msttcorefonts/Arial_Italic.ttf")
showtext_auto()
source("scripts_vis/kaiju_beta_nr_euk_species.R")
source("scripts_vis/dbrda_kaiju.R")
source("scripts_vis/figure_s1_s2.R")
source("scripts_vis/figure_5.R")
source("scripts_vis/figure_7A.R")

dir.create("figure/blog", showWarnings = FALSE)
# 1 radar + dbrda
q_11 <- (radar_NEC + labs(tag = "a") + radar_Diarrhea)*theme(legend.position = "none") + 
theme(plot.tag.position = "topleft")

p_dbrda1 <- ggplot(db_rda.site, aes(CAP1, CAP2)) +
  geom_point(aes(color = group#,shape=SOW
                 ),size=4) +
  scale_color_manual(values = c(cols[1], cols[2], cols[3])) +
  xlim(-1,1.5) +
  ylim(-1,1.2) +
  theme_classic() +
  mytheme +
  labs(x = rda1_exp, y = rda2_exp) +
  guides(color=guide_legend(title="FMT")) +
  geom_vline(xintercept = 0, color = 'gray', size = 0.5) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
  geom_segment(data = db_rda.spe.core, aes(x = 0, y = 0, xend = CAP1*0.4, yend = CAP2*0.4), 
               arrow = arrow(length = unit(0.2, 'cm')), size = 0.3, color = 'blue') +
  geom_text(data = db_rda.spe.core, aes(CAP1*0.43, CAP2*0.46, label=name), color='black', size=3, family = "Arial Italic")

p_dbrda2 <- p_dbrda1 + labs(tag = "b") + theme(plot.tag.position = "topleft")
layout_man <- c(
  area(t = 1, b = 4, l = 1, r = 6),
  area(t = 1, b = 1, l = 4, r = 6)
)
p_dbrda3 <- p_dbrda2 + p_bar  + plot_layout(design = layout_man)

q_11 / p_dbrda3 +
  plot_layout(guides = "collect")
ggsave("figure/blog/1donor.pdf", height = 9, width = 9)

# 2 lactobacilli abundance/replication/tree
p_11 <- p_16s_lac_box + labs(tag = "a") + p_kaiju_lac_box + theme(plot.tag.position = "topleft") 
p_22 <- p_irep_lac_mean + labs(tag = "b") + theme(plot.tag.position = "topleft")

p_31 <-  (tree_lre + labs(tag = "c") + tree_lcr) * 
theme(plot.title = element_text(family = "Arial Italic")) + theme(plot.tag.position = "topleft") + 
  plot_layout(guides = "collect")

(p_11 + p_22 + plot_layout(guides = "collect")) / p_31
ggsave("figure/blog/2lactobacilli.pdf", height = 6, width=8)
