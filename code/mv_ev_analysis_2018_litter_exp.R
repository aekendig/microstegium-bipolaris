##### info ####

# file: mv_ev_analysis_2018_litter_exp
# author: Amy Kendig
# date last edited: 2/3/20
# goal: evaluate correlations between Mv and Ev in the litter experiment


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(cowplot)

# import data
micro <- read_csv("intermediate-data/mv_germination_covariates_2018_litter_exp.csv")
elymus <- read_csv("./intermediate-data/ev_damage_seeds_covariates_2018_litter_exp.csv")


#### edit data ####

# combine data
dat <- inner_join(micro, elymus)


#### figures ####

# text sizes
sm_txt = 6
lg_txt = 8

# scatterplots
scat_fig <- ggplot(dat, aes(x = mv_germ_planted_bg_jun, y = seeds)) +
  geom_point(aes(fill = site), shape = 21) +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = lg_txt),
        axis.text = element_text(color = "black", size = sm_txt),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(color = "black", size = lg_txt),
        legend.text = element_text(color = "black", size = sm_txt),
        legend.key.height = unit(0.3, "cm"),
        legend.position = "none") +
  scale_fill_manual(values = c("black", "blue", "purple", "red"), name = "Site")

# June germination vs. seeds
fig_germ_jn_seeds <- scat_fig +
  xlab("Background and planted Mv germinants in June") +
  ylab("Number of Ev seeds")

# July germination vs. seeds
fig_germ_jl_seeds <- scat_fig %+%
  aes(x = mv_germ_planted_bg_jul) +
  xlab("Background and planted Mv germinants in July") +
  ylab("Number of Ev seeds")

# June germination vs. damage
fig_germ_jn_dam <- scat_fig %+%
  aes(y = prop_dam) +
  xlab("Background and planted Mv germinants in June") +
  ylab("Proportion of leaf area\nwith lesions")

# July germination vs. damage
fig_germ_jl_dam <- scat_fig %+%
  aes(x = mv_germ_planted_bg_jul, y = prop_dam) +
  xlab("Background and planted Mv germinants in July") +
  ylab("Proportion of leaf area\nwith lesions")

pdf("./output/mv_germination_ev_raw_2018_litter_exp.pdf", width = 6, height = 4)
plot_grid(fig_germ_jn_seeds, fig_germ_jl_seeds, fig_germ_jn_dam, fig_germ_jl_dam,
          ncol = 2,
          labels = letters[1:4],
          label_size = lg_txt)
dev.off()