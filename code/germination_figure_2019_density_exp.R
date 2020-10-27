##### info ####

# file: germination_figure_2019_density_exp
# author: Amy Kendig
# date last edited: 10/25/20
# goal: figure of regression fits


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(cowplot)

# import figures
load("output/elymus_germination_figure_2019_density_exp.rda")
load("output/microstegium_germination_figure_2018_density_exp.rda")


#### visualize ####

# general theme
fig_theme <- theme(panel.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   axis.text = element_text(size = 9, color = "black"),
                   axis.title = element_text(size = 10),
                   legend.text = element_text(size = 9),
                   legend.title = element_text(size = 9),
                   strip.background = element_blank(),
                   strip.text = element_text(size = 10),
                   strip.placement = "outside",
                   plot.title = element_text(size = 10, hjust = 0.5, face = "italic"))

# figures
evGermD2Fig <- evGermD2Plot +
  ggtitle("E. virginicus") +
  scale_fill_viridis_d(name = "Disease treatment", direction = -1, end = 0.6) +
  scale_color_viridis_d(name = "Disease treatment", direction = -1, end = 0.6) +
  ylab("Proportion germinated") +
  xlab("Disease treatment") +
  fig_theme +
  theme(legend.position = "none")

mvGermD1Fig <- mvGermD1Plot +
  ggtitle("M. vimineum") +
  scale_fill_viridis_d(name = "Disease treatment", direction = -1, end = 0.6) +
  scale_color_viridis_d(name = "Disease treatment", direction = -1, end = 0.6) +
  ylab("Proportion germinated") +
  xlab("Disease treatment") +
  fig_theme +
  theme(legend.position = "none")

pdf("output/germination_figure_2019_density_exp.pdf", width = 4.5, height = 2.5)
plot_grid(mvGermD1Fig,evGermD2Fig,
          nrow = 1,
          label_size = 10,
          labels = LETTERS[1:2])
dev.off()