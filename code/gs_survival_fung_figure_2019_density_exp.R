##### info ####

# file: gs_survival_fung_figure_2019_density_exp
# author: Amy Kendig
# date last edited: 10/27/20
# goal: figure of regression fits to fungicide treatment


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(cowplot)

# import figures
load("output/elymus_seedling_gs_survival_bh_figure_2019_density_exp.rda")
load("output/elymus_adult_gs_survival_fung_figure_2019_density_exp.rda")
load("output/microstegium_gs_survival_bh_figure_2019_density_exp.rda")


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
                   plot.title = element_text(size = 10, hjust = 0.5))

# figures
mvGsSurvD2BhFig <- mvGsSurvD2BhPlot +
  ggtitle(expression(italic("M. vimineum"))) +
  scale_color_viridis_d(direction = -1, end = 0.6) +
  ylab("Proportion established") +
  xlab("Disease treatment") +
  coord_cartesian(ylim = c(0.5, 1)) +
  fig_theme +
  theme(legend.position = "none")

evSGsSurvD2BhFig <- evSGsSurvD2BhPlot +
  ggtitle(expression(paste(italic("E. virginicus"), " seedling", sep = ""))) +
  scale_color_viridis_d(direction = -1, end = 0.6) +
  ylab("Proportion established") +
  xlab("Disease treatment") +
  coord_cartesian(ylim = c(0.5, 1)) +
  fig_theme +
  theme(legend.position = "none")

evAGsSurvD2FuFig <- evAGsSurvD2FuPlot +
  ggtitle(expression(paste(italic("E. virginicus"), " adult", sep = ""))) +
  scale_color_viridis_d(direction = -1, end = 0.6) +
  ylab("Proportion alive") +
  xlab("Disease treatment") +
  coord_cartesian(ylim = c(0.5, 1)) +
  fig_theme +
  theme(legend.position = "none")

pdf("output/gs_survival_fung_figure_2019_density_exp.pdf", width = 6, height = 2.5)
plot_grid(mvGsSurvD2BhFig, evSGsSurvD2BhFig, evAGsSurvD2FuFig,
          nrow = 1,
          label_size = 10,
          labels = LETTERS[1:3])
dev.off()