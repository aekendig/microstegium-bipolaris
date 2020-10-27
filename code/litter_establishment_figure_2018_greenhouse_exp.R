##### info ####

# file: litter_establishment_figure_2018_greenhouse_exp
# author: Amy Kendig
# date last edited: 10/23/20
# goal: figure of regression fits


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(cowplot)
library(viridis)

# import figures
load("output/microstegium_litter_establishment_figure_2018_greenhouse_exp.rda")
load("output/elymus_litter_establishment_figure_2018_greenhouse_exp.rda")


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
mvLitEstFig <- mvLitEstPlot +
  ggtitle(expression(italic("M. vimineum"))) +
  coord_cartesian(ylim = c(0.67, 0.87)) +
  ylab("Proportion established") +
  xlab(expression(paste("Litter (g ", m^-2, ")", sep = ""))) +
  fig_theme

evLitEstFig <- evLitEstPlot +
  ggtitle(expression(italic("E. virginicus"))) +
  coord_cartesian(ylim = c(0.67, 0.87)) +
  ylab("Proportion established") +
  xlab(expression(paste("Litter (g ", m^-2, ")", sep = ""))) +
  fig_theme

pdf("output/litter_establishment_figure_2018_greenhouse_exp.pdf", width = 4.5, height = 2.5)
plot_grid(mvLitEstFig, evLitEstFig,
          nrow = 1,
          label_size = 10,
          labels = LETTERS[1:2])
dev.off()