##### info ####

# file: biomass_figure_2019_greenhouse_exp
# author: Amy Kendig
# date last edited: 10/27/20
# goal: figure of raw data


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(cowplot)

# import figures
load("output/elymus_biomass_figure_2019_greenhouse_exp.rda")
load("output/microstegium_biomass_figure_2019_greenhouse_exp.rda")


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
mvBioGHFig <- mvFungPlot +
  ggtitle(expression(italic("M. vimineum"))) +
  scale_color_viridis_d(name = "Disease treatment", direction = -1, end = 0.6) +
  ylab("log(Biomass (g))") +
  xlab("Disease treatment") +
  fig_theme +
  theme(legend.position = "none")

evBioGHFig <- evFungPlot +
  ggtitle(expression(italic("E. virginicus"))) +
  scale_color_viridis_d(name = "Disease treatment", direction = -1, end = 0.6) +
  ylab("log(Biomass (g))") +
  fig_theme +
  theme(legend.position = "none")

pdf("output/biomass_figure_2019_greenhouse_exp.pdf", width = 4.5, height = 2.5)
plot_grid(mvBioGHFig, evBioGHFig,
          nrow = 1,
          label_size = 10,
          labels = LETTERS[1:2])
dev.off()