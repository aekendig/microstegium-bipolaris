##### info ####

# file: gs_survival_figure_2019_density_exp
# author: Amy Kendig
# date last edited: 10/23/20
# goal: figure of regression fits


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(cowplot)

# import models
load("output/elymus_seedling_gs_survival_figure_2019_density_exp.rda")
load("output/elymus_adult_gs_survival_figure_2019_density_exp.rda")
load("output/microstegium_gs_survival_figure_2019_density_exp.rda")


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

# plant group labels
plant_groups <- c("E. virginicus adult",
                  "E. virginicus seedling",
                  "M. vimineum")
names(plant_groups) <- c("Ev adult", "Ev seedling", "Mv seedling")

# figures
mvGsSurvD2Fig <- mvGsSurvD2Plot +
  ggtitle(expression(italic("M. vimineum"))) +
  scale_fill_viridis_d(name = "Disease treatment", direction = -1, end = 0.6) +
  scale_color_viridis_d(name = "Disease treatment", direction = -1, end = 0.6) +
  facet_grid(~ background, scales = "free_x", switch = "x", labeller = labeller(background = plant_groups)) +
  ylab("Proportion surviving") +
  xlab(expression(paste("Density (plants ", m^-1, ")", sep = ""))) +
  fig_theme +
  theme(legend.position = c(0.87, 0.27),
        strip.text = element_blank(),
        axis.title.x = element_blank())

evSGsSurvD2Fig <- evSGsSurvD2Plot +
  ggtitle(expression(paste(italic("E. virginicus"), " seedling", sep = ""))) +
  scale_fill_viridis_d(name = "Disease treatment", direction = -1, end = 0.6) +
  scale_color_viridis_d(name = "Disease treatment", direction = -1, end = 0.6) +
  facet_grid(~ background, scales = "free_x", switch = "x", labeller = labeller(background = plant_groups)) +
  ylab("Proportion surviving") +
  xlab(expression(paste("Density (plants ", m^-1, ")", sep = ""))) +
  fig_theme +
  theme(legend.position = "none",
        strip.text = element_blank(),
        axis.title.x = element_blank())

evAGsSurvD2Fig <- evAGsSurvD2Plot +
  ggtitle(expression(paste(italic("E. virginicus"), " adult", sep = ""))) +
  scale_fill_viridis_d(name = "Disease treatment", direction = -1, end = 0.6) +
  scale_color_viridis_d(name = "Disease treatment", direction = -1, end = 0.6) +
  facet_grid(~ background, scales = "free_x", switch = "x", labeller = labeller(background = plant_groups)) +
  ylab("Proportion surviving") +
  xlab(expression(paste("Density (plants ", m^-1, ")", sep = ""))) +
  fig_theme +
  theme(legend.position = "none")

pdf("output/gs_survival_figure_2019_density_exp.pdf", width = 6, height = 7)
plot_grid(mvGsSurvD2Fig, evSGsSurvD2Fig, evAGsSurvD2Fig,
          rel_heights = c(0.8, 0.8, 1),
          nrow = 3,
          label_size = 10,
          labels = LETTERS[1:3])
dev.off()