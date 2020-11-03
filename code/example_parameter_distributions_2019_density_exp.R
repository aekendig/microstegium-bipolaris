##### info ####

# file: example_parameter_distributions_2019_density_exp
# author: Amy Kendig
# date last edited: 10/28/20
# goal: simulate populations with parameters derived from data


#### set-up ####

# clear workspace
rm(list = ls())

# number of samples
n_samps <- 300

# load packages
library(tidyverse)
library(brms)
library(tidybayes)

# load scripts
source("code/microstegium_seed_production_parameter_2019_density_exp.R")


#### simulations ####

# create data frame
dat <- tibble(iter = 1:n_samps) %>%
  expand_grid(disease = c(0, 1)) %>%
  rowwise() %>%
  mutate(seeds = Y_A_fun(disease, 0, 0, 0, iter)) %>%
  ungroup()


#### figures ####
pdf("output/example_establishment_distributions_fungicide_2019_density_experiment.pdf", width = 2.5, height = 2.5)
dat %>%
  filter(disease == 0) %>%
  ggplot(aes(y = disease, x = seeds)) +
  stat_halfeyeh(fill = "#22A884", point_color = NA, interval_color = NA) +
  xlab("Seeds produced") +

  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11, color = "black"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 13),
        axis.ticks.y = element_blank())
dev.off()

pdf("output/example_establishment_distributions_water_2019_density_experiment.pdf", width = 2.5, height = 2.5)
dat %>%
  filter(disease == 1) %>%
  ggplot(aes(y = disease, x = seeds)) +
  stat_halfeyeh(fill = "#430A54", point_color = NA, interval_color = NA) +
  xlab("Seeds produced") +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 11, color = "black"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 13),
        axis.ticks.y = element_blank())
dev.off()