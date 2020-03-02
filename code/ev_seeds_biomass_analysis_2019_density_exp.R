##### info ####

# file: ev_seeds_biomass_analysis_2019_density_exp
# author: Amy Kendig
# date last edited: 2/13/20
# goal: evaluate the effects of density treatments and environmental covariates on the biomass and seed production of Elymus


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(cowplot)

# import data
spike <- read_csv("./data/ev_spikelets_2019_density_litter_exp.csv")
plots <- read_csv("./data/plot_treatments_for_figures_2018_2019_density_exp.csv")


#### edit data ####

# use conversion from 2018 for now (ev_seeds_data_processing_2018.R)
# seeds per mg spikelet weight
conv <- 0.07431734

# add location columns
# estimate seed number
# add plot data
elymus <- spike %>%
  mutate(treatment = str_extract(plot, "[aA-zZ]+") %>% recode("W" = "water", "F" = "fungicide"),
         plot = str_extract(plot, "[0-9]+") %>% as.numeric(),
         ID = str_remove(plant, "EV"),
         age = case_when(ID == "A" ~ "adult",
                         TRUE ~ "seedling"),
         seeds = spikelet_weight.g * 1000 * conv) %>%
  left_join(plots)  %>%
  mutate(background_density_tot = case_when(background == "Ev adult" ~ background_density + 1,
                                            TRUE ~ background_density + 3))


#### output intermediate data ####
write_csv(elymus, "intermediate-data/ev_seeds_biomass_covariates_2019_density_exp.csv")

