##### info ####

# file: covariate_data_processing_2018_litter_exp
# author: Amy Kendig
# date last edited: 1/7/20
# goal: create a dataset of covariates


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(GGally)

# import data
soil <- read_csv("./data/soil_moisture_jun_2018_litter_exp.csv")
trees <- read_csv("./data/canopy_cover_trees_jun_2018_density_litter_exp.csv")
canopy <- read_csv("./data/canopy_cover_oct_2018_litter_exp.csv")
plots <- read_csv("./data/plot_treatments_2018_litter_exp.csv")


#### edit data ####

# use proportion
soil2 <- soil %>%
  mutate(soil_moiture.prop = soil_moisture.vwc / 100,
         soil_moisture.centered = soil_moiture.prop - mean(soil_moiture.prop)) %>%
  select(site, plot, soil_moiture.prop, soil_moisture.centered)

# look at canopy notes
unique(canopy$processing_notes)

#  proportion canopy cover
canopy2 <- canopy %>%
  mutate(canopy_cover.prop = case_when(type == "o" ~ 1 - ((count * 1.04) / 100),
                                       type == "c" ~ (count * 1.04) / 100),
         canopy_cover.centered = canopy_cover.prop - mean(canopy_cover.prop)) %>%
  select(site, plot, canopy_cover.prop, canopy_cover.centered)
  
# site level proportion canopy cover, 
canopy_site <- trees %>%
  filter(experiment == "litter" & !is.na(count)) %>%
  mutate(canopy_cover.prop = case_when(type == "o" ~ 1 - ((count * 1.04) / 100),
                                       type == "c" ~ (count * 1.04) / 100))  %>%
  group_by(site) %>%
  summarise(canopy_cover.prop = mean(canopy_cover.prop)) %>%
  ungroup() %>%
  mutate(canopy_cover.centered = canopy_cover.prop - mean(canopy_cover.prop))

# site level basal tree stand area
stand_site <- trees %>%
  filter(experiment == "litter" & !is.na(trees)) %>%
  mutate(stand_area.m2ha = trees * 10)  %>%
  group_by(site) %>%
  summarise(stand_area.m2ha = mean(stand_area.m2ha)) %>%
  ungroup() %>%
  mutate(stand_area.scaled = (stand_area.m2ha - mean(stand_area.m2ha)) / sd(stand_area.m2ha))

# merge plot data
plot_dat <- full_join(plots, soil2) %>%
  full_join(canopy2)

# merge site data
site_dat <- full_join(canopy_site, stand_site)


#### visualize data ####

# plot data
plot_dat %>%
  select(litter, litter_weight.g, soil_moiture.prop, canopy_cover.prop) %>%
  ggpairs()
# weak relationships among the variables

# site data
cor.test(site_dat$canopy_cover.prop, site_dat$stand_area.m2ha)
# uncorrelated


#### save data ####
write_csv(plot_dat, "plot_covariates_2018_litter_exp.csv")
write_csv(site_dat, "site_covariates_2018_litter_exp.csv")
