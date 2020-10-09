##### info ####

# file: elymus_adult_seeds_per_biomass_model_2019_density_exp
# author: Amy Kendig
# date last edited: 10/8/20
# goal: analyze Elymus adult seeds per unit biomass


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)

# import data
bioD2Dat <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")
seedD2Dat <- read_csv("intermediate-data/ev_processed_seeds_both_year_conversion_2019_density_exp.csv")


#### edit data ####

# add columns
# remove missing data
# select Elymus adults
evASeedD2Dat <- seedD2Dat %>%
  group_by(site, plot, treatment, sp, ID) %>%
  summarise(seeds = sum(seeds)) %>%
  ungroup() %>%
  full_join(bioD2Dat %>%
              rename(bio.g = weight)) %>%
  mutate(seeds = replace_na(seeds, 0),
         seeds_per_bio = seeds / bio.g,
         seeds_prod = ifelse(seeds > 0, 1, 0),
         log_seeds_per_bio = log(seeds_per_bio),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         log_bio.g = log(bio.g),
         treatment = recode(treatment, water = "control"),
         plotr = ifelse(treatment == "fungicide", plot + 10, plot)) %>%
  filter(!is.na(seeds_per_bio) & ID == "A")


#### initial visualizations ####

# figure
ggplot(evASeedD2Dat, aes(log_bio.g, log_seeds_per_bio, color = treatment)) +
  geom_point(size = 2, alpha = 0.5) +
  theme_bw()

# histogram
ggplot(evASeedD2Dat, aes(seeds_per_bio)) +
  geom_histogram() + 
  theme_bw()
# 7 zeros

# yes/no seeds
ggplot(evASeedD2Dat, aes(log_bio.g, seeds_prod, color = treatment)) +
  geom_point(size = 2, alpha = 0.5) +
  theme_bw()
# the range looks very similar to Elymus seedlings -- combine the data
