##### info ####

# file: mv_seeds_per_biomass_analysis_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 7/21/20
# goal: analyze Mv seeds production per unit biomass


# see if there's an adjustment for fungicide
# analyze the raw and adjusted data (adjusted for some missing values)

#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(cowplot)
library(glmmTMB)

# import data
seeds18 <- read_csv("./intermediate-data/mv_processed_seeds_2018_density_exp.csv")
seeds19 <- read_csv("intermediate-data/mv_plant_level_seeds_2019_density_exp.csv")
bio18 <- read_csv("./data/mv_biomass_oct_2018_density_exp.csv")
bio19 <- read_csv("./data/mv_biomass_seeds_2019_density_exp.csv")
plots <- read_csv("data/plot_treatments_for_analyses_2018_2019_density_exp.csv")
treat <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")
sev18 <- read_csv("intermediate-data/focal_leaf_scans_2018_density_exp.csv")
sev19 <- read_csv("intermediate-data/all_leaf_scans_2019_density_exp.csv")


#### edit data ####

# both years of data
dat <- seeds19 %>%
  select(site, plot, treatment, sp, plant, seeds) %>%
  left_join(bio19 %>%
              select(site, plot, treatment, sp, plant, biomass_weight.g)) %>%
  mutate(yearf = "Year 2") %>%
  full_join(seeds18 %>%
              select(site, plot, treatment, total_seeds) %>%
              rename(seeds = total_seeds) %>%
              left_join(bio18 %>%
                          mutate(site = case_when(processing_notes == "duplicate - one is D2, probably this one because the other was sorted with other D1 bags; Quynh guessed on which seeds go with which biomass (see email)" ~ "D2", 
                                                  TRUE ~ site),
                                 treatment = case_when(processing_notes == "duplicate - one is 7F, ID = A (not dated), probably this one" ~ "fungicide",
                                                       TRUE ~ treatment)) %>%
                          select(site, plot, treatment, bio.g) %>%
                          rename(biomass_weight.g = bio.g)) %>%
              mutate(sp = "Mv",
                     yearf = "Year 1")) %>%
  mutate(seeds_per_biomass = seeds / biomass_weight.g,
         fungicide = case_when(treatment == "water" ~ 0,
                               TRUE ~ 1))

# combine with plot data
dat_plot <- dat %>%
  left_join(plots)


#### visualize ####

# seeds and biomass relationships
ggplot(dat, aes(biomass_weight.g, seeds, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ yearf, scales = "free")
# only related on the individual-scale

# seeds per biomass across density
ggplot(dat_plot, aes(background_density, seeds_per_biomass)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, aes(group = treatment)) +
  stat_summary(geom = "point", fun = "mean", shape = 21, size = 2, aes(fill = treatment)) +
  facet_grid(yearf ~ background, scales = "free_x")
# looks like a beverton-holt function for year 2
# more evidence for stress in no background plots (acting differently, allocating more to seeds)

# biomass across density
ggplot(dat_plot, aes(background_density, biomass_weight.g)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, aes(group = treatment)) +
  stat_summary(geom = "point", fun = "mean", shape = 21, size = 2, aes(fill = treatment)) +
  facet_grid(yearf ~ background, scales = "free_x")

# seeds across density
ggplot(dat_plot, aes(background_density, seeds)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, aes(group = treatment)) +
  stat_summary(geom = "point", fun = "mean", shape = 21, size = 2, aes(fill = treatment)) +
  facet_grid(yearf ~ background, scales = "free_x")


#### seeds per biomass regression ####

mv_seed_bio_mod <- glmmTMB(seeds ~ yearf * fungicide * biomass_weight.g + (1|site), data = dat, family = gaussian)
summary(mv_seed_bio_mod)
stepAIC(mv_seed_bio_mod)
# keep full model


#### output ####
save(mv_seed_bio_mod, file = "output/mv_seeds_per_biomass_model_2018_2019_density_exp.rda")
