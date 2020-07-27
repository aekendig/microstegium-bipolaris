##### info ####

# file: ev_seeds_per_biomass_analysis_2019_density_exp
# author: Amy Kendig
# date last edited: 7/21/20
# goal: analyze Ev seeds production per unit biomass


#### set-up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(cowplot)
library(glmmTMB)

# import data
bio <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")
seed <- read_csv("intermediate-data/ev_processed_seeds_both_year_conversion_2019_density_exp.csv")
plots <- read_csv("data/plot_treatments_for_analyses_2018_2019_density_exp.csv")
treat <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")
sev <- read_csv("intermediate-data/all_leaf_scans_2019_density_exp.csv")


#### edit data ####

dat <- bio %>%
  select(-c(seeds_collected:processing_notes)) %>%
  rename(biomass_weight.g = weight) %>%
  full_join(seed %>%
              group_by(site, plot, treatment, sp, ID) %>%
              summarise(seeds = sum(seeds)) %>%
              ungroup) %>%
              mutate(age = case_when(ID == "A" ~ "adult",
                                     TRUE ~ "seedling"),
                     seeds = replace_na(seeds, 0),
                     seeds_per_biomass = seeds / biomass_weight.g,
                     fungicide = case_when(treatment == "water" ~ 0,
                                           TRUE ~ 1),
                     log_seeds = log(seeds + 1),
                     log_biomass = log(biomass_weight.g),
                     plot_trt = paste(plot, substr(treatment, 1, 1) %>% toupper(), sep = ""))

# combine with plot data
dat_plot <- dat %>%
  filter(!is.na(seeds_per_biomass)) %>%
  left_join(plots)



#### visualize ####

# seeds and biomass relationships
ggplot(dat, aes(biomass_weight.g, seeds, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ age, scales = "free")

# log-transformed seeds and biomass relationships
ggplot(dat, aes(log_biomass, log_seeds, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ age, scales = "free")

# seeds per biomass across density
ggplot(dat_plot, aes(background_density, seeds_per_biomass)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, aes(group = treatment)) +
  stat_summary(geom = "point", fun = "mean", shape = 21, size = 2, aes(fill = treatment)) +
  facet_grid(age ~ background, scales = "free")
# does not look like beverton-holt

# biomass across density
ggplot(dat_plot, aes(background_density, biomass_weight.g)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, aes(group = treatment)) +
  stat_summary(geom = "point", fun = "mean", shape = 21, size = 2, aes(fill = treatment)) +
  facet_grid(age ~ background, scales = "free")

# seeds across density
ggplot(dat_plot, aes(background_density, seeds)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0, aes(group = treatment)) +
  stat_summary(geom = "point", fun = "mean", shape = 21, size = 2, aes(fill = treatment)) +
  facet_grid(age ~ background, scales = "free")


#### seeds per biomass regression ####

ev_seed_bio_mod <- glmmTMB(log_seeds ~ age * fungicide * log_biomass + (1|site/plot), data = dat, family = gaussian)
summary(ev_seed_bio_mod)
stepAIC(ev_seed_bio_mod)
# keep full model

# remove treatment effect because mean values cause fluctuations in Elymus seedling population that make simulations difficult to interpret
ev_seed_bio_no_treat_mod <- glmmTMB(log_seeds ~ age * log_biomass + (1|site/plot_trt), data = dat, family = gaussian)
summary(ev_seed_bio_no_treat_mod)
stepAIC(ev_seed_bio_no_treat_mod)

#### output ####
save(ev_seed_bio_mod, file = "output/ev_seeds_per_biomass_model_2019_density_exp.rda")
save(ev_seed_bio_no_treat_mod, file = "output/ev_seeds_per_biomass_no_treatment_model_2019_density_exp.rda")
