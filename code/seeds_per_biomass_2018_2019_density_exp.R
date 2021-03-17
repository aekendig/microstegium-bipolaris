##### info ####

# file: seeds_per_biomass_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 3/17/21
# goal: analyses of seeds per g biomass


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(tidybayes) # for mean_hdi

# import biomass data
mvBioD1Dat <- read_csv("./intermediate-data/mv_processed_biomass_oct_2018_density_exp.csv")
# mv_biomass_data_processing_2018_density_exp.R
mvBioD2Dat <- read_csv("data/mv_biomass_seeds_2019_density_exp.csv")
evBioD2Dat <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")

# import seed data
mvSeedD1Dat <- read_csv("intermediate-data/mv_processed_seeds_2018_density_exp.csv")
# mv_seeds_data_processing_2018_density_exp.R and mv_biomass_data_processing_2018_density_exp.R
mvSeedD2Dat <- read_csv("intermediate-data/mv_plant_level_seeds_2019_density_exp.csv") 
# mv_seeds_data_processing_2019_density_exp.R
# adj values substitute averages of other two plants in the plot when a plant is missing data
evSeedD2Dat <- read_csv("intermediate-data/ev_processed_seeds_both_year_conversion_2019_density_exp.csv") 
# ev_seeds_data_processing_2019.R and ev_seeds_data_processing_2018.R

# model functions
source("code/brms_model_fitting_functions.R")

# fungicide effect function
post_pred_fun <- function(mod){
  
  posterior_samples(mod) %>%
    rename(b_fung_bio = "b_fungicide:log_bio") %>%
    transmute(waterSlope = b_log_bio,
              fungSlope = b_log_bio + b_fung_bio,
              slope = 100 * (fungSlope - waterSlope) / waterSlope) %>%
    pivot_longer(cols = everything(),
                 names_to = "parameter",
                 values_to = "effect") %>%
    group_by(parameter) %>% mean_hdi(effect)
  
}


#### edit data ####

# Ev dat
evD2Dat <- evSeedD2Dat %>%
  group_by(site, plot, treatment, sp, ID) %>%
  summarise(seeds = sum(seeds)) %>%
  ungroup() %>%
  full_join(evBioD2Dat %>%
              select(site, plot, treatment, sp, ID, weight)) %>%
  mutate(seeds = replace_na(seeds, 0),
         age = ifelse(ID == "A", "adult", "seedling")) %>%
  rename(biomass_weight.g = weight)

# 2018 data
unique(mvBioD1Dat$processing_notes)

seedsBioD1Dat <- mvSeedD1Dat %>%
  mutate(seeds = seeds_bio + seeds_soil) %>%
  full_join(mvBioD1Dat %>%
              select(site, plot, treatment, bio.g)) %>%
  filter(!is.na(bio.g) & !is.na(seeds)) %>%
  mutate(log_seeds = log(seeds),
         log_bio = log(bio.g),
         fungicide = ifelse(treatment == "fungicide", 1, 0))
# seeds and biomass are all at the scale of one quadrat (0.25 x 0.49 m)

# 2019 data
seedsBioD2Dat <- mvSeedD2Dat %>%
  select(site, plot, treatment, sp, plant, seeds) %>%
  full_join(mvBioD2Dat %>%
              select(site, plot, treatment, sp, plant, biomass_weight.g)) %>%
  mutate(ID = as.character(plant),
         age = "seedling")  %>%
  select(-plant) %>%
  full_join(evD2Dat) %>%
  filter(!is.na(biomass_weight.g) & !is.na(seeds)) %>%
  mutate(log_seeds = log(seeds + 1),
         log_bio = log(biomass_weight.g),
         plant_group = paste(sp, age, sep = "_") %>%
           fct_rev(),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotf = paste(site, plot, substr(treatment, 1, 1), sep = ""))


#### 2018 Mv model ####

# initial visualization
ggplot(seedsBioD1Dat, aes(log_bio, log_seeds, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~x)

# model
seedsBioD1Mod <- brm(data = seedsBioD1Dat, family = gaussian,
                     log_seeds ~ fungicide * log_bio + (1|site),
                     prior <- c(prior(normal(6, 1), class = "Intercept"),
                                prior(normal(0, 1), class = "b")), # use default for sigma
                     iter = 6000, warmup = 1000, chains = 3,
                     control = list(adapt_delta = 0.99))
mod_check_fun(seedsBioD1Mod)

# simulated data
seedsBioD1Sim <- mod_fit_fun(dat = seedsBioD1Dat, mod = seedsBioD1Mod, treatCol = fungicide,
                         xCol = log_bio, minX = min(seedsBioD1Dat$log_bio), maxX = max(seedsBioD1Dat$log_bio), 
                         yCol = log_seeds, f2t = T)
seedsBioD1Sim[[2]]

# posterior means
post_pred_fun(seedsBioD1Mod)

# save
save(seedsBioD1Mod, file = "output/seeds_per_biomass_model_2018_density_exp.rda")


#### 2019 Mv model ####

# initial visualization
ggplot(seedsBioD2Dat, aes(log_bio, log_seeds, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~x) +
  facet_wrap(~ plant_group, scales = "free")

# remove 0 seeds
seedsBioD2Dat2 <- seedsBioD2Dat %>%
  filter(seeds > 0) %>%
  mutate(log_seeds = log(seeds))

ggplot(seedsBioD2Dat2, aes(log_bio, log_seeds, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~x) +
  facet_wrap(~ plant_group, scales = "free")

# model
seedsBioD2Mod <- brm(data = seedsBioD2Dat, family = gaussian,
                     log_seeds ~ plant_group * fungicide * log_bio + (1|site/plotf),
                     prior <- c(prior(normal(4, 1), class = "Intercept"),
                                prior(normal(0, 1), class = "b")), # use default for sigma
                     iter = 6000, warmup = 1000, chains = 3,
                     control = list(adapt_delta = 0.999))
mod_check_fun(seedsBioD2Mod)

# simulated data
seedsBioD2Sim <- mod_fit_fun(dat = seedsBioD2Dat, mod = seedsBioD2Mod, treatCol = fungicide,
                             xCol = log_bio, minX = min(seedsBioD2Dat$log_bio), maxX = max(seedsBioD2Dat$log_bio), 
                             yCol = log_seeds, f2t = T)
seedsBioD2Sim[[2]]
# posterior predictive check suggests individual models for each plant group are needed

# divide data
mvSeedsBioD2Dat <- seedsBioD2Dat2 %>%
  filter(plant_group == "Mv_seedling")

evSSeedsBioD2Dat <- seedsBioD2Dat2 %>%
  filter(plant_group == "Ev_seedling")

evASeedsBioD2Dat <- seedsBioD2Dat2 %>%
  filter(plant_group == "Ev_adult")

# fit models
mvSeedsBioD2Mod <- update(seedsBioD2Mod, formula. = log_seeds ~ fungicide * log_bio + (1|plotf),
                          newdata = mvSeedsBioD2Dat)
# lots of warnings and very slow when I included site in the random effects
mod_check_fun(mvSeedsBioD2Mod)

evSSeedsBioD2Mod <- update(mvSeedsBioD2Mod, newdata = evSSeedsBioD2Dat)
mod_check_fun(evSSeedsBioD2Mod)

evASeedsBioD2Mod <- update(mvSeedsBioD2Mod, formula. = log_seeds ~ fungicide * log_bio + (1|site), 
                           newdata = evASeedsBioD2Dat)
mod_check_fun(evASeedsBioD2Mod)

# simulated data
mvSeedsBioD2Sim <- mod_fit_fun(dat = mvSeedsBioD2Dat, mod = mvSeedsBioD2Mod, treatCol = fungicide,
                             xCol = log_bio, minX = min(mvSeedsBioD2Dat$log_bio), maxX = max(mvSeedsBioD2Dat$log_bio), 
                             yCol = log_seeds, f2t = T)
mvSeedsBioD2Sim[[2]]

evSSeedsBioD2Sim <- mod_fit_fun(dat = evSSeedsBioD2Dat, mod = evSSeedsBioD2Mod, treatCol = fungicide,
                               xCol = log_bio, minX = min(evSSeedsBioD2Dat$log_bio), maxX = max(evSSeedsBioD2Dat$log_bio), 
                               yCol = log_seeds, f2t = T)
evSSeedsBioD2Sim[[2]]

evASeedsBioD2Sim <- mod_fit_fun(dat = evASeedsBioD2Dat, mod = evASeedsBioD2Mod, treatCol = fungicide,
                                xCol = log_bio, minX = min(evASeedsBioD2Dat$log_bio), maxX = max(evASeedsBioD2Dat$log_bio), 
                                yCol = log_seeds, f2t = T)
evASeedsBioD2Sim[[2]]

# posterior means
post_pred_fun(mvSeedsBioD2Mod)
post_pred_fun(evSSeedsBioD2Mod)
post_pred_fun(evASeedsBioD2Mod)

# save
save(mvSeedsBioD2Mod, file = "output/mv_seeds_per_biomass_model_2019_density_exp.rda")
save(evSSeedsBioD2Mod, file = "output/evS_seeds_per_biomass_model_2019_density_exp.rda")
save(evASeedsBioD2Mod, file = "output/evA_seeds_per_biomass_model_2019_density_exp.rda")