##### info ####

# file: cleistogamous_seeds_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 3/17/21
# goal: analyses of ratio of cleistogamous to total seeds


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(tidybayes) # for mean_hdi

# import seed data
mvSeedD1Dat <- read_csv("intermediate-data/mv_processed_seeds_2018_density_exp.csv")
# mv_seeds_data_processing_2018_density_exp.R and mv_biomass_data_processing_2018_density_exp.R
mvSeedD2Dat <- read_csv("intermediate-data/mv_plant_level_seeds_2019_density_exp.csv") 
# mv_seeds_data_processing_2019_density_exp.R
# adj values substitute averages of other two plants in the plot when a plant is missing data

# model functions
source("code/brms_model_fitting_functions.R")

# fungicide effect function
post_pred_fun <- function(mod){
  
  posterior_samples(mod) %>%
    rename(b_fung_seeds = "b_fungicide:log_seeds") %>%
    transmute(waterSlope = b_log_seeds,
              fungSlope = b_log_seeds + b_fung_seeds,
              slope = 100 * (fungSlope - waterSlope) / waterSlope) %>%
    pivot_longer(cols = everything(),
                 names_to = "parameter",
                 values_to = "effect") %>%
    group_by(parameter) %>% mean_hdi(effect)
  
}


#### edit data ####

# 2018 data
unique(mvBioD1Dat$processing_notes)

csD1Dat <- mvSeedD1Dat %>%
  filter(!is.na(seeds_bio) & !is.na(seeds_soil)) %>%
  mutate(seeds = seeds_bio + seeds_soil,
         log_cl_seeds = log(seeds_bio),
         log_seeds = log(seeds),
         fungicide = ifelse(treatment == "fungicide", 1, 0))
# seeds are all at the scale of one quadrat (0.25 x 0.49 m)

# 2019 data
csD2Dat <- mvSeedD2Dat %>%
  filter(!is.na(seeds)) %>%
  mutate(log_cl_seeds = log(stem_seeds),
         log_seeds = log(seeds),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotf = paste(site, plot, substr(treatment, 1, 1), sep = ""))
# seeds are flower_seeds + stem_seeds


#### 2018 Mv model ####

# initial visualization
ggplot(csD1Dat, aes(log_seeds, log_cl_seeds, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~x)

# model
csD1Mod <- brm(data = csD1Dat, family = gaussian,
               log_cl_seeds ~ fungicide * log_seeds + (1|site),
               prior <- c(prior(normal(0, 1), class = "Intercept"),
                          prior(normal(0, 1), class = "b")), # use default for sigma
               iter = 6000, warmup = 1000, chains = 3,
               control = list(adapt_delta = 0.99, max_treedepth = 15))
mod_check_fun(csD1Mod)

# simulated data
csD1Sim <- mod_fit_fun(dat = csD1Dat, mod = csD1Mod, treatCol = fungicide,
                       xCol = log_seeds, minX = min(csD1Dat$log_seeds), maxX = max(csD1Dat$log_seeds), 
                       yCol = log_cl_seeds, f2t = T)
csD1Sim[[2]]

# posterior means
post_pred_fun(csD1Mod)

# save
save(csD1Mod, file = "output/cleistogamous_seeds_model_2018_density_exp.rda")


#### 2019 Mv model ####

# initial visualization
ggplot(csD2Dat, aes(log_seeds, log_cl_seeds, color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~x)

# model
csD2Mod <- brm(data = csD2Dat, family = gaussian,
               log_cl_seeds ~ fungicide * log_seeds + (1|plotf),
               prior <- c(prior(normal(0, 1), class = "Intercept"),
                          prior(normal(0, 1), class = "b")), # use default for sigma
               iter = 6000, warmup = 1000, chains = 3,
               control = list(adapt_delta = 0.99))
mod_check_fun(csD2Mod)

# simulated data
csD2Sim <- mod_fit_fun(dat = csD2Dat, mod = csD2Mod, treatCol = fungicide,
                       xCol = log_seeds, minX = min(csD2Dat$log_seeds), maxX = max(csD2Dat$log_seeds), 
                       yCol = log_cl_seeds, f2t = T)
csD2Sim[[2]]

# posterior means
post_pred_fun(csD2Mod)

# save
save(csD2Mod, file = "output/cleistogamous_seeds_model_2019_density_exp.rda")