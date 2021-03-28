##### info ####

# file: focal_growth_density_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 3/17/21
# goal: analyses of plant growth as a function of density


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(tidybayes) # for mean_hdi

# import plot information
plotsD <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")

# import growth data
growthD1Dat <- read_csv("intermediate-data/focal_processed_growth_2018_density_exp.csv")
# focal_growth_data_processing_2018_density_exp
mvBioD2Dat <- read_csv("data/mv_biomass_seeds_2019_density_exp.csv")
evBioD2Dat <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")

# model functions
source("code/brms_model_fitting_functions.R")

# fungicide effect function
post_pred_fun <- function(mod){
  
  posterior_samples(mod) %>%
    transmute(b0_FE = 100 * (b_b0_treatmentfungicide - b_b0_treatmentwater) / b_b0_treatmentwater,
              alpha_FE = 100 * (b_alpha_treatmentfungicide - b_alpha_treatmentwater) / b_alpha_treatmentwater) %>%
    pivot_longer(cols = everything(),
                 names_to = "parameter",
                 values_to = "effect") %>%
    group_by(parameter) %>% mean_hdi(effect)
  
}

post_pred_fun2 <- function(mod){
  
  posterior_samples(mod) %>%
    transmute(b0_FE = 100 * (b_b0_treatmentfungicide - b_b0_treatmentwater) / b_b0_treatmentwater,
              alpha_FE = 100 * (b_alpha_treatmentfungicide - b_alpha_treatmentwater) / b_alpha_treatmentwater) %>%
    pivot_longer(cols = everything(),
                 names_to = "parameter",
                 values_to = "effect") %>%
    group_by(parameter) %>% median_hdi(effect)
  
}



#### edit data ####

# plant group densities
plotDens <- plotsD %>%
  mutate(Mv_seedling_density = case_when(background == "Mv seedling" ~ background_density + 3, 
                                         TRUE ~ 3),
         Ev_seedling_density = case_when(background == "Ev seedling" ~ background_density + 3, 
                                         TRUE ~ 3),
         Ev_adult_density = case_when(background == "Ev adult" ~ background_density + 1, 
                                      TRUE ~ 1)) %>%
  select(plot, treatment, Mv_seedling_density, Ev_seedling_density, Ev_adult_density)

# missing data
filter(growthD1Dat, is.na(tillers_jul))
filter(mvBioD2Dat, is.na(biomass_weight.g)) # 3
filter(evBioD2Dat, is.na(weight)) # 1 seedling

# combine and separate data
mvD1Dat <- plotDens %>%
  filter(plot %in% 1:4) %>%
  left_join(growthD1Dat) %>%
  mutate(age = ifelse(ID == "A", "adult", "seedling"),
         plant_group = paste(sp, age) %>%
           fct_rev()) %>%
  rename(density = Mv_seedling_density)

mvD2Dat <- plotDens %>%
  filter(plot %in% 1:4) %>%
  left_join(mvBioD2Dat %>%
              rename(ID = plant) %>%
              mutate(ID = as.character(ID)) %>%
              full_join(evBioD2Dat %>%
                          rename(biomass_weight.g = weight))) %>%
  mutate(age = ifelse(ID == "A", "adult", "seedling"),
         plant_group = paste(sp, age) %>%
           fct_rev()) %>%
  rename(density = Mv_seedling_density) %>%
  filter(!is.na(biomass_weight.g))

evSD1Dat <- plotDens %>%
  filter(plot %in% c(1, 5:7)) %>%
  left_join(growthD1Dat) %>%
  mutate(age = ifelse(ID == "A", "adult", "seedling"),
         plant_group = paste(sp, age) %>%
           fct_rev()) %>%
  rename(density = Ev_seedling_density)

evSD2Dat <- plotDens %>%
  filter(plot %in% c(1, 5:7)) %>%
  left_join(mvBioD2Dat %>%
              rename(ID = plant) %>%
              mutate(ID = as.character(ID)) %>%
              full_join(evBioD2Dat %>%
                          rename(biomass_weight.g = weight))) %>%
  mutate(age = ifelse(ID == "A", "adult", "seedling"),
         plant_group = paste(sp, age) %>%
           fct_rev()) %>%
  rename(density = Ev_seedling_density) %>%
  filter(!is.na(biomass_weight.g))

evAD1Dat <- plotDens %>%
  filter(plot %in% c(1, 8:10)) %>%
  left_join(growthD1Dat) %>%
  mutate(age = ifelse(ID == "A", "adult", "seedling"),
         plant_group = paste(sp, age) %>%
           fct_rev()) %>%
  rename(density = Ev_adult_density)

evAD2Dat <- plotDens %>%
  filter(plot %in% c(1, 8:10)) %>%
  left_join(mvBioD2Dat %>%
              rename(ID = plant) %>%
              mutate(ID = as.character(ID)) %>%
              full_join(evBioD2Dat %>%
                          rename(biomass_weight.g = weight))) %>%
  mutate(age = ifelse(ID == "A", "adult", "seedling"),
         plant_group = paste(sp, age) %>%
           fct_rev()) %>%
  rename(density = Ev_adult_density) %>%
  filter(!is.na(biomass_weight.g))

# remove plot 1 (large environmental effect)
mvD1Dat2 <- mvD1Dat %>%
  filter(plot != 1)
evSD1Dat2 <- evSD1Dat %>%
  filter(plot != 1)
evAD1Dat2 <- evAD1Dat %>%
  filter(plot != 1)

mvD2Dat2 <- mvD2Dat %>%
  filter(plot != 1)
evSD2Dat2 <- evSD2Dat %>%
  filter(plot != 1)
evAD2Dat2 <- evAD2Dat %>%
  filter(plot != 1)

# split by species (couldn't specify appropriate priors with them all together)
mvMvD1Dat <- mvD1Dat2 %>%
  filter(plant_group == "Mv seedling")
evSMvD1Dat <- mvD1Dat2 %>%
  filter(plant_group == "Ev seedling")
evAMvD1Dat <- mvD1Dat2 %>%
  filter(plant_group == "Ev adult")

mvEvSD1Dat <- evSD1Dat2 %>%
  filter(plant_group == "Mv seedling")
evSEvSD1Dat <- evSD1Dat2 %>%
  filter(plant_group == "Ev seedling")
evAEvSD1Dat <- evSD1Dat2 %>%
  filter(plant_group == "Ev adult")

mvEvAD1Dat <- evAD1Dat2 %>%
  filter(plant_group == "Mv seedling")
evSEvAD1Dat <- evAD1Dat2 %>%
  filter(plant_group == "Ev seedling")
evAEvAD1Dat <- evAD1Dat2 %>%
  filter(plant_group == "Ev adult")

mvMvD2Dat <- mvD2Dat2 %>%
  filter(plant_group == "Mv seedling")
evSMvD2Dat <- mvD2Dat2 %>%
  filter(plant_group == "Ev seedling")
evAMvD2Dat <- mvD2Dat2 %>%
  filter(plant_group == "Ev adult")

mvEvSD2Dat <- evSD2Dat2 %>%
  filter(plant_group == "Mv seedling")
evSEvSD2Dat <- evSD2Dat2 %>%
  filter(plant_group == "Ev seedling")
evAEvSD2Dat <- evSD2Dat2 %>%
  filter(plant_group == "Ev adult")

mvEvAD2Dat <- evAD2Dat2 %>%
  filter(plant_group == "Mv seedling")
evSEvAD2Dat <- evAD2Dat2 %>%
  filter(plant_group == "Ev seedling")
evAEvAD2Dat <- evAD2Dat2 %>%
  filter(plant_group == "Ev adult")


#### 2019 Mv density models ####

# initial visualization
ggplot(mvD2Dat, aes(density, biomass_weight.g, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ plant_group, scales = "free")

# intercept prior
mvD2Dat2 %>%
  filter(density == 7) %>%
  group_by(plant_group, treatment) %>%
  summarise(bio = mean(biomass_weight.g))

prior_fun(maxX = 40, shape = 20, dist = "gamma")
prior_fun(maxX = 10, shape = 2.5, dist = "gamma")
prior_fun(maxX = 20, shape = 8, dist = "gamma")

# alpha prior
prior_fun(maxX = 10, shape = 0.5, dist = "exponential")

# fit models
mvMvD2Mod <- brm(data = mvMvD2Dat, family = gaussian,
               bf(biomass_weight.g ~ b0/(1 + alpha * density), 
                  b0 ~ 0 + treatment + (1|site), 
                  alpha ~ 0 + treatment, 
                  nl = T),
               prior <- c(prior(gamma(20, 1), nlpar = "b0", lb = 0),
                          prior(exponential(0.5), nlpar = "alpha", lb = 0)), # use default for sigma
               iter = 6000, warmup = 1000, chains = 3) 
mod_check_fun(mvMvD2Mod)

evSMvD2Mod <- update(mvMvD2Mod, newdata = evSMvD2Dat,
                     prior = set_prior("gamma(2.5, 1)", nlpar = "b0", lb = 0),
                     control = list(adapt_delta = 0.99))
mod_check_fun(evSMvD2Mod)

evAMvD2Mod <- update(mvMvD2Mod, newdata = evAMvD2Dat,
                     control = list(adapt_delta = 0.99))
# initially ran with prior shape = 8, but the predictions were much lower than the data
mod_check_fun(evAMvD2Mod)

# simulated data
mvMvD2Sim <- mod_fit_fun(dat = mvMvD2Dat, mod = mvMvD2Mod, treatCol = treatment,
                         xCol = density, minX = 0, maxX = 67, yCol = biomass_weight.g)
mvMvD2Sim[[2]]

evSMvD2Sim <- mod_fit_fun(dat = evSMvD2Dat, mod = evSMvD2Mod, treatCol = treatment,
                         xCol = density, minX = 0, maxX = 67, yCol = biomass_weight.g)
evSMvD2Sim[[2]]

evAMvD2Sim <- mod_fit_fun(dat = evAMvD2Dat, mod = evAMvD2Mod, treatCol = treatment,
                          xCol = density, minX = 0, maxX = 67, yCol = biomass_weight.g)
evAMvD2Sim[[2]]

# posterior means
post_pred_fun(mvMvD2Mod)
post_pred_fun(evSMvD2Mod)
post_pred_fun(evAMvD2Mod)

# save models
save(mvMvD2Mod, file = "output/mv_biomass_mv_density_model_2019_density_exp.rda")
save(evSMvD2Mod, file = "output/evS_biomass_mv_density_model_2019_density_exp.rda")
save(evAMvD2Mod, file = "output/evA_biomass_mv_density_model_2019_density_exp.rda")


#### 2019 Ev seedling density models ####

# initial visualization
ggplot(evSD2Dat, aes(density, biomass_weight.g, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ plant_group, scales = "free")

# intercept prior
evSD2Dat2 %>%
  filter(density == 7) %>%
  group_by(plant_group, treatment) %>%
  summarise(bio = mean(biomass_weight.g))

prior_fun(maxX = 40, shape = 22, dist = "gamma")
prior_fun(maxX = 10, shape = 2.5, dist = "gamma")
prior_fun(maxX = 20, shape = 10, dist = "gamma")

# alpha prior
prior_fun(maxX = 10, shape = 0.5, dist = "exponential")

# fit models
mvEvSD2Mod <- update(mvMvD2Mod, newdata = mvEvSD2Dat,
                     prior = set_prior("gamma(22, 1)", nlpar = "b0", lb = 0))
mod_check_fun(mvEvSD2Mod)

evSEvSD2Mod <- update(mvEvSD2Mod, newdata = evSEvSD2Dat,
                     prior = set_prior("gamma(2.5, 1)", nlpar = "b0", lb = 0),
                     control = list(adapt_delta = 0.99))
mod_check_fun(evSEvSD2Mod)

evAEvSD2Mod <- update(mvEvSD2Mod, newdata = evAEvSD2Dat,
                      prior = set_prior("gamma(10, 1)", nlpar = "b0", lb = 0),
                     control = list(adapt_delta = 0.99))
mod_check_fun(evAEvSD2Mod)

# simulated data
mvEvSD2Sim <- mod_fit_fun(dat = mvEvSD2Dat, mod = mvEvSD2Mod, treatCol = treatment,
                         xCol = density, minX = 0, maxX = 19, yCol = biomass_weight.g)
mvEvSD2Sim[[2]]

evSEvSD2Sim <- mod_fit_fun(dat = evSEvSD2Dat, mod = evSEvSD2Mod, treatCol = treatment,
                          xCol = density, minX = 0, maxX = 19, yCol = biomass_weight.g)
evSEvSD2Sim[[2]]

evAEvSD2Sim <- mod_fit_fun(dat = evAEvSD2Dat, mod = evAEvSD2Mod, treatCol = treatment,
                          xCol = density, minX = 0, maxX = 19, yCol = biomass_weight.g)
evAEvSD2Sim[[2]]

# posterior means
post_pred_fun(mvEvSD2Mod)
post_pred_fun2(mvEvSD2Mod)
post_pred_fun(evSEvSD2Mod)
post_pred_fun(evAEvSD2Mod)

# save models
save(mvEvSD2Mod, file = "output/mv_biomass_evS_density_model_2019_density_exp.rda")
save(evSEvSD2Mod, file = "output/evS_biomass_evS_density_model_2019_density_exp.rda")
save(evAEvSD2Mod, file = "output/evA_biomass_evS_density_model_2019_density_exp.rda")


#### 2019 Ev adult density models ####

# initial visualization
ggplot(evAD2Dat, aes(density, biomass_weight.g, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ plant_group, scales = "free")

# intercept prior
evAD2Dat2 %>%
  filter(density == 3) %>%
  group_by(plant_group, treatment) %>%
  summarise(bio = mean(biomass_weight.g))

prior_fun(maxX = 50, shape = 35, dist = "gamma")
prior_fun(maxX = 10, shape = 3.5, dist = "gamma")
prior_fun(maxX = 20, shape = 12, dist = "gamma")

# alpha prior
prior_fun(maxX = 10, shape = 0.5, dist = "exponential")

# fit models
mvEvAD2Mod <- update(mvMvD2Mod, newdata = mvEvAD2Dat,
                     prior = set_prior("gamma(35, 1)", nlpar = "b0", lb = 0))
mod_check_fun(mvEvAD2Mod)

evSEvAD2Mod <- update(mvEvAD2Mod, newdata = evSEvAD2Dat,
                      prior = set_prior("gamma(3.5, 1)", nlpar = "b0", lb = 0),
                      control = list(adapt_delta = 0.99))
mod_check_fun(evSEvAD2Mod)

evAEvAD2Mod <- update(mvEvAD2Mod, newdata = evAEvAD2Dat,
                      prior = set_prior("gamma(12, 1)", nlpar = "b0", lb = 0),
                      control = list(adapt_delta = 0.99))
mod_check_fun(evAEvAD2Mod)

# simulated data
mvEvAD2Sim <- mod_fit_fun(dat = mvEvAD2Dat, mod = mvEvAD2Mod, treatCol = treatment,
                         xCol = density, minX = 0, maxX = 9, yCol = biomass_weight.g)
mvEvAD2Sim[[2]]

evSEvAD2Sim <- mod_fit_fun(dat = evSEvAD2Dat, mod = evSEvAD2Mod, treatCol = treatment,
                          xCol = density, minX = 0, maxX = 9, yCol = biomass_weight.g)
evSEvAD2Sim[[2]]

evAEvAD2Sim <- mod_fit_fun(dat = evAEvAD2Dat, mod = evAEvAD2Mod, treatCol = treatment,
                          xCol = density, minX = 0, maxX = 9, yCol = biomass_weight.g)
evAEvAD2Sim[[2]]

# posterior means
post_pred_fun(mvEvAD2Mod)
post_pred_fun(evSEvAD2Mod)
post_pred_fun(evAEvAD2Mod)

# save models
save(mvEvAD2Mod, file = "output/mv_biomass_evA_density_model_2019_density_exp.rda")
save(evSEvAD2Mod, file = "output/evS_biomass_evA_density_model_2019_density_exp.rda")
save(evAEvAD2Mod, file = "output/evA_biomass_evA_density_model_2019_density_exp.rda")


#### 2018 Mv density models ####

# initial visualization
ggplot(mvD1Dat, aes(density, tillers_jul, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ plant_group, scales = "free")

# intercept prior
mvD1Dat2 %>%
  filter(density == 7) %>%
  group_by(plant_group, treatment) %>%
  summarise(bio = mean(tillers_jul))

prior_fun(maxX = 50, shape = 35, dist = "gamma")
prior_fun(maxX = 15, shape = 7, dist = "gamma")
prior_fun(maxX = 40, shape = 25, dist = "gamma")

# alpha prior
prior_fun(maxX = 10, shape = 0.5, dist = "exponential")

# fit models
mvMvD1Mod <- brm(data = mvMvD1Dat, family = gaussian,
                 bf(tillers_jul ~ b0/(1 + alpha * density), 
                    b0 ~ 0 + treatment + (1|site), 
                    alpha ~ 0 + treatment, 
                    nl = T),
                 prior <- c(prior(gamma(50, 1), nlpar = "b0", lb = 0),
                            prior(exponential(0.5), nlpar = "alpha", lb = 0)), # use default for sigma
                 iter = 6000, warmup = 1000, chains = 3) 
mod_check_fun(mvMvD1Mod)

evSMvD1Mod <- update(mvMvD1Mod, newdata = evSMvD1Dat,
                     prior = set_prior("gamma(7, 1)", nlpar = "b0", lb = 0),
                     control = list(adapt_delta = 0.99))
mod_check_fun(evSMvD1Mod)

evAMvD1Mod <- update(mvMvD1Mod, newdata = evAMvD1Dat,
                     prior = set_prior("gamma(25, 1)", nlpar = "b0", lb = 0),
                     control = list(adapt_delta = 0.99))
mod_check_fun(evAMvD1Mod)

# simulated data
mvMvD1Sim <- mod_fit_fun(dat = mvMvD1Dat, mod = mvMvD1Mod, treatCol = treatment,
                         xCol = density, minX = 0, maxX = 67, yCol = tillers_jul)
mvMvD1Sim[[2]]

evSMvD1Sim <- mod_fit_fun(dat = evSMvD1Dat, mod = evSMvD1Mod, treatCol = treatment,
                          xCol = density, minX = 0, maxX = 67, yCol = tillers_jul)
evSMvD1Sim[[2]]

evAMvD1Sim <- mod_fit_fun(dat = evAMvD1Dat, mod = evAMvD1Mod, treatCol = treatment,
                          xCol = density, minX = 0, maxX = 67, yCol = tillers_jul)
evAMvD1Sim[[2]]

# posterior means
post_pred_fun(mvMvD1Mod)
post_pred_fun(evSMvD1Mod)
post_pred_fun(evAMvD1Mod)

# save models
save(mvMvD1Mod, file = "output/mv_tillers_mv_density_model_2018_density_exp.rda")
save(evSMvD1Mod, file = "output/evS_tillers_mv_density_model_2018_density_exp.rda")
save(evAMvD1Mod, file = "output/evA_tillers_mv_density_model_2018_density_exp.rda")


#### 2018 Ev seedling density models ####

# initial visualization
ggplot(evSD1Dat, aes(density, tillers_jul, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ plant_group, scales = "free")

# intercept prior
evSD1Dat2 %>%
  filter(density == 7) %>%
  group_by(plant_group, treatment) %>%
  summarise(bio = mean(tillers_jul))

prior_fun(maxX = 50, shape = 37, dist = "gamma")
prior_fun(maxX = 15, shape = 3.5, dist = "gamma")
prior_fun(maxX = 40, shape = 17, dist = "gamma")

# alpha prior
prior_fun(maxX = 10, shape = 0.5, dist = "exponential")

# fit models
mvEvSD1Mod <- update(mvMvD1Mod, newdata = mvEvSD1Dat,
                     prior = set_prior("gamma(37, 1)", nlpar = "b0", lb = 0))
mod_check_fun(mvEvSD1Mod)

evSEvSD1Mod <- update(mvEvSD1Mod, newdata = evSEvSD1Dat,
                      prior = set_prior("gamma(3.5, 1)", nlpar = "b0", lb = 0),
                      control = list(adapt_delta = 0.99))
mod_check_fun(evSEvSD1Mod)

evAEvSD1Mod <- update(mvEvSD1Mod, newdata = evAEvSD1Dat,
                      prior = set_prior("gamma(17, 1)", nlpar = "b0", lb = 0),
                      control = list(adapt_delta = 0.99))
mod_check_fun(evAEvSD1Mod)

# simulated data
mvEvSD1Sim <- mod_fit_fun(dat = mvEvSD1Dat, mod = mvEvSD1Mod, treatCol = treatment,
                          xCol = density, minX = 0, maxX = 19, yCol = tillers_jul)
mvEvSD1Sim[[2]]

evSEvSD1Sim <- mod_fit_fun(dat = evSEvSD1Dat, mod = evSEvSD1Mod, treatCol = treatment,
                           xCol = density, minX = 0, maxX = 19, yCol = tillers_jul)
evSEvSD1Sim[[2]]

evAEvSD1Sim <- mod_fit_fun(dat = evAEvSD1Dat, mod = evAEvSD1Mod, treatCol = treatment,
                           xCol = density, minX = 0, maxX = 19, yCol = tillers_jul)
evAEvSD1Sim[[2]]

# posterior means
post_pred_fun(mvEvSD1Mod)
post_pred_fun(evSEvSD1Mod)
post_pred_fun(evAEvSD1Mod)

# save models
save(mvEvSD1Mod, file = "output/mv_tillers_evS_density_model_2018_density_exp.rda")
save(evSEvSD1Mod, file = "output/evS_tillers_evS_density_model_2018_density_exp.rda")
save(evAEvSD1Mod, file = "output/evA_tillers_evS_density_model_2018_density_exp.rda")


#### 2018 Ev adult density models ####

# initial visualization
ggplot(evAD1Dat, aes(density, tillers_jul, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ plant_group, scales = "free")

# intercept prior
evAD1Dat2 %>%
  filter(density == 3) %>%
  group_by(plant_group, treatment) %>%
  summarise(bio = mean(tillers_jul))

prior_fun(maxX = 50, shape = 35, dist = "gamma")
prior_fun(maxX = 15, shape = 5, dist = "gamma")
prior_fun(maxX = 40, shape = 15, dist = "gamma")

# alpha prior
prior_fun(maxX = 10, shape = 0.5, dist = "exponential")

# fit models
mvEvAD1Mod <- update(mvMvD1Mod, newdata = mvEvAD1Dat,
                     prior = set_prior("gamma(35, 1)", nlpar = "b0", lb = 0))
mod_check_fun(mvEvAD1Mod)

evSEvAD1Mod <- update(mvEvAD1Mod, newdata = evSEvAD1Dat,
                      prior = set_prior("gamma(5, 1)", nlpar = "b0", lb = 0),
                      control = list(adapt_delta = 0.99))
mod_check_fun(evSEvAD1Mod)

evAEvAD1Mod <- update(mvEvAD1Mod, newdata = evAEvAD1Dat,
                      prior = set_prior("gamma(15, 1)", nlpar = "b0", lb = 0),
                      control = list(adapt_delta = 0.99))
mod_check_fun(evAEvAD1Mod)

# simulated data
mvEvAD1Sim <- mod_fit_fun(dat = mvEvAD1Dat, mod = mvEvAD1Mod, treatCol = treatment,
                          xCol = density, minX = 0, maxX = 9, yCol = tillers_jul)
mvEvAD1Sim[[2]]

evSEvAD1Sim <- mod_fit_fun(dat = evSEvAD1Dat, mod = evSEvAD1Mod, treatCol = treatment,
                           xCol = density, minX = 0, maxX = 9, yCol = tillers_jul)
evSEvAD1Sim[[2]]

evAEvAD1Sim <- mod_fit_fun(dat = evAEvAD1Dat, mod = evAEvAD1Mod, treatCol = treatment,
                           xCol = density, minX = 0, maxX = 9, yCol = tillers_jul)
evAEvAD1Sim[[2]]

# posterior means
post_pred_fun(mvEvAD1Mod)
post_pred_fun(evSEvAD1Mod)
post_pred_fun(evAEvAD1Mod)

# save models
save(mvEvAD1Mod, file = "output/mv_tillers_evA_density_model_2018_density_exp.rda")
save(evSEvAD1Mod, file = "output/evS_tillers_evA_density_model_2018_density_exp.rda")
save(evAEvAD1Mod, file = "output/evA_tillers_evA_density_model_2018_density_exp.rda")



#### presentation figure (2019 Mv density) ####

# combine predicted datasets
mvD2Sim <- mvMvD2Sim[[1]] %>%
  mutate(plant_group = "Mv seedling") %>%
  full_join(evSMvD2Sim[[1]] %>%
              mutate(plant_group = "Ev seedling")) %>%
  full_join(evAMvD2Sim[[1]] %>%
              mutate(plant_group = "Ev adult")) %>%
  mutate(plant_group = fct_rev(plant_group),
         treatment = dplyr::recode(treatment, "water" = "water (control)") %>%
           fct_rev())

# labels
labDat <- mvD2Sim %>%
  group_by(plant_group) %>%
  summarise(biomass_weight.g = max(upper)) %>%
  ungroup() %>%
  mutate(treatment = "fungicide")

# rename factors
mvD2Dat3 <- mvD2Dat %>%
  mutate(treatment = dplyr::recode(treatment, "water" = "water (control)") %>%
           fct_rev())

pdf("output/plant_growth_density_2019_density_exp.pdf", width = 10.5, height = 5)
ggplot(mvD2Dat3, aes(x = density, y = biomass_weight.g, color = treatment, fill = treatment)) +
  geom_ribbon(data = mvD2Sim, aes(y = pred, ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
  geom_line(data = mvD2Sim, aes(y = pred)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  geom_text(data = labDat, x = 0, hjust = 0, color = "black", aes(label = plant_group), size = 5) +
  facet_wrap(~plant_group, scales = "free") +
  xlab("Mv seedling density") +
  ylab("Plant biomass (g)") +
  scale_color_viridis_d(end = 0.6, name = "Treatment") +
  scale_fill_viridis_d(end = 0.6, name = "Treatment") +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = c(0.22, 0.85),
        strip.background = element_blank(),
        strip.text = element_blank())
dev.off()


#### values for text ####

mvD2Sim %>%
  filter(density == 67) %>%
  mutate(treatment = dplyr::recode(treatment, "water (control)" = "water")) %>%
  select(-c(lower, upper)) %>%
  pivot_wider(names_from = treatment,
              values_from = pred) %>%
  mutate(diff = (fungicide - water)/water)

mvD2Dat %>%
  filter(density == 67) %>%
  select(site, plot, treatment, ID, plant_group, biomass_weight.g) %>%
  pivot_wider(names_from = treatment,
              values_from = biomass_weight.g) %>%
  filter(!is.na(fungicide) & !(is.na(water))) %>%
  group_by(plant_group) %>%
  summarise(diff = (mean(fungicide) - mean(water)) / mean(water))


#### output ####
save(mvMvD2Mod, file = "output/mv_biomass_mv_density_model_2019_density_exp.rda")
save(mvEvSD2Mod, file = "output/evS_biomass_mv_density_model_2019_density_exp.rda")
save(mvEvAD2Mod, file = "output/evA_biomass_mv_density_model_2019_density_exp.rda")