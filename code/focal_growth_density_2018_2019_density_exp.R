##### info ####

# file: focal_growth_density_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 6/3/21
# goal: analyses of plant growth as a function of density

#### approaches ####

# first approach is Beverton-Holt, a discrete time competition model
# I would like the competition coefficients to translate to a continuous time model
# Ricker model:
# log(Nt+1/Nt) = r - alphaNN x Nt - alphaNM x Mt


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(tidybayes) # for mean_hdi
library(cowplot)

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
  select(plot, treatment, background, density_level, Mv_seedling_density, Ev_seedling_density, Ev_adult_density)

# missing data
filter(growthD1Dat, is.na(tillers_jul))
filter(mvBioD2Dat, is.na(biomass_weight.g)) # 3
filter(evBioD2Dat, is.na(weight)) # 1 seedling

# combine and separate data
d1dat <- growthD1Dat %>%
  left_join(plotDens) %>%
  mutate(density = case_when(plot %in% 2:4 ~ Mv_seedling_density,
                             plot %in% 5:7 ~ Ev_seedling_density,
                             plot%in% 8:10 ~ Ev_adult_density,
                             plot == 1 ~ NA_real_),
         plant_growth = log(tillers_jul/tillers_jun),
         age = ifelse(ID == "A", "adult", "seedling"),
         plant_group = paste(sp, age, sep = "_") %>%
           fct_rev(),
         background = str_replace(background, " ", "_") %>%
           fct_rev(),
         fungicide = ifelse(treatment == "fungicide", 1, 0))

d2dat <- mvBioD2Dat %>%
  rename(ID = plant) %>%
  mutate(ID = as.character(ID)) %>%
  full_join(evBioD2Dat %>%
              rename(biomass_weight.g = weight)) %>%
  left_join(plotDens) %>%
  mutate(density = case_when(plot %in% 2:4 ~ Mv_seedling_density,
                             plot %in% 5:7 ~ Ev_seedling_density,
                             plot%in% 8:10 ~ Ev_adult_density,
                             plot == 1 ~ NA_real_),
         plant_growth = log(biomass_weight.g),
         age = ifelse(ID == "A", "adult", "seedling"),
         plant_group = paste(sp, age, sep = "_") %>%
           fct_rev(),
         background = str_replace(background, " ", "_") %>%
           fct_rev(),
         fungicide = ifelse(treatment == "fungicide", 1, 0)) %>%
  filter(!is.na(biomass_weight.g))

# split by focal
# remove plot 1 (large environmental effect)
mvD1Dat <- d1dat %>%
  filter(plant_group == "Mv_seedling" & plot != 1) %>%
  mutate(background = fct_drop(background))
evSD1Dat <- d1dat %>%
  filter(plant_group == "Ev_seedling" & plot != 1) %>%
  mutate(background = fct_drop(background))
evAD1Dat <- d1dat %>%
  filter(plant_group == "Ev_adult" & plot != 1) %>%
  mutate(background = fct_drop(background))

mvD2Dat <- d2dat %>%
  filter(plant_group == "Mv_seedling" & plot != 1) %>%
  mutate(background = fct_drop(background))
evSD2Dat <- d2dat %>%
  filter(plant_group == "Ev_seedling" & plot != 1) %>%
  mutate(background = fct_drop(background))
evAD2Dat <- d2dat %>%
  filter(plant_group == "Ev_adult" & plot != 1) %>%
  mutate(background = fct_drop(background))

# split by species pairs
mvMvD1Dat <- mvD1Dat %>%
  filter(plot %in% 2:4)
evSMvD1Dat <- evSD1Dat %>%
  filter(plot %in% 2:4)
evAMvD1Dat <- evAD1Dat %>%
  filter(plot %in% 2:4)

mvEvSD1Dat <- mvD1Dat %>%
  filter(plot %in% 5:7)
evSEvSD1Dat <- evSD1Dat %>%
  filter(plot %in% 5:7)
evAEvSD1Dat <- evAD1Dat %>%
  filter(plot %in% 5:7)

mvEvAD1Dat <- mvD1Dat %>%
  filter(plot %in% 8:10)
evSEvAD1Dat <- evSD1Dat %>%
  filter(plot %in% 8:10)
evAEvAD1Dat <- evAD1Dat %>%
  filter(plot %in% 8:10)

mvMvD2Dat <- mvD2Dat %>%
  filter(plot %in% 2:4)
evSMvD2Dat <- evSD2Dat %>%
  filter(plot %in% 2:4)
evAMvD2Dat <- evAD2Dat %>%
  filter(plot %in% 2:4)

mvEvSD2Dat <- mvD2Dat %>%
  filter(plot %in% 5:7)
evSEvSD2Dat <- evSD2Dat %>%
  filter(plot %in% 5:7)
evAEvSD2Dat <- evAD2Dat %>%
  filter(plot %in% 5:7)

mvEvAD2Dat <- mvD2Dat %>%
  filter(plot %in% 8:10)
evSEvAD2Dat <- evSD2Dat %>%
  filter(plot %in% 8:10)
evAEvAD2Dat <- evAD2Dat %>%
  filter(plot %in% 8:10)


#### 2019 Mv models ####

# initial visualization
ggplot(mvD2Dat, aes(density, biomass_weight.g, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ background, scales = "free")

ggplot(mvD2Dat, aes(density, plant_growth, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ background, scales = "free")

# intercept prior
mvD2Dat %>%
  filter(density_level == "low") %>%
  group_by(background, treatment) %>%
  summarise(bio = mean(biomass_weight.g))

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

mvEvSD2Mod <- update(mvMvD2Mod, newdata = mvEvSD2Dat,
                     prior = set_prior("gamma(22, 1)", nlpar = "b0", lb = 0))
mod_check_fun(mvEvSD2Mod)

mvEvAD2Mod <- update(mvMvD2Mod, newdata = mvEvAD2Dat,
                     prior = set_prior("gamma(35, 1)", nlpar = "b0", lb = 0))
mod_check_fun(mvEvAD2Mod)

mvRickD2Mod <- brm(data = mvD2Dat, family = gaussian,
                   plant_growth ~ fungicide * density * background + (1|site),
                   prior <- c(prior(normal(3, 1), class = "Intercept"),
                              prior(normal(0, 1), class = "b")), # use default for sigma
                   iter = 6000, warmup = 1000, chains = 3, 
                   control = list(adapt_delta = 0.9999)) 
mod_check_fun(mvRickD2Mod)

# simulated data
mvMvD2Sim <- mod_fit_fun(dat = mvMvD2Dat, mod = mvMvD2Mod, treatCol = treatment,
                         xCol = density, minX = 0, maxX = 67, yCol = biomass_weight.g)
mvMvD2Sim[[2]]
mvEvSD2Sim <- mod_fit_fun(dat = mvEvSD2Dat, mod = mvEvSD2Mod, treatCol = treatment,
                          xCol = density, minX = 0, maxX = 19, yCol = biomass_weight.g)
mvEvSD2Sim[[2]]
mvEvAD2Sim <- mod_fit_fun(dat = mvEvAD2Dat, mod = mvEvAD2Mod, treatCol = treatment,
                          xCol = density, minX = 0, maxX = 9, yCol = biomass_weight.g)
mvEvAD2Sim[[2]]
mvRickD2Sim <- mod_fit_fun2(dat = mvD2Dat, mod = mvRickD2Mod, treatCol1 = fungicide, treatCol2 = background,
                         xCol = density, minX = 0, maxX = 67, yCol = plant_growth, f2t = T)
mvRickD2Sim[[2]]

# posterior means
post_pred_fun(mvMvD2Mod)
post_pred_fun(mvEvSD2Mod)
post_pred_fun2(mvEvSD2Mod)
post_pred_fun(mvEvAD2Mod)

# save models
save(mvMvD2Mod, file = "output/mv_biomass_mv_density_model_2019_density_exp.rda")
save(mvEvSD2Mod, file = "output/mv_biomass_evS_density_model_2019_density_exp.rda")
save(mvEvAD2Mod, file = "output/mv_biomass_evA_density_model_2019_density_exp.rda")
save(mvRickD2Mod, file = "output/mv_biomass_Ricker_model_2019_density_exp.rda")


#### 2019 Ev seedling models ####

# initial visualization
ggplot(evSD2Dat, aes(density, biomass_weight.g, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ background, scales = "free")

ggplot(evSD2Dat, aes(density, plant_growth, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ background, scales = "free")

# intercept prior
evSD2Dat %>%
  filter(density_level == "low") %>%
  group_by(background, treatment) %>%
  summarise(bio = mean(biomass_weight.g))

prior_fun(maxX = 10, shape = 2.5, dist = "gamma")

# alpha prior
prior_fun(maxX = 10, shape = 0.5, dist = "exponential")

# fit model
evSEvSD2Mod <- update(mvEvSD2Mod, newdata = evSEvSD2Dat,
                      prior = set_prior("gamma(2.5, 1)", nlpar = "b0", lb = 0),
                      control = list(adapt_delta = 0.99))
mod_check_fun(evSEvSD2Mod)

evSMvD2Mod <- update(mvMvD2Mod, newdata = evSMvD2Dat,
                     prior = set_prior("gamma(2.5, 1)", nlpar = "b0", lb = 0),
                     control = list(adapt_delta = 0.99))
mod_check_fun(evSMvD2Mod)

evSEvAD2Mod <- update(mvEvAD2Mod, newdata = evSEvAD2Dat,
                      prior = set_prior("gamma(3.5, 1)", nlpar = "b0", lb = 0),
                      control = list(adapt_delta = 0.99))
mod_check_fun(evSEvAD2Mod)

evSRickD2Mod <- update(mvRickD2Mod, newdata = evSD2Dat,
                     prior = set_prior("normal(0.8, 1)", class = "Intercept"))
mod_check_fun(evSRickD2Mod)

# simulated data
evSEvSD2Sim <- mod_fit_fun(dat = evSEvSD2Dat, mod = evSEvSD2Mod, treatCol = treatment,
                           xCol = density, minX = 0, maxX = 19, yCol = biomass_weight.g)
evSEvSD2Sim[[2]]
evSMvD2Sim <- mod_fit_fun(dat = evSMvD2Dat, mod = evSMvD2Mod, treatCol = treatment,
                          xCol = density, minX = 0, maxX = 67, yCol = biomass_weight.g)
evSMvD2Sim[[2]]
evSEvAD2Sim <- mod_fit_fun(dat = evSEvAD2Dat, mod = evSEvAD2Mod, treatCol = treatment,
                           xCol = density, minX = 0, maxX = 9, yCol = biomass_weight.g)
evSEvAD2Sim[[2]]
evSRickD2Sim <- mod_fit_fun2(dat = evSD2Dat, mod = evSRickD2Mod, treatCol1 = fungicide, treatCol2 = background,
                            xCol = density, minX = 0, maxX = 67, yCol = plant_growth, f2t = T)
evSRickD2Sim[[2]]

# posterior means
post_pred_fun(evSEvSD2Mod)
post_pred_fun(evSMvD2Mod)
post_pred_fun(evSEvAD2Mod)

# save models
save(evSEvSD2Mod, file = "output/evS_biomass_evS_density_model_2019_density_exp.rda")
save(evSMvD2Mod, file = "output/evS_biomass_mv_density_model_2019_density_exp.rda")
save(evSEvAD2Mod, file = "output/evS_biomass_evA_density_model_2019_density_exp.rda")
save(evSRickD2Mod, file = "output/evS_biomass_Ricker_model_2019_density_exp.rda")


#### 2019 Ev adult models ####

# initial visualization
ggplot(evAD2Dat, aes(density, biomass_weight.g, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ background, scales = "free")

ggplot(evAD2Dat, aes(density, plant_growth, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ background, scales = "free")

# intercept prior
evAD2Dat %>%
  filter(density_level == "low") %>%
  group_by(background, treatment) %>%
  summarise(bio = mean(biomass_weight.g))

prior_fun(maxX = 50, shape = 35, dist = "gamma")
prior_fun(maxX = 10, shape = 3.5, dist = "gamma")
prior_fun(maxX = 20, shape = 12, dist = "gamma")

# alpha prior
prior_fun(maxX = 10, shape = 0.5, dist = "exponential")

# fit models
evAEvAD2Mod <- update(mvEvAD2Mod, newdata = evAEvAD2Dat,
                      prior = set_prior("gamma(12, 1)", nlpar = "b0", lb = 0),
                      control = list(adapt_delta = 0.99))
mod_check_fun(evAEvAD2Mod)

evAMvD2Mod <- update(mvMvD2Mod, newdata = evAMvD2Dat,
                     control = list(adapt_delta = 0.99))
# initially ran with prior shape = 8, but the predictions were much lower than the data
mod_check_fun(evAMvD2Mod)

evAEvSD2Mod <- update(mvEvSD2Mod, newdata = evAEvSD2Dat,
                      prior = set_prior("gamma(10, 1)", nlpar = "b0", lb = 0),
                      control = list(adapt_delta = 0.99))
mod_check_fun(evAEvSD2Mod)

evARickD2Mod <- update(mvRickD2Mod, newdata = evAD2Dat,
                       prior = set_prior("normal(1.5, 1)", class = "Intercept"))
mod_check_fun(evARickD2Mod)


# simulated data
evAEvAD2Sim <- mod_fit_fun(dat = evAEvAD2Dat, mod = evAEvAD2Mod, treatCol = treatment,
                           xCol = density, minX = 0, maxX = 9, yCol = biomass_weight.g)
evAEvAD2Sim[[2]]
evAMvD2Sim <- mod_fit_fun(dat = evAMvD2Dat, mod = evAMvD2Mod, treatCol = treatment,
                          xCol = density, minX = 0, maxX = 67, yCol = biomass_weight.g)
evAMvD2Sim[[2]]
evAEvSD2Sim <- mod_fit_fun(dat = evAEvSD2Dat, mod = evAEvSD2Mod, treatCol = treatment,
                           xCol = density, minX = 0, maxX = 19, yCol = biomass_weight.g)
evAEvSD2Sim[[2]]
evARickD2Sim <- mod_fit_fun2(dat = evAD2Dat, mod = evARickD2Mod, treatCol1 = fungicide, treatCol2 = background,
                             xCol = density, minX = 0, maxX = 67, yCol = plant_growth, f2t = T)
evARickD2Sim[[2]]

# posterior means
post_pred_fun(evAEvAD2Mod)
post_pred_fun(evAMvD2Mod)
post_pred_fun(evAEvSD2Mod)

# save models
save(evAEvAD2Mod, file = "output/evA_biomass_evA_density_model_2019_density_exp.rda")
save(evAMvD2Mod, file = "output/evA_biomass_mv_density_model_2019_density_exp.rda")
save(evAEvSD2Mod, file = "output/evA_biomass_evS_density_model_2019_density_exp.rda")
save(evARickD2Mod, file = "output/evA_biomass_Ricker_model_2019_density_exp.rda")


#### 2018 Mv models ####

# initial visualization
ggplot(mvD1Dat, aes(density, tillers_jul, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ background, scales = "free")

ggplot(mvD1Dat, aes(density, plant_growth, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ background, scales = "free")

# intercept prior
mvD1Dat %>%
  filter(density_level == "low") %>%
  group_by(background, treatment) %>%
  summarise(bio = mean(tillers_jul))

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

mvEvSD1Mod <- update(mvMvD1Mod, newdata = mvEvSD1Dat,
                     prior = set_prior("gamma(37, 1)", nlpar = "b0", lb = 0))
mod_check_fun(mvEvSD1Mod)

mvEvAD1Mod <- update(mvMvD1Mod, newdata = mvEvAD1Dat,
                     prior = set_prior("gamma(35, 1)", nlpar = "b0", lb = 0))
mod_check_fun(mvEvAD1Mod)

mvD1Dat2 <- mvD1Dat %>% filter(!is.na(plant_growth) & plant_growth != Inf & plant_growth != -Inf) # remove zero tiller counts
mvRickD1Mod <- update(mvRickD2Mod, newdata = mvD1Dat2,
                       prior = set_prior("normal(1, 1)", class = "Intercept"))
mod_check_fun(mvRickD1Mod)

# simulated data
mvMvD1Sim <- mod_fit_fun(dat = mvMvD1Dat, mod = mvMvD1Mod, treatCol = treatment,
                         xCol = density, minX = 0, maxX = 67, yCol = tillers_jul)
mvMvD1Sim[[2]]
mvEvSD1Sim <- mod_fit_fun(dat = mvEvSD1Dat, mod = mvEvSD1Mod, treatCol = treatment,
                          xCol = density, minX = 0, maxX = 19, yCol = tillers_jul)
mvEvSD1Sim[[2]]
mvEvAD1Sim <- mod_fit_fun(dat = mvEvAD1Dat, mod = mvEvAD1Mod, treatCol = treatment,
                          xCol = density, minX = 0, maxX = 9, yCol = tillers_jul)
mvEvAD1Sim[[2]]
mvRickD1Sim <- mod_fit_fun2(dat = mvD1Dat2, mod = mvRickD1Mod, treatCol1 = fungicide, treatCol2 = background,
                             xCol = density, minX = 0, maxX = 67, yCol = plant_growth, f2t = T)
mvRickD1Sim[[2]]

# posterior means
post_pred_fun(mvMvD1Mod)
post_pred_fun(mvEvSD1Mod)
post_pred_fun(mvEvAD1Mod)

# save models
save(mvMvD1Mod, file = "output/mv_tillers_mv_density_model_2018_density_exp.rda")
save(mvEvSD1Mod, file = "output/mv_tillers_evS_density_model_2018_density_exp.rda")
save(mvEvAD1Mod, file = "output/mv_tillers_evA_density_model_2018_density_exp.rda")
save(mvRickD1Mod, file = "output/mv_tillers_Ricker_model_2018_density_exp.rda")


#### 2018 Ev seedling models ####

# initial visualization
ggplot(evSD1Dat, aes(density, tillers_jul, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ background, scales = "free")

ggplot(evSD1Dat, aes(density, plant_growth, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ background, scales = "free")

# intercept prior
evSD1Dat %>%
  filter(density_level == "low") %>%
  group_by(background, treatment) %>%
  summarise(bio = mean(tillers_jul))

prior_fun(maxX = 15, shape = 3.5, dist = "gamma")

# alpha prior
prior_fun(maxX = 10, shape = 0.5, dist = "exponential")

# fit models
evSEvSD1Mod <- update(mvEvSD1Mod, newdata = evSEvSD1Dat,
                      prior = set_prior("gamma(3.5, 1)", nlpar = "b0", lb = 0),
                      control = list(adapt_delta = 0.99))
mod_check_fun(evSEvSD1Mod)

evSMvD1Mod <- update(mvMvD1Mod, newdata = evSMvD1Dat,
                     prior = set_prior("gamma(7, 1)", nlpar = "b0", lb = 0),
                     control = list(adapt_delta = 0.99))
mod_check_fun(evSMvD1Mod)

evSEvAD1Mod <- update(mvEvAD1Mod, newdata = evSEvAD1Dat,
                      prior = set_prior("gamma(5, 1)", nlpar = "b0", lb = 0),
                      control = list(adapt_delta = 0.99))
mod_check_fun(evSEvAD1Mod)

evSRickD1Mod <- update(mvRickD2Mod, newdata = evSD1Dat,
                      prior = set_prior("normal(1, 1)", class = "Intercept"))
mod_check_fun(evSRickD1Mod)

evSD1Dat2 <- evSD1Dat %>% filter(!is.na(plant_growth) & plant_growth != Inf & plant_growth != -Inf) # remove zero tiller counts
evSRickD1Mod <- update(mvRickD2Mod, newdata = evSD1Dat2,
                      prior = set_prior("normal(-0.3, 1)", class = "Intercept"),
                      control = list(adapt_delta = 0.99999))
mod_check_fun(evSRickD1Mod)

# simulated data
evSEvSD1Sim <- mod_fit_fun(dat = evSEvSD1Dat, mod = evSEvSD1Mod, treatCol = treatment,
                           xCol = density, minX = 0, maxX = 19, yCol = tillers_jul)
evSEvSD1Sim[[2]]
evSMvD1Sim <- mod_fit_fun(dat = evSMvD1Dat, mod = evSMvD1Mod, treatCol = treatment,
                          xCol = density, minX = 0, maxX = 67, yCol = tillers_jul)
evSMvD1Sim[[2]]
evSEvAD1Sim <- mod_fit_fun(dat = evSEvAD1Dat, mod = evSEvAD1Mod, treatCol = treatment,
                           xCol = density, minX = 0, maxX = 9, yCol = tillers_jul)
evSEvAD1Sim[[2]]
evSRickD1Sim <- mod_fit_fun2(dat = evSD1Dat2, mod = evSRickD1Mod, treatCol1 = fungicide, treatCol2 = background,
                            xCol = density, minX = 0, maxX = 67, yCol = plant_growth, f2t = T)
evSRickD1Sim[[2]]

# posterior means
post_pred_fun(evSEvSD1Mod)
post_pred_fun(evSMvD1Mod)
post_pred_fun(evSEvAD1Mod)

# save models
save(evSEvSD1Mod, file = "output/evS_tillers_evS_density_model_2018_density_exp.rda")
save(evSMvD1Mod, file = "output/evS_tillers_mv_density_model_2018_density_exp.rda")
save(evSEvAD1Mod, file = "output/evS_tillers_evA_density_model_2018_density_exp.rda")
save(evSRickD1Mod, file = "output/evS_tillers_Ricker_model_2018_density_exp.rda")


#### 2018 Ev adult models ####

# initial visualization
ggplot(evAD1Dat, aes(density, tillers_jul, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ background, scales = "free")

ggplot(evAD1Dat, aes(density, plant_growth, color = treatment)) +
  stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0) +
  stat_summary(geom = "line", fun = "mean") +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ background, scales = "free")

# intercept prior
evAD1Dat %>%
  filter(density_level == "low") %>%
  group_by(background, treatment) %>%
  summarise(bio = mean(tillers_jul))

prior_fun(maxX = 15, shape = 5, dist = "gamma")

# alpha prior
prior_fun(maxX = 10, shape = 0.5, dist = "exponential")

# fit models
evAEvAD1Mod <- update(mvEvAD1Mod, newdata = evAEvAD1Dat,
                      prior = set_prior("gamma(15, 1)", nlpar = "b0", lb = 0),
                      control = list(adapt_delta = 0.99))
mod_check_fun(evAEvAD1Mod)

evAMvD1Mod <- update(mvMvD1Mod, newdata = evAMvD1Dat,
                     prior = set_prior("gamma(25, 1)", nlpar = "b0", lb = 0),
                     control = list(adapt_delta = 0.99))
mod_check_fun(evAMvD1Mod)

evAEvSD1Mod <- update(mvEvSD1Mod, newdata = evAEvSD1Dat,
                      prior = set_prior("gamma(17, 1)", nlpar = "b0", lb = 0),
                      control = list(adapt_delta = 0.99))
mod_check_fun(evAEvSD1Mod)

evAD1Dat2 <- evAD1Dat %>% filter(!is.na(plant_growth) & plant_growth != Inf & plant_growth != -Inf) # remove zero tiller counts
evARickD1Mod <- update(mvRickD2Mod, newdata = evAD1Dat2,
                       prior = set_prior("normal(0.2, 1)", class = "Intercept"),
                       control = list(adapt_delta = 0.999999, max_treedepth = 15))
mod_check_fun(evARickD1Mod)

# simulated data
evAEvAD1Sim <- mod_fit_fun(dat = evAEvAD1Dat, mod = evAEvAD1Mod, treatCol = treatment,
                           xCol = density, minX = 0, maxX = 9, yCol = tillers_jul)
evAEvAD1Sim[[2]]
evAMvD1Sim <- mod_fit_fun(dat = evAMvD1Dat, mod = evAMvD1Mod, treatCol = treatment,
                          xCol = density, minX = 0, maxX = 67, yCol = tillers_jul)
evAMvD1Sim[[2]]
evAEvSD1Sim <- mod_fit_fun(dat = evAEvSD1Dat, mod = evAEvSD1Mod, treatCol = treatment,
                           xCol = density, minX = 0, maxX = 19, yCol = tillers_jul)
evAEvSD1Sim[[2]]
evARickD1Sim <- mod_fit_fun2(dat = evAD1Dat2, mod = evARickD1Mod, treatCol1 = fungicide, treatCol2 = background,
                             xCol = density, minX = 0, maxX = 67, yCol = plant_growth, f2t = T)
evARickD1Sim[[2]]

# posterior means
post_pred_fun(evAEvAD1Mod)
post_pred_fun(evAMvD1Mod)
post_pred_fun(evAEvSD1Mod)

# save models
save(evAEvAD1Mod, file = "output/evA_tillers_evA_density_model_2018_density_exp.rda")
save(evAMvD1Mod, file = "output/evA_tillers_mv_density_model_2018_density_exp.rda")
save(evAEvSD1Mod, file = "output/evA_tillers_evS_density_model_2018_density_exp.rda")
save(evARickD1Mod, file = "output/evA_tillers_Ricker_model_2018_density_exp.rda")


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


#### Ricker fit figure ####

# prediction function
pair_pred_fun <- function(year, foc, bg, trt){
  
  # data
  if(year == "2018"){
    dat <- d1dat %>% filter(plant_group == foc & background == bg & treatment == trt)
  }else{
    dat <- d2dat %>% filter(plant_group == foc & background == bg & treatment == trt)
  }
  
  # model
  if(year == "2018" & foc == "Mv_seedling"){
    mod <- mvRickD1Mod
  }else if(year == "2018" & foc == "Ev_seedling"){
    mod <- evSRickD1Mod
  }else if(year == "2018" & foc == "Ev_adult"){
    mod <- evARickD1Mod
  }else if(year == "2019" & foc == "Mv_seedling"){
    mod <- mvRickD2Mod
  }else if(year == "2019" & foc == "Ev_seedling"){
    mod <- evSRickD2Mod
  }else if(year == "2019" & foc == "Ev_adult"){
    mod <- evARickD2Mod
  }
  
  # sequence of density values
  # other variables
  # predicted values
  simDat <- tibble(density = seq(min(dat$density), max(dat$density), length.out = 100)) %>%
    mutate(background = bg,
           plotf = "A",
           fungicide = case_when(trt == "water" ~ 0,
                                 trt == "fungicide" ~ 1)) %>%
    mutate(plant_growth = fitted(mod, newdata = ., allow_new_levels = T)[, "Estimate"],
           lower = fitted(mod, newdata = ., allow_new_levels = T)[, "Q2.5"],
           upper = fitted(mod, newdata = ., allow_new_levels = T)[, "Q97.5"]) %>%
    select(-background)
  
  # output
  return(simDat)
  
}

# predicted data
predDat <- tibble(focal = c("Ev_adult", "Ev_seedling", "Mv_seedling")) %>%
  expand_grid(tibble(background = c("Ev_adult", "Ev_seedling", "Mv_seedling"))) %>%
  expand_grid(tibble(treatment = c("water", "fungicide"))) %>%
  expand_grid(tibble(year = c("2018", "2019"))) %>%
  mutate(pred = pmap(list(year, focal, background, treatment), pair_pred_fun)) %>%
  unnest(pred) %>%
  mutate(focal = fct_recode(focal, "Mv" = "Mv_seedling") %>%
           str_replace("_", " "),
         background = fct_recode(background, "Mv" = "Mv_seedling") %>%
           str_replace("_", " "),
         treatment = fct_recode(treatment, "control (water)" = "water") %>%
           fct_relevel("control (water)"))

# combine alphas
alphasDat <- posterior_samples(mvRickD1Mod) %>%
  mutate(focal = "a", year = "2018") %>%
  full_join(posterior_samples(evSRickD1Mod) %>%
              mutate(focal = "s", year = "2018")) %>%
  full_join(posterior_samples(evARickD1Mod) %>%
              mutate(focal = "p", year = "2018")) %>%
  full_join(posterior_samples(mvRickD2Mod) %>%
              mutate(focal = "a", year = "2019")) %>%
  full_join(posterior_samples(evSRickD2Mod) %>%
              mutate(focal = "s", year = "2019")) %>%
  full_join(posterior_samples(evARickD2Mod) %>%
              mutate(focal = "p", year = "2019")) %>%
  rename_with(str_replace, pattern = ":", replacement = "_") %>%
  rename_with(str_replace, pattern = ":", replacement = "_") %>%
  mutate(a_water = b_density,
         s_water = b_density + b_density_backgroundEv_seedling,
         p_water = b_density + b_density_backgroundEv_adult,
         a_fungicide = b_density + b_fungicide_density,
         s_fungicide = s_water + b_fungicide_density + b_fungicide_density_backgroundEv_seedling,
         p_fungicide = p_water + b_fungicide_density + b_fungicide_density_backgroundEv_adult) %>%
  select(year, focal, a_water, s_water, p_water, a_fungicide, s_fungicide, p_fungicide) %>%
  as_tibble() %>%
  pivot_longer(cols = -c(year, focal),
               names_to = c(".value", "treatment"),
               names_pattern = "(.)_(.*)") %>%
  pivot_longer(cols = -c(year, focal, treatment),
               names_to = "background",
               values_to = "alpha") %>%
  group_by(year, treatment, focal, background) %>%
  mean_hdi(alpha) %>%
  ungroup() %>%
  mutate(sig = case_when((.lower < 0 & .upper < 0) | (.lower > 0 & .upper > 0) ~ "omits 0",
                         TRUE ~ "includes 0"),
         focal = fct_recode(focal, "Ev adult" = "p", "Ev seedling" = "s", "Mv" = "a"),
         background = fct_recode(background, "Ev adult" = "p", "Ev seedling" = "s", "Mv" = "a"),
         treatment = fct_recode(treatment, "control (water)" = "water") %>%
           fct_relevel("control (water)"))

# add to predicted dataset
predDat2 <- predDat %>%
  inner_join(alphasDat %>%
              select(year, treatment, focal, background, sig))

# raw data
figDat <- d1dat %>%
  select(site, plot, treatment, sp, age, ID, plant_group, background, density, plant_growth) %>%
  mutate(year = "2018") %>%
  full_join(d2dat %>%
              select(site, plot, treatment, sp, age, ID, plant_group, background, density, plant_growth) %>%
              mutate(year = "2019"))%>%
  mutate(focal = str_replace(plant_group, "_", " ") %>%
           fct_recode(Mv = "Mv seedling"),
         background = str_replace(background, "_", " ") %>%
           fct_recode(Mv = "Mv seedling"),
         treatment = fct_recode(treatment, "control (water)" = "water") %>%
           fct_relevel("control (water)")) %>%
  filter(background != "none" & !is.na(plant_growth) & plant_growth != Inf & plant_growth != -Inf)

# sig alpha values
alphasDat2 <- alphasDat %>%
  filter(sig == "omits 0") %>%
  left_join(predDat2 %>%
              group_by(year, background) %>%
              summarise(density = max(density))) %>%
  left_join(predDat2 %>%
              group_by(year, focal) %>%
              summarise(upper = max(upper)) %>%
              ungroup() %>%
              full_join(figDat %>%
                          group_by(year, focal) %>%
                          summarise(raw = max(plant_growth))) %>%
              rowwise() %>%
              mutate(plant_growth = max(c(raw, upper))) %>%
              ungroup() %>%
              select(-c(raw, upper))) %>%
  rename(param = alpha) %>%
  mutate(param = round(param, 2),
         plant_growth = case_when(treatment == "fungicide" & background == "Mv" & focal == "Ev seedling" ~ plant_growth - 0.5,
                                  treatment == "fungicide" & background == "Mv" & focal == "Mv" ~ plant_growth - 0.6,
                                  TRUE ~ plant_growth + 0.1))

# figure settings
fig_theme <- theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.background = element_blank(),
        legend.position = "none",
        legend.margin = margin(-0.1, 0, 0.2, 2, unit = "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        strip.placement = "outside")

col_pal = c("black", "#238A8DFF")

yearText <- tibble(year = c("2018", "2019"),
                   density = c(3.2, 3),
                   plant_growth = c(0.87, 2.7),
                   background = "Ev adult",
                   focal = "Ev adult",
                   treatment = "fungicide")

textSize = 3

# water figure
pairD1Fig <- ggplot(filter(predDat2, year == "2018"), aes(x = density, y = plant_growth)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment), alpha = 0.3) +
  geom_line(aes(color = treatment, linetype = sig)) +
  geom_point(data = filter(figDat, year == "2018"), aes(color = treatment), alpha = 0.5, size = 0.5) +
  geom_text(data = filter(alphasDat2, year == "2018"), 
            aes(label = paste("alpha", " == ", param, sep = ""), 
                color = treatment), parse = T, hjust = 1, vjust = 1, show.legend = F, size = textSize) +
  geom_text(data = filter(yearText, year == "2018"), aes(label = year), color = "black", size = textSize, hjust = 0) +
  facet_grid(rows = vars(focal),
             cols = vars(background),
             scales = "free",
             switch = "both") +
  scale_linetype_manual(values = c("dashed", "solid"), guide = F) +
  scale_color_manual(values = col_pal, name = "Treatment") +
  scale_fill_manual(values = col_pal, name = "Treatment") +
  xlab(expression(paste("Competitor density (", m^-2, ")", sep = ""))) +
  ylab("Plant growth") +
  fig_theme +
  theme(legend.position = "bottom",
        legend.direction = "horizontal")

# fungicide figure
pairD2Fig <- ggplot(filter(predDat2, year == "2019"), aes(x = density, y = plant_growth)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treatment), alpha = 0.3) +
  geom_line(aes(color = treatment, linetype = sig)) +
  geom_point(data = filter(figDat, year == "2019"), aes(color = treatment), alpha = 0.5, size = 0.5) +
  geom_text(data = filter(alphasDat2, year == "2019"), 
            aes(label = paste("alpha", " == ", param, sep = ""), 
                color = treatment), parse = T, hjust = 1, vjust = 1, show.legend = F, size = textSize) +
  geom_text(data = filter(yearText, year == "2019"), aes(label = year), color = "black", size = textSize, hjust = 0) +
  facet_grid(rows = vars(focal),
             cols = vars(background),
             scales = "free",
             switch = "both") +
  scale_linetype_manual(values = c("dashed", "solid"), guide = F) +
  scale_color_manual(values = col_pal, name = "Treatment") +
  scale_fill_manual(values = col_pal, name = "Treatment") +
  xlab(expression(paste("Competitor density (", m^-2, ")", sep = ""))) +
  ylab("Plant growth") +
  fig_theme +
  theme(axis.title.y = element_blank(),
        strip.text.y = element_blank())

# legend
leg <- get_legend(pairD1Fig)

# combine plots
combFig <- plot_grid(pairD1Fig + theme(legend.position = "none"), pairD2Fig,
                     nrow = 1,
                     labels = LETTERS[1:2],
                     rel_widths = c(1, 0.9),
                     label_x = c(0, -0.01))

# combine
pdf("output/focal_growth_pairwise_figure_2018_2019_density_exp.pdf", width = 7, height = 3.5)
plot_grid(combFig, leg,
          nrow = 2,
          rel_heights = c(0.8, 0.06))
dev.off()

# save alpha data
write_csv(alphasDat, "output/focal_growth_pairwise_coefficients_2018_2019_density_exp.csv")
