##### info ####

# file: size_no_background_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 3/16/21
# goal: plot-level biomass and density relationships


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)
library(tidybayes) # for mean_hdi

# import data
growthD1Dat <- read_csv("intermediate-data/focal_processed_growth_2018_density_exp.csv")
# focal_growth_data_processing_2018_density_exp
mvBioD2Dat <- read_csv("data/mv_biomass_seeds_2019_density_exp.csv")
evBioD2Dat <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")


#### edit data ####

# year 1
growthD1Dat2 <- growthD1Dat %>%
  filter(plot == 1) %>%
  mutate(age = ifelse(ID == "A", "adult", "seedling"),
         plant_group = paste(sp, age, sep = "_") %>%
           fct_rev(),
         fungicide = ifelse(treatment == "fungicide", 1, 0))

# units
unique(growthD1Dat2$days)
# per month

# year 2
growthD2Dat <- mvBioD2Dat %>%
  mutate(ID = as.character(plant)) %>%
  full_join(evBioD2Dat %>%
              rename(biomass_weight.g = weight)) %>%
  mutate(age = ifelse(ID == "A", "adult", "seedling"),
         plant_group = paste(sp, age, sep = "_") %>%
           fct_rev(),
         log_bio = log(biomass_weight.g),
         fungicide = ifelse(treatment == "fungicide", 1, 0)) %>%
  filter(!is.na(biomass_weight.g))


#### initial visualizations ####

ggplot(growthD1Dat2, aes(x = tiller_growth)) +
  geom_density() +
  facet_wrap(~ plant_group)

ggplot(growthD1Dat2, aes(x = height_growth)) +
  geom_density() +
  facet_wrap(~ plant_group)

ggplot(growthD2Dat, aes(x = log_bio)) +
  geom_density() +
  facet_wrap(~ plant_group)


#### year 1 model ####

growthD1Dat2 %>%
  filter(plant_group == "Mv_seedling" & fungicide == 0) %>%
  summarise(growth = mean(tiller_growth))

growthD1Mod <- brm(data = growthD1Dat2, family = gaussian,
                   tiller_growth ~ fungicide * plant_group + (1|site),
                   prior <- c(prior(normal(1.6, 1), class = "Intercept"),
                              prior(normal(0, 1), class = "b")), # use default for sigma
                   iter = 6000, warmup = 1000, chains = 3,
                   control = list(adapt_delta = 0.999))

summary(growthD1Mod)
pp_check(growthD1Mod, nsamples = 100)
plot(growthD1Mod)


#### year 1 model ####

growthD2Dat %>%
  filter(plant_group == "Mv_seedling" & fungicide == 0) %>%
  summarise(growth = mean(log_bio))

growthD2Mod <- brm(data = growthD2Dat, family = gaussian,
                   log_bio ~ fungicide * plant_group + (1|site),
                   prior <- c(prior(normal(2.5, 1), class = "Intercept"),
                              prior(normal(0, 1), class = "b")), # use default for sigma
                   iter = 6000, warmup = 1000, chains = 3,
                   control = list(adapt_delta = 0.99999))

summary(growthD2Mod)
pp_check(growthD2Mod, nsamples = 100)
plot(growthD2Mod)


#### values for text ####

posterior_samples(growthD1Mod) %>%
  rename(b_fungicide_Ev_adult = "b_fungicide:plant_groupEv_adult",
         b_fungicide_Ev_seedling = "b_fungicide:plant_groupEv_seedling") %>%
  transmute(mvWater = b_Intercept,
            mvFung = b_Intercept + b_fungicide,
            evSWater = b_Intercept + b_plant_groupEv_seedling,
            evSFung = b_Intercept + b_plant_groupEv_seedling + b_fungicide + b_fungicide_Ev_seedling,
            evAWater = b_Intercept + b_plant_groupEv_adult,
            evAFung = b_Intercept + b_plant_groupEv_adult + b_fungicide + b_fungicide_Ev_adult) %>%
  mutate(mvFE = 100 * (mvFung - mvWater) / mvWater,
         evSFE = 100 * (evSFung - evSWater) / evSWater,
         evAFE = 100 * (evAFung - evAWater) / evAWater) %>%
  pivot_longer(cols = everything(),
               names_to = "treatment",
               values_to = "growth") %>%
  group_by(treatment) %>%
  mean_hdi(growth)

posterior_samples(growthD2Mod) %>%
  rename(b_fungicide_Ev_adult = "b_fungicide:plant_groupEv_adult",
         b_fungicide_Ev_seedling = "b_fungicide:plant_groupEv_seedling") %>%
transmute(mvWater = exp(b_Intercept),
          mvFung = exp(b_Intercept + b_fungicide),
          evSWater = exp(b_Intercept + b_plant_groupEv_seedling),
          evSFung = exp(b_Intercept + b_plant_groupEv_seedling + b_fungicide + b_fungicide_Ev_seedling),
          evAWater = exp(b_Intercept + b_plant_groupEv_adult),
          evAFung = exp(b_Intercept + b_plant_groupEv_adult + b_fungicide + b_fungicide_Ev_adult)) %>%
  mutate(mvFE = 100 * (mvFung - mvWater) / mvWater,
         evSFE = 100 * (evSFung - evSWater) / evSWater,
         evAFE = 100 * (evAFung - evAWater) / evAWater) %>%
  pivot_longer(cols = everything(),
               names_to = "treatment",
               values_to = "growth") %>%
  group_by(treatment) %>%
  mean_hdi(growth)


#### output ####

save(growthD1Mod, file = "output/focal_growth_no_background_model_2018_density_exp.rda")
save(growthD2Mod, file = "output/focal_growth_no_background_model_2019_density_exp.rda")
