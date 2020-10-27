##### info ####

# file: elymus_adult_biomass_fung_model_2019_density_exp
# author: Amy Kendig
# date last edited: 10/27/20
# goal: estimate Elymus adult biomass based on fungicide treatment


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)

# import data
bioD2Dat <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")
plotsD <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")


#### edit data ####

# plant group densities
plotDens <- plotsD %>%
  mutate(mv_seedling_density = case_when(background == "Mv seedling" ~ background_density, TRUE ~ 0),
         ev_seedling_density = case_when(background == "Ev seedling" ~ background_density, TRUE ~ 0),
         ev_adult_density = case_when(background == "Ev adult" ~ background_density, TRUE ~ 0))

# biomass
# remove missing data
# select Elymus seedlings
# add columns
evABioD2Dat <- bioD2Dat %>%
  filter(!is.na(weight) & ID == "A") %>%
  rename(bio.g = weight) %>%
  left_join(plotDens) %>%
  mutate(fungicide = ifelse(treatment == "fungicide", 1, 0),
         log_bio.g = log(bio.g),
         treatment = recode(treatment, water = "control"))

#### greenhouse regression ####

# fit in elymus_seedling_biomass_model_2019_density_exp.R
# biomass change
# 0.03


#### fit regression ####

# initial fit
evABioFuMod1 <- brm(log_bio.g ~ fungicide + treatment + (1|site),
                  data = evABioD2Dat, family = gaussian,
                  prior <- c(prior(normal(2, 10), class = "Intercept"),
                             prior(normal(0.03, 0.001), class = "b", coef = "treatmentfungicide"),
                             prior(normal(0, 10), class = "b"),
                             prior(cauchy(0, 1), class = "sd"),
                             prior(cauchy(0, 1), class = "sigma")),
                  iter = 6000, warmup = 1000, chains = 1)
# 2 divergent transition
summary(evABioFuMod1)

# increase chains and adapt delta
evABioFuMod2 <- update(evABioFuMod1, chains = 3,
                     control = list(adapt_delta = 0.999))
summary(evABioFuMod2)
plot(evABioFuMod2)
pp_check(evABioFuMod2, nsamples = 50)


#### visualize ####

# simulate fit
fitDat <- tibble(fungicide = c(0, 1), treatment = c("control", "fungicide")) %>%
  mutate(log_bio.g = fitted(evABioFuMod2, newdata = ., re_formula = NA)[, "Estimate"],
         log_bio.g_lower = fitted(evABioFuMod2, newdata = ., re_formula = NA)[, "Q2.5"],
         log_bio.g_upper = fitted(evABioFuMod2, newdata = ., re_formula = NA)[, "Q97.5"],
         treatment = ifelse(treatment == "control", "water", treatment))

# change levels on raw data
evABioD2Dat2 <- evABioD2Dat %>%
  mutate(treatment = ifelse(treatment == "control", "water", treatment)) %>%
  group_by(site, treatment) %>%
  summarise(log_bio.g = mean(log_bio.g)) %>%
  ungroup()

# fit figure
evABioFuPlot <- ggplot(fitDat, aes(treatment, log_bio.g, color = treatment)) +
  geom_point(data = evABioD2Dat2, shape = 15, size = 2, position = position_jitter(width = 0.1, height = 0), alpha = 0.5) +
  geom_errorbar(aes(ymin = log_bio.g_lower, ymax = log_bio.g_upper), width = 0) +
  geom_point(size = 4, shape = 16) +
  theme_bw()



#### output ####
save(evABioFuMod2, file = "output/elymus_adult_biomass_fung_model_2019_density_exp.rda")
save(evABioFuPlot, file = "output/elymus_adult_biomass_fung_figure_2019_density_exp.rda")
