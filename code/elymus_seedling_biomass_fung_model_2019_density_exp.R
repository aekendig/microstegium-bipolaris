##### info ####

# file: elymus_seedling_biomass_fung_model_2019_density_exp
# author: Amy Kendig
# date last edited: 10/27/20
# goal: estimate Elymus seedling biomass based on and fungicide treatment


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)

# import data
bioD2Dat <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")
fungDat <- read_csv("data/ev_biomass_dec_2019_fungicide_exp.csv")
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
evSBioD2Dat <- bioD2Dat %>%
  filter(!is.na(weight) & ID != "A") %>%
  rename(bio.g = weight) %>%
  left_join(plotDens) %>%
  mutate(fungicide = ifelse(treatment == "fungicide", 1, 0),
         log_bio.g = log(bio.g),
         treatment = recode(treatment, water = "control"),
         plotr = ifelse(treatment == "fungicide", plot + 10, plot))

# fungicide experiment
# combine live and dead weight
# fungicide column
evFungDat <- fungDat %>%
  group_by(treatment, pot, sp) %>%
  summarise(weight.g = sum(weight.g)) %>%
  ungroup() %>%
  mutate(fungicide = ifelse(treatment == "fungicide", 1, 0),
         log_bio.g = log(weight.g))


#### fit greenhouse regression ####

# initial fit
evBioGhMod1 <- brm(log_bio.g ~ fungicide,
                   data = evFungDat, family = gaussian,
                   prior <- c(prior(normal(0, 10), class = Intercept),
                              prior(normal(0, 10), class = b),
                              prior(cauchy(0, 1), class = sigma)),
                   iter = 6000, warmup = 1000, chains = 1)
summary(evBioGhMod1)

# increase chains
evBioGhMod2 <- update(evBioGhMod1, chains = 3)
summary(evBioGhMod2)
plot(evBioGhMod2)
pp_check(evBioGhMod2, nsamples = 50)

# biomass change
evFungDat %>%
  group_by(fungicide) %>%
  summarise(log_mean = mean(log_bio.g),
            mean = mean(weight.g))
# log-scale: add 0.03 with fungicide
# non-transformed scale: multiply by 1.03 (e^0.03) with fungicide


#### fit regression ####

# initial fit
evSBioFuMod1 <- brm(log_bio.g ~ fungicide + treatment + (1|site/plotr),
                    data = evSBioD2Dat, family = gaussian,
                    prior <- c(prior(normal(1.5, 10), class = "Intercept"),
                               prior(normal(0.03, 0.001), class = "b", coef = "treatmentfungicide"),
                               prior(normal(0, 10), class = "b"),
                               prior(cauchy(0, 1), class = "sd"),
                               prior(cauchy(0, 1), class = "sigma")),
                    iter = 6000, warmup = 1000, chains = 1)
# 16 divergent transitions
summary(evSBioFuMod1)

# increase chains and adapt delta
evSBioFuMod2 <- update(evSBioFuMod1, chains = 3,
                     control = list(adapt_delta = 0.999, max_treedepth = 15))
summary(evSBioFuMod2)
plot(evSBioFuMod2)
pp_check(evSBioFuMod2, nsamples = 50)


#### visualize ####

# simulate fit
fitDat <- tibble(fungicide = c(0, 1), treatment = c("control", "fungicide")) %>%
  mutate(log_bio.g = fitted(evSBioFuMod2, newdata = ., re_formula = NA)[, "Estimate"],
         log_bio.g_lower = fitted(evSBioFuMod2, newdata = ., re_formula = NA)[, "Q2.5"],
         log_bio.g_upper = fitted(evSBioFuMod2, newdata = ., re_formula = NA)[, "Q97.5"],
         treatment = ifelse(treatment == "control", "water", treatment))

# change levels on raw data
evSBioD2Dat2 <- evSBioD2Dat %>%
  mutate(treatment = ifelse(treatment == "control", "water", treatment)) %>%
  group_by(site, treatment) %>%
  summarise(log_bio.g = mean(log_bio.g)) %>%
  ungroup()

# fit figure
evSBioFuPlot <- ggplot(fitDat, aes(treatment, log_bio.g, color = treatment)) +
  geom_point(data = evSBioD2Dat2, shape = 15, size = 2, position = position_jitter(width = 0.1, height = 0), alpha = 0.5) +
  geom_errorbar(aes(ymin = log_bio.g_lower, ymax = log_bio.g_upper), width = 0) +
  geom_point(size = 4, shape = 16) +
  theme_bw()

# fungicide figure
evFungPlot <- ggplot(evFungDat, aes(treatment, log_bio.g, color = treatment)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_cl_boot") +
  stat_summary(geom = "point", size = 2, fun = "mean") +
  theme_bw() +
  theme(legend.position = "none")


#### output ####
save(evSBioFuMod2, file = "output/elymus_seedling_biomass_fung_model_2019_density_exp.rda")
save(evSBioFuPlot, file = "output/elymus_seedling_biomass_fung_figure_2019_density_exp.rda")
save(evFungPlot, file = "output/elymus_biomass_figure_2019_greenhouse_exp.rda")
