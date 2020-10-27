##### info ####

# file: microstegium_biomass_fung_model_2019_density_exp
# author: Amy Kendig
# date last edited: 10/27/20
# goal: estimate Microstegium biomass based on fungicide treatment


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)

# import data
bioD2Dat <- read_csv("data/mv_biomass_seeds_2019_density_exp.csv")
fungDat <- read_csv("data/mv_biomass_seeds_height_jun_2019_fungicide_exp.csv")
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
mvBioD2Dat <- bioD2Dat %>%
  filter(!is.na(biomass_weight.g)) %>%
  rename(bio.g = biomass_weight.g) %>%
  left_join(plotDens) %>%
  mutate(fungicide = ifelse(treatment == "fungicide", 1, 0),
         log_bio.g = log(bio.g),
         treatment = recode(treatment, water = "control"),
         plotr = ifelse(treatment == "fungicide", plot + 10, plot))

# fungicide experiment
# remove dead plant
# fungicide column
unique(fungDat$notes)

mvFungDat <- fungDat %>%
  filter(is.na(notes)) %>%
  mutate(fungicide = ifelse(treatment == "fungicide", 1, 0),
         log_bio.g = log(weight.g))


#### fit greenhouse regression ####

# initial fit
mvBioGhMod1 <- brm(log_bio.g ~ fungicide,
                   data = mvFungDat, family = gaussian,
                   prior <- c(prior(normal(0, 10), class = Intercept),
                              prior(normal(0, 10), class = b),
                              prior(cauchy(0, 1), class = sigma)),
                   iter = 6000, warmup = 1000, chains = 1)
summary(mvBioGhMod1)

# increase chains
mvBioGhMod2 <- update(mvBioGhMod1, chains = 3)
summary(mvBioGhMod2)
plot(mvBioGhMod2)
pp_check(mvBioGhMod2, nsamples = 50)

# biomass change
# -0.16


#### fit regression ####

# initial fit
mvBioFuMod1 <- brm(log_bio.g ~ fungicide + treatment + (1|site/plotr),
                    data = mvBioD2Dat, family = gaussian,
                    prior <- c(prior(normal(4, 10), class = "Intercept"),
                               prior(normal(-0.16, 0.001), class = "b", coef = "treatmentfungicide"),
                               prior(normal(0, 10), class = "b"),
                               prior(cauchy(0, 1), class = "sd"),
                               prior(cauchy(0, 1), class = "sigma")),
                    iter = 6000, warmup = 1000, chains = 1)

# 8 divergent transitions
summary(mvBioFuMod1)

# increase chains and adapt delta
mvBioFuMod2 <- update(mvBioFuMod1, chains = 3,
                    control = list(adapt_delta = 0.999999, max_treedepth = 15))
summary(mvBioFuMod2)
plot(mvBioFuMod2)
pp_check(mvBioFuMod2, nsamples = 50)


#### visualize ####

# simulate fit
fitDat <- tibble(fungicide = c(0, 1), treatment = c("control", "fungicide")) %>%
  mutate(log_bio.g = fitted(mvBioFuMod2, newdata = ., re_formula = NA)[, "Estimate"],
         log_bio.g_lower = fitted(mvBioFuMod2, newdata = ., re_formula = NA)[, "Q2.5"],
         log_bio.g_upper = fitted(mvBioFuMod2, newdata = ., re_formula = NA)[, "Q97.5"],
         treatment = ifelse(treatment == "control", "water", treatment))

# change levels on raw data
mvBioD2Dat2 <- mvBioD2Dat %>%
  mutate(treatment = ifelse(treatment == "control", "water", treatment)) %>%
  group_by(site, treatment) %>%
  summarise(log_bio.g = mean(log_bio.g)) %>%
  ungroup()

# fit figure
mvBioFuPlot <- ggplot(fitDat, aes(treatment, log_bio.g, color = treatment)) +
  geom_point(data = mvBioD2Dat2, shape = 15, size = 2, position = position_jitter(width = 0.1, height = 0), alpha = 0.5) +
  geom_errorbar(aes(ymin = log_bio.g_lower, ymax = log_bio.g_upper), width = 0) +
  geom_point(size = 4, shape = 16) +
  theme_bw()

# fungicide figure
mvFungPlot <- ggplot(mvFungDat, aes(treatment, log_bio.g, color = treatment)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_cl_boot") +
  stat_summary(geom = "point", size = 2, fun = "mean") +
  theme_bw() +
  theme(legend.position = "none")


#### output ####
save(mvBioFuMod2, file = "output/microstegium_biomass_fung_model_2019_density_exp.rda")
save(mvBioFuPlot, file = "output/microstegium_biomass_fung_figure_2019_density_exp.rda")
save(mvFungPlot, file = "output/microstegium_biomass_figure_2019_greenhouse_exp.rda")
