##### info ####

# file: microstegium_biomass_fung_model_2018_2019_density_exp
# author: Amy Kendig
# date last edited: 11/11/20
# goal: estimate Microstegium biomass based on fungicide treatment


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)

# import data
bioD1Dat <- read_csv("./data/mv_biomass_oct_2018_density_exp.csv")
bioD2Dat <- read_csv("data/mv_biomass_seeds_2019_density_exp.csv")
plotsD <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")


#### edit data ####

# plant group densities
plotDens <- plotsD %>%
  mutate(mv_seedling_density = case_when(background == "Mv seedling" ~ background_density, TRUE ~ 0),
         ev_seedling_density = case_when(background == "Ev seedling" ~ background_density, TRUE ~ 0),
         ev_adult_density = case_when(background == "Ev adult" ~ background_density, TRUE ~ 0))

# update errors
# make biomass per m^2
mvBioD1Dat <- bioD1Dat %>%
  mutate(site = case_when(processing_notes == "duplicate - one is D2, probably this one because the other was sorted with other D1 bags; Quynh guessed on which seeds go with which biomass (see email)" ~ "D2", 
                          TRUE ~ site),
         treatment = case_when(processing_notes == "duplicate - one is 7F, ID = A (not dated), probably this one" ~ "fungicide",
                               TRUE ~ treatment)) %>%
  left_join(plotDens) %>%
  mutate(fungicide = ifelse(treatment == "fungicide", 1, 0),
         log_bio.g = log(bio.g),
         treatment = recode(treatment, water = "control"),
         plotr = ifelse(treatment == "fungicide", plot + 10, plot),
         yearf = "year 1")

# add columns
mvBioD2Dat <- bioD2Dat %>%
  filter(!is.na(biomass_weight.g)) %>%
  rename(bio.g = biomass_weight.g) %>%
  left_join(plotDens) %>%
  mutate(fungicide = ifelse(treatment == "fungicide", 1, 0),
         log_bio.g = log(bio.g),
         treatment = recode(treatment, water = "control"),
         plotr = ifelse(treatment == "fungicide", plot + 10, plot),
         yearf = "year 2")

# combine data
mvBioDat <- full_join(mvBioD1Dat, mvBioD2Dat)


#### visualize ####

ggplot(mvBioDat, aes(treatment, log_bio.g)) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  facet_wrap(~ yearf)


#### fit greenhouse regression in microstegium_biomass_fung_model_2019_density_exp.R ####


#### fit regression ####

# initial fit
mvBioFuMod1 <- brm(log_bio.g ~ fungicide * yearf + treatment + (1|site/plotr),
                    data = mvBioDat, family = gaussian,
                    prior <- c(prior(normal(3, 10), class = "Intercept"),
                               prior(normal(-0.16, 0.001), class = "b", coef = "treatmentfungicide"),
                               prior(normal(0, 10), class = "b"),
                               prior(cauchy(0, 1), class = "sd"),
                               prior(cauchy(0, 1), class = "sigma")),
                    iter = 6000, warmup = 1000, chains = 1)

# 9 divergent transitions
summary(mvBioFuMod1)

# increase chains and adapt delta
mvBioFuMod2 <- update(mvBioFuMod1, chains = 3,
                    control = list(adapt_delta = 0.99999999, max_treedepth = 15))
summary(mvBioFuMod2)
plot(mvBioFuMod2)
pp_check(mvBioFuMod2, nsamples = 50)


#### visualize ####

# simulate fit
fitDat <- tibble(fungicide = rep(c(0, 1), 2), treatment = rep(c("control", "fungicide"), 2), yearf = rep(c("year 1", "year 2"), each = 2)) %>%
  mutate(log_bio.g = fitted(mvBioFuMod2, newdata = ., re_formula = NA)[, "Estimate"],
         log_bio.g_lower = fitted(mvBioFuMod2, newdata = ., re_formula = NA)[, "Q2.5"],
         log_bio.g_upper = fitted(mvBioFuMod2, newdata = ., re_formula = NA)[, "Q97.5"],
         treatment = ifelse(treatment == "control", "water", treatment))

# change levels on raw data
mvBioDat2 <- mvBioDat %>%
  mutate(treatment = ifelse(treatment == "control", "water", treatment)) %>%
  group_by(yearf, site, treatment) %>%
  summarise(log_bio.g = mean(log_bio.g)) %>%
  ungroup()

# fit figure
(mvBioFuPlot <- ggplot(fitDat, aes(treatment, log_bio.g, color = treatment)) +
  geom_point(data = mvBioDat2, size = 2, position = position_jitter(width = 0.1, height = 0), alpha = 0.5, aes(shape = site)) +
  geom_errorbar(aes(ymin = log_bio.g_lower, ymax = log_bio.g_upper), width = 0) +
  geom_point(size = 4, shape = 16) +
  facet_wrap(~ yearf) +
  theme_bw())


#### output ####
save(mvBioFuMod2, file = "output/microstegium_biomass_fung_model_2018_2019_density_exp.rda")
save(mvBioFuPlot, file = "output/microstegium_biomass_fung_figure_2018_2019_density_exp.rda")
