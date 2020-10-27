##### info ####

# file: microstegium_biomass_model_2019_density_exp
# author: Amy Kendig
# date last edited: 10/26/20
# goal: estimate Microstegium biomass based on density and fungicide treatments


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
plotsD2 <- read_csv("data/plot_treatments_for_analyses_2018_2019_density_exp.csv")


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


#### initial visualizations ####

# modify dataset so that zero plots are repeated
vizDat <- bioD2Dat %>%
  filter(!is.na(biomass_weight.g)) %>%
  rename(bio.g = biomass_weight.g) %>%
  left_join(plotsD2) %>%
  mutate(log_bio.g = log(bio.g))

# log-transformed biomass
ggplot(vizDat, aes(background_density, log_bio.g, color = treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  facet_wrap(~background, scales = "free_x") +
  theme_bw()


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

# remove plots with no background (negative effect on growth)
mvBioD2Dat2 <- mvBioD2Dat %>%
  filter(background != "none")

# initial fit
mvBioMod1 <- brm(bf(log_bio.g ~ logv - log(1 + alphaA * mv_seedling_density + alphaS * ev_seedling_density + alphaP * ev_adult_density),
                     logv ~ fungicide + treatment + (1|site/plotr),
                     alphaA ~ 0 + treatment,
                     alphaS ~ 0 + treatment,
                     alphaP ~ 0 + treatment,
                     nl = T),
                  data = mvBioD2Dat2, family = gaussian,
                  prior <- c(prior(normal(4, 10), nlpar = "logv", class = "b", coef = "Intercept"),
                             prior(normal(-0.16, 0.001), nlpar = "logv", class = "b", coef = "treatmentfungicide"),
                             prior(normal(0, 10), nlpar = "logv", class = "b"),
                             prior(exponential(0.5), nlpar = "alphaA", lb = 0),
                             prior(exponential(0.5), nlpar = "alphaS", lb = 0),
                             prior(exponential(0.5), nlpar = "alphaP", lb = 0),
                             prior(cauchy(0, 1), nlpar = "logv", class = "sd"),
                             prior(cauchy(0, 1), class = "sigma")),
                  iter = 6000, warmup = 1000, chains = 1)
# 6 divergent transitions
summary(mvBioMod1)

# increase chains and adapt delta
mvBioMod2 <- update(mvBioMod1, chains = 3,
                    control = list(adapt_delta = 0.9999, max_treedepth = 15))
summary(mvBioMod2)
plot(mvBioMod2)
pp_check(mvBioMod2, nsamples = 50)


#### visualize ####

# simulation data
simDat <- tibble(mv_seedling_density = c(seq(0, 64, length.out = 100), rep(0, 200)),
                 ev_seedling_density = c(rep(0, 100), seq(0, 16, length.out = 100), rep(0, 100)),
                 ev_adult_density = c(rep(0, 200), seq(0, 8, length.out = 100)),
                 background = rep(c("Mv seedling", "Ev seedling", "Ev adult"), each = 100)) %>%
  expand_grid(fungicide = c(0, 1)) %>%
  mutate(treatment = ifelse(fungicide == 1, "fungicide", "control"),
         site = NA,
         plotr = NA,
         background_density = case_when(background == "Mv seedling" ~ mv_seedling_density,
                                        background == "Ev seedling" ~ ev_seedling_density,
                                        TRUE ~ ev_adult_density))

# simulate fit
fitDat <- simDat %>%
  mutate(log_bio.g = fitted(mvBioMod2, newdata = ., re_formula = NA)[, "Estimate"],
         log_bio.g_lower = fitted(mvBioMod2, newdata = ., re_formula = NA)[, "Q2.5"],
         log_bio.g_upper = fitted(mvBioMod2, newdata = ., re_formula = NA)[, "Q97.5"],
         treatment = ifelse(treatment == "control", "water", treatment))

# change levels on raw data
mvBioD2Dat3 <- mvBioD2Dat2 %>%
  mutate(treatment = ifelse(treatment == "control", "water", treatment))

# fit figure
mvBioPlot <- ggplot(fitDat, aes(background_density, log_bio.g, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = log_bio.g_lower, ymax = log_bio.g_upper), alpha = 0.5, color = NA) +
  geom_line(size = 1.5) +
  stat_summary(data = mvBioD2Dat3, geom = "errorbar", width = 0, fun.data = "mean_cl_boot") +
  stat_summary(data = mvBioD2Dat3, geom = "point", size = 2, fun = "mean") +
  facet_wrap(~ background, scales = "free_x") +
  theme_bw()

# fungicide figure
mvFungPlot <- ggplot(mvFungDat, aes(treatment, log_bio.g, color = treatment)) +
  stat_summary(geom = "errorbar", width = 0, fun.data = "mean_cl_boot") +
  stat_summary(geom = "point", size = 2, fun = "mean") +
  theme_bw() +
  theme(legend.position = "none")


#### output ####
save(mvBioMod2, file = "output/microstegium_biomass_model_2019_density_exp.rda")
save(mvBioPlot, file = "output/microstegium_biomass_figure_2019_density_exp.rda")
save(mvFungPlot, file = "output/microstegium_biomass_figure_2019_greenhouse_exp.rda")
