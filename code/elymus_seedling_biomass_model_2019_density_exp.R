##### info ####

# file: elymus_seedling_biomass_model_2019_density_exp
# author: Amy Kendig
# date last edited: 10/6/20
# goal: estimate Elymus seedling biomass based on density and fungicide treatments


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


#### initial visualizations ####

# modify dataset so that zero plots are repeated
vizDat <- bioD2Dat %>%
  filter(!is.na(weight) & ID != "A") %>%
  rename(bio.g = weight) %>%
  left_join(plotsD2) %>%
  mutate(log_bio.g = log(bio.g),
         treatment = recode(treatment, water = "control"))

# non-transformed biomass
ggplot(vizDat, aes(background_density, bio.g, color = treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  facet_wrap(~background, scales = "free_x") +
  theme_bw()

# log-transformed biomass
ggplot(vizDat, aes(background_density, log_bio.g, color = treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  facet_wrap(~background, scales = "free_x") +
  theme_bw()


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

# remove plots with no background (negative effect on growth)
evSBioD2Dat2 <- evSBioD2Dat %>%
  filter(background != "none")

# initial fit
evSBioMod1 <- brm(bf(log_bio.g ~ logv - log(1 + alphaA * mv_seedling_density + alphaS * ev_seedling_density + alphaP * ev_adult_density),
                     logv ~ fungicide + treatment + (1|site/plotr),
                     alphaA ~ 0 + treatment,
                     alphaS ~ 0 + treatment,
                     alphaP ~ 0 + treatment,
                     nl = T),
                  data = evSBioD2Dat2, family = gaussian,
                  prior <- c(prior(normal(1.5, 10), nlpar = "logv", class = "b", coef = "Intercept"),
                             prior(normal(0.03, 0.001), nlpar = "logv", class = "b", coef = "treatmentfungicide"),
                             prior(normal(0, 10), nlpar = "logv", class = "b"),
                             prior(exponential(0.5), nlpar = "alphaA", lb = 0),
                             prior(exponential(0.5), nlpar = "alphaS", lb = 0),
                             prior(exponential(0.5), nlpar = "alphaP", lb = 0),
                             prior(cauchy(0, 1), nlpar = "logv", class = "sd"),
                             prior(cauchy(0, 1), class = "sigma")),
                  iter = 6000, warmup = 1000, chains = 1)
# 7 divergent transitions
summary(evSBioMod1)

# increase chains and adapt delta
evSBioMod2 <- update(evSBioMod1, chains = 3,
                     control = list(adapt_delta = 0.999, max_treedepth = 15))
summary(evSBioMod2)
plot(evSBioMod2)
pp_check(evSBioMod2, nsamples = 50)


#### visualize ####

# simulation data
simDat <- tibble(mv_seedling_density = c(seq(0, 64, length.out = 100), rep(0, 200)),
                 ev_seedling_density = c(rep(0, 100), seq(0, 16, length.out = 100), rep(0, 100)),
                 ev_adult_density = c(rep(0, 200), seq(0, 8, length.out = 100)),
                 background = rep(c("Mv seedling", "Ev seedling", "Ev adult"), each = 100)) %>%
  expand_grid(fungicide = c(0, 1)) %>%
  mutate(treatment = ifelse(fungicide == 1, "fungicide", "control"),
         site = NA,
         background_density = case_when(background == "Mv seedling" ~ mv_seedling_density,
                                        background == "Ev seedling" ~ ev_seedling_density,
                                        TRUE ~ ev_adult_density))

# simulate fit
fitDat <- simDat %>%
  mutate(bio.g = exp(fitted(evSBioMod2, newdata = ., re_formula = NA)[, "Estimate"]),
         bio.g_lower = exp(fitted(evSBioMod2, newdata = ., re_formula = NA)[, "Q2.5"]),
         bio.g_upper = exp(fitted(evSBioMod2, newdata = ., re_formula = NA)[, "Q97.5"]))

# fit figure
ggplot(fitDat, aes(background_density, bio.g, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = bio.g_lower, ymax = bio.g_upper), alpha = 0.5, color = NA) +
  geom_line(size = 1.5) +
  stat_summary(data = vizDat, geom = "errorbar", width = 0, fun.data = "mean_cl_boot") +
  stat_summary(data = vizDat, geom = "point", size = 2, fun = "mean") +
  facet_wrap(~ background, scales = "free_x") +
  theme_bw()


#### output ####
save(evSBioMod2, file = "output/elymus_seedling_biomass_model_2019_density_exp.rda")
