##### info ####

# file: elymus_seedling_seed_production_model_2019_density_exp
# author: Amy Kendig
# date last edited: 11/12/20
# goal: estimate Elymus seedling seed production based on density and fungicide treatments


#### set up ####

# clear all existing data
rm(list=ls())

# load packages
library(tidyverse)
library(brms)

# import data
seedD2Dat <- read_csv("intermediate-data/ev_processed_seeds_both_year_conversion_2019_density_exp.csv")
bioD2Dat <- read_csv("data/ev_biomass_seeds_oct_2019_density_exp.csv")
plotsD <- read_csv("data/plot_treatments_2018_2019_density_exp.csv")
plotsD2 <- read_csv("data/plot_treatments_for_analyses_2018_2019_density_exp.csv")


#### edit data ####

# plant group densities
plotDens <- plotsD %>%
  mutate(mv_seedling_density = case_when(background == "Mv seedling" ~ background_density, TRUE ~ 0),
         ev_seedling_density = case_when(background == "Ev seedling" ~ background_density, TRUE ~ 0),
         ev_adult_density = case_when(background == "Ev adult" ~ background_density, TRUE ~ 0))

# check for notes about seeds
unique(seedD2Dat$spikelet_notes)
filter(seedD2Dat, is.na(spikelet_weight.g))
# no missing data

# add columns
# remove missing data
evSSeedD2Dat <- seedD2Dat %>%
  group_by(site, plot, treatment, sp, ID) %>%
  summarise(seeds = sum(seeds)) %>%
  ungroup() %>%
  full_join(bioD2Dat %>%
              select(site, plot, treatment, sp, ID)) %>%
  left_join(plotDens) %>%
  mutate(seeds = replace_na(seeds, 0),
         log_seeds = log(seeds + 1),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotr = ifelse(treatment == "fungicide", plot + 10, plot),
         age = ifelse(ID == "A", "adult", "seedling"),
         treatment = ifelse(treatment == "water", "control", treatment)) %>%
  filter(!is.na(seeds) & age == "seedling")


#### initial visualizations ####

# modify dataset so that zero plots are repeated
vizDat <- evSSeedD2Dat %>%
  mutate(treatment = ifelse(treatment == "control", "water", treatment)) %>%
  select(-c(background, background_sp)) %>%
  left_join(plotsD2)

# non-transformed biomass
ggplot(vizDat, aes(background_density, seeds, color = treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  facet_wrap(~background, scales = "free_x") +
  theme_bw()

# log-transformed biomass
ggplot(vizDat, aes(background_density, log_seeds, color = treatment)) +
  stat_summary(geom = "point", fun = "mean", size = 2) +
  stat_summary(geom = "errorbar", fun.data = "mean_cl_boot", width = 0) +
  facet_wrap(~background, scales = "free_x") +
  theme_bw()


#### fit regression ####

# remove plots with no background (negative effect on growth)
evSSeedD2Dat2 <- evSSeedD2Dat %>%
  filter(background != "none")

# initial fit
evSSeedD2Mod1 <- brm(bf(log_seeds ~ logS - log(1 + alphaA * mv_seedling_density + alphaS * ev_seedling_density + alphaP * ev_adult_density),
                        logS ~ treatment + (1|site),
                        alphaA ~ 0 + treatment,
                        alphaS ~ 0 + treatment,
                        alphaP ~ 0 + treatment,
                     nl = T),
                  data = evSSeedD2Dat2, family = gaussian,
                  prior <- c(prior(normal(3, 10), nlpar = "logS", class = "b", coef = "Intercept"),
                             prior(normal(0, 10), nlpar = "logS", class = "b"),
                             prior(exponential(0.5), nlpar = "alphaA", lb = 0),
                             prior(exponential(0.5), nlpar = "alphaS", lb = 0),
                             prior(exponential(0.5), nlpar = "alphaP", lb = 0),
                             prior(cauchy(0, 1), nlpar = "logS", class = "sd"),
                             prior(cauchy(0, 1), class = "sigma")),
                  iter = 6000, warmup = 1000, chains = 1)
# 3 divergent transitions
summary(evSSeedD2Mod1)

# increase chains and adapt delta
evSSeedD2Mod2 <- update(evSSeedD2Mod1, chains = 3,
                     control = list(adapt_delta = 0.999))
summary(evSSeedD2Mod2)
plot(evSSeedD2Mod2)
pp_check(evSSeedD2Mod2, nsamples = 50)


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
  mutate(log_seeds = fitted(evSSeedD2Mod2, newdata = ., re_formula = NA)[, "Estimate"],
         log_seeds_lower = fitted(evSSeedD2Mod2, newdata = ., re_formula = NA)[, "Q2.5"],
         log_seeds_upper = fitted(evSSeedD2Mod2, newdata = ., re_formula = NA)[, "Q97.5"],
         treatment = ifelse(treatment == "control", "water", treatment))

# change levels on raw data
evSSeedD2Dat3 <- evSSeedD2Dat2 %>%
  mutate(treatment = ifelse(treatment == "control", "water", treatment))

# fit figure
(evSSeedD2Plot <- ggplot(fitDat, aes(background_density, log_seeds, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = log_seeds_lower, ymax = log_seeds_upper), alpha = 0.5, color = NA) +
  geom_line(size = 1.5) +
  stat_summary(data = evSSeedD2Dat3, geom = "errorbar", width = 0, fun.data = "mean_cl_boot") +
  stat_summary(data = evSSeedD2Dat3, geom = "point", size = 2, fun = "mean") +
  facet_wrap(~ background, scales = "free_x") +
  theme_bw())


#### output ####
save(evSSeedD2Mod2, file = "output/elymus_seedling_seed_model_2019_density_exp.rda")
save(evSSeedD2Plot, file = "output/elymus_seedling_seed_figure_2019_density_exp.rda")
