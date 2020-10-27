##### info ####

# file: elymus_adult_seed_production_model_2019_density_exp
# author: Amy Kendig
# date last edited: 10/26/20
# goal: estimate Elymus adult seed production based on density and fungicide treatments


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
unique(bioD2Dat$processing_notes)

# add columns
# remove missing data
evASeedD2Dat <- seedD2Dat %>%
  group_by(site, plot, treatment, sp, ID) %>%
  summarise(seeds = sum(seeds)) %>%
  ungroup() %>%
  full_join(bioD2Dat %>%
              select(site, plot, treatment, sp, ID)) %>%
  left_join(plotDens) %>%
  mutate(seeds = replace_na(seeds, 0),
         log_seeds = log(seeds),
         fungicide = ifelse(treatment == "fungicide", 1, 0),
         plotr = ifelse(treatment == "fungicide", plot + 10, plot),
         age = ifelse(ID == "A", "adult", "seedling"),
         treatment = ifelse(treatment == "water", "control", treatment)) %>%
  filter(!is.na(seeds) & age == "adult")


#### initial visualizations ####

# modify dataset so that zero plots are repeated
vizDat <- evASeedD2Dat %>%
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
evASeedD2Dat2 <- evASeedD2Dat %>%
  filter(background != "none")

# initial fit
evASeedD2Mod1 <- brm(bf(seeds ~ maxS/(1 + gammaA * mv_seedling_density + gammaS * ev_seedling_density + gammaP * ev_adult_density),
                        maxS ~ treatment + (1|site),
                        gammaA ~ 0 + treatment,
                        gammaS ~ 0 + treatment,
                        gammaP ~ 0 + treatment,
                     nl = T),
                  data = evASeedD2Dat2, family = gaussian,
                  prior <- c(prior(normal(200, 50), nlpar = "maxS", class = "b", coef = "Intercept"),
                             prior(normal(0, 10), nlpar = "maxS", class = "b"),
                             prior(exponential(0.5), nlpar = "gammaA", lb = 0),
                             prior(exponential(0.5), nlpar = "gammaS", lb = 0),
                             prior(exponential(0.5), nlpar = "gammaP", lb = 0),
                             prior(cauchy(0, 1), nlpar = "maxS", class = "sd"),
                             prior(cauchy(0, 1), class = "sigma")),
                  iter = 6000, warmup = 1000, chains = 1)
# 108 divergent transitions
summary(evASeedD2Mod1)

# increase chains and adapt delta
evASeedD2Mod2 <- update(evASeedD2Mod1, chains = 3,
                     control = list(adapt_delta = 0.99999))
summary(evASeedD2Mod2)
plot(evASeedD2Mod2)
pp_check(evASeedD2Mod2, nsamples = 50)


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
  mutate(seeds = fitted(evASeedD2Mod2, newdata = ., re_formula = NA)[, "Estimate"],
         seeds_lower = fitted(evASeedD2Mod2, newdata = ., re_formula = NA)[, "Q2.5"],
         seeds_upper = fitted(evASeedD2Mod2, newdata = ., re_formula = NA)[, "Q97.5"],
         treatment = ifelse(treatment == "control", "water", treatment))

# change levels on raw data
evASeedD2Dat3 <- evASeedD2Dat2 %>%
  mutate(treatment = ifelse(treatment == "control", "water", treatment))

# fit figure
evASeedD2Plot <- ggplot(fitDat, aes(background_density, seeds, color = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = seeds_lower, ymax = seeds_upper), alpha = 0.5, color = NA) +
  geom_line(size = 1.5) +
  stat_summary(data = evASeedD2Dat3, geom = "errorbar", width = 0, fun.data = "mean_cl_boot") +
  stat_summary(data = evASeedD2Dat3, geom = "point", size = 2, fun = "mean") +
  facet_wrap(~ background, scales = "free_x") +
  theme_bw()


#### output ####
save(evASeedD2Mod2, file = "output/elymus_adult_seed_model_2019_density_exp.rda")
save(evASeedD2Plot, file = "output/elymus_adult_seed_figure_2019_density_exp.rda")
